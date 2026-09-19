#!/usr/bin/env Rscript
args<-commandArgs(TRUE);stopifnot(length(args)==3)
input<-args[1];out<-args[2];config<-args[3]
options(mc.cores=1)
suppressPackageStartupMessages({library(data.table);library(edgeR);library(fgsea);library(BiocParallel)})
register(SerialParam())
wr<-function(x,name)fwrite(x,file.path(out,name),sep='\t')
settings<-fread(file.path(out,'settings.tsv'));s<-setNames(settings$value,settings$key)
patients<-fread(file.path(out,'patients.tsv'))
m<-patients[as.logical(mast_de_eligible)]
stopifnot(nrow(m)>0,!anyDuplicated(m$sample_id),all(m$n_mast_cells>s['mast_min_cells_exclusive']))
m[,pfs6:=factor(pfs6,levels=c('Long','Short'))];m[,batch:=factor(batch)]
stopifnot(all(table(m$pfs6)>=2),nlevels(m$batch)>=2)
ct<-fread(file.path(input,'mast_counts.tsv.gz'))
counts<-as.matrix(ct[,m$sample_id,with=FALSE]);rownames(counts)<-ct$gene
design<-model.matrix(~batch+pfs6,m)
stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design))
wr(data.table(sample_id=m$sample_id,design),'mast_design.tsv')
y<-DGEList(counts);keep<-filterByExpr(y,design)
wr(data.table(gene=rownames(y),retained=keep),'mast_filter_audit.tsv')
y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
test<-glmQLFTest(fit,coef='pfs6Short')
deg<-as.data.table(topTags(test,n=Inf,sort.by='none')$table,keep.rownames='gene')
wr(deg,'mast_DEG.tsv')
ranks<-sort(setNames(sign(deg$logFC)*sqrt(deg$F),deg$gene),decreasing=TRUE)
wr(data.table(gene=names(ranks),ranking_metric=unname(ranks)),'mast_ranked_genes.tsv')
defs<-fread(file.path(config,'mast_gene_sets.tsv'));sets<-split(defs$gene,defs$module)
stopifnot(!anyDuplicated(defs[,.(module,gene)]))
set.seed(as.integer(s['gsea_seed']))
fg<-as.data.table(fgseaMultilevel(sets,ranks,minSize=as.integer(s['gsea_min_size']),maxSize=as.integer(s['gsea_max_size']),eps=as.numeric(s['gsea_eps']),BPPARAM=SerialParam()))
fg[,leadingEdge:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))]
fg<-merge(fg,unique(defs[,.(module,display_label)]),by.x='pathway',by.y='module',all.x=TRUE,sort=FALSE)
wr(fg,'mast_GSEA.tsv')
wr(data.table(module=names(sets),n_defined=lengths(sets),n_ranked=sapply(sets,function(x)sum(x %in% names(ranks))),tested=names(sets) %in% fg$pathway),'mast_GSEA_coverage.tsv')
# Outcome-informed leading-edge score. Export the actual membership on every run.
le<-fg[pathway=='Th2_myeloid_immunoregulatory_output',leadingEdge]
stopifnot(length(le)==1,nchar(le)>0)
genes<-strsplit(le,';',fixed=TRUE)[[1]]
score_y<-calcNormFactors(DGEList(counts),method='TMM')
lc<-cpm(score_y,log=TRUE,prior.count=2)
z<-t(scale(t(lc[genes,,drop=FALSE])));stopifnot(all(is.finite(z)))
sc<-m[,.(sample_id,pfs6)];sc[,leading_edge_score:=colMeans(z)]
wr(sc,'leading_edge_patient_scores.tsv');wr(data.table(gene=genes),'leading_edge_genes.tsv')
short<-sc[pfs6=='Short',leading_edge_score];long<-sc[pfs6=='Long',leading_edge_score]
wr(data.table(n_short=length(short),n_long=length(long),U=as.numeric(wilcox.test(short,long,exact=FALSE)$statistic),
 P_one_sided=wilcox.test(short,long,alternative='greater',exact=FALSE,correct=TRUE)$p.value,
 P_two_sided=wilcox.test(short,long,alternative='two.sided',exact=FALSE,correct=TRUE)$p.value,
 genes=paste(genes,collapse=';')),'leading_edge_statistics.tsv')
capture.output(sessionInfo(),file=file.path(out,'mast_R_sessionInfo.txt'))
cat('Patients:',nrow(m),' genes:',nrow(deg),' sets:',nrow(fg),'\n')
quit(save='no', status=0)
