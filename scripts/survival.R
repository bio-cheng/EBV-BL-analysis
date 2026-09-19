#!/usr/bin/env Rscript
args<-commandArgs(TRUE);stopifnot(length(args)==3);out<-args[2]
suppressPackageStartupMessages({library(data.table);library(survival)})
wr<-function(x,name)fwrite(x,file.path(out,name),sep='\t')
d<-fread(file.path(out,'patients.tsv'));d<-d[as.logical(advanced)]
cps_cuts<-fread(file.path(out,'cps_thresholds.tsv'))$cutoff
coefficients<-list();diagnostics<-list();fits<-list();km<-list();kmstats<-list();risk<-list()
get_cox<-function(q,variables,ep,name,ties='breslow',scale=TRUE){
 q<-copy(q);needed<-c(variables,paste0(ep,c('_time_months','_event')))
 q<-q[complete.cases(q[,..needed])]
 if(nrow(q)<8 || sum(q[[paste0(ep,'_event')]])<3)return(NULL)
 for(v in variables){
  if(scale){iqr<-IQR(q[[v]]);if(!is.finite(iqr)||iqr<=0)return(NULL);q[[v]]<-q[[v]]/iqr}
 }
 form<-as.formula(paste0('Surv(',ep,'_time_months,',ep,'_event)~',paste(variables,collapse='+')))
 f<-coxph(form,data=q,ties=ties,x=TRUE);ss<-summary(f);ci<-ss$conf.int
 co<-data.table(model=name,endpoint=ep,term=rownames(ci),n=f$n,events=f$nevent,HR=ci[,1],CI_low=ci[,3],CI_high=ci[,4],P=ss$coefficients[,5],ties=ties)
 coefficients[[length(coefficients)+1]]<<-co
 ph<-tryCatch(cox.zph(f)$table,error=function(e)NULL)
 if(!is.null(ph))diagnostics[[length(diagnostics)+1]]<<-data.table(model=name,endpoint=ep,term=rownames(ph),PH_P=ph[,'p'])
 fits[[name]]<<-f
 return(f)
}
joint<-d[complete.cases(d[,.(cps,Mast_proportion,PFS_time_months,PFS_event)])]
f1<-get_cox(joint,'cps','PFS','CPS_only')
f2<-get_cox(joint,c('cps','Mast_proportion'),'PFS','CPS_plus_Mast')
if(!is.null(f1)&&!is.null(f2)){
 lr<-2*(as.numeric(logLik(f2))-as.numeric(logLik(f1)))
 wr(data.table(LR=lr,df=1,P=pchisq(lr,1,lower.tail=FALSE),delta_AIC=AIC(f2)-AIC(f1)),'CPS_Mast_nested_model.tsv')
}
for(ep in c('PFS','OS'))for(cut in cps_cuts){
 q<-copy(d);q[,CPS_high:=as.numeric(cps>=cut)]
 if(uniqueN(q$CPS_high,na.rm=TRUE)==2)get_cox(q,'CPS_high',ep,paste(ep,'CPS',cut,sep='_'),scale=FALSE)
}
make_km<-function(q,value,ep,name,il5=FALSE,cutpoint=NULL){
 q<-copy(q);q[,time:=get(paste0(ep,'_time_months'))];q[,event:=get(paste0(ep,'_event'))]
 q<-q[is.finite(time)&event %in% c(0,1)&is.finite(get(value))]
 if(nrow(q)<2)return(NULL)
 cutoff<-if(is.null(cutpoint))median(q[[value]])else cutpoint
 high<-if(il5)q[[value]]>cutoff else q[[value]]>=cutoff
 rule<-if(il5)'high > cutoff' else 'high >= cutoff'
 if(length(unique(high))<2){high<-q[[value]]>cutoff;rule<-'high > cutoff (tie fallback)'}
 if(length(unique(high))<2)return(NULL)
 q[,group:=ifelse(high,'High','Low')]
 lr<-survdiff(Surv(time,event)~group,data=q)
 kmstats[[length(kmstats)+1]]<<-data.table(model=name,endpoint=ep,n=nrow(q),cutoff=cutoff,rule=rule,n_low=sum(!high),n_high=sum(high),logrank_P=pchisq(lr$chisq,1,lower.tail=FALSE))
 for(g in c('Low','High')){
  z<-q[group==g];f<-survfit(Surv(time,event)~1,data=z)
  km[[length(km)+1]]<<-data.table(model=name,endpoint=ep,group=g,time=c(0,f$time),survival=c(1,f$surv),n_censor=c(0,f$n.censor))
  ticks<-seq(0,ceiling(max(q$time)/6)*6,6)
  risk[[length(risk)+1]]<<-data.table(model=name,endpoint=ep,group=g,time=ticks,n_risk=sapply(ticks,function(t)sum(z$time>=t)))
 }
}
il<-d[as.logical(mast_survival_eligible)]
get_cox(il,'IL5','PFS','Mast_IL5',ties='efron');make_km(il,'IL5','PFS','Mast_IL5',il5=TRUE)
raw<-fread(file.path(out,'raw_proportions.tsv'))
for(feature in intersect(c('Myeloid_Mast_TPSAB1','B_Plasma_IgA','T_Cycling_GINS2'),names(raw))){
 q<-merge(d,raw[,c('sample_id',feature),with=FALSE],by='sample_id')
 get_cox(q,feature,'PFS',feature);make_km(q,feature,'PFS',feature)
 if(feature=='Myeloid_Mast_TPSAB1')for(prob in seq(.3,.7,.05))make_km(q,feature,'PFS',paste0('Mast_quantile_',prob),cutpoint=as.numeric(quantile(q[[feature]],prob)))
}
for(ep in c('PFS','OS'))for(cut in cps_cuts)make_km(d,'cps',ep,paste(ep,'CPS',cut,sep='_'),cutpoint=cut)
wr(rbindlist(coefficients,fill=TRUE),'clinical_Cox.tsv');wr(rbindlist(diagnostics,fill=TRUE),'clinical_PH_tests.tsv')
wr(rbindlist(km,fill=TRUE),'KM_coordinates.tsv');wr(rbindlist(kmstats,fill=TRUE),'KM_statistics.tsv');wr(rbindlist(risk,fill=TRUE),'KM_risk_tables.tsv')
capture.output(sessionInfo(),file=file.path(out,'survival_R_sessionInfo.txt'))
quit(save='no', status=0)
