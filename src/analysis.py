"""Patient-level composition, signatures and association analyses.

No absolute paths, no network operations, no automatic patient deletion.
"""
from pathlib import Path
import json
import warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from statsmodels.duration.hazard_regression import PHReg

def read(path):return pd.read_csv(path,sep='\t')
def write(df,out,name,index=False):df.to_csv(Path(out)/name,sep='\t',index=index)
def require(ok,message):
    if not ok:raise ValueError(message)
def bh(p):
    x=np.asarray(p,dtype=float);q=np.full(x.shape,np.nan);m=np.isfinite(x)
    if m.any():q[m]=multipletests(x[m],method='fdr_bh')[1]
    return q
def boolean(s):
    t=s.astype(str).str.lower()
    require(t.isin(['true','false','1','0']).all(),'Boolean columns must contain true/false or 0/1')
    return t.isin(['true','1'])
def pfs_group(time,event,threshold=6):
    # A later-censored patient is classifiable; early censoring is not progression.
    if pd.isna(time):return 'Missing'
    if time==threshold:return 'Boundary'
    if time>threshold:return 'Long'
    if event==1:return 'Short'
    return 'Early'
def load_counts(path):
    x=read(path).set_index('gene')
    require(x.index.is_unique,'Duplicate gene identifiers')
    a=x.to_numpy(dtype=float)
    require(np.isfinite(a).all() and (a>=0).all() and np.equal(a,np.floor(a)).all(),'Counts must be finite nonnegative integers')
    require(x.columns.is_unique,'Duplicate sample identifiers')
    return x
def logcpm(counts):
    total=counts.sum(axis=0)
    require((total>0).all(),'Zero pseudobulk library')
    return np.log2(counts.div(total,axis=1)*1e6+1)
def prepare(data,out,cfg,config_dir):
    require(cfg['pfs_threshold_months']==6,'This publication renderer supports the six-month threshold only')
    c=read(data/'clinical.tsv')
    required=['sample_id','include_baseline','treatment_scene','batch','response_group','cps','PFS_time_months','PFS_event','OS_time_months','OS_event']
    require(set(required)<=set(c),'clinical.tsv missing required columns')
    require(c.sample_id.notna().all() and c.sample_id.is_unique,'One unique sample per patient required')
    c['include_baseline']=boolean(c.include_baseline)
    require(c.treatment_scene.isin(['Neoadjuvant','Line_1','Line_2plus','Unknown']).all(),'Unexpected treatment scene')
    for ep in ['PFS','OS']:
        require(c[ep+'_event'].dropna().isin([0,1]).all(),'Invalid event indicator')
        require((c[ep+'_time_months'].dropna()>=0).all(),'Negative survival time')
    a=read(data/'annotation.tsv')
    require(set(['feature','retained','non_epithelial','major','label'])<=set(a),'annotation.tsv missing columns')
    require(a.feature.is_unique,'Duplicate annotation feature')
    a['retained']=boolean(a.retained);a['non_epithelial']=boolean(a.non_epithelial)
    n=read(data/'cell_counts.tsv').set_index('sample_id')
    require(n.index.is_unique and set(a.feature)==set(n.columns),'Cell count/annotation mismatch')
    arr=n.to_numpy(dtype=float)
    require(np.isfinite(arr).all() and (arr>=0).all() and np.equal(arr,np.floor(arr)).all(),'Invalid cell counts')
    require(set(c.loc[c.include_baseline,'sample_id'])<=set(n.index),'Missing baseline composition')
    c['composition_available']=c.sample_id.isin(n.index)
    c=c[c.include_baseline & c.composition_available].copy().set_index('sample_id')
    require(len(c)>0,'Empty cohort')
    n=n.loc[c.index,a.feature]
    require((n.sum(axis=1)>0).all(),'Empty cell sample')
    raw=n.div(n.sum(axis=1),axis=0)
    keep=a.loc[a.retained,'feature'];k=n.loc[:,keep]
    require((k.sum(axis=1)>0).all(),'Empty retained composition')
    clr=np.log(k+cfg['clr_pseudocount']);clr=clr.sub(clr.mean(axis=1),axis=0)
    require(np.allclose(clr.sum(axis=1),0,atol=1e-8),'CLR closure check failed')
    c['n_mast_cells']=n[cfg['mast_feature']]
    c['n_tpex_cells']=n[cfg['tpex_feature']]
    c['pfs6']=[pfs_group(t,e,cfg['pfs_threshold_months']) for t,e in zip(c.PFS_time_months,c.PFS_event)]
    c['advanced']=c.treatment_scene.isin(['Line_1','Line_2plus'])
    c['mast_eligible']=c.n_mast_cells>cfg['mast_min_cells_exclusive']
    c['mast_de_eligible']=c.advanced & c.mast_eligible & c.pfs6.isin(['Short','Long'])
    c['mast_survival_eligible']=c.advanced & c.mast_eligible
    c['abundance_corr_eligible']=c.advanced & c.pfs6.isin(['Short','Long'])
    # Preserve cohort order for deterministic analysis, not hand-selected patient IDs.
    mast=load_counts(data/'mast_counts.tsv.gz');tpex=load_counts(data/'tpex_counts.tsv.gz')
    require(set(c.index[c.mast_survival_eligible])<=set(mast.columns),'Eligible Mast libraries missing')
    corr_ids=c.index[c.mast_survival_eligible|c.abundance_corr_eligible]
    require(set(corr_ids)<=set(tpex.columns),'Tpex libraries missing for a required correlation')
    if (c.loc[corr_ids,'n_tpex_cells']==0).any():raise ValueError('No target cells in required expression correlation')
    scaling=read(data/'tpex_scaling.tsv')
    require(set(['signature','gene','mean','sd'])<=set(scaling),'Scaling columns missing')
    require(not scaling.duplicated(['signature','gene']).any() and (scaling.sd>0).all(),'Invalid fixed gene scaling')
    require(set(scaling.gene)<=set(tpex.index),'Signature genes absent from Tpex matrix')
    signatures={}
    expr=logcpm(tpex.loc[:,corr_ids])
    for name,z in scaling.groupby('signature',sort=False):
        z=z.set_index('gene');signatures[name]=expr.loc[z.index].sub(z['mean'],axis=0).div(z.sd,axis=0).mean(axis=0)
    require({'IFN_response','effector_MHCII'}<=set(signatures),'Both Tpex signatures required')
    scores=pd.DataFrame(signatures)
    ml=logcpm(mast.loc[:,c.index[c.mast_survival_eligible]])
    c['IL5']=ml.loc['IL5'];c['CCL2']=ml.loc['CCL2']
    c=c.join(scores)
    for name in ['mast_feature','tpex_feature','trm_feature','cdc1_feature']:
        require(cfg[name] in clr,f'Missing retained feature {cfg[name]}')
        c[name+'_CLR']=clr[cfg[name]]
    c['Mast_proportion']=raw[cfg['mast_feature']]
    write(c.reset_index(),out,'patients.tsv')
    write(a,out,'annotation.tsv');write(raw,out,'raw_proportions.tsv',True);write(clr,out,'clr.tsv',True)
    write(scaling,out,'fixed_Tpex_scaling_used.tsv')
    settings={'mast_min_cells_exclusive':cfg['mast_min_cells_exclusive'],'gsea_seed':cfg['gsea_seed'],
              'gsea_min_size':cfg['gsea_min_size'],'gsea_max_size':cfg['gsea_max_size'],'gsea_eps':cfg['gsea_eps']}
    write(pd.DataFrame(settings.items(),columns=['key','value']),out,'settings.tsv')
    write(pd.DataFrame({'cutoff':cfg['cps_thresholds']}),out,'cps_thresholds.tsv')
    gene_sets=read(config_dir/'mast_gene_sets.tsv')
    require(set(['module','gene','display_label'])<=set(gene_sets),'Missing gene-set definitions')
    write(gene_sets,out,'gene_sets_used.tsv')
    audit=read(data/'clinical.tsv')[['sample_id','include_baseline']]
    audit=audit.merge(c.reset_index()[['sample_id','n_mast_cells','pfs6','advanced','mast_eligible','mast_de_eligible']],on='sample_id',how='left')
    write(audit,out,'eligibility_audit.tsv')
    # Patient-equal major composition (no pseudocount), unlike Cox denominator.
    prop=k.div(k.sum(axis=1),axis=0)
    mapping=a.set_index('feature').major
    major=prop.T.groupby(mapping.loc[prop.columns],sort=False).sum().T
    require(np.allclose(major.sum(axis=1),1),'Major composition closure failed')
    major=major.join(c.treatment_scene)
    write(major.reset_index(),out,'major_patient_proportions.tsv')
    write(major.groupby('treatment_scene').mean(numeric_only=True).reset_index(),out,'major_group_means.tsv')
    return c

def cox_one(df,feature,ep,ties='breslow'):
    q=df[[feature,ep+'_time_months',ep+'_event']].replace([np.inf,-np.inf],np.nan).dropna()
    base={'feature':feature,'endpoint':ep,'n':len(q),'events':int(q[ep+'_event'].sum()),'ties':ties}
    scale=q[feature].quantile(.75)-q[feature].quantile(.25)
    if len(q)<8 or base['events']<3 or q[feature].nunique()<3 or not np.isfinite(scale) or scale<=0:
        return dict(base,status='not_estimable',P=np.nan,HR=np.nan,CI_low=np.nan,CI_high=np.nan,IQR=scale)
    x=(q[feature]-q[feature].median())/scale
    try:
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always')
            fit=PHReg(q[ep+'_time_months'],x.to_numpy()[:,None],status=q[ep+'_event'],ties=ties).fit(disp=False)
        beta=float(fit.params[0]);se=float(fit.bse[0]);p=float(fit.pvalues[0])
        valid=np.isfinite([beta,se,p]).all()
        return dict(base,status='ok' if valid else 'nonfinite',P=p,HR=np.exp(beta),CI_low=np.exp(beta-1.96*se),CI_high=np.exp(beta+1.96*se),IQR=scale,warnings='; '.join(str(w.message) for w in ws))
    except (ValueError,np.linalg.LinAlgError) as e:
        return dict(base,status='fit_failed',P=np.nan,HR=np.nan,CI_low=np.nan,CI_high=np.nan,IQR=scale,warnings=str(e))

def associations(data,out,cfg):
    c=read(out/'patients.tsv').set_index('sample_id');a=read(out/'annotation.tsv')
    raw=read(out/'raw_proportions.tsv').set_index('sample_id');clr=read(out/'clr.tsv').set_index('sample_id')
    advanced=c[boolean(c.advanced)]
    infer=a.loc[boolean(a.retained),'feature']
    rows=[cox_one(advanced.join(raw[[f]]),f,ep) for ep in ['PFS','OS'] for f in infer]
    stats=pd.DataFrame(rows);stats['q_all_retained']=stats.groupby('endpoint').P.transform(bh)
    write(stats,out,'subtype_survival_all.tsv')
    non_epi=set(a.loc[boolean(a.non_epithelial)&boolean(a.retained),'feature'])
    displayed=stats[stats.feature.isin(non_epi)&stats.endpoint.eq('PFS')].sort_values('P').head(12)
    write(displayed,out,'Fig3e_top12_PFS.tsv')
    # All CPS threshold comparisons, with a separate BH family per threshold.
    rows=[]
    for cutoff in cfg['cps_thresholds']:
        ids=c.index[c.cps.notna()]
        for f in infer:
            low=clr.loc[ids[c.loc[ids,'cps']<cutoff],f];high=clr.loc[ids[c.loc[ids,'cps']>=cutoff],f]
            if not len(low) or not len(high):continue
            rows.append({'threshold':cutoff,'feature':f,'n_low':len(low),'n_high':len(high),
                         'effect_high_minus_low':high.median()-low.median(),
                         'P':mannwhitneyu(high,low,alternative='two-sided',method='asymptotic').pvalue})
    cps=pd.DataFrame(rows);cps['q']=cps.groupby('threshold').P.transform(bh);write(cps,out,'CPS_composition.tsv')
    # No minimum-Mast gate in the original 21-person abundance correlations.
    rows=[]
    for cohort,mask,x,features in [
      ('abundance_classifiable',boolean(c.abundance_corr_eligible),'mast_feature_CLR',['tpex_feature_CLR','trm_feature_CLR','cdc1_feature_CLR','effector_MHCII']),
      (f"Mast_gt{cfg['mast_min_cells_exclusive']}_all",boolean(c.mast_survival_eligible),'IL5',['tpex_feature_CLR','IFN_response','effector_MHCII'])]:
        q=c.loc[mask]
        for y in features:
            pairs=q[[x,y]].dropna()
            require(len(pairs)==len(q),'Missing correlation values: no silent patient dropping')
            rho,p=spearmanr(pairs[x],pairs[y])
            rows.append({'cohort':cohort,'x':x,'y':y,'n':len(pairs),'rho':rho,'P_two_sided':p})
    write(pd.DataFrame(rows),out,'correlations.tsv')
    q=c[boolean(c.mast_de_eligible)];rows=[]
    for gene in ['IL5','CCL2']:
        short=q.loc[q.pfs6.eq('Short'),gene];long=q.loc[q.pfs6.eq('Long'),gene]
        res=mannwhitneyu(short,long,alternative='greater',method='asymptotic')
        rows.append({'feature':gene,'n_short':len(short),'n_long':len(long),'U':res.statistic,
                     'P_one_sided':res.pvalue,'P_two_sided':mannwhitneyu(short,long,alternative='two-sided',method='asymptotic').pvalue})
    write(pd.DataFrame(rows),out,'Mast_gene_boxplot_statistics.tsv')
    optional=[]
    if (data/'cd8_program_scores.tsv').exists():
        scores=read(data/'cd8_program_scores.tsv')
        require(set(['sample_id','program','score','n_cells'])<=set(scores),'Invalid CD8 score table')
        require(not scores.duplicated(['sample_id','program']).any(),'Duplicate CD8 patient/program rows')
        rows=[]
        for program,d in scores.groupby('program'):
            q=advanced.join(d[d.n_cells>=20].set_index('sample_id')[['score']],how='inner')
            for ep in ['PFS','OS']:
                row=cox_one(q,'score',ep,ties='efron');row['program']=program;rows.append(row)
        ss=pd.DataFrame(rows);ss['q']=ss.groupby('endpoint').P.transform(bh);write(ss,out,'CD8_program_survival.tsv')
    else:optional.append('FigS3h skipped: cd8_program_scores.tsv not supplied (upstream scoring not claimed).')
    if (data/'external_mast_expression.tsv').exists():
        d=read(data/'external_mast_expression.tsv');genes=['MS4A2','KIT','HDC','LTC4S','IL1RL1','CMA1']
        require(set(genes+['sample_id','response_group'])<=set(d),'External table requires six raw TPM columns')
        require(d.sample_id.is_unique and d.response_group.isin(['R','NR']).all(),'Invalid external patient labels')
        require(np.isfinite(d[genes]).all().all() and (d[genes]>=0).all().all(),'Invalid TPM')
        d['score']=np.log2(d[genes]+1).mean(axis=1)
        r=d.loc[d.response_group.eq('R'),'score'];nr=d.loc[d.response_group.eq('NR'),'score']
        stat=mannwhitneyu(nr,r,alternative='greater',method='asymptotic')
        write(d,out,'external_Mast_source.tsv');write(pd.DataFrame([{'n_R':len(r),'n_NR':len(nr),'U':stat.statistic,'P_one_sided':stat.pvalue}]),out,'external_Mast_statistics.tsv')
    else:optional.append('Fig3j skipped: external_mast_expression.tsv not supplied.')
    (out/'optional_stages.txt').write_text('\n'.join(optional) or 'All optional association stages ran.')
