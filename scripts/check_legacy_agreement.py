#!/usr/bin/env python3
"""Compare private rerun outputs to archived aggregate statistics; no uploads."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
def read(p):return pd.read_csv(p,sep='\t')
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--legacy-root',type=Path,required=True);p.add_argument('--results',type=Path,required=True);a=p.parse_args()
    b=a.legacy_root;o=a.results;checks=[]
    def compare(label,x,y,keys,cols):
        x=x.set_index(keys).sort_index();y=y.set_index(keys).sort_index()
        if not x.index.equals(y.index):raise AssertionError(label+': row identifiers differ')
        for col in cols:
            z=x[col].to_numpy(float);w=y[col].to_numpy(float)
            ok=np.allclose(z,w,rtol=1e-7,atol=1e-10,equal_nan=True)
            checks.append({'check':label+':'+col,'rows':len(x),'passed':bool(ok),'max_absolute_difference':float(np.nanmax(np.abs(z-w)))})
            if not ok:raise AssertionError(label+': '+col)
    source=b/'Mast_gt5_PFS6_DEG_GSEA_20260919/tables'
    compare('DEG',read(o/'mast_DEG.tsv'),read(source/'05_DEG_Batch_adjusted.tsv'),['gene'],['logFC','logCPM','F','PValue','FDR'])
    compare('GSEA',read(o/'mast_GSEA.tsv'),read(source/'07_GSEA_Batch_adjusted.tsv'),['pathway'],['NES','ES','pval','padj','size'])
    old=read(b/'CPS_Mast_joint_PFS_20260910/tables/03_Cox_coefficients.tsv')
    old['term']=old.term.replace({'CPS_IQR':'cps','Mast_IQR':'Mast_proportion'})
    new=read(o/'clinical_Cox.tsv');new=new[new.model.isin(['CPS_only','CPS_plus_Mast'])]
    compare('CPS_Mast',new,old,['model','term'],['n','events','HR','CI_low','CI_high','P'])
    expected={'baseline':29,'advanced':25,'mast_DE':18,'mast_expression':22,'abundance_classifiable':21}
    q=read(o/'patients.tsv')
    actual={'baseline':len(q),'advanced':int(q.advanced.sum()),'mast_DE':int(q.mast_de_eligible.sum()),'mast_expression':int(q.mast_survival_eligible.sum()),'abundance_classifiable':int(q.abundance_corr_eligible.sum())}
    for name,v in expected.items():
        ok=actual[name]==v;checks.append({'check':name,'passed':ok,'expected':v,'actual':actual[name]})
        if not ok:raise AssertionError(name)
    cor=read(o/'correlations.tsv').set_index(['cohort','y'])
    benchmarks=[('abundance_classifiable','trm_feature_CLR',21,-.6233766233766234,.0025346926925616103),
      ('abundance_classifiable','tpex_feature_CLR',21,-.5610389610389611,.00814481615489852),
      ('abundance_classifiable','cdc1_feature_CLR',21,-.6363636363636362,.0019262960071318917),
      ('abundance_classifiable','effector_MHCII',21,-.4753246753246753,.02943042122641138),
      ('Mast_gt5_all','IFN_response',22,-.5013234712904683,.017459146157518948)]
    for cohort,target,n,rho,pval in benchmarks:
        row=cor.loc[(cohort,target)];ok=int(row.n)==n and np.isclose(row.rho,rho,atol=1e-10) and np.isclose(row.P_two_sided,pval,atol=1e-10)
        checks.append({'check':cohort+':'+target,'passed':bool(ok)})
        if not ok:raise AssertionError(target)
    (o/'legacy_agreement.json').write_text(json.dumps(checks,indent=2))
    print(f'PASS: {len(checks)} aggregate regression checks. No individual-level data exported.')
if __name__=='__main__':main()
