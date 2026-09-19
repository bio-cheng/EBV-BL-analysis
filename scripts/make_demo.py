#!/usr/bin/env python3
"""Generate synthetic input tables for software testing."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def generate(dest):
    dest=Path(dest)
    if dest.exists() and any(dest.iterdir()):raise ValueError('Use a new empty demo directory')
    dest.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(104729)
    def wr(x,name,index=False):x.to_csv(dest/name,sep='\t',index=index)
    ids=[f'SYNTHETIC_{i:03d}' for i in range(32)]
    pfs=np.r_[rng.uniform(1,5.5,16),rng.uniform(7,25,16)]
    events=np.ones(32,dtype=int);events[[2,9,23,30]]=0
    c=pd.DataFrame({'sample_id':ids,'include_baseline':True,'treatment_scene':'Line_1',
      'batch':['Batch_A','Batch_B']*16,'response_group':['R','NR']*16,
      'cps':rng.choice([0,2,7,15,40,80],32),'PFS_time_months':pfs,'PFS_event':events,
      'OS_time_months':pfs+rng.uniform(3,16,32),'OS_event':rng.binomial(1,.7,32)})
    c.loc[[0,31],'treatment_scene']='Neoadjuvant';wr(c,'clinical.tsv')
    features=['Myeloid_Mast_TPSAB1','T_CD8_Tpex_CXCL13','T_CD8_TRM_ZNF683','Myeloid_DC_cDC1',
              'B_Plasma_IgA','T_Cycling_GINS2']+[f'Synthetic_feature_{i:02d}' for i in range(16)]
    major=['Mast','T cells','T cells','Myeloid','B cells','T cells']+['NK cells','Stromal','Epithelial','Neutrophils']*4
    ann=pd.DataFrame({'feature':features,'retained':True,'non_epithelial':[x!='Epithelial' for x in major],'major':major,'label':features})
    wr(ann,'annotation.tsv')
    n=pd.DataFrame(rng.integers(10,300,(32,len(features))),index=ids,columns=features);n.index.name='sample_id';wr(n,'cell_counts.tsv',True)
    eff=['CCL5','CCL4','PRF1','CTSW','GNLY','NKG7','HLA-DRA','HLA-DRB1']
    inf=['STAT1','IRF1','ISG15','IFIT1','IFIT3','MX1','OAS1','GBP1']
    defs=pd.read_csv(ROOT/'config/mast_gene_sets.tsv',sep='\t')
    genes=sorted(set(defs.gene)|set(eff+inf)|{f'SYNGENE{i:04d}' for i in range(500)})
    for name in ['mast','tpex']:
        mu=rng.uniform(10,150,(len(genes),1))*rng.lognormal(0,.2,(1,32))
        if name=='mast':
            for gene in ['IL5','IL13','CCL2','TGFB1']:mu[genes.index(gene),:16]*=2
        counts=pd.DataFrame(rng.negative_binomial(10,10/(10+mu)),index=genes,columns=ids);counts.index.name='gene'
        wr(counts,name+'_counts.tsv.gz',True)
        if name=='tpex':expr=np.log2(counts.div(counts.sum(axis=0),axis=1)*1e6+1)
    rows=[]
    for sig,gs in [('effector_MHCII',eff),('IFN_response',inf)]:
        for g in gs:rows.append({'signature':sig,'gene':g,'mean':expr.loc[g].mean(),'sd':expr.loc[g].std(ddof=1)})
    wr(pd.DataFrame(rows),'tpex_scaling.tsv')
    (dest/'SYNTHETIC_DATA_NOTICE.txt').write_text('Entirely simulated data. Not patient data. Not scientific evidence.\n')
    return dest
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out-dir',type=Path,required=True)
    generate(p.parse_args().out_dir)
