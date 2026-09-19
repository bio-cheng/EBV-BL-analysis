#!/usr/bin/env python3
"""Convert archived analysis tables to the pipeline input schema.

Patient tables are written outside this repository. The adapter uses the archived
project layout and performs no network operations.
"""
import argparse
import ast
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def read(p):return pd.read_csv(p,sep='\t')
def assignment(path,key):
    tree=ast.parse(path.read_text())
    return next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id==key for t in n.targets))
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--legacy-root',type=Path,required=True,help='Legacy tme_analysis directory')
    p.add_argument('--dest',type=Path,required=True,help='Private destination OUTSIDE the release repository')
    p.add_argument('--exclude-clinical-sample',action='append',default=[],help='Exact clinical specimen label; may repeat')
    p.add_argument('--build-gene-definitions',action='store_true',help='Write only non-patient gene metadata into config/')
    a=p.parse_args();b=a.legacy_root.resolve();dest=a.dest.resolve()
    if dest==ROOT or ROOT in dest.parents:raise ValueError('Private import destination must be outside the public repository')
    if dest.exists() and any(dest.iterdir()):raise ValueError('Use an empty private destination')
    dest.mkdir(parents=True,exist_ok=True)
    c=read(b/'tables/clinical_manifest_with_analysis_fields.tsv')
    excluded=c.clinical_raw_sample_id.isin(a.exclude_clinical_sample)
    for name in a.exclude_clinical_sample:
        if not c.clinical_raw_sample_id.eq(name).any():raise ValueError('Requested exclusion not found')
    c['include_baseline']=c.composition_match.astype(str).str.lower().eq('true') & ~excluded
    cols=['old_id_sample_id','include_baseline','treatment_scene','batch','response_group','cps','PFS_time_months','PFS_event','OS_time_months','OS_event']
    c[cols].rename(columns={'old_id_sample_id':'sample_id'}).to_csv(dest/'clinical.tsv',sep='\t',index=False)
    counts=read(b/'tables/source_data_all67_cell_counts.tsv').rename(columns={'old_id_sample_id':'sample_id'})
    counts.to_csv(dest/'cell_counts.tsv',sep='\t',index=False)
    labels=read(b/'tables/plot_display_label_dictionary.tsv').set_index('canonical_feature').plot_label
    olddict=read(b.parents[1]/'treatment_line_cell_composition/reannotated_BMS67_annotation_dictionary.tsv')
    rename=assignment(b/'run_tme_clinical_analysis.py','CANONICAL_FEATURE_RENAMES')
    olddict['feature']=olddict.feature.replace(rename)
    kept=set(read(b/'tables/source_data_inferred60_closed_proportions.tsv').columns)-{'old_id_sample_id'}
    mapping={'T_cell':'T cells','B_cell':'B cells','Epi_cell':'Epithelial','Stromal':'Stromal','Neutrophils':'Neutrophils','Mast':'Mast','Myeloid':'Myeloid'}
    olddict['major']=olddict.cell_type_major.map(mapping)
    olddict.loc[olddict.feature.isin(['NK_Cyto_FCGR3A','NK_Cyto_FGFBP2','gdT_NKlike_KLRC2']),'major']='NK cells'
    olddict.loc[olddict.feature.eq('Myeloid_Mast_TPSAB1'),'major']='Mast'
    olddict.loc[olddict.feature.str.startswith('Myeloid_Neutro'),'major']='Neutrophils'
    olddict['retained']=olddict.feature.isin(kept);olddict['non_epithelial']=olddict.cell_type_major.ne('Epi_cell')
    olddict['label']=olddict.feature.map(labels).fillna(olddict.feature)
    assert olddict.feature.is_unique and olddict.major.notna().all()
    olddict[['feature','retained','non_epithelial','major','label']].to_csv(dest/'annotation.tsv',sep='\t',index=False)
    for subtype,name in [('Myeloid_Mast_TPSAB1','mast'),('T_CD8_Tpex_CXCL13','tpex')]:
        shutil.copyfile(b/'pseudobulk'/subtype/'raw_counts.tsv.gz',dest/(name+'_counts.tsv.gz'))
    ef=read(b/'IL5_Tpex_cytotoxicity_20260917/custom8_gene_patient_expression.tsv').pivot(index='sample_id',columns='gene',values='Tpex_log2CPM_plus1')
    genes=['STAT1','IRF1','ISG15','IFIT1','IFIT3','MX1','OAS1','GBP1']
    sub=b/'mast_Tpex_function_support_20260830/tables'
    ref_ids=read(sub/'05_primary_Tpex_patient_function_scores.tsv').sample_id
    interferon=read(sub/'02_Tpex_signature_gene_patient_expression.tsv').pivot(index='sample_id',columns='gene',values='Tpex_log2CPM_plus1').loc[ref_ids,genes]
    ref=[]
    for name,frame in [('effector_MHCII',ef),('IFN_response',interferon)]:
        assert len(frame)==20
        for gene in frame:ref.append({'signature':name,'gene':gene,'mean':frame[gene].mean(),'sd':frame[gene].std(ddof=1)})
    pd.DataFrame(ref).to_csv(dest/'tpex_scaling.tsv',sep='\t',index=False)
    scores=read(b/'immune_program_response_survival_20260904/tables/04_CD8_Mast_immune_program_patient_scores.tsv')
    scores[scores.lineage.eq('CD8 T')][['old_id_sample_id','program_code','score','n_cells']].rename(columns={'old_id_sample_id':'sample_id','program_code':'program'}).to_csv(dest/'cd8_program_scores.tsv',sep='\t',index=False)
    ex=read(b/'figure_replacement_panels_20260918/tables/Fig3i_external_patient_source.tsv')
    # The frozen source contains log2(TPM+1), not raw TPM: explicitly invert it.
    genes=['MS4A2','KIT','HDC','LTC4S','IL1RL1','CMA1']
    ex[genes]=np.exp2(ex[genes])-1
    ex['response_group']=ex.response_group.map({'Responder':'R','Nonresponder':'NR','R':'R','NR':'NR'})
    assert ex.response_group.notna().all()
    ex[['sample_id','response_group']+genes].to_csv(dest/'external_mast_expression.tsv',sep='\t',index=False)
    if a.build_gene_definitions:
        src=b/'fig4_PFS6_replacement_20260907/tables/07_original_module_gene_definitions.tsv'
        defs=read(src)
        labels=assignment(b/'fig4_PFS6_replacement_20260907/scripts/03_enlarge_publication_fonts.py','short')
        defs['display_label']=defs.module.map(labels)
        defs[['module','gene','display_label','evidence_basis','public_source_id','public_source_url']].to_csv(ROOT/'config/mast_gene_sets.tsv',sep='\t',index=False)
        cd8=read(b/'Tcell_function_NR_vs_R_20260903/tables/09_CD8_five_program_definitions.tsv')
        cd8[['module','program_label','genes','reference','reference_url']].to_csv(ROOT/'config/cd8_gene_sets.tsv',sep='\t',index=False)
        (ROOT/'docs/gene_definition_provenance.json').write_text(json.dumps({'mast_source_basename':src.name,'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'note':'Only curated gene definitions and source metadata included; patient-derived scores and gene SDs removed.'},indent=2))
    print('Private inputs written outside repository. Do not commit or publish them.')
if __name__=='__main__':main()
