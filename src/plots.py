"""Render vector analysis panels from pipeline output tables."""
import os
os.environ.setdefault('MPLCONFIGDIR',str(__import__('tempfile').gettempdir()+'/ebv_mast_mpl'))
from pathlib import Path
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import statsmodels.api as sm
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Liberation Sans'],
 'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,
 'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':300,'legend.frameon':False})
BLUE='#2166AC';RED='#B2182B'
def read(o,n):return pd.read_csv(o/n,sep='\t')
def save(fig,f,name):
    fig.savefig(f/(name+'.pdf'),bbox_inches='tight');fig.savefig(f/(name+'.png'),bbox_inches='tight',dpi=300);plt.close(fig)
def forest(d,f,name,label='feature'):
    d=d.dropna(subset=['HR','CI_low','CI_high'])
    fig,ax=plt.subplots(figsize=(7,max(2,len(d)*.38)))
    for i,row in d.reset_index(drop=True).iterrows():
        color=RED if row.HR>1 else BLUE
        ax.plot([row.CI_low,row.CI_high],[i,i],c=color,lw=1)
        ax.plot(row.HR,i,'o',mfc=color if row.P<.05 else 'white',mec=color,ms=5)
        ax.text(1.02,i,f'{row.HR:.2f} [{row.CI_low:.2f}, {row.CI_high:.2f}]  P={row.P:.3g}',transform=ax.get_yaxis_transform(),va='center',fontsize=9)
    ax.set_yticks(range(len(d)),d[label]);ax.invert_yaxis();ax.set_xscale('log');ax.axvline(1,c='#999999',ls='--',lw=.7);ax.set_xlabel('Hazard ratio (95% CI)')
    save(fig,f,name)
def render(o):
    f=o/'figures';f.mkdir(exist_ok=True)
    if (o/'mast_DEG.tsv').exists():
        d=read(o,'mast_DEG.tsv');sig=d.PValue<.05
        fig,ax=plt.subplots(figsize=(4,4))
        ax.scatter(d.logFC,-np.log10(d.PValue.clip(lower=1e-300)),s=5,lw=0,c=np.where(sig,np.where(d.logFC>0,RED,BLUE),'#CCCCCC'),alpha=.7)
        ax.axhline(-np.log10(.05),ls='--',c='#888888',lw=.6)
        for side in [-1,1]:
            for i,(_,row) in enumerate(d[sig & (np.sign(d.logFC)==side)].nsmallest(5,'PValue').iterrows()):
                ax.annotate(row.gene,(row.logFC,-np.log10(row.PValue)),xytext=(side*5,5+i*2),textcoords='offset points',fontsize=7)
        ax.set(xlabel='log2FC (PFS <6 − PFS >6)',ylabel='−log10(P)',title=f'Mast: P<0.05 {sig.sum()}; FDR<0.05 {(d.FDR<.05).sum()}')
        save(fig,f,'Fig4a_volcano')
        g=read(o,'mast_GSEA.tsv').sort_values('NES',ascending=False)
        cfg=__import__('json').loads((o/'run_manifest.json').read_text())['config']
        shown=g[~g.pathway.isin(cfg['overview_omitted_sets'])].reset_index(drop=True)
        fig,ax=plt.subplots(figsize=(8,6.5))
        ax.barh(range(len(shown)),shown.NES,color=[RED if x>0 else BLUE for x in shown.NES],alpha=.6)
        ax.set_yticks(range(len(shown)),shown.display_label);ax.invert_yaxis();ax.axvline(0,c='#888888',lw=.6);ax.set_xlabel('NES (positive: short PFS)')
        for i,row in shown.iterrows():ax.text(1.02,i,f'P={row.pval:.3g}  FDR={row.padj:.3g}',transform=ax.get_yaxis_transform(),va='center',fontsize=9,fontweight='bold' if row.pval<.05 else 'normal')
        save(fig,f,'FigS4a_GSEA_overview')
        ranks=read(o,'mast_ranked_genes.tsv');defs=read(o,'gene_sets_used.tsv');v=ranks.ranking_metric.to_numpy();n=len(v)
        with PdfPages(f/'GSEA_all_curves.pdf') as pdf:
            for _,row in g.iterrows():
                genes=set(defs.loc[defs.module.eq(row.pathway),'gene']);hit=ranks.gene.isin(genes).to_numpy()
                es=np.r_[0,np.cumsum(np.where(hit,abs(v)/abs(v[hit]).sum(),-1/(n-hit.sum())))]
                if not np.isclose(es[np.argmax(abs(es))],row.ES,atol=1e-8):raise ValueError('Running ES mismatch')
                fig,axes=plt.subplots(3,1,figsize=(5,4.2),gridspec_kw={'height_ratios':[4,.4,1]},sharex=True)
                color=RED if row.NES>0 else BLUE
                axes[0].plot(np.arange(n+1),es,c=color);axes[0].axhline(0,c='#999999',ls='--',lw=.6)
                axes[0].set_title(row.display_label);axes[0].set_ylabel('Enrichment score')
                axes[0].text(.98,.97,f'NES={row.NES:.2f}\nP={row.pval:.3g}\nFDR={row.padj:.3g}',transform=axes[0].transAxes,ha='right',va='top')
                axes[1].vlines(np.flatnonzero(hit)+1,0,1,color=color,lw=.6);axes[1].axis('off')
                axes[2].fill_between(np.arange(n),v,0,color='#999999');axes[2].set(xlabel='Rank in ordered gene list',ylabel='Rank metric')
                fig.text(.12,-.05,'\n'.join(textwrap.wrap('Leading edge: '+str(row.leadingEdge).replace(';',', '),75)),fontsize=8)
                fig.tight_layout();pdf.savefig(fig,bbox_inches='tight')
                if row.pathway in ['Antigen_presentation_MHC','Th2_myeloid_immunoregulatory_output']:
                    save(fig,f,'Fig4b_'+row.pathway)
                else:plt.close(fig)
    c=read(o,'patients.tsv')
    if (o/'correlations.tsv').exists():
        rows=read(o,'correlations.tsv')
        with PdfPages(f/'Fig3_Fig4_correlations.pdf') as pdf:
            for _,row in rows.iterrows():
                flag='abundance_corr_eligible' if row.cohort=='abundance_classifiable' else 'mast_survival_eligible'
                d=c[c[flag].astype(str).str.lower().eq('true')]
                fig,ax=plt.subplots(figsize=(3.6,3.6))
                for group,color,marker in [('Long',BLUE,'o'),('Short',RED,'s'),('Early','#AAAAAA','^'),('Boundary','#777777','D'),('Missing','#777777','x')]:
                    if not d.pfs6.eq(group).any():continue
                    q=d[d.pfs6.eq(group)];ax.scatter(q[row.x],q[row.y],s=25,c=color,marker=marker,label=group,edgecolor='white',lw=.4)
                fit=sm.OLS(d[row.y],sm.add_constant(d[row.x])).fit();grid=np.linspace(d[row.x].min(),d[row.x].max(),100)
                pred=fit.get_prediction(sm.add_constant(grid)).summary_frame()
                ax.plot(grid,pred['mean'],c='#777777');ax.fill_between(grid,pred.mean_ci_lower,pred.mean_ci_upper,color='#EEEEEE',zorder=0)
                ax.set(xlabel=row.x,ylabel=row.y,title=f'n={row.n}; ρ={row.rho:.2f}; P={row.P_two_sided:.3g}')
                ax.legend(loc='upper center',bbox_to_anchor=(.5,-.25),ncol=3,fontsize=8)
                pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
    boxdata=c[c.mast_de_eligible.astype(str).str.lower().eq('true')].copy()
    if (o/'leading_edge_patient_scores.tsv').exists():boxdata=boxdata.merge(read(o,'leading_edge_patient_scores.tsv')[['sample_id','leading_edge_score']],on='sample_id')
    if (o/'Mast_gene_boxplot_statistics.tsv').exists():
        stats=read(o,'Mast_gene_boxplot_statistics.tsv').set_index('feature')
        if (o/'leading_edge_statistics.tsv').exists():stats=pd.concat([stats,read(o,'leading_edge_statistics.tsv').assign(feature='leading_edge_score').set_index('feature')])
        with PdfPages(f/'Fig4c_FigS4b_boxplots.pdf') as pdf:
            for gene,row in stats.iterrows():
                arrays=[boxdata.loc[boxdata.pfs6.eq(g),gene].to_numpy() for g in ['Long','Short']]
                fig,ax=plt.subplots(figsize=(3,3.5));bp=ax.boxplot(arrays,positions=[0,1],widths=.5,patch_artist=True,showfliers=False)
                rng=np.random.default_rng(20260919)
                for i,(values,color) in enumerate(zip(arrays,[BLUE,RED])):
                    bp['boxes'][i].set(facecolor=color,alpha=.2);ax.scatter(rng.normal(i,.055,len(values)),values,s=18,c=color,edgecolor='white',lw=.4)
                ax.set_xticks([0,1],[f'>6 mo\nn={len(arrays[0])}',f'<6 mo\nn={len(arrays[1])}'])
                ax.set(xlabel='PFS',ylabel=gene,title=f'One-sided P={row.P_one_sided:.3g}')
                pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
    if (o/'Fig3e_top12_PFS.tsv').exists():forest(read(o,'Fig3e_top12_PFS.tsv'),f,'Fig3e_top12_PFS')
    if (o/'clinical_Cox.tsv').exists():
        d=read(o,'clinical_Cox.tsv');q=d[d.model.isin(['CPS_only','CPS_plus_Mast'])].copy();q['label']=q.model+': '+q.term
        forest(q,f,'Fig3g_CPS_Mast',label='label')
    if (o/'CD8_program_survival.tsv').exists():
        for ep,q in read(o,'CD8_program_survival.tsv').groupby('endpoint'):forest(q,f,'FigS3h_CD8_'+ep,label='program')
    if (o/'KM_coordinates.tsv').exists():
        coords=read(o,'KM_coordinates.tsv');stats=read(o,'KM_statistics.tsv');risks=read(o,'KM_risk_tables.tsv')
        with PdfPages(f/'PFS_OS_KM_panels.pdf') as pdf:
            for (name,ep),d in coords.groupby(['model','endpoint']):
                if name.startswith('Mast_quantile_'):continue
                fig,(ax,ra)=plt.subplots(2,1,figsize=(4,4),gridspec_kw={'height_ratios':[4,1]})
                for i,(g,color) in enumerate([('Low',BLUE),('High',RED)]):
                    q=d[d.group.eq(g)];ax.step(q.time,q.survival,where='post',c=color,label=g)
                    cens=q[q.n_censor>0];ax.plot(cens.time,cens.survival,'+',c=color,ms=5)
                    risk=risks[(risks.model==name)&(risks.endpoint==ep)&(risks.group==g)]
                    for z in risk.itertuples():ra.text(z.time,1-i,str(z.n_risk),ha='center',va='center',fontsize=8)
                st=stats[(stats.model==name)&(stats.endpoint==ep)].iloc[0]
                ax.set(title=f'{name}: {ep}\nLog-rank P={st.logrank_P:.3g}',ylabel='Survival probability',ylim=(-.03,1.05));ax.legend(loc='upper left',bbox_to_anchor=(1,1))
                ra.set(xlim=ax.get_xlim(),ylim=(-.5,1.5),yticks=[0,1],yticklabels=['High','Low'],xlabel=ep+' (months)');ra.set_title('No. at risk',fontsize=9,loc='left')
                fig.tight_layout();pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
    m=read(o,'major_group_means.tsv').set_index('treatment_scene')
    palette=['#F08A4B','#8064A2','#B783C5','#87B7D0','#D95F5F','#A85D32','#6F9DA4','#76A55A']
    fig,ax=plt.subplots(figsize=(4.5,4));bottom=np.zeros(len(m))
    for color,col in zip(palette,m.columns):ax.bar(m.index,m[col],bottom=bottom,label=col,color=color);bottom+=m[col]
    ax.set(ylabel='Patient-equal mean proportion',ylim=(0,1));ax.tick_params(axis='x',rotation=25);ax.legend(loc='upper left',bbox_to_anchor=(1,1));save(fig,f,'FigS3d_major_composition')
