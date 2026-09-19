# Coverage and provenance

Panel numbering is provisional. The table maps each analysis to its implementation
and output. Scripts render individual vector panels; final manuscript page layout
and withdrawn exploratory panels are outside the current scope.

| Content | Implementation / output | Coverage |
|---|---|---|
| Subtype PFS/OS screen; Fig3e top12 | src/analysis.py; subtype_survival_all.tsv, Fig3e_top12_PFS | Recomputed from cell counts |
| Mast/CPS joint Cox; Fig3g | scripts/survival.R; clinical_Cox.tsv | Recomputed |
| Mast / selected subtype / CPS KM | scripts/survival.R; KM tables, PFS_OS_KM_panels | Recomputed; includes sensitivity curves in tables |
| Fig3h Mast correlations with TRM/Tpex/cDC1 | src/analysis.py; correlations.tsv | Recomputed CLR and signatures using frozen reference scaling |
| Major composition by treatment scene | src/analysis.py; major_group_means.tsv | Recomputed patient-equal bars |
| CPS-associated composition | src/analysis.py; CPS_composition.tsv | Statistics only; legacy volcano layout not ported |
| CD8 program PFS/OS | src/analysis.py; CD8_program_survival.tsv | Cox/forest recomputed from optional frozen scores |
| External Mast marker signature | src/analysis.py; external_Mast_statistics.tsv | Score/test recomputed from optional TPM; no external-data download |
| Fig4a Mast PFS6 volcano | scripts/mast.R; mast_DEG.tsv | Recomputed, >5 cells |
| Fig4b / S4a focused GSEA | scripts/mast.R; mast_GSEA.tsv; src/plots.py | Recomputed 16-set tests, 14-set overview |
| Fig4c / S4b IL5/CCL2 boxplots | src/analysis.py, src/plots.py | Recomputed short/long comparison, excludes early censoring |
| Fig4 IL5 PFS KM | scripts/survival.R | Recomputed, includes right censoring |
| Fig4f/g IL5 correlations with Tpex abundance/function | src/analysis.py | Recomputed, >5-cell cohort |
| Th2 leading-edge comparison | scripts/mast.R | Recomputed from current leading edge |
| Atlas UMAP, subtype UMAP, marker/receptor dotplots | Not ported | Requires annotated single-cell objects and verified upstream settings |
| Sample/Harmony preprocessing; tissue/QC | Not ported | Sample factor specified; remaining upstream parameters require verification |
| Final multi-panel page montage, clinical heatmap | Not ported | Layout/source preparation remains outside this package |

Principal archived analysis directories (relative to the original analysis root):

- `Mast_gt5_PFS6_DEG_GSEA_20260919`: current sample-level DEG/GSEA reference.
- `CPS_Mast_joint_PFS_20260910`: CPS joint models.
- `IL5_Tpex_cytotoxicity_20260917`: frozen Tpex composites and correlation cohorts.
- `immune_program_response_survival_20260904`: overall-CD8 program scores.
- `figure_replacement_panels_20260918`: external signature source.
- `fig4_PFS6_replacement_20260907`: original module definitions and labels.

The importer reads these sources without modifying them. Earlier figure versions
use >10 cells, whereas this release defaults to >5 cells. The two thresholds
produce different cohorts and results.
