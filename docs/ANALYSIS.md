# Analysis decisions

## Denominators and patient sets

Raw subtype abundance is subtype count divided by the total of all annotated
cell counts (67 features in the archived input). Subtype Cox screening uses this
raw proportion. Retained composition uses 60 features: CLR is
`ln(count + 0.5) - mean_j ln(count_j + 0.5)` within each patient. Major-lineage
stacked bars renormalize retained counts to one and average patient proportions
within treatment scene, giving patients equal weight. They are not pooled cell
counts and do not have the same denominator as raw-proportion Cox models.

All inferential tests operate at patient level. Survival analyses restrict to
Line_1/Line_2plus; the clinical overview can include neoadjuvant samples. Mast
expression requires >5 Mast cells. Abundance correlations in the classifiable
PFS cohort have no Mast expression-count gate. The eligible cohort and exclusions
are exported for each run. Tpex expression requires available target libraries;
absence is never treated as zero expression.

## Survival

Continuous effects are HR per interquartile-range increase. Subtype screening
uses Breslow ties, n>=8, >=3 events, >=3 distinct values and nonzero IQR.
BH correction is across evaluable retained subtypes within each endpoint.
The displayed top 12 are selected from non-epithelial subtypes by PFS P value;
the complete table is also retained. Selection is exploratory.

CPS-only and CPS-plus-Mast models use identical complete cases. Nested likelihood
ratio, delta AIC and proportional-hazards diagnostics are saved. Adjustment for
CPS alone is not evidence of independence from all clinical confounders.

Most Cox models use Breslow ties. Mast IL5 and overall-CD8 program models retain
the legacy Efron convention. KM uses median abundance (high >= median, with
high > median fallback for all-tied grouping); IL5 uses high > median, matching
its zero-heavy distribution. Cox and log-rank P values must not be interchanged.
Cutpoint sensitivity is exploratory, not a search for a new optimal threshold.

## Expression, GSEA and association tests

Single-gene expression and Tpex scoring use log2(CPM+1) with the full raw library
sum as denominator. Tpex mean gene-z scores use frozen gene-wise reference means
and sample SDs. Overall-CD8 optional input scores use their original reference.

Mast differential expression uses edgeR filterByExpr on the current design, TMM,
robust dispersion estimation and robust quasi-likelihood fitting. Design is
`~ batch + pfs6`; batch values are preserved from metadata, not inferred from
R1 length. GSEA uses signed sqrt(QL F), all tested genes, the versioned 16-set
family and fgseaMultilevel. BH is computed before hiding overview rows. The
two hidden overview sets remain in the tested family. Duplicate rankings can
make Monte Carlo enrichment P values platform/version sensitive.

The Th2 leading-edge score is recomputed from the current GSEA leading edge,
using TMM logCPM (prior.count=2), gene z scores across eligible short/long-PFS
patients and their mean. It is distinct from a complete gene-set score.

Spearman correlations are two-sided. Lines are descriptive ordinary least
squares fits with 95% mean-response confidence bands; their slopes/P values are
not the reported Spearman statistic. Group colors do not define correlation
subgroups. Boxplot rank-sum tests use asymptotic, continuity-corrected inference;
both one-sided (short > long) and two-sided values are exported.

No causal effect, independent validation or treatment-specific predictive
performance is established by these exploratory associations.
