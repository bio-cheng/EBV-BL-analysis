# Validation record — 2026-09-19

## Executed locally

- Complete `--stage all` run on the authorized private inputs: passed, including
  R differential expression/GSEA, survival analysis and PDF rendering.
- Complete `--stage all` run on newly generated synthetic inputs: passed.
- Seven unit tests: passed (censoring/boundary rules, Boolean validation, BH with
  missing values, full-library CPM, integer-count rejection, synthetic preparation
  and strict >5-cell gating).
- Python compilation: passed.
- Twenty-six aggregate regression checks against frozen analysis results: passed
  at rtol=1e-7, atol=1e-10 (exact cohorts; numerical tolerance for statistics).
- Public-file allowlist and private-content-pattern scan: passed. This is a
  limited automated scan, not a guarantee replacing owner review.

## Regression coverage

The 26 checks comprise five DEG columns across 4,398 genes, five GSEA columns
across 16 programs, six joint-Cox columns across three coefficients, five cohort
counts and five correlation benchmarks. Checks can be repeated using
`scripts/check_legacy_agreement.py` with authorized archived source tables.

| Quantity | Reproduced value |
|---|---|
| Baseline analysis patients | 29 |
| Non-neoadjuvant patients | 25 |
| Classifiable PFS abundance-correlation cohort | 21 |
| Non-neoadjuvant Mast >5-cell expression cohort | 22 |
| Mast >5-cell binary PFS DEG cohort | 18 (9 short / 9 long) |
| Tested Mast genes | 4,398 |
| MHC GSEA NES / P / BH q | −1.965420 / 0.00129549 / 0.02072783 |
| CPS-adjusted Mast PFS HR / P | 1.218152 / 0.03274783 |
| Mast CLR–TRM CLR Spearman rho / P | −0.623377 / 0.00253469 |
| Mast CLR–Tpex CLR Spearman rho / P | −0.561039 / 0.00814482 |

The original aggregate reference and all patient-level validation outputs are
stored outside this code package. No individual identifiers are needed in this
record. These results are exploratory and refer to the current analysis version,
not an independent validation cohort.

## Not established by this validation

Fresh installation on another machine, exact final page layout, FASTQ processing,
sample-level Harmony preprocessing, omitted marker/UMAP panels and every legacy
exploratory comparison have not been validated by this package. Optional CD8 and
external-signature analyses ran locally, but are not part of the 26 frozen
regression checks. Upstream optional scores remain processed inputs.
