# Private input contract

Keep actual inputs outside this repository. `sample_id` must identify one baseline
specimen per patient, consistently across tables. Empty numeric fields represent
missing values; do not encode missing values as zero. Gene identifiers are unique
human gene symbols. The loader does not perform gene alias conversion.

| File | Columns / orientation |
|---|---|
| `clinical.tsv` | sample_id, include_baseline (true/false), treatment_scene, batch, response_group, cps, PFS_time_months, PFS_event, OS_time_months, OS_event |
| `cell_counts.tsv` | sample_id followed by integer counts for every annotated subtype |
| `annotation.tsv` | feature, retained, non_epithelial, major, label; booleans true/false |
| `mast_counts.tsv.gz` | gene followed by integer raw pseudobulk counts per patient; **all genes**, not just signature genes |
| `tpex_counts.tsv.gz` | same format for CD8 Tpex CXCL13 |
| `tpex_scaling.tsv` | signature, gene, mean, sd; fixed reference mean and sample SD on log2(CPM+1) |
| `cd8_program_scores.tsv` (optional) | sample_id, program, score, n_cells; precomputed overall-CD8 patient scores |
| `external_mast_expression.tsv` (optional) | sample_id, response_group (R/NR), MS4A2, KIT, HDC, LTC4S, IL1RL1, CMA1; **raw TPM**, not log values |

Treatment scenes: `Neoadjuvant`, `Line_1`, `Line_2plus`, `Unknown`.
Events: 1=event, 0=right censoring. Times are months from the endpoint's stated
origin. `include_baseline=false` is the explicit exclusion mechanism; exclusion
reasons belong in the private clinical audit, not public code.

`annotation.tsv` must list all count columns exactly once, including excluded
features. `retained` defines the CLR and stacked-composition feature universe;
`non_epithelial` defines the subset eligible for the displayed subtype ranking.
Major labels: T cells, NK cells, B cells, Myeloid, Neutrophils, Mast, Stromal,
Epithelial. Mast must not also be counted under Myeloid in the major display.

Tpex signatures must include `effector_MHCII` and `IFN_response`. The former uses
CCL5, CCL4, PRF1, CTSW, GNLY, NKG7, HLA-DRA, HLA-DRB1; the latter STAT1, IRF1,
ISG15, IFIT1, IFIT3, MX1, OAS1, GBP1. The original reference cohort had 20
evaluable patients. Re-estimating means/SDs on a different cohort changes the
composite and does not reproduce the frozen analysis. These derived private
parameters are deliberately not shipped.

Optional CD8 scores preserve the original five-program scaling; the public
pipeline recomputes their Cox models but does not reconstruct those input scores
from single cells. Missing optional files are explicitly reported as skipped.

The legacy importer converts archived source tables to this contract, including
inverting log2(TPM+1) in the external frozen table. Do not apply that inversion
to a raw TPM table. Clinical sample exclusion is an explicit importer argument;
the runner never removes named patients based on a hard-coded list.
