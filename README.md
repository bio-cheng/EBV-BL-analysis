# EBV-BL-analysis

Patient-level analyses for the single-cell translational component of the EBVaGC
study (Figures 3 and 4 and related supplementary analyses).

This release candidate contains code and documentation. Patient data, manuscript
files, credentials and workstation paths are excluded. The workflow starts from
annotated cell counts and patient-level pseudobulk counts. Panel numbering is
provisional; see [panel coverage](docs/PANEL_MAP.md).

## 中文说明

本仓库包含 Fig3/4 及相关补充分析的统计与绘图代码。默认 Mast 纳入阈值
为每位患者 >5 cells。临床表、单细胞表达矩阵、固定打分参数、患者级结果均
需要保存在仓库之外；`.gitignore` 不能替代上传前的人工隐私检查。

分析起点为已注释的细胞计数和患者级表达矩阵。FASTQ处理、细胞注释、最终排版
及部分补充图不在当前范围内。原始研究脚本及结果另行保留。

## Quick start

Use Python 3.12 and R 4.3.3. Recorded Python versions are in `requirements.txt`;
R package versions are in `config/R-packages.tsv`. These describe the tested
environment; a cross-platform lockfile is not provided. See [environment notes](docs/ENVIRONMENT.md).

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v

# Private, authorized inputs outside this repository:
python run.py --data-dir ../private-inputs --out-dir ../private-results --rscript Rscript
```

Input and output directories must be disjoint. Existing output is refused unless
`--allow-existing-output` is explicitly supplied. No download or upload occurs.
The runner saves input/code SHA256 hashes, effective parameters, runtime versions,
R session information, statistical tables and vector PDFs. Outputs contain
sensitive patient-level information and require review before publication.

To check the interface without real data:

```bash
python scripts/make_demo.py --out-dir ../synthetic-inputs
python run.py --data-dir ../synthetic-inputs --out-dir ../synthetic-results --stage python
```

Synthetic data are for software testing only; no biological inference is valid.
The demo supports the full runner when R dependencies are available.

## Inputs and stages

See [input schema](data/README.md). All tables are tab-delimited; pseudobulk counts
are gene-by-patient gzip TSVs. Cell counts and raw RNA counts must be integers.

| Stage | Function |
|---|---|
| `validate` | Input checks, cohort eligibility, proportions, CLR, frozen-reference Tpex scores |
| `python` | Above plus subtype Cox screening, correlations, rank-sum tests, optional CD8/external analyses |
| `mast` | Above plus sample-level edgeR QL differential expression and focused GSEA |
| `survival` | Above plus CPS/Mast Cox models and KM coordinates/risk tables |
| `plot` | Replot existing complete analysis outputs; requires `--allow-existing-output` |
| `all` | All statistical stages, followed by plots |

`--config` accepts an edited copy of `config/analysis.json`. Gene definitions are
the versioned tables in `config/`. The optional legacy importer supports the
archived project directory structure. Other inputs must follow the documented
schema. Imported data must be stored outside this repository.

## Statistical scope and caveats

- The patient, not the cell, is the inferential unit.
- PFS <6 months requires an observed progression/death event before six months;
  PFS >6 includes later-censored patients. Early censoring and exactly six months
  are excluded from this strict binary comparison, but not automatically from Cox.
- Mast differential expression: raw pseudobulk counts, edgeR TMM and robust QL,
  `~ batch + pfs6`; positive logFC means short PFS. GSEA uses the entire filtered
  signed `sqrt(F)` ranking, not only significant genes.
- Curated Mast programs are literature-informed gene sets, not all exact public
  database pathways. Their definitions and inherited evidence metadata are supplied.
- Directional boxplots use the specified one-sided short-PFS > long-PFS test;
  two-sided results are also saved. Test directions were selected during
  exploratory analysis rather than prespecified prospectively.
- Leading-edge genes selected from this cohort are not an independent validation.
- Continuous Cox effects, dichotomized KM results and compositional correlations
  are different estimands. Sample sizes and tie handling are documented in
  [analysis decisions](docs/ANALYSIS.md).

## Before publication

1. Confirm manuscript version and replace old >10-cell figures where appropriate.
2. Verify batch definitions, data-sharing approval and controlled-access metadata.
3. Complete upstream preprocessing/QC documentation. The specified Harmony factor
   is sample; implementation awaits verification of the remaining parameters.
4. Select a license with the code owners and add authors/manuscript citation.
5. Run `python scripts/package_release.py --output ../ebv-mast-code.zip` and inspect
   the archive. This command packages a strict code/documentation allowlist.
6. Upload the reviewed code package to the project GitHub repository.

The reported data accession is HRA020437. Access conditions require confirmation.
The earlier study (DOI `10.1016/j.scib.2026.07.049`) includes overlapping patients;
the overlap count requires confirmation. This comparison is not independent
cohort validation.

Repository: [EBV-BL-analysis](https://github.com/bio-cheng/EBV-BL-analysis).
