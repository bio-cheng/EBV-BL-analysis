# Environment

Recorded host: Linux x86_64, Python 3.12.2, R 4.3.3. See requirements.txt and
config/R-packages.tsv for versions actually used. R also requires the transitive
dependencies of edgeR, fgsea and survival, including statmod for robust edgeR.

R installation commands (run separately; network access required):

```r
install.packages(c('BiocManager', 'data.table', 'survival', 'statmod'))
BiocManager::install(version = '3.18') # R 4.3 release family
BiocManager::install(c('edgeR', 'fgsea', 'BiocParallel'))
```

This installs repository-available releases and may not resolve the exact
recorded patch versions. For exact numeric replication, provision the versions
listed in config/R-packages.tsv and inspect the saved sessionInfo. A fresh
installation/container build has not been validated, and a fully resolved
lockfile is not provided. Binary libraries must match the R version.

The local Conda `bin/Rscript` launcher returned 255 even for a successful trivial
program. Its underlying `lib/R/bin/Rscript` returned 0; validation uses that
working executable. The pipeline stops on nonzero exit codes. Use `--rscript` to point
to a tested Rscript executable on your machine.

Arial is preferred; Liberation Sans is the font fallback. Different installed
fonts can change label geometry without changing statistics.
