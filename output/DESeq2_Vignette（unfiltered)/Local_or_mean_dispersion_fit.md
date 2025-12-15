# Local or mean dispersion fit

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

A local smoothed dispersion fit is automatically substitited in the case that the parametric curve doesnât fit the observed dispersion mean relationship. This can be prespecified by providing the argument `fitType="local"` to either *DESeq* or *estimateDispersions*. Additionally, using the mean of gene-wise disperion estimates as the fitted value can be specified by providing the argument `fitType="mean"`.