# How can I include a continuous covariate in the design formula?

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Continuous covariates can be included in the design formula in exactly the same manner as factorial covariates, and then *results* for the continuous covariate can be extracted by specifying `name`. Continuous covariates might make sense in certain experiments, where a constant fold change might be expected for each unit of the covariate. However, in some cases, more meaningful results may be obtained by cutting continuous covariates into a factor defined over a small number of bins (e.g.Â 3-5). In this way, the average effect of each group is controlled for, regardless of the trend over the continuous covariates. In R, *numeric* vectors can be converted into *factors* using the function *cut*.