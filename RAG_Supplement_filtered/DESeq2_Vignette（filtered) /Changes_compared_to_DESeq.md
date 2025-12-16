# Changes compared to DESeq

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

The main changes in the package *DESeq2*, compared to the (older) version *DESeq*, are as follows:

* *RangedSummarizedExperiment* is used as the superclass for storage of input data, intermediate calculations and results.
* Optional, maximum *a posteriori* estimation of GLM coefficients incorporating a zero-centered Normal prior with variance estimated from data (equivalent to Tikhonov/ridge regularization). This adjustment has little effect on genes with high counts, yet it helps to moderate the otherwise large variance in log2 fold change estimates for genes with low counts or highly variable counts. These estimates are now provided by the *lfcShrink* function.
* Maximum *a posteriori* estimation of dispersion replaces the `sharingMode` options `fit-only` or `maximum` of the previous version of the package. This is similar to the dispersion estimation methods of DSS (Wu, Wang, and Wu 2012).
* All estimation and inference is based on the generalized linear model, which includes the two condition case (previously the *exact test* was used).
* The Wald test for significance of GLM coefficients is provided as the default inference method, with the likelihood ratio test of the previous version still available.
* It is possible to provide a matrix of sample-/gene-dependent normalization factors.
* Automatic independent filtering on the mean of normalized counts.
* Automatic outlier detection and handling.