# Expanded model matrices

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

For the specific combination of `lfcShrink` with the type `normal` and using `contrast`, DESeq2 uses *expanded model matrices* to produce shrunken log2 fold change estimates where the shrinkage is independent of the choice of reference level. In all other cases, DESeq2 uses standard model matrices, as produced by `model.matrix`. The expanded model matrices differ from the standard model matrices, in that they have an indicator column (and therefore a coefficient) for each level of factors in the design formula in addition to an intercept. This is described in the DESeq2 paper. Using type `normal` with `coef` uses standard model matrices, as does the `apeglm` shrinkage estimator.