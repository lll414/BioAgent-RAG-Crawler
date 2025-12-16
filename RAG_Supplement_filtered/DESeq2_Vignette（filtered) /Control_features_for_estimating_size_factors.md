# Control features for estimating size factors

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

In some experiments, it may not be appropriate to assume that a minority of features (genes) are affected greatly by the condition, such that the standard median-ratio method for estimating the size factors will not provide correct inference (the log fold changes for features that were truly un-changing will not centered on zero). This is a difficult inference problem for any method, but there is an important feature that can be used: the `controlGenes` argument of `estimateSizeFactors`. If there is any prior information about features (genes) that should not be changing with respect to the condition, providing this set of features to `controlGenes` will ensure that the log fold changes for these features will be centered around 0. The paradigm then becomes:

```
dds <- estimateSizeFactors(dds, controlGenes=ctrlGenes)
dds <- DESeq(dds)
```