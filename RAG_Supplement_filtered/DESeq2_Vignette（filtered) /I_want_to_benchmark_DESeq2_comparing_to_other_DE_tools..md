# I want to benchmark DESeq2 comparing to other DE tools.

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

One aspect which can cause problems for comparison is that, by default, DESeq2 outputs `NA` values for adjusted *p* values based on independent filtering of genes which have low counts. This is a way for the DESeq2 to give extra information on why the adjusted *p* value for this gene is not small. Additionally, *p* values can be set to `NA` based on extreme count outlier detection. These `NA` values should be considered *negatives* for purposes of estimating sensitivity and specificity. The easiest way to work with the adjusted *p* values in a benchmarking context is probably to convert these `NA` values to 1:

```
res$padj <- ifelse(is.na(res$padj), 1, res$padj)
```