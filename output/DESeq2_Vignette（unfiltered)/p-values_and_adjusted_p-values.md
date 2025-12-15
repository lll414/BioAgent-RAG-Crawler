# p-values and adjusted p-values

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

We can order our results table by the smallest *p* value:

```
resOrdered <- res[order(res$pvalue),]
```

We can summarize some basic tallies using the *summary* function.

```
summary(res)
```

```
## 
## out of 8148 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)       : 533, 6.5%
## LFC < 0 (down)     : 536, 6.6%
## outliers [1]       : 0, 0%
## low counts [2]     : 0, 0%
## (mean count < 5)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

How many adjusted p-values were less than 0.1?

```
sum(res$padj < 0.1, na.rm=TRUE)
```

```
## [1] 1069
```

The *results* function contains a number of arguments to customize the results table which is generated. You can read about these arguments by looking up `?results`. Note that the *results* function automatically performs independent filtering based on the mean of normalized counts for each gene, optimizing the number of genes which will have an adjusted *p* value below a given FDR cutoff, `alpha`. Independent filtering is further discussed [below](#indfilt). By default the argument `alpha` is set to \(0.1\). If the adjusted *p* value cutoff will be a value other than \(0.1\), `alpha` should be set to that value:

```
res05 <- results(dds, alpha=0.05)
summary(res05)
```

```
## 
## out of 8148 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 416, 5.1%
## LFC < 0 (down)     : 437, 5.4%
## outliers [1]       : 0, 0%
## low counts [2]     : 0, 0%
## (mean count < 5)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```
sum(res05$padj < 0.05, na.rm=TRUE)
```

```
## [1] 853
```