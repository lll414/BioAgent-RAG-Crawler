# Independent filtering of results

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

The *results* function of the DESeq2 package performs independent filtering by default using the mean of normalized counts as a filter statistic. A threshold on the filter statistic is found which optimizes the number of adjusted *p* values lower than a significance level `alpha` (we use the standard variable name for significance level, though it is unrelated to the dispersion parameter \(\alpha\)). The theory behind independent filtering is discussed in greater detail [below](#indfilttheory). The adjusted *p* values for the genes which do not pass the filter threshold are set to `NA`.

The default independent filtering is performed using the *filtered\_p* function of the [genefilter](http://bioconductor.org/packages/genefilter) package, and all of the arguments of *filtered\_p* can be passed to the *results* function. The filter threshold value and the number of rejections at each quantile of the filter statistic are available as metadata of the object returned by *results*.

For example, we can visualize the optimization by plotting the `filterNumRej` attribute of the results object. The *results* function maximizes the number of rejections (adjusted *p* value less than a significance level), over the quantiles of a filter statistic (the mean of normalized counts). The threshold chosen (vertical line) is the lowest quantile of the filter for which the number of rejections is within 1 residual standard deviation to the peak of a curve fit to the number of rejections over the filter quantiles:

```
metadata(res)$alpha
```

```
## [1] 0.1
```

```
metadata(res)$filterThreshold
```

```
##       0% 
## 5.109586
```

```
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
```


Independent filtering can be turned off by setting `independentFiltering` to `FALSE`.

```
resNoFilt <- results(dds, independentFiltering=FALSE)
addmargins(table(filtering=(res$padj < .1),
                 noFiltering=(resNoFilt$padj < .1)))
```

```
##          noFiltering
## filtering FALSE TRUE  Sum
##     FALSE  7079    0 7079
##     TRUE      0 1069 1069
##     Sum    7079 1069 8148
```