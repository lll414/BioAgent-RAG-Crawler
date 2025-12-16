# More information on results columns

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Information about which variables and tests were used can be found by calling the function *mcols* on the results object.

```
mcols(res)$description
```

```
## [1] "mean of normalized counts for all samples"             
## [2] "log2 fold change (MLE): condition treated vs untreated"
## [3] "standard error: condition treated vs untreated"        
## [4] "Wald statistic: condition treated vs untreated"        
## [5] "Wald test p-value: condition treated vs untreated"     
## [6] "BH adjusted p-values"
```

For a particular gene, a log2 fold change of -1 for `condition treated vs untreated` means that the treatment induces a multiplicative change in observed gene expression level of \(2^{-1} = 0.5\) compared to the untreated condition. If the variable of interest is continuous-valued, then the reported log2 fold change is per unit of change of that variable.

**Note on p-values set to NA**: some values in the results table can be set to `NA` for one of the following reasons:

* If within a row, all samples have zero counts, the `baseMean` column will be zero, and the log2 fold change estimates, *p* value and adjusted *p* value will all be set to `NA`.
* If a row contains a sample with an extreme count outlier then the *p* value and adjusted *p* value will be set to `NA`. These outlier counts are detected by Cookâs distance. Customization of this outlier filtering and description of functionality for replacement of outlier counts and refitting is described [below](#outlier)
* If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted *p* value will be set to `NA`. Description and customization of independent filtering is described [below](#indfilt)