# Exporting results to CSV files

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

A plain-text file of the results can be exported using the base R functions *write.csv* or *write.delim*. We suggest using a descriptive file name indicating the variable and levels which were tested.

```
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_results.csv")
```

Exporting only the results which pass an adjusted *p* value threshold can be accomplished with the *subset* function, followed by the *write.csv* function.

```
resSig <- subset(resOrdered, padj < 0.1)
resSig
```

```
## log2 fold change (MLE): condition treated vs untreated 
## Wald test p-value: condition treated vs untreated 
## DataFrame with 1069 rows and 6 columns
##              baseMean log2FoldChange     lfcSE      stat       pvalue
##             <numeric>      <numeric> <numeric> <numeric>    <numeric>
## FBgn0039155   730.992       -4.61695 0.1667827  -27.6824 1.13645e-168
## FBgn0025111  1504.272        2.90205 0.1261989   22.9958 5.13139e-117
## FBgn0029167  3709.741       -2.19491 0.0957474  -22.9240 2.68013e-116
## FBgn0003360  4344.597       -3.17683 0.1416166  -22.4326 1.89322e-111
## FBgn0035085   638.757       -2.55819 0.1359972  -18.8106  6.18406e-79
## ...               ...            ...       ...       ...          ...
## FBgn0003890 1644.2002      -0.495185  0.199269  -2.48501    0.0129549
## FBgn0010053  178.3259      -0.516894  0.208034  -2.48466    0.0129675
## FBgn0034997 5125.8312       0.250348  0.100801   2.48360    0.0130063
## FBgn0004359   84.1178       0.646359  0.260308   2.48306    0.0130260
## FBgn0027604  283.2465      -0.578675  0.233077  -2.48277    0.0130366
##                     padj
##                <numeric>
## FBgn0039155 9.25983e-165
## FBgn0025111 2.09053e-113
## FBgn0029167 7.27922e-113
## FBgn0003360 3.85649e-108
## FBgn0035085  1.00775e-75
## ...                  ...
## FBgn0003890    0.0991143
## FBgn0010053    0.0991176
## FBgn0034997    0.0993210
## FBgn0004359    0.0993659
## FBgn0027604    0.0993659
```