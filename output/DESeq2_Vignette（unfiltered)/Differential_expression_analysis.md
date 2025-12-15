# Differential expression analysis

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

The standard differential expression analysis steps are wrapped into a single function, *DESeq*. The estimation steps performed by this function are described [below](#theory), in the manual page for `?DESeq` and in the Methods section of the DESeq2 publication (Love, Huber, and Anders 2014).

Results tables are generated using the function *results*, which extracts a results table with log2 fold changes, *p* values and adjusted *p* values. With no additional arguments to *results*, the log2 fold change and Wald test *p* value will be for the **last variable** in the design formula, and if this is a factor, the comparison will be the **last level** of this variable over the **reference level** (see previous [note on factor levels](#factorlevels)). However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the `name` or `contrast` arguments of *results*.

Details about the comparison are printed to the console, directly above the results table. The text, `condition treated vs untreated`, tells you that the estimates are of the logarithmic fold change log2(treated/untreated).

```
dds <- DESeq(dds)
res <- results(dds)
res
```

```
## log2 fold change (MLE): condition treated vs untreated 
## Wald test p-value: condition treated vs untreated 
## DataFrame with 8148 rows and 6 columns
##               baseMean log2FoldChange     lfcSE       stat    pvalue      padj
##              <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
## FBgn0000008   95.28865     0.00399148  0.225010  0.0177391 0.9858470  0.996699
## FBgn0000017 4359.09632    -0.23842494  0.127094 -1.8759764 0.0606585  0.289604
## FBgn0000018  419.06811    -0.10185506  0.146568 -0.6949338 0.4870968  0.822681
## FBgn0000024    6.41105     0.21429657  0.691557  0.3098756 0.7566555  0.939146
## FBgn0000032  990.79225    -0.08896298  0.146253 -0.6082822 0.5430003  0.848881
## ...                ...            ...       ...        ...       ...       ...
## FBgn0261564   1160.028     -0.0857255  0.108354 -0.7911643 0.4288481  0.789246
## FBgn0261565    620.388     -0.2943294  0.140496 -2.0949303 0.0361772  0.206423
## FBgn0261570   3212.969      0.2971841  0.126742  2.3447877 0.0190379  0.133380
## FBgn0261573   2243.936      0.0146611  0.111365  0.1316493 0.8952617  0.977565
## FBgn0261574   4863.807      0.0179729  0.194137  0.0925784 0.9262385  0.986726
```

Note that we could have specified the coefficient or contrast we want to build a results table for, using either of the following equivalent commands:

```
res <- results(dds, name="condition_treated_vs_untreated")
res <- results(dds, contrast=c("condition","treated","untreated"))
```

One exception to the equivalence of these two commands, is that, using `contrast` will additionally set to 0 the estimated LFC in a comparison of two groups, where all of the counts in the two groups are equal to 0 (while other groups have positive counts). As this may be a desired feature to have the LFC in these cases set to 0, one can use `contrast` to build these results tables. More information about extracting specific coefficients from a fitted *DESeqDataSet* object can be found in the help page `?results`. The use of the `contrast` argument is also further discussed [below](#contrasts).

### Log fold change shrinkage for visualization and ranking

Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. To shrink the LFC, we pass the `dds` object to the function `lfcShrink`. Below we specify to use the *apeglm* method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.

We provide the `dds` object and the name or number of the coefficient we want to shrink, where the number refers to the order of the coefficient as it appears in `resultsNames(dds)`.

```
resultsNames(dds)
```

```
## [1] "Intercept"                      "condition_treated_vs_untreated"
```

```
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC
```

```
## log2 fold change (MAP): condition treated vs untreated 
## Wald test p-value: condition treated vs untreated 
## DataFrame with 8148 rows and 5 columns
##               baseMean log2FoldChange     lfcSE    pvalue      padj
##              <numeric>      <numeric> <numeric> <numeric> <numeric>
## FBgn0000008   95.28865     0.00195376  0.152654 0.9858470  0.996699
## FBgn0000017 4359.09632    -0.18810628  0.120870 0.0606585  0.289604
## FBgn0000018  419.06811    -0.06893831  0.122805 0.4870968  0.822681
## FBgn0000024    6.41105     0.01786546  0.199499 0.7566555  0.939146
## FBgn0000032  990.79225    -0.06001511  0.121962 0.5430003  0.848881
## ...                ...            ...       ...       ...       ...
## FBgn0261564   1160.028     -0.0669829 0.0976567 0.4288481  0.789246
## FBgn0261565    620.388     -0.2284564 0.1362122 0.0361772  0.206423
## FBgn0261570   3212.969      0.2395981 0.1237304 0.0190379  0.133380
## FBgn0261573   2243.936      0.0115395 0.0981689 0.8952617  0.977565
## FBgn0261574   4863.807      0.0101618 0.1417667 0.9262385  0.986726
```

Shrinkage estimation is discussed more in a [later section](#moreshrink).

### Speed-up and parallelization thoughts

The above steps should take less than 30 seconds for most analyses. For experiments with complex designs and many samples (e.g.Â dozens of coefficients, ~100s of samples), one may want to have faster computation than provided by the default run of `DESeq`. We have two recommendations:

1. By using the argument `fitType="glmGamPoi"`, one can leverage the faster NB GLM engine written by Constantin Ahlmann-Eltze. Note that glmGamPoiâs interface in DESeq2 requires use of `test="LRT"` and specification of a `reduced` design.
2. One can take advantage of parallelized computation. Parallelizing `DESeq`, `results`, and `lfcShrink` can be easily accomplished by loading the BiocParallel package, and then setting the following arguments: `parallel=TRUE` and `BPPARAM=MulticoreParam(4)`, for example, splitting the job over 4 cores. However, some words of advice on parallelization: first, it is recommend to filter genes where all samples have low counts, to avoid sending data unnecessarily to child processes, when those genes have low power and will be independently filtered anyway; secondly, there is often diminishing returns for adding more cores due to overhead of sending data to child processes, therefore I recommend first starting with small number of additional cores. Note that obtaining `results` for coefficients or contrasts listed in `resultsNames(dds)` is fast and will not need parallelization. As an alternative to `BPPARAM`, one can `register` cores at the beginning of an analysis, and then just specify `parallel=TRUE` to the functions when called.

```
library("BiocParallel")
register(MulticoreParam(4))
```

### p-values and adjusted p-values

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

### Independent hypothesis weighting

A generalization of the idea of *p* value filtering is to *weight* hypotheses to optimize power. A Bioconductor package, [IHW](http://bioconductor.org/packages/IHW), is available that implements the method of *Independent Hypothesis Weighting* (Ignatiadis et al. 2016). Here we show the use of *IHW* for *p* value adjustment of DESeq2 results. For more details, please see the vignette of the [IHW](http://bioconductor.org/packages/IHW) package. The *IHW* result object is stored in the metadata.

**Note:** If the results of independent hypothesis weighting are used in published research, please cite:

> Ignatiadis, N., Klaus, B., Zaugg, J.B., Huber, W. (2016) Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nature Methods*, **13**:7. [10.1038/nmeth.3885](http://dx.doi.org/10.1038/nmeth.3885)

```
# (unevaluated code chunk)
library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult
```

For advanced users, note that all the values calculated by the DESeq2 package are stored in the *DESeqDataSet* object or the *DESeqResults* object, and access to these values is discussed [below](#access).