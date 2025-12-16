# Exploring and exporting results

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

### MA-plot

In DESeq2, the function *plotMA* shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the *DESeqDataSet*. Points will be colored blue if the adjusted *p* value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

```
plotMA(res, ylim=c(-2,2))
```


It is more useful to visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.

```
plotMA(resLFC, ylim=c(-2,2))
```


After calling *plotMA*, one can use the function *identify* to interactively detect the row number of individual genes by clicking on the plot. One can then recover the gene identifiers by saving the resulting indices:

```
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]
```

### Alternative shrinkage estimators

The moderated log fold changes proposed by Love, Huber, and Anders (2014) use a normal prior distribution, centered on zero and with a scale that is fit to the data. The shrunken log fold changes are useful for ranking and visualization, without the need for arbitrary filters on low count genes. The normal prior can sometimes produce too strong of shrinkage for certain datasets. In DESeq2 version 1.18, we include two additional adaptive shrinkage estimators, available via the `type` argument of `lfcShrink`. For more details, see `?lfcShrink`

The options for `type` are:

* `apeglm` is the adaptive t prior shrinkage estimator from the [apeglm](http://bioconductor.org/packages/apeglm) package (Zhu, Ibrahim, and Love 2018). As of version 1.28.0, it is the default estimator.
* `ashr` is the adaptive shrinkage estimator from the [ashr](https://github.com/stephens999/ashr) package (Stephens 2016). Here DESeq2 uses the ashr option to fit a mixture of Normal distributions to form the prior, with `method="shrinkage"`.
* `normal` is the the original DESeq2 shrinkage estimator, an adaptive Normal distribution as prior.

If the shrinkage estimator `apeglm` is used in published research, please cite:

> Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics*. [10.1093/bioinformatics/bty895](https://doi.org/10.1093/bioinformatics/bty895)

If the shrinkage estimator `ashr` is used in published research, please cite:

> Stephens, M. (2016) False discovery rates: a new deal. *Biostatistics*, **18**:2. [10.1093/biostatistics/kxw041](https://doi.org/10.1093/biostatistics/kxw041)

In the LFC shrinkage code above, we specified `coef="condition_treated_vs_untreated"`. We can also just specify the coefficient by the order that it appears in `resultsNames(dds)`, in this case `coef=2`. For more details explaining how the shrinkage estimators differ, and what kinds of designs, contrasts and output is provided by each, see the [extended section on shrinkage estimators](#moreshrink).

```
resultsNames(dds)
```

```
## [1] "Intercept"                      "condition_treated_vs_untreated"
```

```
# because we are interested in treated vs untreated, we set 'coef=2'
resNorm <- lfcShrink(dds, coef=2, type="normal")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
```

```
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
```


**Note:** We have sped up the `apeglm` method so it takes roughly about the same amount of time as `normal`, e.g.Â ~5 seconds for the `pasilla` dataset of ~10,000 genes and 7 samples. If fast shrinkage estimation of LFC is needed, *but the posterior standard deviation is not needed*, setting `apeMethod="nbinomC"` will produce a ~10x speedup, but the `lfcSE` column will be returned with `NA`. A variant of this fast method, `apeMethod="nbinomC*"` includes random starts.

**Note:** If there is unwanted variation present in the data (e.g.Â batch effects) it is always recommend to correct for this, which can be accommodated in DESeq2 by including in the design any known batch variables or by using functions/packages such as `svaseq` in [sva](http://bioconductor.org/packages/sva) (Leek 2014) or the `RUV` functions in [RUVSeq](http://bioconductor.org/packages/RUVSeq) (Risso et al. 2014) to estimate variables that capture the unwanted variation. In addition, the ashr developers have a [specific method](https://github.com/dcgerard/vicar) for accounting for unwanted variation in combination with ashr (Gerard and Stephens 2017).

### Plot counts

It can also be useful to examine the counts of reads for a single gene across the groups. A simple function for making this plot is *plotCounts*, which normalizes counts by the estimated size factors (or normalization factors if these were used) and adds a pseudocount of 1/2 to allow for log scale plotting. The counts are grouped by the variables in `intgroup`, where more than one variable can be specified. Here we specify the gene which had the smallest *p* value from the results table created above. You can select the gene to plot by rowname or by numeric index.

```
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
```


For customized plotting, an argument `returnData` specifies that the function should only return a *data.frame* for plotting with *ggplot*.

```
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
```


### More information on results columns

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

### Rich visualization and reporting of results

**regionReport** An HTML and PDF summary of the results with plots can also be generated using the [regionReport](http://bioconductor.org/packages/regionReport) package. The *DESeq2Report* function should be run on a *DESeqDataSet* that has been processed by the *DESeq* function. For more details see the manual page for *DESeq2Report* and an example vignette in the [regionReport](http://bioconductor.org/packages/regionReport) package.

**Glimma** Interactive visualization of DESeq2 output, including MA-plots (also called MD-plots) can be generated using the [Glimma](http://bioconductor.org/packages/Glimma) package. See the manual page for *glMDPlot.DESeqResults*.

**pcaExplorer** Interactive visualization of DESeq2 output, including PCA plots, boxplots of counts and other useful summaries can be generated using the [pcaExplorer](http://bioconductor.org/packages/pcaExplorer) package. See the *Launching the application* section of the package vignette.

**iSEE** Provides functions for creating an interactive Shiny-based graphical user interface for exploring data stored in SummarizedExperiment objects, including row- and column-level metadata. Particular attention is given to single-cell data in a SingleCellExperiment object with visualization of dimensionality reduction results. [iSEE](https://bioconductor.org/packages/iSEE) is on Bioconductor. An example wrapper function for converting a *DESeqDataSet* to a SingleCellExperiment object for use with *iSEE* can be found at the following gist, written by Federico Marini:

* <https://gist.github.com/federicomarini/4a543eebc7e7091d9169111f76d59de1>

The [iSEEde](https://bioconductor.org/packages/iSEEde) package provides additional panels that facilitate the interactive visualisation of differential expression results in iSEE applications.

**DEvis** DEvis is a powerful, integrated solution for the analysis of differential expression data. This package includes an array of tools for manipulating and aggregating data, as well as a wide range of customizable visualizations, and project management functionality that simplify RNA-Seq analysis and provide a variety of ways of exploring and analyzing data. *DEvis* can be found on [CRAN](https://cran.r-project.org/package=DEVis) and [GitHub](https://github.com/price0416/DEvis).

### Exporting results to CSV files

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