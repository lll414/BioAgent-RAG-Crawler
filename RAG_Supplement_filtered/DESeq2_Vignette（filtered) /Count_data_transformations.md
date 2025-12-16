# Count data transformations

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

In order to test for differential expression, we operate on raw counts and use discrete distributions as described in the previous section on differential expression. However for other downstream analyses â e.g.Â for visualization or clustering â it might be useful to work with transformed versions of the count data.

Maybe the most obvious choice of transformation is the logarithm. Since count values for a gene can be zero in some conditions (and non-zero in others), some advocate the use of *pseudocounts*, i.e.Â transformations of the form:

\[ y = \log\_2(n + n\_0) \]

where *n* represents the count values and \(n\_0\) is a positive constant.

In this section, we discuss two alternative approaches that offer more theoretical justification and a rational way of choosing parameters equivalent to \(n\_0\) above. One makes use of the concept of variance stabilizing transformations (VST) (Tibshirani 1988; Huber et al. 2003; Anders and Huber 2010), and the other is the *regularized logarithm* or *rlog*, which incorporates a prior on the sample differences (Love, Huber, and Anders 2014). Both transformations produce transformed data on the log2 scale which has been normalized with respect to library size or other normalization factors.

The point of these two transformations, the VST and the *rlog*, is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low. Both VST and *rlog* use the experiment-wide trend of variance over mean, in order to transform the data to remove the experiment-wide trend. Note that we do not require or desire that all the genes have *exactly* the same variance after transformation. Indeed, in a figure below, you will see that after the transformations the genes with the same mean do not have exactly the same standard deviations, but that the experiment-wide trend has flattened. It is those genes with row variance above the trend which will allow us to cluster samples into interesting groups.

**Note on running time:** if you have many samples (e.g.Â 100s), the *rlog* function might take too long, and so the *vst* function will be a faster choice. The rlog and VST have similar properties, but the rlog requires fitting a shrinkage term for each sample and each gene which takes time. See the DESeq2 paper for more discussion on the differences (Love, Huber, and Anders 2014).

### Blind dispersion estimation

The two functions, *vst* and *rlog* have an argument `blind`, for whether the transformation should be blind to the sample information specified by the design formula. When `blind` equals `TRUE` (the default), the functions will re-estimate the dispersions using only an intercept. This setting should be used in order to compare samples in a manner wholly unbiased by the information about experimental groups, for example to perform sample QA (quality assurance) as demonstrated below.

However, blind dispersion estimation is not the appropriate choice if one expects that many or the majority of genes (rows) will have large differences in counts which are explainable by the experimental design, and one wishes to transform the data for downstream analysis. In this case, using blind dispersion estimation will lead to large estimates of dispersion, as it attributes differences due to experimental design as unwanted *noise*, and will result in overly shrinking the transformed values towards each other. By setting `blind` to `FALSE`, the dispersions already estimated will be used to perform transformations, or if not present, they will be estimated using the current design formula. Note that only the fitted dispersion estimates from mean-dispersion trend line are used in the transformation (the global dependence of dispersion on mean for the entire experiment). So setting `blind` to `FALSE` is still for the most part not using the information about which samples were in which experimental group in applying the transformation.

### Extracting transformed values

These transformation functions return an object of class *DESeqTransform* which is a subclass of *RangedSummarizedExperiment*. For ~20 samples, running on a newly created `DESeqDataSet`, *rlog* may take 30 seconds, while *vst* takes less than 1 second. The running times are shorter when using `blind=FALSE` and if the function *DESeq* has already been run, because then it is not necessary to re-estimate the dispersion values. The *assay* function is used to extract the matrix of normalized values.

```
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
```

```
##              treated1  treated2  treated3 untreated1 untreated2 untreated3
## FBgn0000008  7.777746  7.984491  7.765448   7.735026   7.807720   7.997378
## FBgn0000017 11.954934 12.035700 12.028610  12.049838  12.295662  12.471723
## FBgn0000018  9.219162  9.089390  9.028791   9.374577   9.173538   9.052143
##             untreated4
## FBgn0000008   7.832458
## FBgn0000017  12.089675
## FBgn0000018   9.142848
```

### Variance stabilizing transformation

Above, we used a parametric fit for the dispersion. In this case, the closed-form expression for the variance stabilizing transformation is used by the *vst* function. If a local fit is used (option `fitType="locfit"` to *estimateDispersions*) a numerical integration is used instead. The transformed data should be approximated variance stabilized and also includes correction for size factors or normalization factors. The transformed data is on the log2 scale for large counts.

### Regularized log transformation

The function *rlog*, stands for *regularized log*, transforming the original count data to the log2 scale by fitting a model with a term for each sample and a prior distribution on the coefficients which is estimated from the data. This is the same kind of shrinkage (sometimes referred to as regularization, or moderation) of log fold changes used by *DESeq* and *nbinomWaldTest*. The resulting data contains elements defined as:

\[ \log\_2(q\_{ij}) = \beta\_{i0} + \beta\_{ij} \]

where \(q\_{ij}\) is a parameter proportional to the expected true concentration of fragments for gene *i* and sample *j* (see formula [below](#theory)), \(\beta\_{i0}\) is an intercept which does not undergo shrinkage, and \(\beta\_{ij}\) is the sample-specific effect which is shrunk toward zero based on the dispersion-mean trend over the entire dataset. The trend typically captures high dispersions for low counts, and therefore these genes exhibit higher shrinkage from the *rlog*.

Note that, as \(q\_{ij}\) represents the part of the mean value \(\mu\_{ij}\) after the size factor \(s\_j\) has been divided out, it is clear that the rlog transformation inherently accounts for differences in sequencing depth. Without priors, this design matrix would lead to a non-unique solution, however the addition of a prior on non-intercept betas allows for a unique solution to be found.

### Effects of transformations on the variance

The figure below plots the standard deviation of the transformed data, across samples, against the mean, using the shifted logarithm transformation, the regularized log transformation and the variance stabilizing transformation. The shifted logarithm has elevated standard deviation in the lower count range, and the regularized log to a lesser extent, while for the variance stabilized data the standard deviation is roughly constant along the whole dynamic range.

Note that the vertical axis in such plots is the square root of the variance over all samples, so including the variance due to the experimental conditions. While a flat curve of the square root of variance over the mean may seem like the goal of such transformations, this may be unreasonable in the case of datasets with many true differences due to the experimental conditions.

```
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
```

```
## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
## â¹ Please use tidy evaluation idioms with `aes()`.
## â¹ See also `vignette("ggplot2-in-packages")` for more information.
## â¹ The deprecated feature was likely used in the vsn package.
##   Please report the issue to the authors.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```


```
meanSdPlot(assay(vsd))
```


```
meanSdPlot(assay(rld))
```
