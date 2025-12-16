# Sample-/gene-dependent normalization factors

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

In some experiments, there might be gene-dependent dependencies which vary across samples. For instance, GC-content bias or length bias might vary across samples coming from different labs or processed at different times. We use the terms *normalization factors* for a gene x sample matrix, and *size factors* for a single number per sample. Incorporating normalization factors, the mean parameter \(\mu\_{ij}\) becomes:

\[ \mu\_{ij} = NF\_{ij} q\_{ij} \]

with normalization factor matrix *NF* having the same dimensions as the counts matrix *K*. This matrix can be incorporated as shown below. We recommend providing a matrix with row-wise geometric means of 1, so that the mean of normalized counts for a gene is close to the mean of the unnormalized counts. This can be accomplished by dividing out the current row geometric means.

```
normFactors <- normFactors / exp(rowMeans(log(normFactors)))
normalizationFactors(dds) <- normFactors
```

These steps then replace *estimateSizeFactors* which occurs within the *DESeq* function. The *DESeq* function will look for pre-existing normalization factors and use these in the place of size factors (and a message will be printed confirming this).

The methods provided by the [cqn](http://bioconductor.org/packages/cqn) or [EDASeq](http://bioconductor.org/packages/EDASeq) packages can help correct for GC or length biases. They both describe in their vignettes how to create matrices which can be used by DESeq2. From the formula above, we see that normalization factors should be on the scale of the counts, like size factors, and unlike offsets which are typically on the scale of the predictors (i.e.Â the logarithmic scale for the negative binomial GLM). At the time of writing, the transformation from the matrices provided by these packages should be:

```
cqnOffset <- cqnObject$glm.offset
cqnNormFactors <- exp(cqnOffset)
EDASeqNormFactors <- exp(-1 * EDASeqOffset)
```