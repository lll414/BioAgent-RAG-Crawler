# Approach to count outliers

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

RNA-seq data sometimes contain isolated instances of very large counts that are apparently unrelated to the experimental or study design, and which may be considered outliers. There are many reasons why outliers can arise, including rare technical or experimental artifacts, read mapping problems in the case of genetically differing samples, and genuine, but rare biological events. In many cases, users appear primarily interested in genes that show a consistent behavior, and this is the reason why by default, genes that are affected by such outliers are set aside by DESeq2, or if there are sufficient samples, outlier counts are replaced for model fitting. These two behaviors are described below.

The *DESeq* function calculates, for every gene and for every sample, a diagnostic test for outliers called *Cookâs distance*. Cookâs distance is a measure of how much a single sample is influencing the fitted coefficients for a gene, and a large value of Cookâs distance is intended to indicate an outlier count. The Cookâs distances are stored as a matrix available in `assays(dds)[["cooks"]]`.

The *results* function automatically flags genes which contain a Cookâs distance above a cutoff for samples which have 3 or more replicates. The *p* values and adjusted *p* values for these genes are set to `NA`. At least 3 replicates are required for flagging, as it is difficult to judge which sample might be an outlier with only 2 replicates. This filtering can be turned off with `results(dds, cooksCutoff=FALSE)`.

With many degrees of freedom â i.,e., many more samples than number of parameters to be estimated â it is undesirable to remove entire genes from the analysis just because their data include a single count outlier. When there are 7 or more replicates for a given sample, the *DESeq* function will automatically replace counts with large Cookâs distance with the trimmed mean over all samples, scaled up by the size factor or normalization factor for that sample. This approach is conservative, it will not lead to false positives, as it replaces the outlier value with the value predicted by the null hypothesis. This outlier replacement only occurs when there are 7 or more replicates, and can be turned off with `DESeq(dds, minReplicatesForReplace=Inf)`.

The default Cookâs distance cutoff for the two behaviors described above depends on the sample size and number of parameters to be estimated. The default is to use the 99% quantile of the F(p,m-p) distribution (with *p* the number of parameters including the intercept and *m* number of samples). The default for gene flagging can be modified using the `cooksCutoff` argument to the *results* function. For outlier replacement, *DESeq* preserves the original counts in `counts(dds)` saving the replacement counts as a matrix named `replaceCounts` in `assays(dds)`. Note that with continuous variables in the design, outlier detection and replacement is not automatically performed, as our current methods involve a robust estimation of within-group variance which does not extend easily to continuous covariates. However, users can examine the Cookâs distances in `assays(dds)[["cooks"]]`, in order to perform manual visualization and filtering if necessary.

**Note on many outliers:** if there are very many outliers (e.g.Â many hundreds or thousands) reported by `summary(res)`, one might consider further exploration to see if a single sample or a few samples should be removed due to low quality. The automatic outlier filtering/replacement is most useful in situations which the number of outliers is limited. When there are thousands of reported outliers, it might make more sense to turn off the outlier filtering/replacement (*DESeq* with `minReplicatesForReplace=Inf` and *results* with `cooksCutoff=FALSE`) and perform manual inspection: First it would be advantageous to make a PCA plot as described above to spot individual sample outliers; Second, one can make a boxplot of the Cookâs distances to see if one sample is consistently higher than others (here this is not the case):

```
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
```