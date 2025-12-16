# Filtering criteria

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

The goal of independent filtering is to filter out those tests from the procedure that have no, or little chance of showing significant evidence, without even looking at their test statistic. Typically, this results in increased detection power at the same experiment-wide type I error. Here, we measure experiment-wide type I error in terms of the false discovery rate.

A good choice for a filtering criterion is one that

1. is statistically independent from the test statistic under the null hypothesis,
2. is correlated with the test statistic under the alternative, and
3. does not notably change the dependence structure â if there is any â between the tests that pass the filter, compared to the dependence structure between the tests before filtering.

The benefit from filtering relies on property (2), and we will explore it further below. Its statistical validity relies on property (1) â which is simple to formally prove for many combinations of filter criteria with test statistics â and (3), which is less easy to theoretically imply from first principles, but rarely a problem in practice. We refer to (Bourgon, Gentleman, and Huber 2010) for further discussion of this topic.

A simple filtering criterion readily available in the results object is the mean of normalized counts irrespective of biological condition, and so this is the criterion which is used automatically by the *results* function to perform independent filtering. Genes with very low counts are not likely to see significant differences typically due to high dispersion. For example, we can plot the \(-\log\_{10}\) *p* values from all genes over the normalized mean counts:

```
plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
```
