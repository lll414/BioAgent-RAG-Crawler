# Independent filtering and multiple testing

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

### Filtering criteria

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


### Why does it work?

Consider the *p* value histogram below It shows how the filtering ameliorates the multiple testing problem â and thus the severity of a multiple testing adjustment â by removing a background set of hypotheses whose *p* values are distributed more or less uniformly in [0,1].

```
dds_demo <- makeExampleDESeqDataSet(
  n=5000, m=6,
  interceptMean=3,
  interceptSD=1,
  betaSD=1
)
dds_demo <- DESeq(dds_demo)
res_demo <- results(dds_demo)
```

```
use <- res_demo$baseMean > metadata(res_demo)$filterThreshold
h1 <- hist(res_demo$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res_demo$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
```

Histogram of p values for all tests. The area shaded in blue indicates the subset of those that pass the filtering, the area in khaki those that do not pass:

```
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))
```
