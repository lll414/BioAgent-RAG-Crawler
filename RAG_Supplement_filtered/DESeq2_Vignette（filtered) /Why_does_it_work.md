# Why does it work?

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

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
