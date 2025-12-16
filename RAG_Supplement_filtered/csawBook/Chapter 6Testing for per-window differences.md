Source URL: https://bioconductor.org/books/release/csawBook/chap-stats.html
Date Scraped: 2025-12-15

---

# Chapter 6 Testing for per-window differences



## 6.1 Overview

Low counts per window are typically observed in ChIP-seq datasets, even for genuine binding sites.
Any statistical analysis to identify DB sites must be able to handle discreteness in the data.
Count-based models are ideal for this purpose.
In this guide, the quasi-likelihood (QL) framework in the *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* package is used (Lund et al. [2012](#ref-lund2012)).
Counts are modelled using NB distributions that account for overdispersion between biological replicates (Robinson and Smyth [2008](#ref-robinson2008)).
Each window can then be tested for significant DB between conditions.

Of course, any statistical method can be used if it is able to accept a count matrix and a vector of normalization factors (or more generally, a matrix of offsets).
The choice of *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* is primarily motivated by its performance relative to some published alternatives (Law et al. [2014](#ref-law2014))[2](#fn2).

## 6.2 Setting up for *[edgeR](https://bioconductor.org/packages/3.22/edgeR)*

A `DGEList` object is first constructed from the `SummarizedExperiment` contaiing our filtered count matrix.
If normalization factors or offsets are present in the `RangedSummarizedExperiment` object – see Chapter~[5](chap-norm.html#chap-norm) – they will automatically be inserted into the `DGEList`.
Otherwise, they can be manually passed to the `asDGEList()` function.
If offsets are available, they will generally override the normalization factors in the downstream *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* analysis.

View set-up code

```
library(csaw)
y <- asDGEList(filtered.data)
```

The experimental design is described by a design matrix.
In this case, the only relevant factor is the cell type of each sample.
A generalized linear model (GLM) will be fitted to the counts for each window using the specified design (McCarthy, Chen, and Smyth [2012](#ref-mccarthy2012)).
This provides a general framework for the analysis of complex experiments with multiple factors.
In this case, our design matrix contains an intercept representing the average abundance in the ESC group,
plus a `cell.type` coefficient representing the log-fold change of the TN group over the ESC group.

```
cell.type <- sub("NF-YA ([^ ]+) .*", "\\1", head(tf.data$Description, -1))
cell.type
```

```
## [1] "ESC" "ESC" "TN"  "TN"
```

```
design <- model.matrix(~factor(cell.type))
colnames(design) <- c("intercept", "cell.type")
design
```

```
##   intercept cell.type
## 1         1         0
## 2         1         0
## 3         1         1
## 4         1         1
## attr(,"assign")
## [1] 0 1
## attr(,"contrasts")
## attr(,"contrasts")$`factor(cell.type)`
## [1] "contr.treatment"
```

Readers are referred to the user’s guide in *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* for more details on parametrization of the design matrix.

## 6.3 Estimating the dispersions

### 6.3.1 Stabilising estimates with empirical Bayes

Under the QL framework, both the QL and NB dispersions are used to model biological variability in the data (Lund et al. [2012](#ref-lund2012)).
The former ensures that the NB mean-variance relationship is properly specified with appropriate contributions from the Poisson and Gamma components.
The latter accounts for variability and uncertainty in the dispersion estimate.
However, limited replication in most ChIP-seq experiments means that each window does not contain enough information for precise estimation of either dispersion.

This problem is overcome in *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* by sharing information across windows.
For the NB dispersions, a mean-dispersion trend is fitted across all windows to model the mean-variance relationship (McCarthy, Chen, and Smyth [2012](#ref-mccarthy2012)).
The raw QL dispersion for each window is estimated after fitting a GLM with the trended NB dispersion.
Another mean-dependent trend is fitted to the raw QL estimates.  
An empirical Bayes (EB) strategy is then used to stabilize the raw QL dispersion estimates by shrinking them towards the second trend (Lund et al. [2012](#ref-lund2012)).
The ideal amount of shrinkage is determined from the variability of the dispersions.

```
library(edgeR)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0592  0.1017  0.1047  0.1043  0.1076  0.1087
```

```
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.post)
```

```
## Length  Class   Mode 
##      0   NULL   NULL
```

The effect of EB stabilisation can be visualized by examining the biological coefficient of variation (for the NB dispersion) and the quarter-root deviance (for the QL dispersion).
Plot such as those in Figure [6.1](chap-stats.html#fig:nfya-disp-plot) can also be used to decide whether the fitted trend is appropriate.
Sudden irregulaties may be indicative of an underlying structure in the data which cannot be modelled with the mean-dispersion trend.
Discrete patterns in the raw dispersions are indicative of low counts and suggest that more aggressive filtering is required.

```
par(mfrow=c(1,2))
o <- order(y$AveLogCPM)
plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological coefficient of variation"))
plotQLDisp(fit)
```

![Fitted trend in the NB dispersion (left) or QL dispersion (right) as a function of the average abundance for each window. For the NB dispersion, the square root is shown as the biological coefficient of variation. For the QL dispersion, the shrunken estimate is also shown for each window.](modelling_files/figure-html/nfya-disp-plot-1.png)

Figure 6.1: Fitted trend in the NB dispersion (left) or QL dispersion (right) as a function of the average abundance for each window. For the NB dispersion, the square root is shown as the biological coefficient of variation. For the QL dispersion, the shrunken estimate is also shown for each window.

For most sequencing count data, we expect to see a decreasing trend that plateaus with increasing average abundance.
This reflects the greater reliability of large counts, where the effects of stochasticity and technical artifacts (e.g., mapping errors, PCR duplicates) are averaged out.
In some cases, a strong trend may also be observed where the NB dispersion drops sharply with increasing average abundance.
It is difficult to accurately fit an empirical curve to these strong trends, and as a consequence, the dispersions at high abundances may be overestimated.
Filtering of low-abundance regions (as described in Chapter [4](chap-filter.html#chap-filter)) provides some protection by removing the strongest part of the trend.
This has an additional benefit of removing those tests that have low power due to the magnitude of the dispersions.

```
relevant <- rowSums(assay(data)) >= 20 # weaker filtering than 'filtered.data'
yo <- asDGEList(data[relevant,], norm.factors=filtered.data$norm.factors)
yo <- estimateDisp(yo, design)
oo <- order(yo$AveLogCPM)

plot(yo$AveLogCPM[oo], sqrt(yo$trended.dispersion[oo]), type="l", lwd=2,
     ylim=c(0, max(sqrt(yo$trended))), xlab=expression("Ave."~Log[2]~"CPM"), 
     ylab=("Biological coefficient of variation"))
lines(y$AveLogCPM[o], sqrt(y$trended[o]), lwd=2, col="grey")
legend("topright", c("raw", "filtered"), col=c("black", "grey"), lwd=2)
```

![Fitted trend in the NB dispersions before (black) and after (grey) removing low-abundance windows.](modelling_files/figure-html/nfya-disp-rmcheck-1.png)

Figure 6.2: Fitted trend in the NB dispersions before (black) and after (grey) removing low-abundance windows.

Note that only the trended dispersion will be used in the downstream steps – the common and tagwise values are only shown in Figure [6.1](chap-stats.html#fig:nfya-disp-plot) for diagnostic purposes.
Specifically, the common BCV provides an overall measure of the variability in the data set, averaged across all windows.
The tagwise BCVs should also be dispersed above and below the fitted trend, indicating that the fit was successful.

### 6.3.2 Modelling variable dispersions between windows

Any variability in the dispersions across windows is modelled in *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* by the prior degrees of freedom (d.f.).
A large value for the prior d.f. indicates that the variability is low.
This means that more EB shrinkage can be performed to reduce uncertainty and maximize power.
However, strong shrinkage is not appropriate if the dispersions are highly variable.
Fewer prior degrees of freedom (and less shrinkage) are required to maintain type I error control.

```
summary(fit$df.prior)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.386  53.732  53.732  53.677  53.732  53.732
```

On occasion, the estimated prior degrees of freedom will be infinite.
This is indicative of a strong batch effect where the dispersions are consistently large.
A typical example involves uncorrected differences in IP efficiency across replicates.
In severe cases, the trend may fail to pass through the bulk of points as the variability is too low to be properly modelled in the QL framework.
This problem is usually resolved with appropriate normalization.

Note that the prior degrees of freedom should be robustly estimated (Phipson et al. [2016](#ref-phipson2016)).
Obviously, this protects against large positive outliers (e.g., highly variable windows) but it also protects against near-zero dispersions at low counts.
These will manifest as large negative outliers after a log transformation step during estimation (Smyth [2004](#ref-smyth2004)).
Without robustness, incorporation of these outliers will inflate the observed variability in the dispersions.
This results in a lower estimated prior d.f. and reduced DB detection power.

## 6.4 Testing for DB windows

We identify windows with significant differential binding with respect to specific factors of interest in our design matrix.
In the QL framework, pp-values are computed using the QL F-test (Lund et al. [2012](#ref-lund2012)).
This is more appropriate than using the likelihood ratio test as the F-test accounts for uncertainty in the dispersion estimates.
Associated statistics such as log-fold changes and log-counts per million are also computed for each window.

```
results <- glmQLFTest(fit, contrast=c(0, 1))
head(results$table)
```

```
##   logFC  logCPM     F   PValue
## 1 1.238  0.5410 5.102 0.027837
## 2 1.408  1.3780 7.872 0.006901
## 3 1.574  1.3513 9.335 0.003446
## 4 1.091  0.8299 3.842 0.055004
## 5 1.174 -0.3423 3.498 0.066704
## 6 1.095 -0.2668 2.818 0.098842
```

The null hypothesis here is that the cell type has no effect.
The `contrast` argument in the `glmQLFTest()` function specifies which factors are of interest[3](#fn3).
In this case, a contrast of `c(0, 1)` defines the null hypothesis as `0*intercept + 1*cell.type = 0`, i.e., that the log-fold change between cell types is zero.
DB windows can then be identified by rejecting the null.
Users may find it more intuitive to express the contrast through *[limma](https://bioconductor.org/packages/3.22/limma)*’s `makeContrasts()` function,
or by directly supplying the names of the relevant coefficients to `glmQLFTest()`, as shown below.

```
colnames(design)
```

```
## [1] "intercept" "cell.type"
```

```
# Same as above.
results2 <- glmQLFTest(fit, coef="cell.type")
```

The log-fold change for each window is similarly defined from the contrast as `0*intercept + 1*cell.type`, i.e., the value of the `cell.type` coefficient.
Recall that this coefficient represents the log-fold change of the TN group over the ESC group.
Thus, in our analysis, positive log-fold changes represent increase in binding of TN over ESC, and vice versa for negative log-fold changes.
One could also define the contrast as `c(0, -1)`, in which case the interpretation of the log-fold changes would be reversed.

Once the significance statistics have been calculated, they can be stored in row metadata of the `RangedSummarizedExperiment` object.
This ensures that the statistics and coordinates are processed together, e.g., when subsetting to select certain windows.

```
rowData(filtered.data) <- cbind(rowData(filtered.data), results$table)
```

## 6.5 What to do without replicates

Designing a ChIP-seq experiment without any replicates is strongly discouraged.
Without replication, the reproducibility of any findings cannot be determined.
Nonetheless, it may be helpful to salvage some information from datasets that lack replicates.
This is done by supplying a “reasonable” value for the NB dispersion during GLM fitting (e.g., 0.05 - 0.1, based on past experience).
DB windows are then identified using the likelihood ratio test.

```
fit.norep <- glmFit(y, design, dispersion=0.05)
results.norep <- glmLRT(fit.norep, contrast=c(0, 1))
head(results.norep$table)
```

```
##   logFC  logCPM     LR    PValue
## 1 1.237  0.5410  8.407 3.738e-03
## 2 1.407  1.3780 13.170 2.845e-04
## 3 1.573  1.3513 16.102 6.002e-05
## 4 1.089  0.8299  7.159 7.458e-03
## 5 1.173 -0.3423  5.506 1.896e-02
## 6 1.093 -0.2668  4.993 2.546e-02
```

Obviously, this approach has a number of pitfalls.
The lack of replicates means that the biological variability in the data cannot be modelled.
Thus, it becomes impossible to gauge the sensibility of the supplied NB dispersions in the analysis.
Another problem is spurious DB due to inconsistent PCR duplication between libraries.
Normally, inconsistent duplication results in a large QL dispersion for the affected window, such that significance is downweighted.
However, estimation of the QL dispersion is not possible without replicates.
This means that duplicates may need to be removed to protect against false positives.

## 6.6 Examining replicate similarity with MDS plots

As a quality control measure, the window counts can be used to examine the similarity of replicates through multi-dimensional scaling (MDS) plots (Figure [6.3](chap-stats.html#fig:nfya-mds-qc)).
The distance between each pair of libraries is computed as the square root of the mean squared log-fold change across the top set of bins with the highest absolute log-fold changes.
A small top set visualizes the most extreme differences whereas a large set visualizes overall differences.
Checking a range of `top` values may be useful when the scope of DB is unknown.
Again, counting with large bins is recommended as fold changes will be undefined in the presence of zero counts.

```
par(mfrow=c(2,2), mar=c(5,4,2,2))
adj.counts <- cpm(y, log=TRUE)
for (top in c(100, 500, 1000, 5000)) {
    plotMDS(adj.counts, main=top, col=c("blue", "blue", "red", "red"),
        labels=c("es.1", "es.2", "tn.1", "tn.2"), top=top)
}
```

![MDS plots computed with varying numbers of top windows with the strongest log-fold changes between libaries. In each plot, each library is marked with its name and colored according to its cell type.](modelling_files/figure-html/nfya-mds-qc-1.png)

Figure 6.3: MDS plots computed with varying numbers of top windows with the strongest log-fold changes between libaries. In each plot, each library is marked with its name and colored according to its cell type.

Replicates from different groups should form separate clusters in the MDS plot, as observed above.
This indicates that the results are reproducible and that the effect sizes are large.
Mixing between replicates of different conditions indicates that the biological difference has no effect on protein binding, or that the data is too variable for any effect to manifest.
Any outliers should also be noted as their presence may confound the downstream analysis.
In the worst case, outlier samples may need to be removed to obtain sensible results.

## Session information

View session info

### Bibliography

Chen, Y., A. T. L. Lun, and G. K. Smyth. 2014. “Differential Expression Analysis of Complex RNA-Seq Experiments Using EdgeR.” In *Statistical Analysis of Next Generation Sequence Data*, edited by S. Datta and D. S. Nettleton. New York: Springer.

Law, C. W., Y. Chen, W. Shi, and G. K. Smyth. 2014. “Voom: precision weights unlock linear model analysis tools for RNA-seq read counts.” *Genome Biol.* 15 (2): R29.

Lund, S. P., D. Nettleton, D. J. McCarthy, and G. K. Smyth. 2012. “Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates.” *Stat. Appl. Genet. Mol. Biol.* 11 (5).

McCarthy, D. J., Y. Chen, and G. K. Smyth. 2012. “Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.” *Nucleic Acids Res.* 40 (10): 4288–97.

Phipson, B., S. Lee, I. J. Majewski, W. S. Alexander, and G. K. Smyth. 2016. “Robust Hyperparameter Estimation Protects Against Hypervariable Genes and Improves Power to Detect Differential Expression.” *Ann. Appl. Stat.* 10 (2): 946–63.

Robinson, M. D., and G. K. Smyth. 2008. “Small-sample estimation of negative binomial dispersion, with applications to SAGE data.” *Biostatistics* 9 (2): 321–32.

Smyth, G. K. 2004. “Linear models and empirical bayes methods for assessing differential expression in microarray experiments.” *Stat. Appl. Genet. Mol. Biol.* 3: Article 3.

---

2. This author’s desire to increase his h-index may also be a factor (Chen, Lun, and Smyth [2014](#ref-chen2014))![↩︎](chap-stats.html#fnref2)
3. Specification of the contrast is explained in greater depth in the *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* user’s manual.[↩︎](chap-stats.html#fnref3)