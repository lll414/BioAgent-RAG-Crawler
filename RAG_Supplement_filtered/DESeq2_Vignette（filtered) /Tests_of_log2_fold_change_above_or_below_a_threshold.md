# Tests of log2 fold change above or below a threshold

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

It is also possible to provide thresholds for constructing Wald tests of significance. Two arguments to the *results* function allow for threshold-based Wald tests: `lfcThreshold`, which takes a numeric of a non-negative threshold value, and `altHypothesis`, which specifies the kind of test. Note that the *alternative hypothesis* is specified by the user, i.e.Â those genes which the user is interested in finding, and the test provides *p* values for the null hypothesis, the complement of the set defined by the alternative. The `altHypothesis` argument can take one of the following four values, where \(\beta\) is the log2 fold change specified by the `name` argument, and \(x\) is the `lfcThreshold`.

* `greaterAbs` - \(|\beta| > x\) - tests are two-tailed
* `lessAbs` - \(|\beta| < x\) - *p* values are the maximum of the upper and lower tests
* `greater` - \(\beta > x\)
* `less` - \(\beta < -x\)

The four possible values of `altHypothesis` are demonstrated in the following code and visually by MA-plots in the following figures.

```
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
```
