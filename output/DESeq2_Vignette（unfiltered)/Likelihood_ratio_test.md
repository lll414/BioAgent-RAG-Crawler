# Likelihood ratio test

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

DESeq2 offers two kinds of hypothesis tests: the Wald test, where we use the estimated standard error of a log2 fold change to test if it is equal to zero, and the likelihood ratio test (LRT). The LRT examines two models for the counts, a *full* model with a certain number of terms and a *reduced* model, in which some of the terms of the *full* model are removed. The test determines if the increased likelihood of the data using the extra terms in the *full* model is more than expected if those extra terms are truly zero.

The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the *deviance* captures the difference in likelihood between a full and a reduced model.

The likelihood ratio test can be performed by specifying `test="LRT"` when using the *DESeq* function, and providing a reduced design formula, e.g.Â one in which a number of terms from `design(dds)` are removed. The degrees of freedom for the test is obtained from the difference between the number of parameters in the two models. A simple likelihood ratio test, if the full design was `~condition` would look like:

```
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)
```

If the full design contained other variables, such as a batch variable, e.g.Â `~batch + condition` then the likelihood ratio test would look like:

```
dds <- DESeq(dds, test="LRT", reduced=~batch)
res <- results(dds)
```