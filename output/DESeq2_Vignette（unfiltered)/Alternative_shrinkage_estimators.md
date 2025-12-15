# Alternative shrinkage estimators

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

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