# Independent hypothesis weighting

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

A generalization of the idea of *p* value filtering is to *weight* hypotheses to optimize power. A Bioconductor package, [IHW](http://bioconductor.org/packages/IHW), is available that implements the method of *Independent Hypothesis Weighting* (Ignatiadis et al. 2016). Here we show the use of *IHW* for *p* value adjustment of DESeq2 results. For more details, please see the vignette of the [IHW](http://bioconductor.org/packages/IHW) package. The *IHW* result object is stored in the metadata.

**Note:** If the results of independent hypothesis weighting are used in published research, please cite:

> Ignatiadis, N., Klaus, B., Zaugg, J.B., Huber, W. (2016) Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nature Methods*, **13**:7. [10.1038/nmeth.3885](http://dx.doi.org/10.1038/nmeth.3885)

```
# (unevaluated code chunk)
library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult
```

For advanced users, note that all the values calculated by the DESeq2 package are stored in the *DESeqDataSet* object or the *DESeqResults* object, and access to these values is discussed [below](#access).