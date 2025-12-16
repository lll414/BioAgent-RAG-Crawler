# Quick start

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Here we show the most basic steps for a differential expression analysis. There are a variety of steps upstream of DESeq2 that result in the generation of counts or estimated counts for each sample, which we will discuss in the sections below. This code chunk assumes that you have a count matrix called `cts` and a table of sample information called `coldata`. The `design` indicates how to model the samples, here, that we want to measure the effect of the condition, controlling for batch differences. The two factor variables `batch` and `condition` should be columns of `coldata`.

```
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ batch + condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_trt_vs_untrt")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
```

The following starting functions will be explained below:

* If you have performed transcript quantification (with *Salmon*, *kallisto*, *RSEM*, etc.) you could import the data with *tximport*, which produces a list, and then you can use `DESeqDataSetFromTximport()`.
* If you imported quantification data with *tximeta*, which produces a *SummarizedExperiment* with additional metadata, you can then use `DESeqDataSet()`.
* If you have *htseq-count* files, you can use `DESeqDataSetFromHTSeq()`.