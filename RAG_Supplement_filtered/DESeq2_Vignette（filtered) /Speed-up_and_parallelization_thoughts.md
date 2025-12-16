# Speed-up and parallelization thoughts

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

The above steps should take less than 30 seconds for most analyses. For experiments with complex designs and many samples (e.g.Â dozens of coefficients, ~100s of samples), one may want to have faster computation than provided by the default run of `DESeq`. We have two recommendations:

1. By using the argument `fitType="glmGamPoi"`, one can leverage the faster NB GLM engine written by Constantin Ahlmann-Eltze. Note that glmGamPoiâs interface in DESeq2 requires use of `test="LRT"` and specification of a `reduced` design.
2. One can take advantage of parallelized computation. Parallelizing `DESeq`, `results`, and `lfcShrink` can be easily accomplished by loading the BiocParallel package, and then setting the following arguments: `parallel=TRUE` and `BPPARAM=MulticoreParam(4)`, for example, splitting the job over 4 cores. However, some words of advice on parallelization: first, it is recommend to filter genes where all samples have low counts, to avoid sending data unnecessarily to child processes, when those genes have low power and will be independently filtered anyway; secondly, there is often diminishing returns for adding more cores due to overhead of sending data to child processes, therefore I recommend first starting with small number of additional cores. Note that obtaining `results` for coefficients or contrasts listed in `resultsNames(dds)` is fast and will not need parallelization. As an alternative to `BPPARAM`, one can `register` cores at the beginning of an analysis, and then just specify `parallel=TRUE` to the functions when called.

```
library("BiocParallel")
register(MulticoreParam(4))
```