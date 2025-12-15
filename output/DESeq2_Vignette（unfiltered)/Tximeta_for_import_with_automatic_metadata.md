# Tximeta for import with automatic metadata

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Another Bioconductor package, [tximeta](https://bioconductor.org/packages/tximeta) (Love et al. 2020), extends *tximport*, offering the same functionality, plus the additional benefit of automatic addition of annotation metadata for commonly used transcriptomes (GENCODE, Ensembl, RefSeq for human and mouse). See the [tximeta](https://bioconductor.org/packages/tximeta) package vignette for more details. *tximeta* produces a *SummarizedExperiment* that can be loaded easily into *DESeq2* using the `DESeqDataSet` function, with an example in the *tximeta* package vignette, and below:

```
coldata <- samples
coldata$files <- files
coldata$names <- coldata$run
```

```
library("tximeta")
se <- tximeta(coldata)
ddsTxi <- DESeqDataSet(se, design = ~ condition)
```

The `ddsTxi` object here can then be used as `dds` in the following analysis steps. If *tximeta* recognized the reference transcriptome as one of those with a pre-computed hashed checksum, the `rowRanges` of the `dds` object will be pre-populated. Again, see the *tximeta* vignette for full details.