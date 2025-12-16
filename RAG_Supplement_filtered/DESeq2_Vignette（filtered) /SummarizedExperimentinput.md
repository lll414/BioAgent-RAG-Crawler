# SummarizedExperimentinput

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

If one has already created or obtained a *SummarizedExperiment*, it can be easily input into DESeq2 as follows. First we load the package containing the `airway` dataset.

```
library("airway")
data("airway")
se <- airway
```

The constructor function below shows the generation of a *DESeqDataSet* from a *RangedSummarizedExperiment* `se`.

```
library("DESeq2")
ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
ddsSE
```

```
## class: DESeqDataSet 
## dim: 63677 8 
## metadata(2): '' version
## assays(1): counts
## rownames(63677): ENSG00000000003 ENSG00000000005 ... ENSG00000273492
##   ENSG00000273493
## rowData names(10): gene_id gene_name ... seq_coord_system symbol
## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
## colData names(9): SampleName cell ... Sample BioSample
```