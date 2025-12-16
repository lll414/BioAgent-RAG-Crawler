Source URL: https://bioconductor.org/books/release/OSTA/pages/mult-diff-spatial-patterns.html
Date Scraped: 2025-12-15

---

1. [Multiple-sample analyses](https://bioconductor.org/books/release/OSTA/pages/mult-diff-spatial-patterns.html)
2. [33  Differential spatial patterns](https://bioconductor.org/books/release/OSTA/pages/mult-diff-spatial-patterns.html)

# 33  Differential spatial patterns

## 33.1 Preamble

### 33.1.1 Introduction

Recent advances in spatial transcriptomics (ST) techniques have facilitated the generation of multi-sample datasets across various experimental conditions (e.g. healthy vs. diseased, or multiple treatments). When multi-sample, multi-condition ST data are available, differential analysis can be performed to identify genes with distinct spatial expression patterns that change between conditions either across domains or within any domain; we call these **differential spatial pattern (DSP)** genes.

A DSP gene may be spatially variable in one condition and display spatially uniform expression in another one, or it could be spatially variable in all conditions, but with distinct spatial expression patterns. Conversely, genes that show spatially uniform abundance, or that are spatially variable with the same spatial structure in all conditions, are not considered DSP, because their spatial distribution is unchanged. Another way to think about DSP is that it is the spatial analog of ‘differential state analysis’ for single cell analyses. That is, here we focus on differences in expression within any domain (or different across domains) between conditions.

Identifying DSP genes can provide valuable insights into how gene expression spatial structures change across experimental conditions, helping to understand the underlying biological processes.

### 33.1.2 Dependencies

Code

```
library(edgeR)
library(DESpace)
library(ggplot2)
library(ggspavis)
library(muSpaData)
library(patchwork)
```

### 33.1.3 Load data

In this demo, we will analyze a Stereo-seq dataset from axolotl brain tissues collected at various stages of regeneration ([Wei et al. 2022](https://bioconductor.org/books/release/OSTA/pages/mult-diff-spatial-patterns.html#ref-Wei2022-Stereo-seq-axolotl)), which includes multiple samples (i.e. serial sections) measured under multiple experimental conditions (regeneration stages).

In the original dataset, which is publicly available through the [Spatial Transcript Omics DataBase](https://db.cngb.org/stomics/datasets/STDS0000056/data), there are 5 stages with multiple sections, totaling 16 samples (3 to 4 sections per stage).

For computational reasons, we use 3 stages in this chapter – 2, 10, and 20 days post injury (DPI) – and only two sections per stage. The analysis framework uses pre-computed spatial domains generated with *[Banksy](https://bioconductor.org/packages/3.22/Banksy)* ([Singhal et al. 2024](https://bioconductor.org/books/release/OSTA/pages/mult-diff-spatial-patterns.html#ref-Singhal2024-BANKSY)) (see [Chapter 28](https://bioconductor.org/books/release/OSTA/pages/ind-clustering.html)). For preprocessing details of this dataset, see the *[muSpaData](https://bioconductor.org/packages/3.22/muSpaData)* Bioconductor package.

Code

```
spe <- Wei22_example()
spatialCoordsNames(spe) <- c("x", "y")
spe
```

```
##  class: SpatialExperiment 
##  dim: 5000 55660 
##  metadata(0):
##  assays(2): counts logcounts
##  rownames(5000): AMEX60DD009830 AMEX60DD009962 ... AMEX60DD004094
##    AMEX60DD054542
##  rowData names(2): gene_name gene_id
##  colnames(55660): CELL.17879.10DPI_1 CELL.17922.10DPI_1 ...
##    CELL.9123.2DPI_2 CELL.9124.2DPI_2
##  colData names(5): sample_id condition Banksy_smooth sdimx sdimy
##  reducedDimNames(0):
##  mainExpName: NULL
##  altExpNames(0):
##  spatialCoords names(2) : x y
##  imgData names(1): sample_id
```

Code

```
plotCoords(spe, 
    annotate="Banksy_smooth", 
    in_tissue=NULL,
    x_coord="sdimx", 
    y_coord="sdimy", 
    y_reverse=FALSE,
    sample_id="sample_id") + 
    theme(legend.key.size=unit(0, "lines")) + 
    scale_color_manual(values=unname(pals::trubetskoy()))
```

![](mult-diff-spatial-patterns_files/figure-html/plt-clu-xy-1.png)

## 33.2 DESpace

Here, we briefly demonstrate how to identify DSP genes using the *[DESpace](https://bioconductor.org/packages/3.22/DESpace)* ([Cai, Robinson, and Tiberi 2024](https://bioconductor.org/books/release/OSTA/pages/mult-diff-spatial-patterns.html#ref-Cai2024-DESpace)) Bioconductor package. For more details, see the [package vignette](https://peicai.github.io/DESpace/articles/DSP.html). In brief, spot- or cell-level counts are aggregated into pseudobulk counts for each sample and domain, and these are compared across domains and conditions. DSP focuses on testing whether (condition x domain) interaction terms in the model are different from zero. In other words, does the expression of a gene change between conditions, but change differently for some domains compared to others?

Code

```
# run 'DESpace'
dsp <- dsp_test(spe, 
    sample_col="sample_id", 
    condition_col="condition", 
    cluster_col="Banksy_smooth",
    # return full statistics 
    # from 'edgeR::glmFit/LRT'
    verbose=TRUE)
```

Code

```
# extract gene-level results
names(dsp_global <- dsp$gene_results)
```

```
##   [1] "gene_id"                          "logFC.condition20DPI.cluster_id1"
##   [3] "logFC.condition2DPI.cluster_id1"  "logFC.condition20DPI.cluster_id2"
##   [5] "logFC.condition2DPI.cluster_id2"  "logFC.condition20DPI.cluster_id3"
##   [7] "logFC.condition2DPI.cluster_id3"  "logFC.condition20DPI.cluster_id4"
##   [9] "logFC.condition2DPI.cluster_id4"  "logCPM"                          
##  [11] "LR"                               "PValue"                          
##  [13] "FDR"
```

Code

```
# count significant DSP genes 
# (at 5% FDR significance level)
table(dsp_global$FDR <= 0.05)
```

```
##  
##  FALSE  TRUE 
##   4799   201
```

In order to identify the key spatial cluster (domain) where expression changes across conditions, we use the `individual_dsp()` function, which focuses on domain-specific DSP.

Code

```
dsp_clu <- individual_dsp(spe, 
    sample_col="sample_id", 
    condition_col="condition",
    cluster_col="Banksy_smooth")
```

Code

```
# results for cluster 2
dsp_clu2 <- dsp_clu$`2`
head(dsp_clu2, n=3)
```

```
##                        gene_id logFC.condition20DPI.cluster_id2
##  AMEX60DD014721 AMEX60DD014721                       -0.5801491
##  AMEX60DD014991 AMEX60DD014991                        2.6942759
##  AMEX60DD055246 AMEX60DD055246                       -0.2661447
##                 logFC.condition2DPI.cluster_id2   logCPM       LR
##  AMEX60DD014721                        1.270049 9.582333 76.16657
##  AMEX60DD014991                        3.015124 7.330327 73.02536
##  AMEX60DD055246                       -2.298603 5.661946 60.91208
##                       PValue          FDR
##  AMEX60DD014721 2.888278e-17 1.444139e-13
##  AMEX60DD014991 1.389135e-16 3.472837e-13
##  AMEX60DD055246 5.930751e-14 9.884585e-11
```

Code

```
# top DSP gene for cluster 2
top_dsp <- dsp_clu2$gene_id[1]

# get gene symbols from Ensembl identifiers
idx <- match(top_dsp, rowData(spe)$gene_name)
(.gs <- rowData(spe)$gene_id[idx])
```

```
##  [1] "ECM1"
```

## 33.3 Visualization

DSP genes can be further investigated, for example, by plotting the spatial expression of genes across conditions, or the log2-transformed counts per million of genes across different clusters and conditions.

Code

```
.ids <- levels(spe$sample_id)
lapply(seq_along(.ids), \(.) {
    .spe <- spe[, spe$sample_id == .ids[.]]
    p <- FeaturePlot(.spe, 
        feature=top_dsp,
        platform="Stereo-seq",
        coordinates=c("sdimx", "sdimy"),
        cluster_col="Banksy_smooth", cluster="2",
        diverging=TRUE, low="gray95", high="blue")
    # ignore alignment
    free(p[[1]] + ggtitle(.ids[.])) 
}) |> wrap_plots(ncol=3)
```

![](mult-diff-spatial-patterns_files/figure-html/plt-xy-1.png)

The boxplots below show the average log-CPM for cluster 2 and for all other clusters (excluding cluster 2) across different time points. The largest difference between cluster 2 and the remaining clusters is observed at 2 DPI, with the difference decreasing over time. At 20 DPI, no significant log-CPM differences are observed between the clusters. This suggests that spatial patterns change across conditions (i.e. time stages).

Code

```
# calculate log-CPM
log_cpm <- edgeR::cpm(dsp$estimated_y, log=TRUE)
nms <- colnames(log_cpm)
# wrangling
dat <- data.frame(
    log_cpm=log_cpm[top_dsp, ],
    Banksy=factor(sub(".*_", "", nms)), 
    sample_id=sub("(_[0-9]+)$", "", nms),
    day=as.numeric(sub("([0-9]+)DPI.*", "\\1", nms)))
# visualization
ggplot(dat, aes(factor(day), log_cpm)) + 
    geom_jitter(aes(col=Banksy), size=2, width=0.1) + 
    geom_boxplot(aes(fill=ifelse(Banksy == "2", "cluster 2", "other"))) + 
    scale_x_discrete("Days post injury", breaks=unique(dat$day)) + 
    scale_fill_manual(NULL, values=c("limegreen", "gray")) + 
    labs(title=.gs, y="log2 counts per million (logCPM)") + 
    theme(legend.position="right")
```

![](mult-diff-spatial-patterns_files/figure-html/plt-cpm-1.png)

## 33.4 Appendix

Back to top