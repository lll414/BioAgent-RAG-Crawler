Source URL: https://bioconductor.org/books/release/OSTA/pages/img-intermediate-processing.html
Date Scraped: 2025-12-15

---

1. [Imaging-based platforms](https://bioconductor.org/books/release/OSTA/pages/img-introduction.html)
2. [20  Intermediate processing](https://bioconductor.org/books/release/OSTA/pages/img-intermediate-processing.html)

# 20  Intermediate processing

## 20.1 Preamble

### 20.1.1 Introduction

In this demo, we will continue processing the Xenium dataset from a human breast cancer biopsy section collected by Janesick et al. ([2023](https://bioconductor.org/books/release/OSTA/pages/img-intermediate-processing.html#ref-Janesick2023-high-res)), which has already undergone quality control and minimal filtering; see [Chapter 19](https://bioconductor.org/books/release/OSTA/pages/img-quality-control.html). The processing steps carried out here – namely, normalization, dimension reduction, and clustering – lay the foundating for a variety of downstream analysis tasks that will be covered in the next chapters. Here, we will run standard non-spatial and spatially-aware approaches; for a more comprehensive overview of methodology and tools, see [Chapter 27](https://bioconductor.org/books/release/OSTA/pages/ind-dimensionality-reduction.html) on dimension reduction and [Chapter 28](https://bioconductor.org/books/release/OSTA/pages/ind-clustering.html) on clustering.

### 20.1.2 Dependencies

Code

```
library(scran)
library(scater)
library(igraph)
library(Banksy)
library(ggplot2)
library(ggspavis)
library(patchwork)
library(OSTA.data)
library(SpatialExperiment)
# set seed for random number generation
# in order to make results reproducible
set.seed(20000229)
# load data from preceding 
# chapter (post quality control)
(spe <- readRDS("img-spe_qc.rds"))
```

```
##  class: SpatialExperiment 
##  dim: 313 140268 
##  metadata(0):
##  assays(1): counts
##  rownames(313): ABCC11 ACTA2 ... ZEB2 ZNF562
##  rowData names(3): ID Symbol Type
##  colnames(140268): 2 3 ... 167779 167780
##  colData names(12): cell_id transcript_counts ... detected total
##  reducedDimNames(0):
##  mainExpName: NULL
##  altExpNames(0):
##  spatialCoords names(2) : x_centroid y_centroid
##  imgData names(1): sample_id
```

Before getting started, we will retrieve the authors’ cell type labels, which were obtained by transferring scFFPE-seq annotations (supervised); these comprise 20 subpopulations.

Code

```
# get annotations from 'BiocFileCache'
# (data has been retrieved already)
id <- "Xenium_HumanBreast1_Janesick"
pa <- OSTA.data_load(id, mol=FALSE)
dir.create(td <- tempfile())
unzip(pa, "annotation.csv", exdir=td)
df <- read.csv(list.files(td, full.names=TRUE))
# add annotations as cell metadata
cs <- match(spe$cell_id, df$Barcode)
spe$Label <- df$Annotation[cs]
```

## 20.2 Normalization

Library size-based normalization, as typically used for scRNA-seq data, has been shown to be problematic for ST data, especially so for targeted panels underlying current commercial imaging-based ST platforms ([Atta et al. 2024](https://bioconductor.org/books/release/OSTA/pages/img-intermediate-processing.html#ref-Atta2024-gene-count); [Bhuva et al. 2024](https://bioconductor.org/books/release/OSTA/pages/img-intermediate-processing.html#ref-Bhuva2024-library-size)). For lack of a better approach, we here use standard log-library size normalization. We caution readers, however, to keep an eye out in the literature for attempts to provide a better strategy.

See also [Chapter 26](https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html).

Code

```
spe <- logNormCounts(spe)
```

## 20.3 Feature selection

At this stage, we would typically perform selection of (e.g., highly variable) features; see [Chapter 11](https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html). The dataset at hand, however, is targeted and relatively low-plex, so that features have already been selected by design (e.g., different targets will be included in immuno-oncology and neuroscience panels).

## 20.4 Dimension Reduction

As a baseline, we will perform principal component analysis (PCA), which underlies many standard scRNA-seq analysis pipelines, such as (spatially unaware) graph-based clustering based on a shared nearest neighbor (SNN) graph and the Leiden or Louvain algorithm for community detection.

Code

```
spe <- runPCA(spe, ncomponents=20)
```

For comparison, we will also perform spatially-aware dimension reduction with *[BANKSY](https://bioconductor.org/packages/3.22/BANKSY)* ([Singhal et al. 2024](https://bioconductor.org/books/release/OSTA/pages/img-intermediate-processing.html#ref-Singhal2024-BANKSY)); see [Chapter 27](https://bioconductor.org/books/release/OSTA/pages/ind-dimensionality-reduction.html).

Code

```
spe <- computeBanksy(spe, assay_name="logcounts")
spe <- runBanksyPCA(spe, npcs=20, lambda=0.2)
```

To not confuse different types of PCs, we rename `reducedDims` to end in `_sp` and `_tx` for spatially aware and unaware results, respectively.

Code

```
reducedDimNames(spe) <- c("PCA_tx", "PCA_sp")
```

## 20.5 Clustering

Here, we perform standard graph-based clustering by (i) constructing a shared nearest neighbor (SNN) graph and (ii) using the Leiden algorithm for community detection. By basing the SNN graph on standard and `BANKSY` PCs, respectively, we can obtain non-spatial as well as spatially aware assignments:

Code

```
pcs <- c(Leiden="PCA_tx", Banksy="PCA_sp")
for (. in names(pcs)) {
    # build cellular shared nearest-neighbor (SNN) graph
    g <- buildSNNGraph(spe, use.dimred=pcs[.], type="jaccard", k=20)
    # cluster using Leiden community detection algorithm
    k <- cluster_leiden(g, objective_function="modularity", resolution=0.8)
    spe[[.]] <- factor(k$membership)
}
```

## 20.6 Visualization

Let’s visualize the assignment obtains from non-spatial and spatially aware clustering:

Code

```
spe$in_tissue <- 1; spe$x_centroid <- spe$y_centroid <- NULL
lapply(c("Label", "Leiden", "Banksy"), \(.) {
    plotCoords(spe, annotate=., point_size = 0.1)
}) |>
    wrap_plots(nrow=1) &
    theme(legend.key.size=unit(0, "lines")) &
    scale_color_manual(values=unname(pals::trubetskoy()))
```

![](img-intermediate-processing_files/figure-html/plt-clu-xy-1.png)

## 20.7 Appendix

### Save data

Code

```
colLabels(spe) <- spe$Banksy
saveRDS(spe, "img-spe_cl.rds")
```

Back to top