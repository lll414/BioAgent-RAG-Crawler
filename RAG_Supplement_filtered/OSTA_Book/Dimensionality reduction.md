Source URL: https://bioconductor.org/books/release/OSTA/pages/ind-dimensionality-reduction.html
Date Scraped: 2025-12-15

---

1. [Platform-independent analyses](https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html)
2. [27  Dimensionality reduction](https://bioconductor.org/books/release/OSTA/pages/ind-dimensionality-reduction.html)

# 27  Dimensionality reduction

## 27.1 Preamble

### 27.1.1 Introduction

In single-cell omics data analysis, dimensionality reduction (DR) techniques are often categorized as linear (e.g., multi-dimensional scaling (MDS), linear discriminant analysis (LDA), principal component analysis (PCA)), or non-linear (e.g., t-distributed stochastic neighbor embedding (t-SNE), uniform manifold approximation and projection (UMAP)); see [OSCA](https://bioconductor.org/books/release/OSCA.basic/dimensionality-reduction.html).

In the context of ST data, we would instead like to distinguish between spatially aware and non-spatial methods, i.e., whether or not DR incorporates physical locations or is based on molecular profiles only. Spatially aware DR may be used as input for standard single-cell clustering approaches (say, SNN-graph based Leiden), while non-spatial DR may be used for spatially aware clustering. While these approach may differ in methodological detail, they are arguably similar conceptually. We will try and demonstrate as much here.

### 27.1.2 Dependencies

Code

```
library(dplyr)
library(Banksy)
library(scater)
library(ggplot2)
library(OSTA.data)
library(patchwork)
library(SpatialExperiment)
# set seed for random number generation
# in order to make results reproducible
set.seed(20000229)
# load Xenium data post quality control
spe <- readRDS("img-spe_qc.rds")
```

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

## 27.2 Principal component analysis (PCA)

### 27.2.1 Non-spatial

As a baseline, we will perform principal component analysis (PCA), which underlies many standard scRNA-seq analysis pipelines, such as (spatially unaware) graph-based clustering based on a shared nearest neighbor (SNN) graph and the Leiden or Louvain algorithm for community detection.

A standard approach is to apply PCA to the set of top HVGs, and retain a subset of PCs for subsequent steps. This is done for two main reasons: (i) to reduce noise due to random variation in expression of biologically uninformative genes, which are assumed to have expression patterns independent of each other, and (ii) to improve computational efficiency. Because our data are relatively low-plex in this example, we use all features instead.

Here, we use the implementation of PCA provided in the *[scater](https://bioconductor.org/packages/3.22/scater)* package ([McCarthy et al. 2017](https://bioconductor.org/books/release/OSTA/pages/ind-dimensionality-reduction.html#ref-McCarthy2017-scater)):

For large-scale datasets (100,000s of cells), argument `BSPARAM=RandomParam()` can be used to decrease runtime by approximating the singular value decomposition (SVD). Because this implementation uses randomization, and seed for random number generation should be set in order to make results reproducible.

Code

```
spe <- logNormCounts(spe)
spe <- runPCA(spe, ncomponents=20)
```

### 27.2.2 Spatially aware

*[BANKSY](https://bioconductor.org/packages/3.22/BANKSY)* ([Singhal et al. 2024](https://bioconductor.org/books/release/OSTA/pages/ind-dimensionality-reduction.html#ref-Singhal2024-BANKSY)) computes PCs on a spatial-neighborhood augmented matrix, thereby embedding cells in a product space of their own and their local neighborhood’s (average) transcriptome, representing cell state and microenvironment. A key parameter for this method is λ∈[0,1] (argument `lambda` in `runBanksyPCA()`), which controls the spatial component’s weight; notably, **when λ=0, `BANKSY` reduces to non-spatial clustering**. Secondly, `k_geom` determines the number of neighbors to use for computing local transcriptomic neighborhoods.

In Visium, k=6 would correspond to first-order, k=18 to first- and second-order neighbors.

Code

```
spe <- computeBanksy(spe, assay_name="logcounts", k_geom=20)
spe <- runBanksyPCA(spe, npcs=20, lambda=0.4)
```

To not confuse different types of PCs, we rename the corresponding `reducedDims` to end in `_sp` and `_tx` for spatially aware and unaware results, respectively.

Code

```
reducedDimNames(spe) <- c("PCA_tx", "PCA_sp")
```

## 27.3 Uniform manifold approximation and projection (UMAP)

We can also perform non-linear dimensionality reduction using the UMAP algorithm, applied to the set of top PCs (default 50). The first two UMAP dimensions are often used for visualization purposes (although interpretation of the UMAP dimensions is difficult and can be misleading). In addition, the UMAP algorithm is stochastic, so the visual layout of the UMAP plot will change from one run to the next; setting a seed for random number generation ensures that the visualization is reproducible.

Note that UMAP plots should be interpreted with caution ([Chari and Pachter 2023](https://bioconductor.org/books/release/OSTA/pages/ind-dimensionality-reduction.html#ref-Chari2023-specious)). Specifically, the UMAP algorithm can introduce misleading structure in a 2D plot, and the relative distances between points in the 2D plot are not meaningful. For this reason, we recommend using UMAP only for generating visualizations, not clustering – any cluster labels shown on a UMAP plot should be generated using a separate clustering algorithm. Similar arguments also apply to the t-SNE algorithm, which may be used as an alternative to UMAP.

Here, we run UMAP on both non-spatial and spatially aware PCs, retaining the first two UMAP components to be used for visualization:

Code

```
spe <- runUMAP(spe, dimred="PCA_tx", name="UMAP_tx")
spe <- runUMAP(spe, dimred="PCA_sp", name="UMAP_sp")
```

## 27.4 Exploratory

A useful visualization is to color cells by their PCs in physical space. This will help highlight key drivers of transcriptional variability, e.g., between major biological compartments such as epithelia, immune and stromal cells. Exemplary plots of PCs 1-3 are rendered below:

Code

```
df <- data.frame(colData(spe), spatialCoords(spe), 
                 reducedDim(spe, "PCA_tx"))
lapply(paste0("PC", seq(3)), \(.) {
    fd <- df[order(abs(df[[.]])), ]
    ggplot(fd, aes(x_centroid, y_centroid, col=.data[[.]]))
}) |>
    wrap_plots(nrow=1) &
    geom_point(shape=16, stroke=0, size=0.2) &
    scale_color_gradientn(colors=pals::jet()) &
    coord_equal() & theme_void() & theme(
        legend.position="bottom",
        legend.key.width=unit(0.8, "lines"),
        legend.key.height=unit(0.4, "lines"))
```

![](ind-dimensionality-reduction_files/figure-html/plot-xy-pcs-1.png)

PCs as RGB

A neat trick in this context is to rescale PC 1-3 values between 0 and 1, and interpret them as RGB values for visualization:

Code

```
pcs <- reducedDim(spe, "PCA_tx")[, seq(3)]
pcs <- sweep(pcs, 2, colMins(pcs), `-`)
pcs <- sweep(pcs, 2, colMaxs(pcs), `/`)
rgb <- apply(pcs, 1, \(.) rgb(.[1], .[2], .[3]))
df <- data.frame(colData(spe), spatialCoords(spe), rgb)
```

E.g., we see that colors corresponding to opposing PC values are clearly separated in space, indicating distinct drivers of transcriptional variability. In some regions, we see color mixtures, indicating cells of ambiguous transcriptional profiles (with respect to what is being captured by the 3 PCs considered here).

Code

```
ggplot(df, aes(x_centroid, y_centroid, col=rgb)) +
    geom_point(shape=16, stroke=0, size=0.4) +
    scale_color_identity() + coord_equal() + theme_void()
```

![](ind-dimensionality-reduction_files/figure-html/plot-xy-rgb-1.png)

Code

```
.plot_umap <- \(dr) plotReducedDim(
    object=spe, dimred=dr, colour_by="Label", 
    point_shape=16, point_size=0) + ggtitle(dr)
.plot_umap("UMAP_tx") +
.plot_umap("UMAP_sp") +
    plot_layout(guides="collect") &
    guides(col=guide_legend(override.aes=list(alpha=1, size=2))) &
    theme_void() & theme(aspect.ratio=1, 
        legend.key.size=unit(0, "lines"),
        plot.title=element_text(hjust=0.5))
```

![](ind-dimensionality-reduction_files/figure-html/plt-umap-1.png)

## 27.5 Appendix

Back to top