Source URL: https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html
Date Scraped: 2025-12-15

---

1. [Sequencing-based platforms](https://bioconductor.org/books/release/OSTA/pages/seq-introduction.html)
2. [12  Deconvolution](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html)

# 12  Deconvolution

## 12.1 Introduction

Sequencing-based ST data can contain zero to multiple cells per spot, which might be fully or only partially covered by cells, depending on the spatial resolution of the platform and the tissue cell density (see also [Chapter 8](https://bioconductor.org/books/release/OSTA/pages/seq-introduction.html) and the schematic Figure below). This aspect of the data implies that there may be a mixture of cell types in a spot and thus a mixture of transcriptional programs.

![](../images/seq-deconvolution.png)

Schematic illustrating multiple cells in a 10x Genomics Visium spot (diameter of 55µm; purple line) overlaying the hematoxylin and eosin (H&E) stained image. The center-to-center distance between spots is 100µm (yellow line), and a typical immune cell’s diameter is around 10µm (cyan line).

To help understand these mixtures, at least 20 deconvolution techniques have been proposed for spot-level ST data. Some methods require borrowing insights from a scRNA-seq reference dataset, while others can be reference-free. Based on their underlying algorithms, Li et al. ([2023](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Li2023-benchmark-deconvolution)) have grouped methods into five categories:

* **probabilistic-based**: use Bayesian inference, likelihood estimation, or probabilistic modeling to estimate cell type compositions while incorporating uncertainty. Available tools include
  + *[CARDspa](https://bioconductor.org/packages/3.22/CARDspa)* ([Ma and Zhou 2022](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Ma2022-CARD))
  + *[spacexr](https://bioconductor.org/packages/3.22/spacexr)*’s `RCTD` ([Dylan M. Cable et al. 2022b](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Cable2022-RCDT))
  + *[SpatialDecon](https://bioconductor.org/packages/3.22/SpatialDecon)* ([Danaher et al. 2022](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Danaher2022-SpatialDecon))
  + *[STdeconvolve](https://github.com/JEFworks-Lab/STdeconvolve)* ([Miller et al. 2022](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Miller2022-STdeconvolve)) in R, and
  + *[cell2location](https://github.com/BayraktarLab/cell2location)* ([Kleshchevnikov et al. 2022](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Kleshchevnikov2022-cell2location))
  + *[scvi-tools](https://github.com/scverse/scvi-tools)*’s `DestVI` ([Lopez et al. 2022](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Lopez2022-DestVI))
  + *[std-poisson](https://github.com/SpatialTranscriptomicsResearch/std-poisson)* ([Berglund et al. 2018](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Berglund2018-std-poisson))
  + *[stereoscope](https://github.com/almaan/stereoscope)* ([Andersson et al. 2020](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Andersson2020-stereoscope))
  + *[STRIDE](https://github.com/DongqingSun96/STRIDE)* ([Sun et al. 2022](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Sun2022-STRIDE)) in Python.
* **non-negative matrix factorization (NMF)-based**: decompose gene expression data into latent components representing different cell types. Available tools include *[Giotto](https://github.com/drieslab/Giotto)*’s `SpatialDWLS` ([Chen et al. 2025](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Chen2025-Giotto-Suite)) and *[SPOTlight](https://bioconductor.org/packages/3.22/SPOTlight)* ([Elosua-Bayes et al. 2021](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-ElosuaBayes2021-SPOTlight)) in R, and *[NMFreg\_tutorial](https://github.com/tudaga/NMFreg_tutorial)* in Python.
* **graph-based**: use graph neural networks or graph-based optimization to model spatial relationships. Available tools include *[SD2](https://github.com/leihouyeung/SD2)* ([Li et al. 2022](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Li2022-SD2)) in R, and *[DSTG](https://github.com/Su-informatics-lab/DSTG)* ([Song and Su 2021](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Song2021-DSTG)) and *[SpiceMIx](https://github.com/ma-compbio/SpiceMIx)* ([Chidester et al. 2023](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Chidester2023-SpiceMix)) in Python.
* **optimal transport (OT)-based**: infer spatial gene expression distributions by mapping scRNA-seq and ST data. Available tools include *[SpaOTsc](https://github.com/zcang/SpaOTsc)* ([Cang and Nie 2020](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Cang2020-SpaOTsc)) and *[novosparc](https://github.com/rajewsky-lab/novosparc)* ([Nitzan et al. 2019](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Nitzan2019-novoSpaRc)) in Python.
* **deep learning-based**: align and integrate single-cell and spatial transcriptomics data with neural networks. For example, *[Tangram](https://github.com/broadinstitute/Tangram)* ([Biancalani et al. 2021](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Biancalani2021-Tangram)) in Python.

Among these, `std-poisson`, `STdeconvolve`, and `SpiceMix` are reference-free methods. Methods that incorporate spatial location information are `CARD`, `DSTG`, `SD2`, `Tangram`, `cell2location`, `DestVI`, `std-poisson`, and `SpiceMix`.

In this section, we will demonstrate deconvolution of cell types per spot, using `RCTD` on Visium and Visium HD datasets.

## 12.2 Dependencies

Code

```
library(BiocParallel)
library(CARDspa)
library(DropletUtils)
library(ggplot2)
library(ggspavis)
library(OSTA.data)
library(patchwork)
library(pheatmap)
library(scran)
library(scater)
library(spacexr)
library(SpatialExperiment)
library(VisiumIO)
```

In this example of Visium breast cancer data ([Janesick et al. 2023](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Janesick2023-high-res)), we perform cell type deconvolution without a single-cell (Chromium) reference and compare the concordance with the provided Visium annotation provided by 10x Genomics.

Code

```
# retrieve dataset from OSF repository
id <- "Visium_HumanBreast_Janesick"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)

# read into 'SpatialExperiment'
vis <- TENxVisium(
    spacerangerOut=file.path(td, "outs"), 
    processing="filtered", 
    format="h5", 
    images="lowres") |> 
    import()

# retrieve spot annotations & add as metadata
df <- read.csv(file.path(td, "annotation.csv"))
cs <- match(colnames(vis), df$Barcode)
vis$anno <- factor(df$Annotation[cs])

# set gene symbols as feature names
rownames(vis) <- make.unique(rowData(vis)$Symbol)
vis
```

```
##  class: SpatialExperiment 
##  dim: 18085 4992 
##  metadata(2): resources spatialList
##  assays(1): counts
##  rownames(18085): SAMD11 NOC2L ... MT-ND6 MT-CYB
##  rowData names(3): ID Symbol Type
##  colnames(4992): AACACCTACTATCGAA-1 AACACGTGCATCGCAC-1 ...
##    TGTTGGCCAGACCTAC-1 TGTTGGCCTACACGTG-1
##  colData names(5): in_tissue array_row array_col sample_id anno
##  reducedDimNames(0):
##  mainExpName: Gene Expression
##  altExpNames(0):
##  spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
##  imgData names(4): sample_id image_id data scaleFactor
```

Code

```
xy <- spatialCoords(vis) * scaleFactors(vis)
ys <- nrow(imgRaster(vis)) - range(xy[, 2])
xs <- range(xy[, 1])
box <- geom_rect(
    xmin=xs[1], xmax=xs[2], ymin=ys[1], ymax=ys[2], 
    col="black", fill=NA, linetype=2, linewidth=2/3)
plotVisium(vis, spots=FALSE, point_size=1) + box + 
    plotVisium(vis, point_size=1, zoom=TRUE) + 
    plot_layout(nrow=1) & facet_null()
```

![](seq-deconvolution_files/figure-html/plt-hne-1.png)

Deconvolution is performed after quality control, as detailed in [Chapter 10](https://bioconductor.org/books/release/OSTA/pages/seq-quality-control.html), and is usually performed on unnormalized and untransformed (i.e. raw) counts. Here, we quickly check some typically spot-level metrics.

Code

```
sub <- list(mt=grep("^MT-", rownames(vis)))
vis <- addPerCellQCMetrics(vis, subsets=sub)
```

Code

```
vis$log_sum <- log1p(vis$sum)
plotCoords(vis, 
    annotate="log_sum") + 
    ggtitle("log library size") + 
plotCoords(vis, 
    annotate="subsets_mt_percent") + 
    ggtitle("% mitochondrial") + 
ggplot(
    data.frame(colData(vis)), 
    aes(x=sum, y=subsets_mt_percent)) + 
    geom_point() + geom_density_2d() +
    scale_x_log10() + scale_y_sqrt() +
    theme(aspect.ratio=2/3) +
plot_layout(nrow=1) & theme(
    legend.key.width=unit(0.5, "lines"), 
    legend.key.height=unit(1, "lines")) & 
    scale_color_gradientn(colors=pals::jet())
```

![](seq-deconvolution_files/figure-html/plt-xy-qc-1.png)

A few spots have low library sizes, and can be removed.

Code

```
vis <- vis[, vis$sum > 1000]
```

We first visualize the spot-level cell type annotation provided by 10x Genomics.

Code

```
plotCoords(vis, 
    annotate="anno", point_size=1, 
    pal=unname(pals::trubetskoy())) + 
    theme(legend.key.size=unit(0, "lines"))
```

![](seq-deconvolution_files/figure-html/visanno-1.png)

Now, we load the single-cell (Chromium) reference data for the Visium dataset. To streamline the demonstration, we consolidate some of the cell type annotations provided by 10x Genomics (i.e. `Annotation`) into more generalized categories (i.e. `Annogrp`).

Code

```
# retrieve dataset from OSF repo
id <- "Chromium_HumanBreast_Janesick"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)

# read into 'SingleCellExperiment'
fs <- list.files(td, full.names=TRUE)
h5 <- grep("h5$", fs, value=TRUE)
sce <- read10xCounts(h5, col.names=TRUE)

# use gene symbols as feature names
rownames(sce) <- make.unique(rowData(sce)$Symbol)

# retrieve cell type labels
csv <- grep("csv$", fs, value=TRUE)
cd <- read.csv(csv, row.names=1)

# ignore mixtures
lab <- cd$Annotation
lab[grepl("Hyb", lab)] <- NA 

# simplify annotations
pat <- c(
    "B Cell"="B", "T Cell"="T", "Mac"="macro", "Mast"="mast", 
    "DCs"="dendritic", "Peri"="perivas", "End"="endo", 
    "Str"="stromal", "Inv"="tumor", "Myo"="myoepi")
for (. in names(pat)) 
    lab[grep(., lab)] <- pat[.]
lab <- gsub("\\s", "", lab)

# add as cell metadata
table(cd$Annogrp <- lab)
```

```
##  
##          B     DCIS1     DCIS2         T dendritic      endo     macro 
##       1463      1863      2159      4742       313      1055      3724 
##       mast    myoepi   perivas   stromal     tumor 
##         92      1839       285      2611      5897
```

Code

```
colData(sce)[names(cd)] <- cd[colnames(sce), ]
```

We only keep the Chromium data with an annotation and are not labeled as “Hybrid”, as these correspond to mixed subpopulations.

Code

```
sce <- sce[, !is.na(sce$Annogrp)]
dim(sce)
```

```
##  [1] 18082 26031
```

## 12.3 RCTD

Next, we perform deconvolution with *[spacexr](https://bioconductor.org/packages/3.22/spacexr)* (also known as RCTD)([Dylan M. Cable et al. 2022a](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Cable2022-RCTD)). By default, `runRctd()`’s `rctd_mode = "doublet"` specifies at most two subpopulations coexist in a data unit (i.e. within a spot); here, we set `rctd_mode = "full"` in order to allow for an arbitrary number of subpopulations to be fit instead.

Note that RCTD can also be adapted to Visium HD data with `rctd_mode = "doublet"`, as demonstrated by ([de Oliveira et al. 2025](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-deOliveira2025-high-def)) and [Chapter 16](https://bioconductor.org/books/release/OSTA/pages/seq-workflow-visium-hd-seg.html).

Code

```
rctd_data <- createRctd(vis, sce, cell_type_col="Annogrp")
(res <- runRctd(rctd_data, max_cores=4, rctd_mode="full"))
```

```
##  class: SpatialExperiment 
##  dim: 12 4962 
##  metadata(4): spatial_rna config cell_type_info internal_vars
##  assays(1): weights
##  rownames(12): B DCIS1 ... stromal tumor
##  rowData names(0):
##  colnames(4962): AACACCTACTATCGAA-1 AACACGTGCATCGCAC-1 ...
##    TGTTGGCCAGACCTAC-1 TGTTGGCCTACACGTG-1
##  colData names(1): sample_id
##  reducedDimNames(0):
##  mainExpName: NULL
##  altExpNames(0):
##  spatialCoords names(2) : x y
##  imgData names(0):
```

Weights inferred by `RCTD` should be normalized such that proportions of cell types sum to 1 for each spot:

Code

```
# scale weights such that they sum to 1
ws <- assay(res)
ws <- sweep(ws, 2, colSums(ws), `/`)

ws_rctd <- data.frame(t(as.matrix(ws)))
round(ws_rctd[1:5, 1:5], 2)
```

```
##                        B DCIS1 DCIS2    T dendritic
##  AACACCTACTATCGAA-1 0.00  0.00  0.00 0.00      0.01
##  AACACGTGCATCGCAC-1 0.04  0.01  0.00 0.04      0.00
##  AACACTTGGCAAGGAA-1 0.00  0.00  0.03 0.01      0.01
##  AACAGGAAGAGCATAG-1 0.03  0.00  0.00 0.08      0.03
##  AACAGGATTCATAGTT-1 0.00  0.00  0.01 0.00      0.00
```

Code

```
# add proportion estimates to colData
colData(vis)[names(ws_rctd)] <- ws_rctd[colnames(vis), ]
```

## 12.4 CARD

Another method that can be used is `CARD`. First, we rename the columns of spatial coordinates for `CARD`.

Code

```
# realize delayed matrices, as CARD does 
# not yet support delayed matrix handling
counts(sce) <- as(counts(sce), "sparseMatrix")
counts(vis) <- as(counts(vis), "sparseMatrix")
colnames(spatialCoords(vis)) <- c("x", "y")
```

Next, we perform the `CARD` deconvolution. Here, we demonstrate `CARD`’s interoperability with `SingleCellExperiment` and `SpatialExperiment`. The deconvolution result matrix is already normalized such that the sum of cell type proportions for each spot is equal to 1.

Note: `CARD` can also take a reference matrix, a reference cell type annotation column, a spatial count matrix, and a spatial coordinates data frame as separate items in `sc_count`, `sc_meta`, `spatial_count`, and `spatial_location`, respectively. However, we encourage simplifying the process by using existing Bioconductor classes.

Code

```
set.seed(2025)
CARD_obj <- CARD_deconvolution(
    spe=vis,
    sce=sce,
    sc_count=NULL,
    sc_meta=NULL,
    spatial_count=NULL,
    spatial_location=NULL,
    ct_varname="Annogrp",
    ct_select=NULL,      # use all 'sce$Annogrp' cell types
    sample_varname=NULL, # use all 'sce' as one 'ref' sample 
    mincountgene=100,
    mincountspot=5)
ws_card <- CARD_obj$Proportion_CARD

# order cell type names alphabetically, as for RCTD
ws_card <- data.frame(ws_card[, colnames(ws_rctd)])
round(ws_card[1:5, 1:5], 2)
```

```
##                     B DCIS1 DCIS2    T dendritic
##  AACACCTACTATCGAA-1 0     0     0 0.00         0
##  AACACGTGCATCGCAC-1 0     0     0 0.00         0
##  AACACTTGGCAAGGAA-1 0     0     0 0.01         0
##  AACAGGAAGAGCATAG-1 0     0     0 0.02         0
##  AACAGGATTCATAGTT-1 0     0     0 0.00         0
```

## 12.5 Visualization

First, we define a couple accessory functions.

Code

```
.plt_xy <- \(ws, vis, col, point_size) {
    xy <- spatialCoords(vis)[rownames(ws), ]
    colnames(xy) <- c("x", "y")
    df <- cbind(ws, xy)
    ggplot(df, aes(x, y, col=.data[[col]])) + 
        coord_equal() + theme_void() + 
        geom_point(size=point_size)
}

.plt_decon <- \(ws, vis) {
    ps <- lapply(names(ws), \(.) .plt_xy(ws, vis, col=., point_size=0.3))
    ps |> wrap_plots(nrow=3) & theme(
        legend.key.width=unit(0.5, "lines"),
        legend.key.height=unit(1, "lines")) &
        scale_color_gradientn(colors=pals::jet())
}
```

We can visualize deconvolution weights in x-y space, i.e., coloring by the proportion of a given cell type estimated to fall within a given spot:

* RCDT
* CARD

Code

```
.plt_decon(ws=ws_rctd, vis)
```

![](seq-deconvolution_files/figure-html/plt-dec-xy-rctd-1.png)

Code

```
.plt_decon(ws=ws_card, vis)
```

![](seq-deconvolution_files/figure-html/plt-dec-xy-card-1.png)

The deconvolution results can also be viewed as a heatmap, where rows = cells and columns = clusters:

Code

```
plot_heat_ws <- \(ws, string){
    p <- pheatmap(ws, 
        show_rownames=FALSE, show_colnames=TRUE, main=string,
        cellwidth=12, treeheight_row=5, treeheight_col=5)
    return(p)
}
plot_heat_ws(ws_rctd, string="RCTD") 
plot_heat_ws(ws_card, string="CARD")
```

![](seq-deconvolution_files/figure-html/heatmap-1.png)

![](seq-deconvolution_files/figure-html/heatmap-2.png)

In both methods, we see that more than half of the spots are estimated to have a stromal proportion of more than 50%. Few spots have an intense and distinct signal for cancerous subpopulations, DCIS1 and DCIS2. For the following analysis, we focus on `RCTD` as an example.

For comparison with spot annotations provided by 10x Genomics, we include majority voted cell type from deconvolution by `RCTD`. Note that, because stromal cells show broad signals across the entire tissue, to better investigate immune cell signals, we remove stromal from the majority vote calculation for an alternative label: `RCTD_no_stroma`.

Code

```
ws <- ws_rctd
# derive majority vote label
ids <- names(ws)[apply(ws, 1, which.max)]
names(ids) <- rownames(ws)
vis$RCTD <- factor(ids[colnames(vis)])

# derive majority vote excluding stromal cells
ws_no_stroma <- ws[, colnames(ws) != "stromal"]
ids_no_stroma <- names(ws_no_stroma)[apply(ws_no_stroma, 1, which.max)]
names(ids_no_stroma) <- rownames(ws)
vis$RCTD_no_stroma <- factor(ids_no_stroma[colnames(vis)])
```

We can visualize these three annotations spatially:

Code

```
lapply(
    c("anno", "RCTD", "RCTD_no_stroma"), 
    \(.) plotCoords(vis, annotate=.)) |>
    wrap_plots(nrow=1) &
    theme(legend.key.size=unit(0, "lines")) &
    scale_color_manual(values=unname(pals::trubetskoy()))
```

![](seq-deconvolution_files/figure-html/plt-clu-xy-1.png)

Note the strong stromal signals and macrophages being the second most common cell type for stromal cells. To help characterize subpopulations from deconvolution, we can view their distribution against the provided annotation:

Code

```
cd <- data.frame(colData(vis))
df <- as.data.frame(with(cd, table(RCTD, anno)))
fd <- as.data.frame(with(cd, table(RCTD_no_stroma, anno)))
ggplot(df, 
    aes(Freq, RCTD, fill=anno)) + 
    ggtitle("RCTD") +
ggplot(fd, 
    aes(Freq, RCTD_no_stroma, fill=anno)) + 
    ggtitle("RCTD_no_stroma") +
plot_layout(nrow=1, guides="collect") &
    labs(x="Proportion", y=NULL) &
    coord_cartesian(expand=FALSE) &
    geom_col(width=1, col="white", position="fill") &
    scale_fill_manual(values=unname(pals::trubetskoy())) &
    theme_minimal() & theme(
        aspect.ratio=1,
        legend.key.size=unit(2/3, "lines"),
        plot.title=element_text(hjust=0.5))
```

![](seq-deconvolution_files/figure-html/plt-clu-fq-1.png)

Next, we can investigate the agreement between the provided annotation against the two deconvolution majority vote labels.

Code

```
hm <- \(mat, string) pheatmap(
    mat, show_rownames=TRUE, show_colnames=TRUE, main=string,
    cellwidth=10, cellheight=10, treeheight_row=5, treeheight_col=5)
hm(prop.table(table(vis$anno, vis$RCTD), 2), string="RCTD")
hm(prop.table(table(vis$anno, vis$RCTD_no_stroma), 2), string="RCTD_no_stroma")
```

![](seq-deconvolution_files/figure-html/dec-clu-nostroma-1.png)

![](seq-deconvolution_files/figure-html/dec-clu-nostroma-2.png)

Overall, we observe agreement between the provided spot labels and the `RCTD` deconvolution derived annotations. Before cleaning up stromal, some immune cell types, such as dendritic and mast, never had a chance to have the highest cell type proportion. On the left panel, among among all the spots annotated by `RCTD` as T cells, nearly all of them are from the “immune” type in the provided annotation. Strong agreements are also observed for spots with cell type of “DCIS1”, “DCIS2”, and “Invasive tumor”.

PC regression

Next, we prepare the principal components (PCs) needed to perform PC regression:

Code

```
# log-library size normalization
vis <- logNormCounts(vis)
# feature selection 
dec <- modelGeneVar(vis)
hvg <- getTopHVGs(dec, prop=0.1)
# dimension reduction 
set.seed(1234)
vis <- runPCA(vis, ncomponents=20, subset_row=hvg)
```

We fit the deconvolution result of each cell type against the first 10 PCs to obtain 10 regressions.

Code

```
idx <- rownames(ws)
ids <- colnames(ws)
pcs <- reducedDim(vis, "PCA")
pcs <- pcs[idx, seq_len(10)]
pcr <- lapply(ids, \(id) {
    fit <- summary(lm(pcs ~ ws[[id]]))
    r2 <- sapply(fit, \(.) .$adj.r.squared)
    data.frame(id, pc=seq_along(r2), r2)
}) |> do.call(what=rbind)
```

Here we plot the coefficient of determination of the first 10 PCs for each cell type.

Code

```
pcr$id <- factor(pcr$id, ids)
pal <- pals::trubetskoy()
ggplot(pcr, aes(pc, r2, col=id)) +
    geom_line(show.legend=FALSE) + geom_point() +
    scale_color_manual("predictor", values=unname(pal)) +
    scale_x_continuous(breaks=c(1, seq(5, 20, 5))) +
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
    labs(x="principal component", y="coeff. of determination") +
    guides(col=guide_legend(override.aes=list(size=2))) +
    coord_cartesian(xlim=c(1, 10)) +
    theme_minimal() + theme(
      panel.grid.minor=element_blank(),
      legend.key.size=unit(0, "lines"))
```

![](seq-deconvolution_files/figure-html/pcr-plot-1.png)

Let’s inspect the key drivers of (expression) variability in terms of PCs. Considering deconvolution results from above, we can see that, e.g.:

* PC1 distinguishes stromal, tumor, macrophage from the rest of the tissue
* PC3, PC4 and PC5 separate DCIS1, T and endothelial cells, respectively

Code

```
# retrieve top-10 PCs
pcs <- reducedDim(vis, "PCA")
pcs <- pcs[rownames(ws), seq_len(10)]
# specify subpopulations & PCs to visualize
var <- c("DCIS1", "T", "endo")
var <- c(var, colnames(pcs)[3:5])
# visualize deconvolution weights alongside PCs
lapply(var, \(.) {
    .plt_xy(
        cbind(ws, pcs), vis, col=., point_size=0.3) +
        scale_color_gradientn(., colors=pals::jet())
}) |>
    wrap_plots(nrow=2) & theme(
        plot.title=element_blank(),
        legend.key.width=unit(0.5, "lines"),
        legend.key.height=unit(1, "lines"))
```

![](seq-deconvolution_files/figure-html/plt-pcs-xy-1.png)

Note that the direction of each PC is irrelevant from how much variation it explains.

In conclusion, deconvolution-based cell type proportion estimates are able to recapitulate PCs and, in turn, expression variability.

Apart from being a tool for spot deconvolution, `RCTD` can be used as a label transfer tool to annotate imaging-based ST data, such as for Xenium and MERSCOPE. For this, the default `doublet_mode = "doublet"` should be used, and a certainty score would be returned to indicate doublets with two predicted cell types.

## 12.6 Appendix

### Benchmarks

Benchmarking studies of deconvolution methods often require generating synthetic spots to establish a ground truth for cell type proportions. This process involves either simulating artificial tissue patterns or aggregating counts from scRNA-seq or imaging-based spatial transcriptomics data into spots. However, it is equally important to evaluate how these methods perform in tissues with highly spatially-heterogeneous cell type compositions, such as real cancer samples. Below are three comprehensive benchmarking studies, two of which that incorporate both artificial and real datasets.

* Sang-aram et al. ([2023](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Sang-aram2023-Spotless)) developed a pipeline to benchmark 11 deconvolution methods, including `RCTD`, across 63 synthetic, 3 binned, and 2 real datasets. `RCTD` and `cell2location` were the most recommended methods. Figure 2 gives an overview of the benchmarking results.
* Li et al. ([2023](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Li2023-benchmark-deconvolution)) benchmarked 18 deconvolution methods, including `RCTD`, across 50 simulated and real datasets. Among these methods, `CARD`, `cell2location`, and `Tangram` are highly recommended. Figure 1 summarizes the method performance, and Figure 4 gives a flowchart of how to decide on which method to use.
* Gaspard-Boulinc et al. ([2025](https://bioconductor.org/books/release/OSTA/pages/seq-deconvolution.html#ref-Gaspard-Boulinc2025-deconvolution)) review and compare available cell-type deconvolution methods, and provide a continuously updated web-based summary table.

Back to top