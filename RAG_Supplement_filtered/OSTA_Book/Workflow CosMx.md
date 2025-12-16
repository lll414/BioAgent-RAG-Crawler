Source URL: https://bioconductor.org/books/release/OSTA/pages/img-workflow-cosmx.html
Date Scraped: 2025-12-15

---

1. [Imaging-based platforms](https://bioconductor.org/books/release/OSTA/pages/img-introduction.html)
2. [25  Workflow: CosMx](https://bioconductor.org/books/release/OSTA/pages/img-workflow-cosmx.html)

# 25  Workflow: CosMx

## 25.1 Preamble

### 25.1.1 Introduction

In this demo, we will analyze a 1k-plex CosMx data (Bruker) of a mouse brain sample; see [Chapter 6](https://bioconductor.org/books/release/OSTA/pages/bkg-example-datasets.html). Following quality control and preprocessing, we will perform cell type annotation, identification of marker genes, and cell-cell interaction analysis.

### 25.1.2 Dependencies

Code

```
library(DESpace)
library(ggplot2)
library(OSTA.data)
library(patchwork)
library(pheatmap)
library(scater)
library(scrapper)
library(scRNAseq)
library(sf)
library(SingleR)
library(SpaceTrooper)
library(SpatialExperiment)
library(SpatialExperimentIO)
library(spatialFDA)
library(SpatialFeatureExperiment)
library(spdep)
library(Voyager)
# set seed for random number generation
# in order to make results reproducible
set.seed(112358)
```

Code

```
# retrieve CosMx dataset from OSF repo
id <- "CosMx1k_MouseBrain2"
pa <- OSTA.data_load(id, mol=FALSE)
dir.create(td <- tempfile())
unzip(pa, exdir=td)
cos <- readCosmxSXE(td, addTx=FALSE)
# prepare data for 'SpaceTrooper'
cos <- updateCosmxSPE(cos, td, sampleName="CosMx")
cos <- readAndAddPolygonsToSPE(cos)
cos$in_tissue <- TRUE
cos
```

```
##  class: SpatialExperiment 
##  dim: 950 48556 
##  metadata(4): fov_positions fov_dim polygons technology
##  assays(1): counts
##  rownames(950): Chrna4 Slc6a1 ... Cck Aqp4
##  rowData names(0):
##  colnames(48556): f1_c1 f1_c10 ... f99_c98 f99_c99
##  colData names(22): fov cellID ... polygons in_tissue
##  reducedDimNames(0):
##  mainExpName: NULL
##  altExpNames(1): NegPrb
##  spatialCoords names(2) : CenterX_global_px CenterY_global_px
##  imgData names(1): sample_id
```

## 25.2 Quality control

We can have an overview of the data, by plotting the spatial distribution of the cells, represented by their centroids. This is useful to have a global bird’s-eye-view of the tissue.

Code

```
plotCentroids(cos)
```

![](img-workflow-cosmx_files/figure-html/plot-xy-1.png)

Unlike, Xenium and MERSCOPE, CosMx does not stitch fields of view (FOVs) in a single image, but rather treats each FOV independently. It is therefore useful to visualize a map of the FOVs to understand their spatial arrangement.

Code

```
plotCellsFovs(cos)
```

![](img-workflow-cosmx_files/figure-html/plot-fov-map-1.png)

[Chapter 19](https://bioconductor.org/books/release/OSTA/pages/img-quality-control.html) gives a thorough overview of quality control steps for imaging-based spatial transcriptomics data, including CosMX-specific metrics. Here, we will simply compute *[SpaceTrooper](https://bioconductor.org/packages/3.22/SpaceTrooper)*’s aggregated QC score and keep the cells that pass the threshold.

Code

```
lys <- c("NegPrb", "Negative", "SystemControl")
cos <- spatialPerCellQC(cos, rmZeros=TRUE, negProbList=lys)
cos <- computeQCScore(cos)
plotCentroids(cos, colourBy="QC_score")
```

![](img-workflow-cosmx_files/figure-html/qc-cosmx-1.png)

Code

```
cos <- computeQCScoreFlags(cos, qsThreshold=0.5)
plotCentroids(cos, colourBy="low_qcscore") + 
    theme(legend.key.size=rel(0)) +
    guides(col=guide_legend(override.aes=list(size=2)))
```

![](img-workflow-cosmx_files/figure-html/filter-qc-cosmx-1.png)

Code

```
cos <- cos[, !cos$low_qcscore]
```

## 25.3 Processing

For the sake of runtime, we will perform downstream analyses only on a small portion of the tissue, defined by the following FOVs:

Code

```
fs <- c(72:74, 90:92, 97:99, 114:116)
plotZoomFovsMap(cos, fovs=fs)
```

![](img-workflow-cosmx_files/figure-html/subset-fovs-cosmx-1.png)

Code

```
sub <- cos[, cos$fov %in% fs]
```

Next, we’ll log-normalize counts by area, and perform principal component analysis (PCA) on all 950 RNA targets:

See [Chapter 26](https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html) for more details on normalization.

Code

```
# cell area-based normalization
logcounts(sub) <- sweep(assay(sub), 2, sub$Area, `/`)
# principal component analysis
sub <- runPCA(sub)
```

Let’s visualize the expression of the first two principal components (PCs) in the spatial context.

Code

```
# add PCs as cell metadata
pcs <- reducedDim(sub, "PCA")
colData(sub) <- cbind(colData(sub), pcs)
# visualize PCs 1 & 2 in space
plotCentroids(sub, colourBy="PC1") +
plotCentroids(sub, colourBy="PC2")
```

![](img-workflow-cosmx_files/figure-html/plt-gs-xy-1.png)

## 25.4 Annotation

Here, we use the Zeisel et al. ([2015](https://bioconductor.org/books/release/OSTA/pages/img-workflow-cosmx.html#ref-zeisel2015cell)) dataset as a reference to annotate the cell types with *[SingleR](https://bioconductor.org/packages/3.22/SingleR)*.

Code

```
sceZ <- ZeiselBrainData()
sceZ <- logNormCounts(sceZ)
```

Code

```
pred <- SingleR(test=sub, ref=sceZ, labels=sceZ$level1class, de.method="wilcox")
sub$SingleR_label <- factor(pred$pruned.labels)

nk <- nlevels(sub$SingleR_label)
pal <- hcl.colors(nk, "Spectral")

plotPolygons(sub, 
    colourBy="SingleR_label") +
    scale_fill_manual(values=pal)
```

![](img-workflow-cosmx_files/figure-html/annot-cosmx-1.png)

Note that the annotation performed here is a transfer label from a reference dataset that only partially matches with the tissue under study and purposely annotated at a coarse level (level 1). Ideally, one can use a more appropriate reference dataset, or better yet, build a custom reference from single-cell RNA-seq data of the same tissue.

A different, often preferable, approach consists in using, in addition to the gene expression profiles, the image data and the spatial context of the cells. *[InSituType](https://github.com/Nanostring-Biostats/InSituType)* ([Danaher et al. 2022](https://bioconductor.org/books/release/OSTA/pages/img-workflow-cosmx.html#ref-danaher2022insitutype)) provides a framework to perform such analysis in unsupervised, semi-supervised, or supervised manners. Please refer to its documentation for more details.

### 25.4.1 Marker genes

Here, we identify cluster-related spatially-variable genes using the *[DESpace](https://bioconductor.org/packages/3.22/DESpace)* package ([Cai, Robinson, and Tiberi 2024](https://bioconductor.org/books/release/OSTA/pages/img-workflow-cosmx.html#ref-Cai2024-DESpace)). In this workflow, we use this as a way to identify marker genes for the cell types annotated above.

Code

```
res <- svg_test(spe=sub, cluster_col="SingleR_label")
head(res$gene_results)
```

```
##          gene_id       LR   logCPM PValue FDR
##  Itm2a     Itm2a 2596.761 11.03859      0   0
##  Flt1       Flt1 2601.543 11.05386      0   0
##  Cldn5     Cldn5 2825.269 11.11885      0   0
##  Gja1       Gja1 4982.193 11.29054      0   0
##  Gpr37l1 Gpr37l1 4989.928 11.30344      0   0
##  Bcan       Bcan 3567.078 11.61163      0   0
```

We focus on the genes with the smallest p-values, and compute the average expression of each gene in each cell type. We can visualize this with a heatmap.

Code

```
df <- res$gene_results
gs <- df$gene_id[rank(df$PValue, ties.method="min") == 1]
mu <- t(apply(counts(sub)[gs, ], 1, tapply, sub$SingleR_label, mean))

pheatmap(t(mu), 
    scale="column", show_colnames=FALSE,
    cluster_rows=TRUE, cluster_cols=TRUE, 
    main="Average expression of top SVGs per cell type")
```

![](img-workflow-cosmx_files/figure-html/plt-svg-1.png)

We can also visualize the spatial distribution of some of these genes, e.g., Mbp that marks oligodendrocytes, and Calm1, Calm2, and Snap25 that mark different subsets of neurons. Note also how some genes show a difference of expression within the pyramidal neuron population, suggesting the presence of subtypes.

Code

```
gs <- c("Mbp", "Calm1", "Calm2", "Snap25")
ps <- lapply(gs, \(g) {
    sub[[g]] <- logcounts(sub)[g, ]
    plotPolygons(sub, colourBy=g) + ggtitle(g)
}) 
dy <- range(logcounts(sub)[gs, ])
wrap_plots(ps, guides="collect") &
    scale_fill_viridis_c(
        "expression", 
        limits=dy, breaks=dy,
        labels=c("low", "high")) &
    theme(plot.title=element_text(hjust=0.5))
```

![](img-workflow-cosmx_files/figure-html/plt-svg-genes-1.png)

## 25.5 Neighborhood-based gene expression analyses

We start by creating a cell neighborhood graph, based on a distance threshold of 50 microns. To do so, we use the *[SpatialFeatureExperiment](https://bioconductor.org/packages/3.22/SpatialFeatureExperiment)* package ([Moses et al. 2023](https://bioconductor.org/books/release/OSTA/pages/img-workflow-cosmx.html#ref-Moses2023-Voyager)).

Note that an alternative approach uses the k nearest neighbors (kNN) method to define the neighborhood graph. Here, we prefer to use a distance-based approach, as it better captures the morphological structure of this tissue area.

Code

```
# convert to SFE & remove cells with NA labels
sfe <- as(sub, "SpatialFeatureExperiment")
sfe <- sfe[, !is.na(sfe$SingleR_label)]

# this is needed for Voyager's plotting functions
colnames(sfe$polygons)[6] <- "geometry"
st_geometry(sfe$polygons) <- "geometry"
colGeometry(sfe, "cellSeg") <- sfe$polygons

colGraph(sfe, "poly2nb") <- 
    findSpatialNeighbors(sfe,
        type="cellSeg", 
        method="poly2nb", 
        style="W", snap=50)

plotColGraph(sfe, 
    colGraphName="poly2nb", 
    colGeometryName="cellSeg") + 
    theme_void()
```

![](img-workflow-cosmx_files/figure-html/cci-graph-cosmx-1.png)

Given the neighborhood graph, we can calculate local measures of spatial autocorrelation for specific genes in the dataset.

Here, we focus on Snap25 and Calm1, neuronal marker genes identified above that show interesting spatial expression distributions (see above). First, we compute the local Moran’s I coefficient.

Code

```
gs <- c("Snap25", "Calm1")
sfe <- runUnivariate(sfe,
    features=gs,
    type="localmoran",
    zero.policy=TRUE,
    colGraphName="poly2nb",
    colGeometryName="cellSeg")

plotLocalResult(sfe, 
    name="localmoran",
    features=gs, ncol=2,
    colGeometryName="cellSeg",
    divergent=TRUE, diverge_center=0)
```

![](img-workflow-cosmx_files/figure-html/cci-moran-snap25-1.png)

From the plots, we can see that Snap25 shows positive spatial autocorrelation in a very localized part of the tissue, while Calm1 shows a more diffuse autocorrelation structure.

We focus on Snap25 and visualize its autocorrelation structure using a Moran scatter plot, which shows the relationship between the expression of the gene in each cell and the average expression in its neighbors.

Code

```
sfe <- runUnivariate(sfe,
    feature=gs[1],
    zero.policy=TRUE,
    type="moran.plot",
    colGraphName="poly2nb",
    colGeometryName="cellSeg")

moranPlot(sfe,
    feature=gs[1],
    graphName="poly2nb",
    swap_rownames="symbol")
```

![](img-workflow-cosmx_files/figure-html/cci-moran-scatter-snap25-1.png)

We can then classify cells based on the significance of their local Moran’s I values, and visualize the resulting clusters in the spatial context.

Code

```
res <- localResults(sfe)[["localmoran"]][[1]]
res$locClust <- with(res, ifelse(
    `-log10p_adj` > -log10(0.05), 
    as.character(mean), "non-siginificant"))
localResults(sfe)[["localmoran"]][[1]] <- res

plotLocalResult(sfe,
    features=gs[1],
    name="localmoran",
    attribute="locClust",
    colGeometryName="cellSeg")
```

![](img-workflow-cosmx_files/figure-html/cci-moran-scatter-sig-1.png)

In this case, we can see a significant positive autocorrelation in the lower left part of the structure, while the rest of the tissue does not show significant autocorrelation.

Note that other local spatial autocorrelation measures, such as Geary’s C coefficient and Getis-Ord statistic, as well as their global versions can be used for similar analyses Moses et al. ([2023](https://bioconductor.org/books/release/OSTA/pages/img-workflow-cosmx.html#ref-Moses2023-Voyager)).

## 25.6 Cell-cell interactions

### 25.6.1 Joint counts

We can also use the neighborhood graph to investigate cell-cell interactions between different cell types annotated above. Here, we use the join count statistic implemented in the *[spdep](https://cran.r-project.org/package=spdep)* package to quantify the tendency of cells of different types to be neighbors.

Code

```
jc <- joincount.multi(sfe$SingleR_label, colGraph(sfe, "poly2nb"), zero.policy=TRUE)
head(jc[order(abs(jc[, "z-value"]), decreasing=TRUE), ])
```

```
##                                      Joincount   Expected  Variance   z-value
##  Jtot                               1132.15707 2494.52543 265.87369 -83.55200
##  pyramidal CA1:pyramidal CA1         735.66944  351.21205  62.92504  48.46592
##  pyramidal CA1:oligodendrocytes       22.39901  303.28546  75.62290 -32.30012
##  oligodendrocytes:oligodendrocytes   196.71667   65.35907  19.34973  29.86192
##  pyramidal CA1:astrocytes_ependymal   97.55301  379.29404  92.67327 -29.26664
##  pyramidal SS:pyramidal CA1          195.38757  445.19310 106.83889 -24.16781
```

We can see that there is a strong tendency for cells of the same type to be neighbors (high positive z-values). There are also some interesting negative interactions between different cell types: e.g., pyramidal CA1 tend to be away from oligodendrocytes and astrocytes. This highlights the compartimentalized structure of the brain tissue.

### 25.6.2 Point-pattern analysis

A similar result can be obtained using point-pattern analysis, such as Ripley’s K function, and similar methods, using for instance the *[spatialFDA](https://bioconductor.org/packages/3.22/spatialFDA)* package.

Here, we use Besag’s L function to quantify the interaction between pyramidal CA1 neurons, and the interaction between them and oligodendrocytes (cross-L function).

Other methods are available, such as Ripley’s K function, and the pair correlation function ([Emons et al. 2025](https://bioconductor.org/books/release/OSTA/pages/img-workflow-cosmx.html#ref-Emons2025-pasta)). Note that here we are making the quite strong assuption of homogeneity of the point pattern, which is irrealistic in tissue biology. In particular, CSR is a naïve null model that does not take into account tissue structure and cell density variations across the tissue. Hence, these plots should be considered as exploratory analyses rather than formal hypothesis tests.

Code

```
res <- calcMetricPerFov(
    spe=sub, fun="Lest",
    subsetby="sample_id", 
    marks="SingleR_label",
    selection="pyramidal CA1",
    rSeq=seq(0, 500, l=100),
    by="sample_id", ncores=1)

plotMetricPerFov(res, 
    theo=TRUE, correction="iso", 
    x="r", imageId="sample_id")
```

![](img-workflow-cosmx_files/figure-html/cci-CA1-l-1.png)

In this type of plots, the dashed line represents the theoretical value of the metric under complete spatial randomness (CSR). Values above this line indicate clustering, while values below indicate dispersion. Here, we can see that pyramidal CA1 neurons tend to cluster together, as expected.

Note that the L-function lays below the CSR line for the first 50 microns, which is a consequence of the cell segmentation that does not allow cells to overlap, indicating the impossibility of finding two cells at a distance smaller than their size.

While in the plot above we focused on a single cell type, we can also investigate the interaction between two different cell types, e.g., pyramidal CA1 neurons and oligodendrocytes, using the cross-L function.

Code

```
res <- calcMetricPerFov(
    spe=sub, fun="Lcross", 
    subsetby="sample_id", 
    marks="SingleR_label",
    selection=c("pyramidal CA1", "oligodendrocytes"),
    rSeq=seq(0, 500, l=100),
    by="sample_id", ncores=1)

plotMetricPerFov(res, 
    theo=TRUE, correction="iso", 
    x="r", imageId="sample_id")
```

![](img-workflow-cosmx_files/figure-html/cci-cross-l-1.png)

In this case, we can see that pyramidal CA1 neurons and oligodendrocytes tend to be spatially segregated, as indicated by the cross-L function being below the CSR line.

We only briefly discussed cell-cell interactions and spatial exploratory analyses here: we refer interested readers to the comprehensive documentation and tutorials available at the [pasta](https://robinsonlabuzh.github.io/pasta/00-home.html) and [Voyager](https://pachterlab.github.io/voyager/) websites.

## 25.7 Appendix

Back to top