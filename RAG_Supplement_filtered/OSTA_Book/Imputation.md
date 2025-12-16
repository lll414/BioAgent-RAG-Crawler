Source URL: https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html
Date Scraped: 2025-12-15

---

1. [Cross-platform analyses](https://bioconductor.org/books/release/OSTA/pages/crs-spatial-registration.html)
2. [37  Imputation](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html)

# 37  Imputation

## 37.1 Preamble

### 37.1.1 Introduction

As was previously discussed, most commonly employed imaging-based ST assays can resolve 100s to 1000s features. Meanwhile, non-spatially resolved single-cell assays are able to capture transcripts across 10,000 gene features. To allow estimating the expression of genes that may not be available in a ST assay, we would like to integrate spatially resolved data with scRNA-seq data using common features, in order to predict missing features using nearest neighbors (i.e., transcriptionally similar) cells.

### 37.1.2 Dependencies

Code

```
library(dplyr)
library(DropletUtils)
library(ggspavis)
library(harmony)
library(igraph)
library(OSTA.data)
library(patchwork)
library(RANN)
library(scran)
library(scater)
library(SpatialExperiment)
library(SpatialExperimentIO)
```

## 37.2 Data import

The analyses demonstrated here will use 313-plex Xenium data (10x Genomics) of a human breast cancer biopsy section and the Chromium Single Cell Gene Expression Flex applied to FFPE tissues (scFFPE-Seq) assay with around 18,000 features ([Janesick et al. 2023](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-Janesick2023-high-res)).

### 37.2.1 Xenium

Code

```
id <- "Xenium_HumanBreast1_Janesick"
pa <- OSTA.data_load(id, mol=FALSE)
dir.create(td <- tempfile())
unzip(pa, exdir=td)
spe <- readXeniumSXE(td, addTx=FALSE)
spe$sample_id <- "Xenium"
dim(spe)
```

```
##  [1]    313 167780
```

### 37.2.2 Chromium

We now load the single cell FFPE RNA-Seq assay. To streamline the demonstration, we consolidate some of the cell type annotations provided by 10x Genomics (i.e., `Annotation`) into the scRNA-seq assay. The single-cell assay has 18,082 features, some of which are available in the Xenium assay as well.

Code

```
id <- "Chromium_HumanBreast_Janesick"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)
h5 <- file.path(td, "filtered_feature_bc_matrix.h5")
sce <- read10xCounts(h5, col.names=TRUE)
sce$sample_id <- "Chromium"
# set gene symbols as feature names
rownames(sce) <- make.unique(rowData(sce)$Symbol)
# retrieve cell type labels
fs <- list.files(td, full.names=TRUE)
csv <- grep("csv$", fs, value=TRUE)
cd <- read.csv(csv, row.names=1)
colData(sce)[names(cd)] <- cd[colnames(sce), ]
dim(sce)
```

```
##  [1] 18082 30365
```

## 37.3 Analysis

### 37.3.1 Integration

principles of single-cell data integration

Considering an omics assay as a matrix of measurement values where rows/columns correspond to features/observations, there are three possible scenarios for integration ([Argelaguet et al. 2021](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-Argelaguet2021)). Specifically, methods applicable to each type of integration strategy depend on the choice of *anchors* used or, otherwise put, the dimension along which to integrate.

Firstly, experiments may involve multiple (biological or technical) replicates, requiring **horizontal** integration across cells that share a common feature space. Secondly, multi-modal assays can simultaneously extract di↵erent types of features from a given cell. Integration must then be **vertical**, i.e., across features measured from the same cells. Thirdly, technological limitations may require measuring different modalities from separate groups of cells, in which case **diagonal** integration across both di↵erent cells and features is necessary

Single-cell transcriptomics data analyses are most commonly faced with horizontal integration. Methods for this task have been benchmarked extensively ([Tran et al. 2020](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-Tran2020); [Chazarra-Gil et al. 2021](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-ChazarraGil2021); [Luecken et al. 2022](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-Luecken2022)). While these span diverse types of approaches, they often share conceptually similar steps. For example, model-based methods explicitly fit the transcriptional effects of batch assignment, while graph-based methods separately link points within and between batches. `harmony` ([Korsunsky et al. 2019](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-Korsunsky2019-Harmony)), which we demonstrate here, first represents batches in PCA space, and uses maximum diversity clustering followed by linear mixture modeling to compute a corrected embedding.

To accomplish data transfer across cells of either two assays (scRNA-seq and imaging-based ST), we first need to find cells that are similar between the two datasets. To this end, we use the common genes of both assays to integrate both datasets into a joint embedding space. Here, the objective is to find the nearest scRNA-seq cell neighbors of each Xenium cell in the joint embedding space, assuming that these cells will have similar transcriptional profiles.

We can integrate them on a transcriptional level; here, using *[harmony](https://bioconductor.org/packages/3.22/harmony)* ([Korsunsky et al. 2019](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-Korsunsky2019-Harmony)). To this end, we first consolidate the scRNA-seq and Xenium data into one object:

Code

```
# concatenate objects
gs <- intersect(rownames(spe), rownames(sce))
cd <- intersect(names(colData(spe)), names(colData(sce)))
obj <- lapply(list(spe, sce), \(obj) {
    obj <- obj[gs, ]
    y <- as(assay(obj), "dgCMatrix")
    SingleCellExperiment(
        assays=list(counts=y),
        colData=colData(obj)[cd])
}) |> do.call(what=cbind)
# add single-cell annotations
lab <- match(colnames(obj), colnames(sce))
obj$Annotation <- sce$Annotation[lab]
```

We now filter out cells with too few counts, and then perform PC-based `harmony` integration:

Code

```
# minimal filtering
obj <- obj[, colSums(counts(obj)) > 20]
# log-library size normalization 
obj <- logNormCounts(obj)
# principal component analysis
# (w/o additional feature selection)
obj <- runPCA(obj, name=".PCA")
# 'harmony' integration
# (keeping uncorrected PCs)
pcs <- RunHarmony(
    data_mat=reducedDim(obj, ".PCA"), 
    meta_data=obj$sample_id, 
    verbose=FALSE)
reducedDim(obj, "PCA") <- pcs
# dimensionality reduction before/after 
# integration (for visualization only)
obj <- runUMAP(obj, dimred=".PCA", name=".UMAP")
obj <- runUMAP(obj, dimred="PCA", name="UMAP")
```

Let us now visualize the embedding of both Xenium and scRNA-seq cells to assess the mixing of observations from both datasets. We also check if scRNA-seq cells are well-partitioned given their annotations:

Code

```
.obj <- obj[, obj$sample_id == "Chromium"]
plotUMAP(obj, colour_by="sample_id", point_size=0, dimred=".UMAP") + ggtitle("uncorrected") +
plotUMAP(obj, colour_by="sample_id", point_size=0) + ggtitle("corrected") +
plotUMAP(.obj, colour_by="Annotation", point_size=0) + ggtitle("Chromium") +
plot_layout(nrow=1, guides="collect") &
    guides(col=guide_legend(override.aes=list(alpha=1, size=2))) &
    theme_void() & theme(
        aspect.ratio=1, 
        legend.key.size=unit(0, "pt"),
        plot.title=element_text(hjust=0.5))
```

![](crs-imputation_files/figure-html/plt-map-1.png)

quantifying batch effects and their removal

UMAPs like above are perhaps pretty, but not quantiative. We strongly encourage more quantitative metrics to quantify batch effects prior and post correction, e.g., entropy, principal component regression, local inverse simpson index, silhouette width (problematic according to Rautenstrauch and Ohler ([2025](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-Rautenstrauch2025))), or similat. A number of such metrics have been proposed in benchmarks (see previous box), and many are implemeted in *[CellMixS](https://bioconductor.org/packages/3.22/CellMixS)* ([Lütge et al. 2021](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-Lutge2021-CellMixS)).

Here, we will perform principal component regression (PCR) to quantify the variance explained by cluster/batch (before and after correction); after corresion, we want most variation to be attributable to clusters.

Code

```
# helper to perform PCR 
# xs = 'colData' to use as predictor(s)
# dr = 'reducedDim' slot to use as response
.pcr <- \(obj, xs, dr) {
    pcs <- reducedDim(obj, dr)
    lapply(xs, \(x) {
        fit <- summary(lm(pcs ~ obj[[x]]))
        r2 <- sapply(fit, \(.) .$adj.r.squared)
        data.frame(x, pc=seq_along(r2), r2)
    }) |> do.call(what=rbind)  
}

# analysis
xs <- c("sample_id", "Annotation")
pcs <- list(before=".PCA", after="PCA")
pcr <- lapply(pcs, \(dr) .pcr(obj, xs, dr))

# wrangling
df <- bind_rows(pcr, .id="id")
df$id <- factor(df$id, names(pcs))
rownames(df) <- NULL
head(df)
```

```
##        id         x pc          r2
##  1 before sample_id  1 0.002068132
##  2 before sample_id  2 0.069938223
##  3 before sample_id  3 0.005267385
##  4 before sample_id  4 0.009725794
##  5 before sample_id  5 0.150020768
##  6 before sample_id  6 0.355010050
```

Visualizing coefficients of determination (R-squared), we can observe that batch effects were not that strong to begin with, and are virtually absent after correction; meanwhile, clusters explain similar amounts of variation before/after. This is interesting as the UMAPs visualized above would have suggested otherwise.

Code

```
# plot R2 before vs. after
ggplot(df, aes(pc, r2, col=x)) + 
    facet_wrap(~id) +
    geom_line(show.legend=FALSE) + geom_point() +
    scale_x_continuous(breaks=c(1, seq(5, 20, 5))) +
    scale_y_continuous(limits=c(NA, 1), breaks=seq(0, 1, 0.2)) +
    labs(x="principal component", y="coeff. of determination") +
    guides(col=guide_legend("predictor", override.aes=list(size=2))) +
    coord_cartesian(xlim=c(1, 20)) +
    theme_bw() + theme(
        panel.grid.minor=element_blank(),
        legend.key.size=unit(0, "lines"))
```

![](crs-imputation_files/figure-html/plt-pcr-1.png)

### 37.3.2 Selection

Here, we choose a small subset of features to transfer onto the Xenium assay. For this purpose, we test for differentially expressed genes (DEGs) between cell types (`Annotation` column), and select the top-20 markers for each:

Code

```
# normalize
sce <- logNormCounts(sce)
# test for differentially expressed genes 
# (DEGs) to identify subpopulation markers
colLabels(sce) <- sce$Annotation
mgs <- findMarkers(sce, test="wilcox", direction="up")
# get top markers per cell type
top <- lapply(names(mgs), \(k) {
    df <- mgs[[k]][, c("FDR", "summary.AUC")]
    g <- head(rownames(df), 20)
    data.frame(k, g, df[g, ])
}) |> do.call(what=rbind)
rownames(top) <- NULL
head(top)
```

```
##          k     g           FDR summary.AUC
##  1 B Cells  CD74  0.000000e+00   0.9477608
##  2 B Cells MS4A1  0.000000e+00   0.8773396
##  3 B Cells CXCR4 9.594950e-304   0.8271647
##  4 B Cells SP140 4.172406e-229   0.7843702
##  5 B Cells BANK1 1.321886e-228   0.7840307
##  6 B Cells CIITA  0.000000e+00   0.8615623
```

### 37.3.3 Aggregation

Now that we have shown there is a correction for batch effects between Xenium and Chromium assays, we can transfer counts from single-cell onto Xenium data by identifying the k-nearest neighbors (kNNs) of each Xenium cell in the Chromium data, and imputing missing genes by aggregating their profile across kNNs (here, using mean expression).

Code

```
# find k-nearest neighbors across embeddings
# (from Xenium to Chromium cells)
idx <- split(colnames(obj), obj$sample_id)
pcs <- reducedDim(obj, "PCA")
ref <- pcs[idx$Chromium, ]
que <- pcs[idx$Xenium, ]
knn <- nn2(data=ref, query=que, k=k <- 20)
# create adjacency matrix
el <- cbind(
    rep(idx$Xenium, each=k), 
    idx$Chromium[c(t(knn$nn.idx))])
g <- graph_from_edgelist(el, directed=TRUE)
A <- as_adjacency_matrix(g)
# average Chromium across Xenium neighbors
gs <- unique(top$g)
cs <- intersect(rownames(A), idx$Chromium)
ws <- A[idx$Xenium, cs]*(1/k)
y <- logcounts(sce)[gs, cs] %*% t(ws)
```

## 37.4 Downstream

### 37.4.1 Visualization

We now create two new `SpatialExperiment` objects by i) subsetting cells that passed our minimal filtering, and ii) replacing observed by imputed counts in one object:

Code

```
spe <- spe[, idx$Xenium]
spe <- logNormCounts(spe)
sqe <- SpatialExperiment(
    assays=list(logcounts=y), 
    spatialCoords=spatialCoords(spe))
# needed for 'ggspavis' visualization
spe$in_tissue <- sqe$in_tissue <- 1
```

We can now compare imputed to observed counts. Here, we visualize ACTA2, which is included in both assays, as well as KRT17, which also marks myoepithelial cells:

Code

```
.p <- \(obj, col) {
    plotCoords(obj, 
        point_size=0,
        annotate=col, 
        assay_name="logcounts")
}
pal <- scale_color_gradient2("obs.", high="red")
qal <- scale_color_gradient2("imp.", high="blue")
.p(spe, "ACTA2") + pal +
.p(sqe, "ACTA2") + qal +
.p(sqe, "KRT17") + qal +
    plot_layout(nrow=1)
```

![](crs-imputation_files/figure-html/plt-xy-one-1.png)

### 37.4.2 Comparison

As a more systematic comparison, we can visualize the correlation of cell-wise expression values (using shared features only) between observed and imputed data:

Code

```
gs <- intersect(rownames(spe), rownames(sqe))
imp <- as.matrix(t(logcounts(sqe[gs, ])))
obs <- as.matrix(t(logcounts(spe[gs, ])))
cm <- cor(imp, obs, method="pearson")
df <- data.frame(g=gs, p=diag(cm))
df$g <- factor(df$g, df$g[order(-df$p)])
ggplot(df, aes(g, p, fill=p > 0.5)) + geom_col() +
    labs(x="correlation (observed vs. imputed)", y=NULL) +
    geom_hline(yintercept=0.5) +
    theme_minimal() + theme(
        legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
```

![](crs-imputation_files/figure-html/plt-cor-1.png)

We observe the highest correlation between observed and imputed data for KRT7; let’s confirm this by visualizing the gene’s expression in space:

Code

```
.p(spe, "KRT7") + pal +
.p(sqe, "KRT7") + qal +
    plot_layout(nrow=1)
```

![](crs-imputation_files/figure-html/plt-xy-two-1.png)

Genes that exhibit the lowest correlation tend to be rarely expressed. This arguably makes it difficult to identify similar Chromium cells and, in turn, accurately transfer information between both assays. By contrast, highly correlated genes are detected in a decent proportion of cells:

Code

```
# quantify proportion of cells for which
# lowly/highly correlated genes are detected
fq <- \(gs) {
    y <- counts(spe[gs, ])
    round(rowMeans(y > 0), 2)
}
gs <- names(sort(diag(cm)))
fq(tail(gs)) # high corr.
```

```
##     FASN TACSTD2   ERBB2   CXCR4    IL7R    KRT7 
##     0.53    0.57    0.89    0.49    0.26    0.68
```

Code

```
fq(head(gs)) # low corr.
```

```
##   HDC CTSG CPA3 GZMB  KIT CD83 
##  0.30 0.01 0.03 0.07 0.06 0.08
```

## 37.5 Other methods

Alternatively, we can use other R/Python frameworks to integrate imaging-based ST assays with scRNA-seq datasets, and impute missing features. Some of these methods could be used through R using packages such as *[reticulate](https://bioconductor.org/packages/3.22/reticulate)* and *[basilisk](https://bioconductor.org/packages/3.22/basilisk)*; see [Chapter 7](https://bioconductor.org/books/release/OSTA/pages/bkg-python-interoperability.html).

* **Tangram**: is available as a Python module through [scvi-tools](https://docs.scvi-tools.org/en/latest/user_guide/models/tangram.html) ([Biancalani et al. 2021](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-Biancalani2021-Tangram)). The tool performs a softmax optimization to generate a probabilistic mapping between cells of scRNA-seq data and voxels (cells or spots) of ST data, which is then used to transfer labels or feature profiles across datasets.
* **LIGER**: is available in R via *[rliger](https://cran.r-project.org/package=rliger)* ([Welch et al. 2019](https://bioconductor.org/books/release/OSTA/pages/crs-imputation.html#ref-Welch2019-LIGER)), and incorporates a similar workflow performed in this chapter using *[harmony](https://bioconductor.org/packages/3.22/harmony)*. The only difference between these approaches is how the joint embedding is obtained; here, LIGER uses non-negative matrix factorization (NMF). LIGER has been originally proposed to perform this task for scRNA-seq and scATAC-seq assays.

## 37.6 Appendix

Back to top