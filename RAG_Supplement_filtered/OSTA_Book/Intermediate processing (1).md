Source URL: https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html
Date Scraped: 2025-12-15

---

1. [Sequencing-based platforms](https://bioconductor.org/books/release/OSTA/pages/seq-introduction.html)
2. [11  Intermediate processing](https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html)

# 11  Intermediate processing

## 11.1 Preamble

### 11.1.1 Introduction

This chapter demonstrates methods for several intermediate processing steps – normalization, feature selection, and dimensionality reduction – that are required before downstream analysis methods can be applied.

We use methods from the *[scater](https://bioconductor.org/packages/3.22/scater)* ([McCarthy et al. 2017](https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html#ref-McCarthy2017-scater)) and *[scran](https://bioconductor.org/packages/3.22/scran)* ([Lun, McCarthy, and Marioni 2016](https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html#ref-Lun2016-scran)) Bioconductor packages, which were originally developed for scRNA-seq data, with the (simplified) assumption that these methods can be applied to spot-based ST data by treating spots as equivalent to cells. In addition, we discuss alternative spatially-aware methods, and provide links to later parts of the book where these methods are described in more detail.

Following the processing steps, this chapter also includes short demonstrations for subsequent steps, including clustering and identification of marker genes. For more details on these steps, see the later analysis chapters or workflow chapters.

### 11.1.2 Dependencies

Code

```
library(scran)
library(scater)
library(ggspavis)
library(pheatmap)
library(patchwork)
library(BayesSpace)
library(SpatialExperiment)
```

Code

```
# set seed for reproducibility
set.seed(100)
```

### 11.1.3 Load data

Here, we load 10x Genomics Visium data from one postmortem human brain tissue section from the dorsolateral prefrontal cortex (DLPFC) brain region ([Maynard et al. 2021](https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html#ref-Maynard2021-DLPFC)), which has previously undergone quality control and filtering (see also [Chapter 10](https://bioconductor.org/books/release/OSTA/pages/seq-quality-control.html)).

Code

```
# load data
(spe <- readRDS("seq-spe_qc.rds"))
```

```
##  class: SpatialExperiment 
##  dim: 21842 3639 
##  metadata(0):
##  assays(1): counts
##  rownames(21842): ENSG00000243485 ENSG00000238009 ... ENSG00000271254
##    ENSG00000268674
##  rowData names(4): gene_id gene_name feature_type subsets_mito
##  colnames(3639): AAACAAGTATCTCCCA-1 AAACAATCTACTAGCA-1 ...
##    TTGTTTGTATTACACG-1 TTGTTTGTGTAAATTC-1
##  colData names(14): barcode_id sample_id ... subsets_mito_percent
##    total
##  reducedDimNames(0):
##  mainExpName: NULL
##  altExpNames(0):
##  spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
##  imgData names(4): sample_id image_id data scaleFactor
```

## 11.2 Normalization

### 11.2.1 Library size normalization

A simple and fast approach for normalization is “library size normalization”, which consists of calculating log-transformed normalized counts (“logcounts”) using library size factors, treating each spot as equivalent to a single cell. This method is simple, fast, and generally provides a good baseline. It can be calculated using the *[scater](https://bioconductor.org/packages/3.22/scater)* ([McCarthy et al. 2017](https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html#ref-McCarthy2017-scater)) and *[scran](https://bioconductor.org/packages/3.22/scran)* ([Lun, McCarthy, and Marioni 2016](https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html#ref-Lun2016-scran)) packages.

See also [Chapter 26](https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html).

However, library size normalization does not make use of any spatial information. In some datasets, the simplified assumption that each spot can be treated as equivalent to a single cell is not appropriate, and can create problems during analysis. (*“Library size confounds biology in spatial transcriptomics data”* ([Bhuva et al. 2024](https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html#ref-Bhuva2024-library-size)).)

Some alternative non-spatial methods from scRNA-seq workflows (e.g. normalization by deconvolution) are also less appropriate for spot-based ST data, since spots can contain multiple cells from one or more cell types, and datasets can contain multiple samples (e.g. multiple tissue sections), which may result in sample-specific clustering.

Code

```
# calculate library size factors
spe <- computeLibraryFactors(spe)

summary(sizeFactors(spe))
```

```
##     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.04576 0.63035 0.89770 1.00000 1.29000 3.79953
```

Code

```
hist(sizeFactors(spe), breaks = 50, main = "Histogram of size factors")
```

![](seq-intermediate-processing_files/figure-html/logcounts-1.png)

Code

```
# calculate logcounts and store in new assay
spe <- logNormCounts(spe)

assayNames(spe)
```

```
##  [1] "counts"    "logcounts"
```

### 11.2.2 Spatially-aware normalization

*[SpaNorm](https://bioconductor.org/packages/3.22/SpaNorm)* ([Salim et al. 2025](https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html#ref-Salim2025-SpaNorm)) was developed for ST data, using a gene-wise model (e.g. negative binomial) in which variation is decomposed into a library size-dependent (technical) and -independent (biological) component.

## 11.3 Feature selection

### 11.3.1 Highly variable genes (HVGs)

Identifying a set of top “highly variable genes” (HVGs) is a standard step for feature selection in many scRNA-seq workflows, which can also be used as a simple and fast baseline approach in spot-based ST data. This makes the simplified assumption that spots can be treated as equivalent to cells.

Here, we take a standard approach for selecting HVGs; see [Section 29.2.1](https://bioconductor.org/books/release/OSTA/pages/ind-feature-selection-testing.html#sec-ind-feature-selection-testing-hvgs) for more details. We first remove mitochondrial genes, since these tend to be very highly expressed and are not of main biological interest.

Code

```
# identify mitochondrial genes
nm <- rowData(spe)$gene_name
mt <- grepl("^MT-", nm, ignore.case = TRUE)
table(mt)
# remove them
spe <- spe[!mt, ]
```

```
##  mt
##  FALSE  TRUE 
##  21829    13
```

Next, apply methods from *[scran](https://bioconductor.org/packages/3.22/scran)*. This gives a list of top HVGs, which can be used as the input for subsequent steps. The parameter `prop` defines the proportion of HVGs to select – for example, `prop = 0.1` returns the top 10%. (Alternatively, argument `n` can be used to selected a fixed number of top HVGs, say 1000 or 2000.)

Code

```
# fit mean-variance relationship
dec <- modelGeneVar(spe)
# plot
plot(dec$mean, dec$total, xlab = "mean (logexpr)", ylab = "variance (logexpr)")
curve(metadata(dec)$trend(x), add = TRUE, col = "dodgerblue")
```

![](seq-intermediate-processing_files/figure-html/hvgs-1.png)

Code

```
# select top HVGs
sel <- getTopHVGs(dec, prop = 0.1)
# number of HVGs selected
length(sel)
```

```
##  [1] 1420
```

### 11.3.2 Spatially variable genes (SVGs)

Alternatively, spatially-aware methods can be used to identify a set of “spatially variable genes” (SVGs). These methods take the spatial coordinates of the measurements into account, which can generate a more biologically informative ranking of genes associated with biological structure within the tissue area. The set of SVGs may then be used either instead of or complementary to HVGs in subsequent steps.

A number of methods have been developed to identify SVGs. For more details and examples, see [Section 29.3.1](https://bioconductor.org/books/release/OSTA/pages/ind-feature-selection-testing.html#sec-ind-feature-selection-testing-svgs).

## 11.4 Dimensionality reduction

A standard next step is to perform dimensionality reduction using principal component analysis (PCA), applied to the set of top HVGs (or SVGs). The set of top principal components (PCs) can then be used as the input for subsequent steps. See [Chapter 27](https://bioconductor.org/books/release/OSTA/pages/ind-dimensionality-reduction.html) for more details.

Note that we have set a random seed at the start of the chapter for reproducibility, since the fast implementation in `scater::runPCA()` uses a random initialization.

Code

```
spe <- runPCA(spe, subset_row = sel)
```

In addition, we can perform nonlinear dimensionality reduction using the UMAP algorithm, applied to the set of top PCs (default 50). The first two UMAP dimensions can be used for visualization purposes by plotting them on the x and y axes.

Code

```
spe <- runUMAP(spe, dimred = "PCA")
```

## 11.5 Clustering

After completing the intermediate processing steps above, we next demonstrate some short examples of downstream analysis steps.

Here, we run a spatially-aware clustering algorithm, *[BayesSpace](https://bioconductor.org/packages/3.22/BayesSpace)* ([Zhao et al. 2021](https://bioconductor.org/books/release/OSTA/pages/seq-intermediate-processing.html#ref-Zhao2021-BayesSpace)), to identify spatial domains. BayesSpace was developed for sequencing-based ST data, and takes the spatial coordinates of the measurements into account.

For more details on clustering, see [Chapter 28](https://bioconductor.org/books/release/OSTA/pages/ind-clustering.html).

Code

```
# run BayesSpace clustering
.spe <- spatialPreprocess(spe, skip.PCA = TRUE)
.spe <- spatialCluster(.spe, nrep = 1000, burn.in = 100, q = 10, d = 20)

# cluster labels
table(spe$BayesSpace <- factor(.spe$spatial.cluster))
```

```
##  
##    1   2   3   4   5   6   7   8   9  10 
##   35 672 250 527 378 492 316 274 494 201
```

## 11.6 Visualizations

We can visualize the cluster labels in x-y space as spatial domains, alongside the manually annotated reference labels (`ground_truth`) available for this dataset.

Code

```
# using plotting functions from ggspavis package
# and formatting using patchwork package
plotCoords(spe, annotate = "ground_truth") + 
plotCoords(spe, annotate = "BayesSpace") + 
  plot_layout() & 
  theme(legend.key.size = unit(0, "lines")) & 
  scale_color_manual(values = unname(pals::trubetskoy()))
```

![](seq-intermediate-processing_files/figure-html/plt-xy-clu-1.png)

Code

```
table(spe$BayesSpace, spe$ground_truth, useNA = "ifany")
```

```
##      
##       Layer1 Layer2 Layer3 Layer4 Layer5 Layer6  WM <NA>
##    1       3      0      0      0      1     30   0    1
##    2       0      0      4     63    562     43   0    0
##    3       0      0      0      0      0     38 212    0
##    4       0      4    500      0      1      2   0   20
##    5       3    209    165      0      0      1   0    0
##    6       0      0    318    155     19      0   0    0
##    7     265     40      2      0      1      1   0    7
##    8       0      0      0      0      0      1 273    0
##    9       0      0      0      0     88    406   0    0
##    10      2      0      0      0      1    170  28    0
```

Inspecting the spatial distribution of the top PCs, we can observe that the main driver of expression variability is the distinction between white matter (WM) and non-WM cortical layers.

Code

```
pcs <- reducedDim(spe, "PCA")
pcs <- pcs[, seq_len(4)]
lapply(colnames(pcs), \(.) {
    spe[[.]] <- pcs[, .]
    plotCoords(spe, annotate = .)
}) |> 
    wrap_plots(nrow = 1) & coord_equal() & 
    geom_point(shape = 16, stroke = 0, size = 0.2) & 
    scale_color_gradientn(colors = pals::jet(), n.breaks = 3) & 
    theme_void() & theme(
        plot.title = element_text(hjust = 0.5), 
        legend.key.width = unit(0.2, "lines"), 
        legend.key.height = unit(0.8, "lines"))
##  Coordinate system already present.
##  ℹ Adding new coordinate system, which will replace the existing one.
##  Coordinate system already present.
##  ℹ Adding new coordinate system, which will replace the existing one.
##  Coordinate system already present.
##  ℹ Adding new coordinate system, which will replace the existing one.
##  Coordinate system already present.
##  ℹ Adding new coordinate system, which will replace the existing one.
##  Scale for colour is already present.
##  Adding another scale for colour, which will replace the existing scale.
##  Scale for colour is already present.
##  Adding another scale for colour, which will replace the existing scale.
##  Scale for colour is already present.
##  Adding another scale for colour, which will replace the existing scale.
##  Scale for colour is already present.
##  Adding another scale for colour, which will replace the existing scale.
```

![](seq-intermediate-processing_files/figure-html/plot-xy-pcs-1.png)

## 11.7 Differential expression

Having clustered our spots to identify spatial domains, we can test for genes that are differentially expressed (DE) between spatial domains, which can be interpreted as marker genes for the spatial domains; see also [Section 29.2.2](https://bioconductor.org/books/release/OSTA/pages/ind-feature-selection-testing.html#sec-ind-feature-selection-testing-degs).

Here, we will use pairwise t-tests and specifically test for upregulation (as opposed to downregulation), i.e. expression should be higher in the cluster for which a gene is reported to be a marker. See [Chapter 28](https://bioconductor.org/books/release/OSTA/pages/ind-clustering.html) for additional details.

Alternatively, we could also test for spatially variable genes (SVGs) within domains to identify genes with spatial expression patterns that vary independently of the spatial arrangement of domains; see [Section 29.3.1](https://bioconductor.org/books/release/OSTA/pages/ind-feature-selection-testing.html#sec-ind-feature-selection-testing-svgs) for details.

Code

```
# using scran package
mgs <- findMarkers(spe, groups = spe$BayesSpace, direction = "up")
top <- lapply(mgs, \(df) rownames(df)[df$Top <= 2])
length(top <- unique(unlist(top)))
```

```
##  [1] 40
```

We can visualize selected markers as a heatmap (of cluster-wise expression means):

Code

```
# compute cluster-wise averages
pbs <- aggregateAcrossCells(spe, 
    ids = spe$BayesSpace, subset.row = top, 
    use.assay.type = "logcounts", statistics = "mean")
# use gene symbols as feature names
mtx <- t(assay(pbs))
colnames(mtx) <- rowData(pbs)$gene_name
# using pheatmap package
pheatmap(mat = mtx, scale = "column")
```

![](seq-intermediate-processing_files/figure-html/plt-mgs-hm-1.png)

Or spatial plots (of spot-level expression values):

Code

```
# gene-wise spatial plots
gs <- c("MBP", "PLP1", "NRGN", "SNAP25", "NEFL", "HPCAL1")
ps <- lapply(gs, \(.) {
    plotCoords(spe, 
        annotate = ., 
        feature_names = "gene_name", 
        assay_name = "logcounts") })
# figure arrangement
wrap_plots(ps, nrow = 2) & 
  theme(legend.key.width = unit(0.4, "lines"), 
        legend.key.height = unit(0.8, "lines")) & 
  scale_color_gradientn(colors = rev(hcl.colors(9, "Rocket")))
```

![](seq-intermediate-processing_files/figure-html/plt-mgs-xy-1.png)

## 11.8 Appendix

### Save data

Save data object for re-use within later chapters.

Code

```
colLabels(spe) <- spe$BayesSpace
saveRDS(spe, "seq-spe_cl.rds")
```

Back to top