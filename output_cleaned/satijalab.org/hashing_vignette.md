---
source_url: https://satijalab.org/seurat/articles/hashing_vignette
download_date: 2025-11-20
---

# Demultiplexing with hashtag oligos (HTOs)

#### Compiled: 2023-11-14

Source: [`vignettes/hashing_vignette.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/hashing_vignette.Rmd)

`hashing_vignette.Rmd`

Developed in collaboration with the Technology Innovation Group at NYGC, Cell Hashing uses oligo-tagged antibodies against ubiquitously expressed surface proteins to place a “sample barcode” on each single cell, enabling different samples to be multiplexed together and run in a single experiment. For more information, please refer to this [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1).

This vignette will give a brief demonstration on how to work with data produced with Cell Hashing in Seurat. Applied to two datasets, we can successfully demultiplex cells to their the original sample-of-origin, and identify cross-sample doublets.

The demultiplexing function `HTODemux()` implements the following procedure:

* We perform a k-medoid clustering on the normalized HTO values, which initially separates cells into K(# of samples)+1 clusters.
* We calculate a ‘negative’ distribution for HTO. For each HTO, we use the cluster with the lowest average value as the negative group.
* For each HTO, we fit a negative binomial distribution to the negative cluster. We use the 0.99 quantile of this distribution as a threshold.
* Based on these thresholds, each cell is classified as positive or negative for each HTO.
* Cells that are positive for more than one HTOs are annotated as doublets.

## 8-HTO dataset from human PBMCs

Dataset description:

* Data represent peripheral blood mononuclear cells (PBMCs) from eight different donors.
* Cells from each donor are uniquely labeled, using CD45 as a hashing antibody.
* Samples were subsequently pooled, and run on a single lane of the the 10X Chromium v2 system.
* You can download the count matrices for RNA and HTO [here](https://www.dropbox.com/sh/ntc33ium7cg1za1/AAD_8XIDmu4F7lJ-5sp-rGFYa?dl=0), or the FASTQ files from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108313)

### Basic setup

Load packages

```
library(Seurat)
```

Read in data

```
# Load in the UMI matrix
pbmc.umis <- readRDS("/brahms/shared/vignette-data/pbmc_umi_mtx.rds")

# For generating a hashtag count matrix from FASTQ files, please refer to
# https://github.com/Hoohm/CITE-seq-Count.  Load in the HTO count matrix
pbmc.htos <- readRDS("/brahms/shared/vignette-data/pbmc_hto_mtx.rds")

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)
```

```
## [1] "HTO_A" "HTO_B" "HTO_C" "HTO_D" "HTO_E" "HTO_F" "HTO_G" "HTO_H"
```

Setup Seurat object and add in the HTO data

```
# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(pbmc.umis), sparse = T))

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))
```

### Adding HTO data as an independent assay

You can read more about working with multi-modal data [here](/seurat/articles/multimodal_vignette)

```
# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
```

### Demultiplex cells based on HTO enrichment

Here we use the Seurat function `HTODemux()` to assign single cells back to their sample origins.

```
# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using
# the default settings
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)
```

### Visualize demultiplexing results

Output from running `HTODemux()` is saved in the object metadata. We can visualize how many cells are classified as singlets, doublets and negative/ambiguous cells.

```
# Global classification results
table(pbmc.hashtag$HTO_classification.global)
```

```
## 
##  Doublet Negative  Singlet 
##     2598      346    13972
```

Visualize enrichment for selected HTOs with ridge plots

```
# Group cells based on the max HTO signal
Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:2], ncol = 2)
```

![](hashing_vignette_files/figure-html/hashtag_ridge-1.png)

Visualize pairs of HTO signals to confirm mutual exclusivity in singlets

```
FeatureScatter(pbmc.hashtag, feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")
```

![](hashing_vignette_files/figure-html/hashtag_scatter1-1.png)

Compare number of UMIs for singlets, doublets and negative cells

```
Idents(pbmc.hashtag) <- "HTO_classification.global"
VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
```

![](hashing_vignette_files/figure-html/hashtag_vln-1.png)

Generate a two dimensional tSNE embedding for HTOs. Here we are grouping cells by singlets and doublets for simplicity.

```
# First, we will remove negative cells from the object
pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(pbmc.hashtag.subset) <- "HTO"
pbmc.hashtag.subset <- ScaleData(pbmc.hashtag.subset, features = rownames(pbmc.hashtag.subset),
    verbose = FALSE)
pbmc.hashtag.subset <- RunPCA(pbmc.hashtag.subset, features = rownames(pbmc.hashtag.subset), approx = FALSE)
pbmc.hashtag.subset <- RunTSNE(pbmc.hashtag.subset, dims = 1:8, perplexity = 100)
DimPlot(pbmc.hashtag.subset)
```

![](hashing_vignette_files/figure-html/hashtag_sub_tsne-1.png)

```
# You can also visualize the more detailed classification result by running Idents(object) <-
# 'HTO_classification' before plotting. Here, you can see that each of the small clouds on the
# tSNE plot corresponds to one of the 28 possible doublet combinations.
```

Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.

```
# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(pbmc.hashtag, assay = "HTO", ncells = 5000)
```

![](hashing_vignette_files/figure-html/hashtag_heatmap-1.png)

Cluster and visualize cells using the usual scRNA-seq workflow, and examine for the potential presence of batch effects.

```
# Extract the singlets
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

# Select the top 1000 most variable features
pbmc.singlet <- FindVariableFeatures(pbmc.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
pbmc.singlet <- ScaleData(pbmc.singlet, features = VariableFeatures(pbmc.singlet))

# Run PCA
pbmc.singlet <- RunPCA(pbmc.singlet, features = VariableFeatures(pbmc.singlet))
```

```
# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
pbmc.singlet <- FindNeighbors(pbmc.singlet, reduction = "pca", dims = 1:10)
pbmc.singlet <- FindClusters(pbmc.singlet, resolution = 0.6, verbose = FALSE)
pbmc.singlet <- RunTSNE(pbmc.singlet, reduction = "pca", dims = 1:10)

# Projecting singlet identities on TSNE visualization
DimPlot(pbmc.singlet, group.by = "HTO_classification")
```

![](hashing_vignette_files/figure-html/hashtag_tsne-1.png)

## 12-HTO dataset from four human cell lines

Dataset description:

* Data represent single cells collected from four cell lines: HEK, K562, KG1 and THP1
* Each cell line was further split into three samples (12 samples in total).
* Each sample was labeled with a hashing antibody mixture (CD29 and CD45), pooled, and run on a single lane of 10X.
* Based on this design, we should be able to detect doublets both across and within cell types
* You can download the count matrices for RNA and HTO [here](https://www.dropbox.com/sh/c5gcjm35nglmvcv/AABGz9VO6gX9bVr5R2qahTZha?dl=0), and are available on GEO [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108313)

### Create Seurat object, add HTO data and perform normalization

```
# Read in UMI count matrix for RNA
hto12.umis <- readRDS("/brahms/shared/vignette-data/hto12_umi_mtx.rds")

# Read in HTO count matrix
hto12.htos <- readRDS("/brahms/shared/vignette-data/hto12_hto_mtx.rds")

# Select cell barcodes detected in both RNA and HTO
cells.use <- intersect(rownames(hto12.htos), colnames(hto12.umis))

# Create Seurat object and add HTO data
hto12 <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(hto12.umis[, cells.use]), sparse = T),
    min.features = 300)
hto12[["HTO"]] <- CreateAssayObject(counts = t(x = hto12.htos[colnames(hto12), 1:12]))

# Normalize data
hto12 <- NormalizeData(hto12)
hto12 <- NormalizeData(hto12, assay = "HTO", normalization.method = "CLR")
```

### Demultiplex data

```
hto12 <- HTODemux(hto12, assay = "HTO", positive.quantile = 0.99)
```

### Visualize demultiplexing results

Distribution of selected HTOs grouped by classification, displayed by ridge plots

```
RidgePlot(hto12, assay = "HTO", features = c("HEK-A", "K562-B", "KG1-A", "THP1-C"), ncol = 2)
```

![](hashing_vignette_files/figure-html/ridgeplot-1.png)

Visualize HTO signals in a heatmap

```
HTOHeatmap(hto12, assay = "HTO")
```

![](hashing_vignette_files/figure-html/heatmap-1.png)

### Visualize RNA clustering

- Below, we cluster the cells using our standard scRNA-seq workflow. As expected we see four major clusters, corresponding to the cell lines
- In addition, we see small clusters in between, representing mixed transcriptomes that are correctly annotated as doublets.
- We also see within-cell type doublets, that are (perhaps unsurprisingly) intermixed with singlets of the same cell type

```
# Remove the negative cells
hto12 <- subset(hto12, idents = "Negative", invert = TRUE)

# Run PCA on most variable features
hto12 <- FindVariableFeatures(hto12, selection.method = "mean.var.plot")
hto12 <- ScaleData(hto12, features = VariableFeatures(hto12))
hto12 <- RunPCA(hto12)
hto12 <- RunTSNE(hto12, dims = 1:5, perplexity = 100)
DimPlot(hto12) + NoLegend()
```

![](hashing_vignette_files/figure-html/hto_sub_tsne-1.png)

**Session Info**

```
sessionInfo()
```

```
## R version 4.2.2 Patched (2022-11-10 r83330)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.6 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] Seurat_5.0.0       SeuratObject_5.0.0 sp_2.1-1          
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.7.1       Rtsne_0.16             colorspace_2.1-0      
##   [4] deldir_1.0-9           ellipsis_0.3.2         ggridges_0.5.4        
##   [7] rprojroot_2.0.4        RcppHNSW_0.5.0         fs_1.6.3              
##  [10] spatstat.data_3.0-3    farver_2.1.1           leiden_0.4.3          
##  [13] listenv_0.9.0          ggrepel_0.9.4          RSpectra_0.16-1       
##  [16] fansi_1.0.5            codetools_0.2-19       splines_4.2.2         
##  [19] cachem_1.0.8           knitr_1.45             polyclip_1.10-6       
##  [22] spam_2.10-0            jsonlite_1.8.7         ica_1.0-3             
##  [25] cluster_2.1.4          png_0.1-8              uwot_0.1.16           
##  [28] spatstat.sparse_3.0-3  sctransform_0.4.1      shiny_1.7.5.1         
##  [31] compiler_4.2.2         httr_1.4.7             Matrix_1.6-1.1        
##  [34] fastmap_1.1.1          lazyeval_0.2.2         cli_3.6.1             
##  [37] later_1.3.1            formatR_1.14           htmltools_0.5.6.1     
##  [40] tools_4.2.2            igraph_1.5.1           dotCall64_1.1-0       
##  [43] gtable_0.3.4           glue_1.6.2             reshape2_1.4.4        
##  [46] RANN_2.6.1             dplyr_1.1.3            Rcpp_1.0.11           
##  [49] scattermore_1.2        jquerylib_0.1.4        pkgdown_2.0.7         
##  [52] vctrs_0.6.4            nlme_3.1-162           spatstat.explore_3.2-5
##  [55] progressr_0.14.0       lmtest_0.9-40          spatstat.random_3.2-1 
##  [58] xfun_0.40              stringr_1.5.0          globals_0.16.2        
##  [61] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.4       
##  [64] irlba_2.3.5.1          goftest_1.2-3          future_1.33.0         
##  [67] MASS_7.3-58.2          zoo_1.8-12             scales_1.2.1          
##  [70] spatstat.utils_3.0-4   ragg_1.2.5             promises_1.2.1        
##  [73] parallel_4.2.2         RColorBrewer_1.1-3     yaml_2.3.7            
##  [76] gridExtra_2.3          memoise_2.0.1          reticulate_1.34.0     
##  [79] pbapply_1.7-2          ggrastr_1.0.1          ggplot2_3.4.4         
##  [82] sass_0.4.7             stringi_1.7.12         highr_0.10            
##  [85] desc_1.4.2             fastDummies_1.7.3      rlang_1.1.1           
##  [88] pkgconfig_2.0.3        systemfonts_1.0.4      matrixStats_1.0.0     
##  [91] evaluate_0.22          lattice_0.21-9         tensor_1.5            
##  [94] ROCR_1.0-11            purrr_1.0.2            labeling_0.4.3        
##  [97] patchwork_1.1.3        htmlwidgets_1.6.2      cowplot_1.1.1         
## [100] tidyselect_1.2.0       parallelly_1.36.0      RcppAnnoy_0.0.21      
## [103] plyr_1.8.9             magrittr_2.0.3         R6_2.5.1              
## [106] generics_0.1.3         withr_2.5.2            pillar_1.9.0          
## [109] fitdistrplus_1.1-11    abind_1.4-5            survival_3.5-7        
## [112] tibble_3.2.1           future.apply_1.11.0    KernSmooth_2.23-22    
## [115] utf8_1.2.4             spatstat.geom_3.2-7    plotly_4.10.3         
## [118] rmarkdown_2.25         grid_4.2.2             data.table_1.14.8     
## [121] digest_0.6.33          xtable_1.8-4           tidyr_1.3.0           
## [124] httpuv_1.6.12          textshaping_0.3.6      munsell_0.5.0         
## [127] beeswarm_0.4.0         viridisLite_0.4.2      vipor_0.4.5           
## [130] bslib_0.5.1
```
