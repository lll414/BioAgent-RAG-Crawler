---
source_url: https://satijalab.org/seurat/articles/seurat5_integration
download_date: 2025-11-20
---

# Integrative analysis in Seurat v5

#### Compiled: October 31, 2023

Source: [`vignettes/seurat5_integration.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/seurat5_integration.Rmd)

`seurat5_integration.Rmd`

```
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)
```

## Introduction

Integration of single-cell sequencing datasets, for example across experimental batches, donors, or conditions, is often an important step in scRNA-seq workflows. Integrative analysis can help to match shared cell types and states across datasets, which can boost statistical power, and most importantly, facilitate accurate comparative analysis across datasets. In previous versions of Seurat we introduced methods for integrative analysis, including our ‘anchor-based’ integration workflow. Many labs have also published powerful and pioneering methods, including [Harmony](https://github.com/immunogenomics/harmony) and [scVI](https://yoseflab.github.io/software/scvi-tools/), for integrative analysis. We recognize that while the goal of matching shared cell types across datasets may be important for many problems, users may also be concerned about which method to use, or that integration could result in a loss of biological resolution. In Seurat v5, we introduce more flexible and streamlined infrastructure to run different integration algorithms with a single line of code. This makes it easier to explore the results of different integration methods, and to compare these results to a workflow that excludes integration steps. For this vignette, we use a [dataset of human PBMC profiled with seven different technologies](https://www.nature.com/articles/s41587-020-0465-8), profiled as part of a systematic comparative analysis (`pbmcsca`). The data is available as part of our [SeuratData](https://github.com/satijalab/seurat-data) package.

## Layers in the Seurat v5 object

Seurat v5 assays store data in layers. These layers can store raw, un-normalized counts (`layer='counts'`), normalized data (`layer='data'`), or z-scored/variance-stabilized data (`layer='scale.data'`). We can load in the data, remove low-quality cells, and obtain predicted cell annotations (which will be useful for assessing integration later), using our [Azimuth pipeline](https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html).

```
# load in the pbmc systematic comparative analysis dataset
obj <- LoadData("pbmcsca")
obj <- subset(obj, nFeature_RNA > 1000)
obj <- RunAzimuth(obj, reference = "pbmcref")
# currently, the object has two layers in the RNA assay: counts, and data
obj
```

```
## An object of class Seurat 
## 33789 features across 10434 samples within 4 assays 
## Active assay: RNA (33694 features, 0 variable features)
##  2 layers present: counts, data
##  3 other assays present: prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3
##  2 dimensional reductions calculated: integrated_dr, ref.umap
```

The object contains data from nine different batches (stored in the `Method` column in the object metadata), representing seven different technologies. We will aim to integrate the different batches together. In previous versions of Seurat, we would require the data to be represented as nine different Seurat objects. When using Seurat v5 assays, we can instead keep all the data in one object, but simply split the layers. After splitting, there are now 18 layers (a `counts` and `data` layer for each batch). We can also run a standard scRNA-seq analysis (i.e. without integration). Note that since the data is split into layers, normalization and variable feature identification is performed for each batch independently (a consensus set of variable features is automatically identified).

```
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
obj
```

```
## An object of class Seurat 
## 33789 features across 10434 samples within 4 assays 
## Active assay: RNA (33694 features, 0 variable features)
##  18 layers present: counts.Smart-seq2, counts.CEL-Seq2, counts.10x_Chromium_v2_A, counts.10x_Chromium_v2_B, counts.10x_Chromium_v3, counts.Drop-seq, counts.Seq-Well, counts.inDrops, counts.10x_Chromium_v2, data.Smart-seq2, data.CEL-Seq2, data.10x_Chromium_v2_A, data.10x_Chromium_v2_B, data.10x_Chromium_v3, data.Drop-seq, data.Seq-Well, data.inDrops, data.10x_Chromium_v2
##  3 other assays present: prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3
##  2 dimensional reductions calculated: integrated_dr, ref.umap
```

```
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
```

We can now visualize the results of a standard analysis without integration. Note that cells are grouping both by cell type and by underlying method. While a UMAP analysis is just a visualization of this, clustering this dataset would return predominantly batch-specific clusters. Especially if previous cell-type annotations were not available, this would make downstream analysis extremely challenging.

```
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 10434
## Number of edges: 412660
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8981
## Number of communities: 48
## Elapsed time: 0 seconds
```

```
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Method", "predicted.celltype.l2"))
```

![](seurat5_integration_files/figure-html/unintegratedUMAP-1.png)

## Perform streamlined (one-line) integrative analysis

Seurat v5 enables streamlined integrative analysis using the `IntegrateLayers` function. The method currently supports five integration methods. Each of these methods performs integration in low-dimensional space, and returns a dimensional reduction (i.e. `integrated.rpca`) that aims to co-embed shared cell types across batches:

* Anchor-based CCA integration (`method=CCAIntegration`)
* Anchor-based RPCA integration (`method=RPCAIntegration`)
* Harmony (`method=HarmonyIntegration`)
* FastMNN (`method= FastMNNIntegration`)
* scVI (`method=scVIIntegration`)

Note that our anchor-based RPCA integration represents a faster and more conservative (less correction) method for integration. For interested users, we discuss this method in more detail in our [previous RPCA vignette](https://satijalab.org/seurat/articles/integration_rpca)

You can find more detail on each method, and any installation prerequisites, in Seurat’s documentation (for example, `?scVIIntegration`). For example, scVI integration requires `reticulate` which can be installed from CRAN (`install.packages("reticulate")`) as well as `scvi-tools` and its dependencies installed in a conda environment. Please see scVI installation instructions [here](https://docs.scvi-tools.org/en/stable/installation.html).

Each of the following lines perform a new integration using a single line of code:

```
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
```

```
obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)
```

```
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
```

```
obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
```

```
obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "../miniconda3/envs/scvi-env", verbose = FALSE
)
```

For any of the methods, we can now visualize and cluster the datasets. We show this for CCA integration and scVI, but you can do this for any method

```
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "cca_clusters")
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 10434
## Number of edges: 617402
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8037
## Number of communities: 25
## Elapsed time: 1 seconds
```

```
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("Method", "predicted.celltype.l2", "cca_clusters"),
  combine = FALSE, label.size = 2
)

obj <- FindNeighbors(obj, reduction = "integrated.scvi", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "scvi_clusters")
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 10434
## Number of edges: 354664
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.7942
## Number of communities: 22
## Elapsed time: 1 seconds
```

```
obj <- RunUMAP(obj, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
p2 <- DimPlot(
  obj,
  reduction = "umap.scvi",
  group.by = c("Method", "predicted.celltype.l2", "scvi_clusters"),
  combine = FALSE, label.size = 2
)

wrap_plots(c(p1, p2), ncol = 2, byrow = F)
```

![](seurat5_integration_files/figure-html/integratedprojections-1.png)

We hope that by simplifying the process of performing integrative analysis, users can more carefully evaluate the biological information retained in the integrated dataset. For example, users can compare the expression of biological markers based on different clustering solutions, or visualize one method’s clustering solution on different UMAP visualizations.

```
p1 <- VlnPlot(
  obj,
  features = "rna_CD8A", group.by = "unintegrated_clusters"
) + NoLegend() + ggtitle("CD8A - Unintegrated Clusters")
p2 <- VlnPlot(
  obj, "rna_CD8A",
  group.by = "cca_clusters"
) + NoLegend() + ggtitle("CD8A - CCA Clusters")
p3 <- VlnPlot(
  obj, "rna_CD8A",
  group.by = "scvi_clusters"
) + NoLegend() + ggtitle("CD8A - scVI Clusters")
p1 | p2 | p3
```

![](seurat5_integration_files/figure-html/vlnplots-1.png)

```
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
p4 <- DimPlot(obj, reduction = "umap.unintegrated", group.by = c("cca_clusters"))
p5 <- DimPlot(obj, reduction = "umap.rpca", group.by = c("cca_clusters"))
p6 <- DimPlot(obj, reduction = "umap.scvi", group.by = c("cca_clusters"))
p4 | p5 | p6
```

![](seurat5_integration_files/figure-html/umaps-1.png)

Once integrative analysis is complete, you can rejoin the layers - which collapses the individual datasets together and recreates the original `counts` and `data` layers. You will need to do this before performing any differential expression analysis. However, you can always resplit the layers in case you would like to reperform integrative analysis.

```
obj <- JoinLayers(obj)
obj
```

```
## An object of class Seurat 
## 35789 features across 10434 samples within 5 assays 
## Active assay: RNA (33694 features, 2000 variable features)
##  3 layers present: data, counts, scale.data
##  4 other assays present: prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3, mnn.reconstructed
##  12 dimensional reductions calculated: integrated_dr, ref.umap, pca, umap.unintegrated, integrated.cca, integrated.rpca, harmony, integrated.mnn, integrated.scvi, umap.cca, umap.scvi, umap.rpca
```

Lastly, users can also perform integration using sctransform-normalized data (see our [SCTransform vignette](https://satijalab.org/seurat/articles/sctransform_vignette) for more information), by first running SCTransform normalization, and then setting the `normalization.method` argument in `IntegrateLayers`.

```
options(future.globals.maxSize = 3e+09)
obj <- SCTransform(obj)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = F
)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "integrated.dr")
obj <- FindClusters(obj, resolution = 2)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 10434
## Number of edges: 449127
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8354
## Number of communities: 26
## Elapsed time: 1 seconds
```

```
obj <- RunUMAP(obj, dims = 1:30, reduction = "integrated.dr")
```

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
##  [1] patchwork_1.1.3               ggplot2_3.4.4                
##  [3] Azimuth_0.4.6.9004            shinyBS_0.61.1               
##  [5] SeuratWrappers_0.3.19         thp1.eccite.SeuratData_3.1.5 
##  [7] stxBrain.SeuratData_0.1.1     ssHippo.SeuratData_3.1.4     
##  [9] pbmcsca.SeuratData_3.0.0      pbmcref.SeuratData_1.0.0     
## [11] pbmcMultiome.SeuratData_0.1.4 pbmc3k.SeuratData_3.1.4      
## [13] panc8.SeuratData_3.0.2        ifnb.SeuratData_3.0.0        
## [15] hcabm40k.SeuratData_3.0.0     cbmc.SeuratData_3.1.4        
## [17] bmcite.SeuratData_0.3.0       SeuratData_0.2.2.9001        
## [19] Seurat_5.0.1                  SeuratObject_5.0.1           
## [21] sp_2.1-1                     
## 
## loaded via a namespace (and not attached):
##   [1] rsvd_1.0.5                        ica_1.0-3                        
##   [3] RcppRoll_0.3.0                    Rsamtools_2.14.0                 
##   [5] lmtest_0.9-40                     rprojroot_2.0.3                  
##   [7] crayon_1.5.2                      MASS_7.3-58.2                    
##   [9] rhdf5filters_1.10.1               nlme_3.1-162                     
##  [11] rlang_1.1.1                       XVector_0.38.0                   
##  [13] ROCR_1.0-11                       irlba_2.3.5.1                    
##  [15] filelock_1.0.2                    BiocParallel_1.32.6              
##  [17] rjson_0.2.21                      CNEr_1.34.0                      
##  [19] bit64_4.0.5                       glue_1.6.2                       
##  [21] harmony_1.0.3                     poweRlaw_0.70.6                  
##  [23] sctransform_0.4.1                 vipor_0.4.5                      
##  [25] parallel_4.2.2                    spatstat.sparse_3.0-3            
##  [27] AnnotationDbi_1.60.2              SeuratDisk_0.0.0.9020            
##  [29] BiocGenerics_0.44.0               shinydashboard_0.7.2             
##  [31] dotCall64_1.1-0                   spatstat.geom_3.2-7              
##  [33] tidyselect_1.2.0                  SummarizedExperiment_1.28.0      
##  [35] fitdistrplus_1.1-11               XML_3.99-0.14                    
##  [37] tidyr_1.3.0                       zoo_1.8-12                       
##  [39] GenomicAlignments_1.34.1          xtable_1.8-4                     
##  [41] RcppHNSW_0.5.0                    magrittr_2.0.3                   
##  [43] evaluate_0.22                     scuttle_1.8.4                    
##  [45] cli_3.6.1                         zlibbioc_1.44.0                  
##  [47] miniUI_0.1.1.1                    bslib_0.5.1                      
##  [49] fastmatch_1.1-4                   ensembldb_2.22.0                 
##  [51] fastDummies_1.7.3                 shiny_1.7.5.1                    
##  [53] BiocSingular_1.14.0               xfun_0.40                        
##  [55] cluster_2.1.4                     caTools_1.18.2                   
##  [57] KEGGREST_1.38.0                   tibble_3.2.1                     
##  [59] ggrepel_0.9.4                     listenv_0.9.0                    
##  [61] TFMPvalue_0.0.9                   Biostrings_2.66.0                
##  [63] png_0.1-8                         future_1.33.0                    
##  [65] withr_2.5.2                       bitops_1.0-7                     
##  [67] plyr_1.8.9                        cellranger_1.1.0                 
##  [69] AnnotationFilter_1.22.0           pracma_2.4.2                     
##  [71] pillar_1.9.0                      cachem_1.0.8                     
##  [73] GenomicFeatures_1.50.4            fs_1.6.3                         
##  [75] hdf5r_1.3.8                       DelayedMatrixStats_1.20.0        
##  [77] vctrs_0.6.4                       ellipsis_0.3.2                   
##  [79] generics_0.1.3                    tools_4.2.2                      
##  [81] beeswarm_0.4.0                    munsell_0.5.0                    
##  [83] DelayedArray_0.24.0               fastmap_1.1.1                    
##  [85] compiler_4.2.2                    abind_1.4-5                      
##  [87] httpuv_1.6.12                     rtracklayer_1.58.0               
##  [89] plotly_4.10.3                     GenomeInfoDbData_1.2.9           
##  [91] gridExtra_2.3                     lattice_0.21-9                   
##  [93] deldir_1.0-9                      utf8_1.2.4                       
##  [95] later_1.3.1                       dplyr_1.1.3                      
##  [97] BiocFileCache_2.6.1               jsonlite_1.8.7                   
##  [99] scales_1.2.1                      ScaledMatrix_1.6.0               
## [101] pbapply_1.7-2                     sparseMatrixStats_1.10.0         
## [103] lazyeval_0.2.2                    promises_1.2.1                   
## [105] R.utils_2.12.2                    goftest_1.2-3                    
## [107] spatstat.utils_3.0-4              reticulate_1.34.0                
## [109] rmarkdown_2.25                    pkgdown_2.0.7                    
## [111] cowplot_1.1.1                     textshaping_0.3.6                
## [113] glmGamPoi_1.10.2                  Rtsne_0.16                       
## [115] BSgenome_1.66.3                   Biobase_2.58.0                   
## [117] uwot_0.1.16                       igraph_1.5.1                     
## [119] survival_3.5-7                    ResidualMatrix_1.8.0             
## [121] yaml_2.3.7                        systemfonts_1.0.4                
## [123] htmltools_0.5.6.1                 memoise_2.0.1                    
## [125] BiocIO_1.8.0                      IRanges_2.32.0                   
## [127] viridisLite_0.4.2                 digest_0.6.33                    
## [129] RhpcBLASctl_0.23-42               mime_0.12                        
## [131] rappdirs_0.3.3                    spam_2.10-0                      
## [133] RSQLite_2.3.1                     future.apply_1.11.0              
## [135] remotes_2.4.2.1                   data.table_1.14.8                
## [137] blob_1.2.4                        JASPAR2020_0.99.10               
## [139] S4Vectors_0.36.2                  R.oo_1.25.0                      
## [141] TFBSTools_1.36.0                  ragg_1.2.5                       
## [143] styler_1.10.2                     splines_4.2.2                    
## [145] labeling_0.4.3                    Rhdf5lib_1.20.0                  
## [147] googledrive_2.1.1                 ProtGenerics_1.30.0              
## [149] RCurl_1.98-1.12                   hms_1.1.3                        
## [151] rhdf5_2.42.1                      colorspace_2.1-0                 
## [153] BiocManager_1.30.22               ggbeeswarm_0.7.1                 
## [155] GenomicRanges_1.50.2              Signac_1.12.9000                 
## [157] ggrastr_1.0.1                     sass_0.4.7                       
## [159] Rcpp_1.0.11                       BSgenome.Hsapiens.UCSC.hg38_1.4.5
## [161] RANN_2.6.1                        fansi_1.0.5                      
## [163] tzdb_0.4.0                        parallelly_1.36.0                
## [165] R6_2.5.1                          grid_4.2.2                       
## [167] ggridges_0.5.4                    lifecycle_1.0.3                  
## [169] curl_5.1.0                        googlesheets4_1.1.1              
## [171] leiden_0.4.3                      jquerylib_0.1.4                  
## [173] Matrix_1.6-1.1                    desc_1.4.2                       
## [175] RcppAnnoy_0.0.21                  RColorBrewer_1.1-3               
## [177] spatstat.explore_3.2-5            stringr_1.5.0                    
## [179] R.cache_0.16.0                    htmlwidgets_1.6.2                
## [181] beachmat_2.14.2                   polyclip_1.10-6                  
## [183] biomaRt_2.54.1                    purrr_1.0.2                      
## [185] seqLogo_1.64.0                    globals_0.16.2                   
## [187] spatstat.random_3.2-1             progressr_0.14.0                 
## [189] batchelor_1.14.1                  codetools_0.2-19                 
## [191] matrixStats_1.0.0                 GO.db_3.16.0                     
## [193] gtools_3.9.4                      prettyunits_1.1.1                
## [195] SingleCellExperiment_1.20.1       dbplyr_2.3.4                     
## [197] RSpectra_0.16-1                   R.methodsS3_1.8.2                
## [199] GenomeInfoDb_1.34.9               gtable_0.3.4                     
## [201] DBI_1.1.3                         stats4_4.2.2                     
## [203] tensor_1.5                        httr_1.4.7                       
## [205] highr_0.10                        KernSmooth_2.23-22               
## [207] stringi_1.7.12                    presto_1.0.0                     
## [209] progress_1.2.2                    reshape2_1.4.4                   
## [211] farver_2.1.1                      annotate_1.76.0                  
## [213] DT_0.30                           xml2_1.3.5                       
## [215] EnsDb.Hsapiens.v86_2.99.0         shinyjs_2.1.0                    
## [217] BiocNeighbors_1.16.0              restfulr_0.0.15                  
## [219] readr_2.1.4                       scattermore_1.2                  
## [221] bit_4.0.5                         MatrixGenerics_1.10.0            
## [223] spatstat.data_3.0-3               pkgconfig_2.0.3                  
## [225] gargle_1.5.2                      DirichletMultinomial_1.40.0      
## [227] knitr_1.45
```