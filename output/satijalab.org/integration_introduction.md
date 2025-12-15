---
source_url: https://satijalab.org/seurat/articles/integration_introduction
download_date: 2025-11-20
---

# Introduction to scRNA-seq integration

#### Compiled: November 16, 2023

Source: [`vignettes/integration_introduction.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/integration_introduction.Rmd)

`integration_introduction.Rmd`

## Introduction to scRNA-seq integration

Integration of single-cell sequencing datasets, for example across experimental batches, donors, or conditions, is often an important step in scRNA-seq workflows. Integrative analysis can help to match shared cell types and states across datasets, which can boost statistical power, and most importantly, facilitate accurate comparative analysis across datasets. In previous versions of Seurat we introduced methods for integrative analysis, including our ‘anchor-based’ integration workflow. Many labs have also published powerful and pioneering methods, including [Harmony](https://portals.broadinstitute.org/harmony/) and [scVI](https://docs.scvi-tools.org/en/stable/index.html), for integrative analysis. Please see our [Integrating scRNA-seq data with multiple tools](https://satijalab.org/seurat/articles/seurat5_integration) vignette.

## Integration goals

The following tutorial is designed to give you an overview of the kinds of comparative analyses on complex cell types that are possible using the Seurat integration procedure. Here, we address a few key goals:

* Identify cell subpopulations that are present in both datasets
* Obtain cell type markers that are conserved in both control and stimulated cells
* Compare the datasets to find cell-type specific responses to stimulation

## Setup the Seurat objects

For convenience, we distribute this dataset through our [SeuratData](https://github.com/satijalab/seurat-data) package.

```
library(Seurat)
library(SeuratData)
library(patchwork)
```

```
# install dataset
InstallData("ifnb")
```

The object contains data from human PBMC from two conditions, interferon-stimulated and control cells (stored in the `stim` column in the object metadata). We will aim to integrate the two conditions together, so that we can jointly identify cell subpopulations across datasets, and then explore how each group differs across conditions

In previous versions of Seurat, we would require the data to be represented as two different Seurat objects. In Seurat v5, we keep all the data in one object, but simply split it into multiple ‘layers’. To learn more about layers, check out our [Seurat object interaction vignette](https://satijalab.org/seurat/articles/interaction_vignette).

```
# load dataset
ifnb <- LoadData("ifnb")
# split the RNA measurements into two layers one for control cells, one for stimulated cells

ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb
```

## Perform analysis without integration

We can first analyze the dataset without integration. The resulting clusters are defined both by cell type and stimulation condition, which creates challenges for downstream analysis.

```
# run standard anlaysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
```

```
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 13999
## Number of edges: 555146
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8153
## Number of communities: 26
## Elapsed time: 23 seconds
```

```
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))
```

![](integration_introduction_files/figure-html/unintegratedUMAP-1.png)

## Perform integration

We now aim to integrate data from the two conditions, so that cells from the same cell type/subpopulation will cluster together.

We often refer to this procedure as intergration/alignment. When aligning two genome sequences together, identification of shared/homologous regions can help to interpret differences between the sequences as well. Similarly for scRNA-seq integration, our goal is not to remove biological differences across conditions, but to learn shared cell types/states in an initial step - specifically because that will enable us to compare control stimulated and control profiles for these individual cell types.

The Seurat v5 integration procedure aims to return a single dimensional reduction that captures the shared sources of variance across multiple layers, so that cells in a similar biological state will cluster. The method returns a dimensional reduction (i.e. `integrated.cca`) which can be used for visualization and unsupervised clustering analysis. For evaluating performance, we can use cell type labels that are pre-loaded in the `seurat_annotations` metadata column.

```
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE)

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 13999
## Number of edges: 588593
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8454
## Number of communities: 17
## Elapsed time: 26 seconds
```

```
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
```

```
# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
```

![](integration_introduction_files/figure-html/viz-1.png)

To visualize the two conditions side-by-side, we can use the `split.by` argument to show each condition colored by cluster.

```
DimPlot(ifnb, reduction = "umap", split.by = "stim")
```

![](integration_introduction_files/figure-html/split.dim-1.png)

## Identify conserved cell type markers

To identify canonical cell type marker genes that are conserved across conditions, we provide the `FindConservedMarkers()` function. This function performs differential gene expression testing for each dataset/group and combines the p-values using meta-analysis methods from the MetaDE R package. For example, we can calculated the genes that are conserved markers irrespective of stimulation condition in cluster 6 (NK cells).

```
Idents(ifnb) <- "seurat_annotations"
nk.markers <- FindConservedMarkers(ifnb, ident.1 = "NK", grouping.var = "stim", verbose = FALSE)
head(nk.markers)
```

```
##       CTRL_p_val CTRL_avg_log2FC CTRL_pct.1 CTRL_pct.2 CTRL_p_val_adj
## GNLY           0        6.854586      0.943      0.046              0
## NKG7           0        5.358881      0.953      0.085              0
## GZMB           0        5.078135      0.839      0.044              0
## CLIC3          0        5.765314      0.601      0.024              0
## CTSW           0        5.307246      0.537      0.030              0
## KLRD1          0        5.261553      0.507      0.019              0
##       STIM_p_val STIM_avg_log2FC STIM_pct.1 STIM_pct.2 STIM_p_val_adj max_pval
## GNLY           0        6.435910      0.956      0.059              0        0
## NKG7           0        4.971397      0.950      0.081              0        0
## GZMB           0        5.151924      0.897      0.060              0        0
## CLIC3          0        5.505208      0.623      0.031              0        0
## CTSW           0        5.240729      0.592      0.035              0        0
## KLRD1          0        4.852457      0.555      0.027              0        0
##       minimump_p_val
## GNLY               0
## NKG7               0
## GZMB               0
## CLIC3              0
## CTSW               0
## KLRD1              0
```

You can perform these same analysis on the unsupervised clustering results (stored in `seurat_clusters`), and use these conserved markers to annotate cell types in your dataset.

The `DotPlot()` function with the `split.by` parameter can be useful for viewing conserved cell type markers across conditions, showing both the expression level and the percentage of cells in a cluster expressing any given gene. Here we plot 2-3 strong marker genes for each of our 14 clusters.

```
# NEEDS TO BE FIXED AND SET ORDER CORRECTLY
Idents(ifnb) <- factor(Idents(ifnb), levels = c("pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono",
    "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))

markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
    "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
DotPlot(ifnb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
    RotatedAxis()
```

![](integration_introduction_files/figure-html/splitdotplot-1.png)

## Identify differential expressed genes across conditions

Now that we’ve aligned the stimulated and control cells, we can start to do comparative analyses and look at the differences induced by stimulation.

We can aggregate cells of a similar type and condition together to create “pseudobulk” profiles using the `AggregateExpression` command. As an initial exploratory analysis, we can compare pseudobulk profiles of two cell types (naive CD4 T cells, and CD14 monocytes), and compare their gene expression profiles before and after stimulation. We highlight genes that exhibit dramatic responses to interferon stimulation. As you can see, many of the same genes are upregulated in both of these cell types and likely represent a conserved interferon response pathway, though CD14 monocytes exhibit a stronger transcriptional response.

```
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

aggregate_ifnb <- AggregateExpression(ifnb, group.by = c("seurat_annotations", "stim"), return.seurat = TRUE)
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")

p1 <- CellScatter(aggregate_ifnb, "CD14 Mono_CTRL", "CD14 Mono_STIM", highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

p3 <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", highlight = genes.to.label)
p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)

p2 + p4
```

![](integration_introduction_files/figure-html/scatterplots-1.png)

We can now ask what genes change in different conditions for cells of the same type. First, we create a column in the meta.data slot to hold both the cell type and stimulation information and switch the current ident to that column. Then we use `FindMarkers()` to find the genes that are different between stimulated and control B cells. Notice that many of the top genes that show up here are the same as the ones we plotted earlier as core interferon response genes. Additionally, genes like CXCL10 which we saw were specific to monocyte and B cell interferon response show up as highly significant in this list as well.

Please note that p-values obtained from this analysis should be interpreted with caution, as these tests treat each cell as an independent replicate, and ignore inherent correlations between cells originating from the same sample. As discussed [here](https://pubmed.ncbi.nlm.nih.gov/33257685/), DE tests across multiple conditions should expressly utilize multiple samples/replicates, and can be performed after aggregating (‘pseudobulking’) cells from the same sample and subpopulation together. We do not perform this analysis here, as there is a single replicate in the data, but please see our [vignette comparing healthy and diabetic samples](https://satijalab.org/seurat/articles/parsebio_sketch_integration) as an example for how to perform DE analysis across conditions.

```
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)
```

```
##                 p_val avg_log2FC pct.1 pct.2     p_val_adj
## ISG15   5.387767e-159  5.0588481 0.998 0.233 7.571429e-155
## IFIT3   1.945114e-154  6.1124940 0.965 0.052 2.733468e-150
## IFI6    2.503565e-152  5.4933132 0.965 0.076 3.518260e-148
## ISG20   6.492570e-150  3.0549593 1.000 0.668 9.124009e-146
## IFIT1   1.951022e-139  6.2320388 0.907 0.029 2.741772e-135
## MX1     6.897626e-123  3.9798482 0.905 0.115 9.693234e-119
## LY6E    2.825649e-120  3.7907800 0.898 0.150 3.970885e-116
## TNFSF10 4.007285e-112  6.5802175 0.786 0.020 5.631437e-108
## IFIT2   2.672552e-108  5.5525558 0.786 0.037 3.755738e-104
## B2M      5.283684e-98  0.6104044 1.000 1.000  7.425161e-94
## PLSCR1   4.634658e-96  3.8010721 0.793 0.113  6.513085e-92
## IRF7     2.411149e-94  3.1992949 0.835 0.187  3.388388e-90
## CXCL10   3.708508e-86  8.0906108 0.651 0.010  5.211566e-82
## UBE2L6   5.547472e-83  2.5167981 0.851 0.297  7.795863e-79
## PSMB9    1.716262e-77  1.7715351 0.937 0.568  2.411863e-73
```

Another useful way to visualize these changes in gene expression is with the `split.by` option to the `FeaturePlot()` or `VlnPlot()` function. This will display FeaturePlots of the list of given genes, split by a grouping variable (stimulation condition here). Genes such as CD3D and GNLY are canonical cell type markers (for T cells and NK/CD8 T cells) that are virtually unaffected by interferon stimulation and display similar gene expression patterns in the control and stimulated group. IFI6 and ISG15, on the other hand, are core interferon response genes and are upregulated accordingly in all cell types. Finally, CD14 and CXCL10 are genes that show a cell type specific interferon response. CD14 expression decreases after stimulation in CD14 monocytes, which could lead to misclassification in a supervised analysis framework, underscoring the value of integrated analysis. CXCL10 shows a distinct upregulation in monocytes and B cells after interferon stimulation but not in other cell types.

```
FeaturePlot(ifnb, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3, cols = c("grey",
    "red"), reduction = "umap")
```

![](integration_introduction_files/figure-html/feature.heatmaps-1.png)

```
plots <- VlnPlot(ifnb, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "seurat_annotations",
    pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
```

![](integration_introduction_files/figure-html/splitvln-1.png)

## Perform integration with SCTransform-normalized datasets

As an alternative to log-normalization, Seurat also includes support for preprocessing of scRNA-seq using the [sctransform workflow](https://satijalab.org/seurat/articles/sctransform_vignette). The `IntegrateLayers` function also supports SCTransform-normalized data, by setting the `normalization.method` parameter, as shown below.

```
ifnb <- LoadData("ifnb")

# split datasets and process without integration
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- RunUMAP(ifnb, dims = 1:30)
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
```

![](integration_introduction_files/figure-html/sct-1.png)

```
# integrate datasets
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, normalization.method = "SCT", verbose = F)
ifnb <- FindNeighbors(ifnb, reduction = "integrated.dr", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 0.6)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 13999
## Number of edges: 528127
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9060
## Number of communities: 20
## Elapsed time: 21 seconds
```

```
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.dr")
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))
```

![](integration_introduction_files/figure-html/sct-2.png)

```
# perform differential expression
ifnb <- PrepSCTFindMarkers(ifnb)
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
```

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
##  [1] cowplot_1.1.1                  ggplot2_3.4.4                 
##  [3] patchwork_1.1.3                thp1.eccite.SeuratData_3.1.5  
##  [5] stxBrain.SeuratData_0.1.1      ssHippo.SeuratData_3.1.4      
##  [7] pbmcsca.SeuratData_3.0.0       pbmcref.SeuratData_1.0.0      
##  [9] pbmcMultiome.SeuratData_0.1.4  pbmc3k.SeuratData_3.1.4       
## [11] panc8.SeuratData_3.0.2         ifnb.SeuratData_3.0.0         
## [13] hcabm40k.SeuratData_3.0.0      cbmc.SeuratData_3.1.4         
## [15] bonemarrowref.SeuratData_1.0.0 bmcite.SeuratData_0.3.0       
## [17] SeuratData_0.2.2.9001          Seurat_5.0.1             
## [19] testthat_3.2.0                 SeuratObject_5.0.1            
## [21] sp_2.1-1                      
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.4                  spatstat.explore_3.2-5     
##   [3] reticulate_1.34.0           tidyselect_1.2.0           
##   [5] htmlwidgets_1.6.2           grid_4.2.2                 
##   [7] Rtsne_0.16                  devtools_2.4.5             
##   [9] munsell_0.5.0               mutoss_0.1-13              
##  [11] codetools_0.2-19            ragg_1.2.5                 
##  [13] ica_1.0-3                   future_1.33.0              
##  [15] miniUI_0.1.1.1              withr_2.5.2                
##  [17] spatstat.random_3.2-1       colorspace_2.1-0           
##  [19] progressr_0.14.0            Biobase_2.58.0             
##  [21] highr_0.10                  knitr_1.45                 
##  [23] rstudioapi_0.14             stats4_4.2.2               
##  [25] ROCR_1.0-11                 tensor_1.5                 
##  [27] listenv_0.9.0               MatrixGenerics_1.10.0      
##  [29] Rdpack_2.5                  labeling_0.4.3             
##  [31] GenomeInfoDbData_1.2.9      mnormt_2.1.1               
##  [33] polyclip_1.10-6             farver_2.1.1               
##  [35] rprojroot_2.0.4             TH.data_1.1-2              
##  [37] parallelly_1.36.0           vctrs_0.6.4                
##  [39] generics_0.1.3              xfun_0.40                  
##  [41] GenomeInfoDb_1.34.9         R6_2.5.1                   
##  [43] ggbeeswarm_0.7.1            DelayedArray_0.24.0        
##  [45] bitops_1.0-7                spatstat.utils_3.0-4       
##  [47] cachem_1.0.8                promises_1.2.1             
##  [49] scales_1.2.1                multcomp_1.4-25            
##  [51] beeswarm_0.4.0              gtable_0.3.4               
##  [53] globals_0.16.2              processx_3.8.2             
##  [55] goftest_1.2-3               spam_2.10-0                
##  [57] sandwich_3.0-2              rlang_1.1.1                
##  [59] systemfonts_1.0.4           splines_4.2.2              
##  [61] lazyeval_0.2.2              spatstat.geom_3.2-7        
##  [63] yaml_2.3.7                  reshape2_1.4.4             
##  [65] abind_1.4-5                 httpuv_1.6.12              
##  [67] tools_4.2.2                 usethis_2.1.6              
##  [69] ellipsis_0.3.2              jquerylib_0.1.4            
##  [71] RColorBrewer_1.1-3          BiocGenerics_0.44.0        
##  [73] sessioninfo_1.2.2           ggridges_0.5.4             
##  [75] TFisher_0.2.0               Rcpp_1.0.11                
##  [77] plyr_1.8.9                  sparseMatrixStats_1.10.0   
##  [79] zlibbioc_1.44.0             RCurl_1.98-1.12            
##  [81] purrr_1.0.2                 ps_1.7.5                   
##  [83] prettyunits_1.2.0           deldir_1.0-9               
##  [85] pbapply_1.7-2               urlchecker_1.0.1           
##  [87] S4Vectors_0.36.2            zoo_1.8-12                 
##  [89] SummarizedExperiment_1.28.0 ggrepel_0.9.4              
##  [91] cluster_2.1.4               fs_1.6.3                   
##  [93] magrittr_2.0.3              glmGamPoi_1.10.2           
##  [95] data.table_1.14.8           RSpectra_0.16-1            
##  [97] scattermore_1.2             lmtest_0.9-40              
##  [99] RANN_2.6.1                  mvtnorm_1.2-3              
## [101] fitdistrplus_1.1-11         matrixStats_1.0.0          
## [103] pkgload_1.3.3               mime_0.12                  
## [105] evaluate_0.22               xtable_1.8-4               
## [107] IRanges_2.32.0              fastDummies_1.7.3          
## [109] gridExtra_2.3               compiler_4.2.2             
## [111] tibble_3.2.1                KernSmooth_2.23-22         
## [113] crayon_1.5.2                htmltools_0.5.6.1          
## [115] later_1.3.1                 tidyr_1.3.0                
## [117] formatR_1.14                MASS_7.3-58.2              
## [119] rappdirs_0.3.3              Matrix_1.6-1.1             
## [121] brio_1.1.3                  cli_3.6.1                  
## [123] rbibutils_2.2.15            qqconf_1.3.2               
## [125] parallel_4.2.2              dotCall64_1.1-0            
## [127] metap_1.9                   igraph_1.5.1               
## [129] GenomicRanges_1.50.2        pkgconfig_2.0.3            
## [131] sn_2.1.1                    pkgdown_2.0.7              
## [133] numDeriv_2016.8-1.1         plotly_4.10.3              
## [135] spatstat.sparse_3.0-3       vipor_0.4.5                
## [137] bslib_0.5.1                 XVector_0.38.0             
## [139] multtest_2.54.0             stringr_1.5.0              
## [141] callr_3.7.3                 digest_0.6.33              
## [143] sctransform_0.4.1           RcppAnnoy_0.0.21           
## [145] spatstat.data_3.0-3         rmarkdown_2.25             
## [147] leiden_0.4.3                uwot_0.1.16                
## [149] DelayedMatrixStats_1.20.0   shiny_1.7.5.1              
## [151] lifecycle_1.0.4             nlme_3.1-162               
## [153] jsonlite_1.8.7              limma_3.54.1               
## [155] desc_1.4.2                  viridisLite_0.4.2          
## [157] fansi_1.0.5                 pillar_1.9.0               
## [159] lattice_0.21-9              ggrastr_1.0.1              
## [161] plotrix_3.8-2               fastmap_1.1.1              
## [163] httr_1.4.7                  pkgbuild_1.4.2             
## [165] survival_3.5-7              glue_1.6.2                 
## [167] remotes_2.4.2.1             png_0.1-8                  
## [169] presto_1.0.0                stringi_1.7.12             
## [171] sass_0.4.7                  profvis_0.3.7              
## [173] textshaping_0.3.6           RcppHNSW_0.5.0             
## [175] memoise_2.0.1               mathjaxr_1.6-0             
## [177] dplyr_1.1.3                 irlba_2.3.5.1              
## [179] future.apply_1.11.0
```