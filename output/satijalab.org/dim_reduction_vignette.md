---
source_url: https://satijalab.org/seurat/articles/dim_reduction_vignette
download_date: 2025-11-20
---

# Seurat - Dimensional Reduction Vignette

#### Compiled: September 14, 2023

Source: [`vignettes/dim_reduction_vignette.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/dim_reduction_vignette.Rmd)

`dim_reduction_vignette.Rmd`

## Load in the data

This vignette demonstrates how to store and interact with dimensional reduction information (such as the output from `RunPCA()`) in Seurat. For demonstration purposes, we will be using the 2,700 PBMC object that is available via the [SeuratData](https://github.com/satijalab/seurat-data) package.

```
library(Seurat)
library(SeuratData)
pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")
```

## Explore the new dimensional reduction structure

In Seurat v3.0, storing and interacting with dimensional reduction information has been generalized and formalized into the `DimReduc` object. Each dimensional reduction procedure is stored as a `DimReduc` object in the `object@reductions` slot as an element of a named list. Accessing these reductions can be done with the `[[` operator, calling the name of the reduction desired. For example, after running a principle component analysis with `RunPCA()`, `object[['pca']]` will contain the results of the PCA. By adding new elements to the list, users can add additional, and custom, dimensional reductions. Each stored dimensional reduction contains the following slots:

1. **cell.embeddings**: stores the coordinates for each cell in low-dimensional space.
2. **feature.loadings**: stores the weight for each feature along each dimension of the embedding
3. **feature.loadings.projected**: Seurat typically calculate the dimensional reduction on a subset of genes (for example, high-variance genes), and then project that structure onto the entire dataset (all genes). The results of that projection (calculated with `ProjectDim()`) are stored in this slot. Note that the cell loadings will remain unchanged after projection but there are now feature loadings for all feature
4. **stdev**: The standard deviations of each dimension. Most often used with PCA (storing the square roots of the eigenvalues of the covariance matrix) and can be useful when looking at the drop off in the amount of variance that is explained by each successive dimension.
5. **key**: Sets the column names for the cell.embeddings and feature.loadings matrices. For example, for PCA, the column names are PC1, PC2, etc., so the key is “PC”.
6. **jackstraw**: Stores the results of the jackstraw procedure run using this dimensional reduction technique. Currently supported only for PCA.
7. **misc**: Bonus slot to store any other information you might want

To access these slots, we provide the `Embeddings()`,`Loadings()`, and `Stdev()` functions

```
pbmc[["pca"]]
```

```
## A dimensional reduction object with key PC_ 
##  Number of dimensions: 50 
##  Number of cells: 2638 
##  Projected dimensional reduction calculated:  FALSE 
##  Jackstraw run: TRUE 
##  Computed using assay: RNA
```

```
head(Embeddings(pbmc, reduction = "pca")[, 1:5])
```

```
##                      PC_1       PC_2       PC_3       PC_4        PC_5
## AAACATACAACCAC -4.7296855 -0.5184265 -0.7623220 -2.3156790 -0.07160006
## AAACATTGAGCTAC -0.5174029  4.5918957  5.9091921  6.9118856 -1.96243034
## AAACATTGATCAGC -3.1891063 -3.4695154 -0.8313710 -2.0019985 -5.10442765
## AAACCGTGCTTCCG 12.7933021  0.1007166  0.6310221 -0.3687338  0.21838204
## AAACCGTGTATGCG -3.1288078 -6.3481412  1.2507776  3.0191026  7.84739502
## AAACGCACTGGTAC -3.1088963  0.9262125 -0.6482331 -2.3244378 -2.00526763
```

```
head(Loadings(pbmc, reduction = "pca")[, 1:5])
```

```
##                PC_1        PC_2        PC_3        PC_4        PC_5
## PPBP    0.010990202  0.01148426 -0.15176092  0.10403737 0.003299077
## LYZ     0.116231706  0.01472515 -0.01280613 -0.04414540 0.049906881
## S100A9  0.115414362  0.01895146 -0.02368853 -0.05787777 0.085382309
## IGLL5  -0.007987473  0.05454239  0.04901533  0.06694722 0.004603231
## GNLY   -0.015238762 -0.13375626  0.04101340  0.06912322 0.104558611
## FTL     0.118292572  0.01871142 -0.00984755 -0.01555269 0.038743505
```

```
head(Stdev(pbmc, reduction = "pca"))
```

```
## [1] 7.098420 4.495493 3.872592 3.748859 3.171755 2.545292
```

Seurat provides `RunPCA()` (pca),  and `RunTSNE()` (tsne), and  representing dimensional reduction techniques commonly applied to scRNA-seq data. When using these functions, all slots are filled automatically.

We also allow users to add the results of a custom dimensional reduction technique (for example, multi-dimensional scaling (MDS), or [zero-inflated factor analysis](https://github.com/epierson9/ZIFA)), that is computed separately. All you need is a matrix with each cell’s coordinates in low-dimensional space, as shown below.

## Storing a custom dimensional reduction calculation

Though not incorporated as part of the Seurat package, its easy to run multidimensional scaling (MDS) in R. If you were interested in running MDS and storing the output in your Seurat object:

```
# Before running MDS, we first calculate a distance matrix between all pairs of cells.  Here
# we use a simple euclidean distance metric on all genes, using scale.data as input
d <- dist(t(GetAssayData(pbmc, slot = "scale.data")))
# Run the MDS procedure, k determines the number of dimensions
mds <- cmdscale(d = d, k = 2)
# cmdscale returns the cell embeddings, we first label the columns to ensure downstream
# consistency
colnames(mds) <- paste0("MDS_", 1:2)
# We will now store this as a custom dimensional reduction called 'mds'
pbmc[["mds"]] <- CreateDimReducObject(embeddings = mds, key = "MDS_", assay = DefaultAssay(pbmc))

# We can now use this as you would any other dimensional reduction in all downstream functions
DimPlot(pbmc, reduction = "mds", pt.size = 0.5)
```

![](dim_reduction_vignette_files/figure-html/mds-1.png)

```
# If you wold like to observe genes that are strongly correlated with the first MDS coordinate
pbmc <- ProjectDim(pbmc, reduction = "mds")
```

```
## MDS_ 1 
## Positive:  MALAT1, RPS27A, RPS27, RPL3, RPL23A, RPL21, RPL13A, RPS6, RPS3A, RPS3 
##     RPL9, LTB, RPSA, CD3D, RPS25, RPS18, PTPRCAP, RPS12, RPL30, RPL31 
## Negative:  CST3, TYROBP, FCER1G, LST1, FTL, AIF1, FTH1, TYMP, FCN1, LYZ 
##     LGALS1, S100A9, CFD, CD68, SERPINA1, CTSS, IFITM3, SPI1, S100A8, LGALS2 
## MDS_ 2 
## Positive:  NKG7, PRF1, CST7, GZMA, GZMB, B2M, FGFBP2, CTSW, GNLY, HLA-C 
##     GZMH, SPON2, CD247, FCGR3A, CCL5, HLA-A, CCL4, GZMM, KLRD1, CLIC3 
## Negative:  RPL32, RPL18A, HLA-DRA, CD79A, RPL13, MS4A1, RPL11, TCL1A, RPS9, RPL12 
##     LINC00926, HLA-DQB1, HLA-DQA1, HLA-DRB1, RPL28, RPS2, S100A8, HLA-DMA, RPL8, RPLP1
```

```
# Display the results as a heatmap
DimHeatmap(pbmc, reduction = "mds", dims = 1, cells = 500, projected = TRUE, balanced = TRUE)
```

![](dim_reduction_vignette_files/figure-html/mds-2.png)

```
# Explore how the first MDS dimension is distributed across clusters
VlnPlot(pbmc, features = "MDS_1")
```

![](dim_reduction_vignette_files/figure-html/mds-3.png)

```
# See how the first MDS dimension is correlated with the first PC dimension
FeatureScatter(pbmc, feature1 = "MDS_1", feature2 = "PC_1")
```

![](dim_reduction_vignette_files/figure-html/mds-4.png)

```
library(ggplot2)
plot <- DimPlot(pbmc, reduction = "mds", pt.size = 0.5)
ggsave(filename = "../output/images/pbmc_mds.jpg", height = 7, width = 12, plot = plot, quality = 50)
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
##  [1] ggplot2_3.4.3                 thp1.eccite.SeuratData_3.1.5 
##  [3] pbmcsca.SeuratData_3.0.0      pbmcref.SeuratData_1.0.0     
##  [5] pbmcMultiome.SeuratData_0.1.4 pbmc3k.SeuratData_3.1.4      
##  [7] panc8.SeuratData_3.0.2        ifnb.SeuratData_3.1.0        
##  [9] hcabm40k.SeuratData_3.0.0     bmcite.SeuratData_0.3.0      
## [11] SeuratData_0.2.2.9001         Seurat_4.9.9.9059            
## [13] SeuratObject_4.9.9.9091       sp_2.0-0                     
## 
## loaded via a namespace (and not attached):
##   [1] spam_2.9-1             systemfonts_1.0.4      plyr_1.8.8            
##   [4] igraph_1.5.1           lazyeval_0.2.2         splines_4.2.2         
##   [7] RcppHNSW_0.4.1         listenv_0.9.0          scattermore_1.2       
##  [10] digest_0.6.33          htmltools_0.5.6        fansi_1.0.4           
##  [13] magrittr_2.0.3         memoise_2.0.1          tensor_1.5            
##  [16] cluster_2.1.4          ROCR_1.0-11            globals_0.16.2        
##  [19] matrixStats_1.0.0      pkgdown_2.0.7          spatstat.sparse_3.0-2 
##  [22] colorspace_2.1-0       rappdirs_0.3.3         ggrepel_0.9.3         
##  [25] textshaping_0.3.6      xfun_0.40              dplyr_1.1.3           
##  [28] crayon_1.5.2           jsonlite_1.8.7         progressr_0.14.0      
##  [31] spatstat.data_3.0-1    survival_3.5-5         zoo_1.8-12            
##  [34] glue_1.6.2             polyclip_1.10-4        gtable_0.3.4          
##  [37] leiden_0.4.3           future.apply_1.11.0    abind_1.4-5           
##  [40] scales_1.2.1           spatstat.random_3.1-5  miniUI_0.1.1.1        
##  [43] Rcpp_1.0.11            viridisLite_0.4.2      xtable_1.8-4          
##  [46] reticulate_1.31        dotCall64_1.0-2        htmlwidgets_1.6.2     
##  [49] httr_1.4.7             RColorBrewer_1.1-3     ellipsis_0.3.2        
##  [52] ica_1.0-3              farver_2.1.1           pkgconfig_2.0.3       
##  [55] sass_0.4.7             uwot_0.1.16            deldir_1.0-9          
##  [58] utf8_1.2.3             labeling_0.4.3         tidyselect_1.2.0      
##  [61] rlang_1.1.1            reshape2_1.4.4         later_1.3.1           
##  [64] munsell_0.5.0          tools_4.2.2            cachem_1.0.8          
##  [67] cli_3.6.1              generics_0.1.3         ggridges_0.5.4        
##  [70] evaluate_0.21          stringr_1.5.0          fastmap_1.1.1         
##  [73] yaml_2.3.7             ragg_1.2.5             goftest_1.2-3         
##  [76] knitr_1.43             fs_1.6.3               fitdistrplus_1.1-11   
##  [79] purrr_1.0.2            RANN_2.6.1             pbapply_1.7-2         
##  [82] future_1.33.0          nlme_3.1-162           mime_0.12             
##  [85] formatR_1.14           ggrastr_1.0.1          compiler_4.2.2        
##  [88] beeswarm_0.4.0         plotly_4.10.2          png_0.1-8             
##  [91] spatstat.utils_3.0-3   tibble_3.2.1           bslib_0.5.1           
##  [94] stringi_1.7.12         highr_0.10             desc_1.4.2            
##  [97] RSpectra_0.16-1        lattice_0.21-8         Matrix_1.5-3          
## [100] vctrs_0.6.3            pillar_1.9.0           lifecycle_1.0.3       
## [103] spatstat.geom_3.2-5    lmtest_0.9-40          jquerylib_0.1.4       
## [106] RcppAnnoy_0.0.21       data.table_1.14.8      cowplot_1.1.1         
## [109] irlba_2.3.5.1          httpuv_1.6.11          patchwork_1.1.3       
## [112] R6_2.5.1               promises_1.2.1         KernSmooth_2.23-22    
## [115] gridExtra_2.3          vipor_0.4.5            parallelly_1.36.0     
## [118] codetools_0.2-19       fastDummies_1.7.3      MASS_7.3-58.2         
## [121] rprojroot_2.0.3        withr_2.5.0            sctransform_0.3.5     
## [124] parallel_4.2.2         grid_4.2.2             tidyr_1.3.0           
## [127] rmarkdown_2.24         Rtsne_0.16             spatstat.explore_3.2-1
## [130] shiny_1.7.5            ggbeeswarm_0.7.1
```