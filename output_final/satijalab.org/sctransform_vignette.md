---
source_url: https://satijalab.org/seurat/articles/sctransform_vignette
download_date: 2025-11-20
---

# Using sctransform in Seurat

#### Saket Choudhary, Christoph Hafemeister & Rahul Satija

#### Compiled: 2023-10-31

Source: [`vignettes/sctransform_vignette.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/sctransform_vignette.Rmd)

`sctransform_vignette.Rmd`

Biological heterogeneity in single-cell RNA-seq data is often confounded by technical factors including sequencing depth. The number of molecules detected in each cell can vary significantly between cells, even within the same celltype. Interpretation of scRNA-seq data requires effective pre-processing and normalization to remove this technical variability.

In [our manuscript](https://genomebiology-biomedcentral-com/articles/10.1186/s13059-021-02584-9) we introduce a modeling framework for the normalization and variance stabilization of molecular count data from scRNA-seq experiments. This procedure omits the need for heuristic steps including pseudocount addition or log-transformation and improves common downstream analytical tasks such as variable gene selection, dimensional reduction, and differential expression. We named this method `sctransform`.

Inspired by important and rigorous work from [Lause et al](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02451-7), we released an [updated manuscript](https://link.springer.com/article/10.1186/s13059-021-02584-9) and updated the sctransform software to a v2 version, which is now the default in Seurat v5.

```
library(Seurat)
library(ggplot2)
library(sctransform)
```

## Load data and create Seurat object

```
pbmc_data <- Read10X(data.dir = "/brahms/shared/vignette-data/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data)
```

## Apply sctransform normalization

* Note that this single command replaces `NormalizeData()`, `ScaleData()`, and `FindVariableFeatures()`.
* Transformed data will be available in the SCT assay, which is set as the default after running sctransform
* During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage
* In Seurat v5, SCT v2 is applied by default. You can revert to v1 by setting `vst.flavor = 'v1'`
* The [glmGamPoi](https://bioconductor.org/packages/release/bioc/html/glmGamPoi.html) package substantially improves speed and is used by default if installed, with instructions [here](/seurat/articles/install)

```
# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
```

## Perform dimensionality reduction by PCA and UMAP embedding

```
# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
DimPlot(pbmc, label = TRUE)
```

![](sctransform_vignette_files/figure-html/pca-1.png)

**Why can we choose more PCs when using sctransform?**

In the [standard Seurat workflow](/seurat/articles/pbmc3k_tutorial) we focus on 10 PCs for this dataset, though we highlight that the results are similar with higher settings for this parameter. Interestingly, we’ve found that when using sctransform, we often benefit by pushing this parameter even higher. We believe this is because the sctransform workflow performs more effective normalization, strongly removing technical effects from the data.

Even after standard log-normalization, variation in sequencing depth is still a confounding factor (see [Figure 1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)), and this effect can subtly influence higher PCs. In sctransform, this effect is substantially mitigated (see [Figure 3](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)). This means that higher PCs are more likely to represent subtle, but biologically relevant, sources of heterogeneity – so including them may improve downstream analysis.

In addition, sctransform returns 3,000 variable features by default, instead of 2,000. The rationale is similar, the additional variable features are less likely to be driven by technical differences across cells, and instead may represent more subtle biological fluctuations. In general, we find that results produced with sctransform are less dependent on these parameters (indeed, we achieve nearly identical results when using all genes in the transcriptome, though this does reduce computational efficiency). This can help users generate more robust results, and in addition, enables the application of standard analysis pipelines with identical parameter settings that can quickly be applied to new datasets:

For example, the following code replicates the full end-to-end workflow, in a single command:

```
pbmc <- CreateSeuratObject(pbmc_data) %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    SCTransform(vars.to.regress = "percent.mt") %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    FindClusters()
```

**Where are normalized values stored for sctransform?**

The results of sctransfrom are stored in the “SCT” assay. You can learn more about multi-assay data and commands in Seurat in our [vignette](/seurat/articles/multimodal_vignette), [command cheat sheet](/seurat/articles/essential_commands#multi-assay-features), or [developer guide](https://github.com/satijalab/seurat/wiki/Assay).

* `pbmc[["SCT"]]$scale.data` contains the residuals (normalized values), and is used directly as input to PCA. Please note that this matrix is non-sparse, and can therefore take up a lot of memory if stored for all genes. To save memory, we store these values only for variable genes, by setting the return.only.var.genes = TRUE by default in the `SCTransform()` function call.
* To assist with visualization and interpretation, we also convert Pearson residuals back to ‘corrected’ UMI counts. You can interpret these as the UMI counts we would expect to observe if all cells were sequenced to the same depth. If you want to see exactly how we do this, please look at the correct function [here](https://github.com/ChristophH/sctransform/blob/master/R/denoise.R).
* The ‘corrected’ UMI counts are stored in `pbmc[["SCT"]]$counts`. We store log-normalized versions of these corrected counts in `pbmc[["SCT"]]$data`, which are very helpful for visualization.

---

Users can individually annotate clusters based on canonical markers. However, the sctransform normalization reveals sharper biological distinctions compared to the [standard Seurat workflow](/seurat/articles/pbmc3k_tutorial), in a few ways:

* Clear separation of at least 3 CD8 T cell populations (naive, memory, effector), based on CD8A, GZMK, CCL5, CCR7 expression
* Clear separation of three CD4 T cell populations (naive, memory, IFN-activated) based on S100A4, CCR7, IL32, and ISG15
* Additional developmental sub-structure in B cell cluster, based on TCL1A, FCER2
* Additional separation of NK cells into CD56dim vs. bright clusters, based on XCL1 and FCGR3A

```
# These are now standard steps in the Seurat workflow for visualization and clustering
# Visualize canonical marker genes as violin plots.
VlnPlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"),
    pt.size = 0.2, ncol = 4)
```

![](sctransform_vignette_files/figure-html/fplot-1.png)

```
# Visualize canonical marker genes on the sctransform embedding.
FeaturePlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), pt.size = 0.2,
    ncol = 3)
```

![](sctransform_vignette_files/figure-html/fplot-2.png)

```
FeaturePlot(pbmc, features = c("CD3D", "ISG15", "TCL1A", "FCER2", "XCL1", "FCGR3A"), pt.size = 0.2,
    ncol = 3)
```

![](sctransform_vignette_files/figure-html/fplot-3.png)

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
## [1] sctransform_0.4.1  ggplot2_3.4.4      Seurat_5.0.0       SeuratObject_5.0.0
## [5] sp_2.1-1          
## 
## loaded via a namespace (and not attached):
##   [1] spam_2.10-0                 systemfonts_1.0.4          
##   [3] plyr_1.8.9                  igraph_1.5.1               
##   [5] lazyeval_0.2.2              splines_4.2.2              
##   [7] RcppHNSW_0.5.0              listenv_0.9.0              
##   [9] scattermore_1.2             GenomeInfoDb_1.34.9        
##  [11] digest_0.6.33               htmltools_0.5.6.1          
##  [13] fansi_1.0.5                 magrittr_2.0.3             
##  [15] memoise_2.0.1               tensor_1.5                 
##  [17] cluster_2.1.4               ROCR_1.0-11                
##  [19] globals_0.16.2              matrixStats_1.0.0          
##  [21] R.utils_2.12.2              pkgdown_2.0.7              
##  [23] spatstat.sparse_3.0-3       colorspace_2.1-0           
##  [25] ggrepel_0.9.4               textshaping_0.3.6          
##  [27] xfun_0.40                   dplyr_1.1.3                
##  [29] RCurl_1.98-1.12             jsonlite_1.8.7             
##  [31] progressr_0.14.0            spatstat.data_3.0-3        
##  [33] survival_3.5-7              zoo_1.8-12                 
##  [35] glue_1.6.2                  polyclip_1.10-6            
##  [37] gtable_0.3.4                zlibbioc_1.44.0            
##  [39] XVector_0.38.0              leiden_0.4.3               
##  [41] DelayedArray_0.24.0         future.apply_1.11.0        
##  [43] BiocGenerics_0.44.0         abind_1.4-5                
##  [45] scales_1.2.1                spatstat.random_3.2-1      
##  [47] miniUI_0.1.1.1              Rcpp_1.0.11                
##  [49] viridisLite_0.4.2           xtable_1.8-4               
##  [51] reticulate_1.34.0           dotCall64_1.1-0            
##  [53] stats4_4.2.2                htmlwidgets_1.6.2          
##  [55] httr_1.4.7                  RColorBrewer_1.1-3         
##  [57] ellipsis_0.3.2              ica_1.0-3                  
##  [59] farver_2.1.1                pkgconfig_2.0.3            
##  [61] R.methodsS3_1.8.2           sass_0.4.7                 
##  [63] uwot_0.1.16                 deldir_1.0-9               
##  [65] utf8_1.2.4                  labeling_0.4.3             
##  [67] tidyselect_1.2.0            rlang_1.1.1                
##  [69] reshape2_1.4.4              later_1.3.1                
##  [71] munsell_0.5.0               tools_4.2.2                
##  [73] cachem_1.0.8                cli_3.6.1                  
##  [75] generics_0.1.3              ggridges_0.5.4             
##  [77] evaluate_0.22               stringr_1.5.0              
##  [79] fastmap_1.1.1               yaml_2.3.7                 
##  [81] ragg_1.2.5                  goftest_1.2-3              
##  [83] knitr_1.45                  fs_1.6.3                   
##  [85] fitdistrplus_1.1-11         purrr_1.0.2                
##  [87] RANN_2.6.1                  sparseMatrixStats_1.10.0   
##  [89] pbapply_1.7-2               future_1.33.0              
##  [91] nlme_3.1-162                mime_0.12                  
##  [93] formatR_1.14                ggrastr_1.0.1              
##  [95] R.oo_1.25.0                 compiler_4.2.2             
##  [97] beeswarm_0.4.0              plotly_4.10.3              
##  [99] png_0.1-8                   spatstat.utils_3.0-4       
## [101] tibble_3.2.1                bslib_0.5.1                
## [103] glmGamPoi_1.10.2            stringi_1.7.12             
## [105] highr_0.10                  desc_1.4.2                 
## [107] RSpectra_0.16-1             lattice_0.21-9             
## [109] Matrix_1.6-1.1              vctrs_0.6.4                
## [111] pillar_1.9.0                lifecycle_1.0.3            
## [113] spatstat.geom_3.2-7         lmtest_0.9-40              
## [115] jquerylib_0.1.4             RcppAnnoy_0.0.21           
## [117] bitops_1.0-7                data.table_1.14.8          
## [119] cowplot_1.1.1               irlba_2.3.5.1              
## [121] GenomicRanges_1.50.2        httpuv_1.6.12              
## [123] patchwork_1.1.3             R6_2.5.1                   
## [125] promises_1.2.1              KernSmooth_2.23-22         
## [127] gridExtra_2.3               vipor_0.4.5                
## [129] IRanges_2.32.0              parallelly_1.36.0          
## [131] codetools_0.2-19            fastDummies_1.7.3          
## [133] MASS_7.3-58.2               SummarizedExperiment_1.28.0
## [135] rprojroot_2.0.3             withr_2.5.2                
## [137] GenomeInfoDbData_1.2.9      S4Vectors_0.36.2           
## [139] parallel_4.2.2              grid_4.2.2                 
## [141] tidyr_1.3.0                 DelayedMatrixStats_1.20.0  
## [143] rmarkdown_2.25              MatrixGenerics_1.10.0      
## [145] Rtsne_0.16                  spatstat.explore_3.2-5     
## [147] Biobase_2.58.0              shiny_1.7.5.1              
## [149] ggbeeswarm_0.7.1
```
