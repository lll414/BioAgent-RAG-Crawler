---
source_url: https://satijalab.org/seurat/articles/seurat5_sketch_analysis
download_date: 2025-11-20
---

# Sketch-based analysis in Seurat v5

#### Compiled: 2023-10-31

Source: [`vignettes/seurat5_sketch_analysis.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/seurat5_sketch_analysis.Rmd)

`seurat5_sketch_analysis.Rmd`

## Intro: Sketch-based analysis in Seurat v5

As single-cell sequencing technologies continue to improve in scalability in throughput, the generation of datasets spanning a million or more cells is becoming increasingly routine. In Seurat v5, we introduce new infrastructure and methods to analyze, interpret, and explore these exciting datasets.

In this vignette, we introduce a sketch-based analysis workflow to analyze a 1.3 million cell dataset of the developing mouse brain, freely available from 10x Genomics. Analyzing datasets of this size with standard workflows can be challenging, slow, and memory-intensive. Here we introduce an alternative workflow that is highly scalable, even to datasets ranging beyond 10 million cells in size. Our ‘sketch-based’ workflow involves three new features in Seurat v5:

* Infrastructure for on-disk storage of large single-cell datasets

Storing expression matrices in memory can be challenging for extremely large scRNA-seq datasets. In Seurat v5, we introduce support for multiple on-disk storage formats.

* ‘Sketching’ methods to subsample cells from large datasets while preserving rare populations

As introduced in [Hie et al, 2019](https://www.sciencedirect.com/science/article/pii/S2405471219301528), cell sketching methods aim to compactly summarize large single-cell datasets in a small number of cells, while preserving the presence of both abundant and rare cell types. In Seurat v5, we leverage this idea to select subsamples (‘sketches’) of cells from large datasets that are stored on-disk. However, after sketching, the subsampled cells can be stored in-memory, allowing for interactive and rapid visualization and exploration. We store sketched cells (in-memory) and the full dataset (on-disk) as two assays in the same Seurat object. Users can then easily switch between the two versions, providing the flexibility to perform quick analyses on a subset of cells in-memory, while retaining access to the full dataset on-disk.

* Support for ‘bit-packing’ compression and infrastructure

We demonstrate the on-disk capabilities in Seurat v5 using the [BPCells package](https://github.com/bnprks/BPCells) developed by Ben Parks in the Greenleaf Lab. This package utilizes bit-packing compression and optimized, streaming-compatible C++ code to substantially improve I/O and computational performance when working with on-disk data. To run this vignette please install Seurat v5, using the installation instructions found [here](/seurat/articles/install). Additionally, you will need to install the `BPcells` package, using the installation instructions found [here](https://bnprks.github.io/BPCells/#installation).

```
library(Seurat)
library(BPCells)
library(ggplot2)
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)
```

## Create a Seurat object with a v5 assay for on-disk storage

We start by loading the 1.3M dataset from 10x Genomics using the `open_matrix_dir` function from `BPCells`. Note that in our [Introduction to on-disk storage vignette](/seurat/articles/seurat5_bpcells_interaction_vignette), we demonstrate how to create this on-disk representation. This function does not load the dataset into memory, but instead, creates a connection to the data stored on-disk. We then store this on-disk representation in the Seurat object, which is loaded using `readRDS` as per usual.

```
# Read the Seurat object, which contains 1.3M cells stored on-disk as part of the 'RNA' assay
obj <- readRDS("/brahms/hartmana/vignette_data/1p3_million_mouse_brain.rds")
obj
```

```
## An object of class Seurat 
## 27282 features across 1306127 samples within 1 assay 
## Active assay: RNA (27282 features, 0 variable features)
##  1 layer present: counts
```

```
# Note that since the data is stored on-disk, the object size easily fits in-memory (<1GB)
format(object.size(obj), units = "Mb")
```

```
## [1] "596.2 Mb"
```

## ‘Sketch’ a subset of cells, and load these into memory

We select a subset (‘sketch’) of 50,000 cells (out of 1.3M). Rather than sampling all cells with uniform probability, we compute and sample based off a ‘leverage score’ for each cell, which reflects the magnitude of its contribution to the gene-covariance matrix, and its importance to the overall dataset. In [Hao et al, 2022](https://www.biorxiv.org/content/10.1101/2022.02.24.481684v1.full), we demonstrate that the leverage score is highest for rare populations in a dataset. Therefore, our sketched set of 50,000 cells will oversample rare populations, retaining the biological complexity of the sample while drastically compressing the dataset.

The function `SketchData` takes a normalized single-cell dataset (stored either on-disk or in-memory), and a set of variable features. It returns a Seurat object with a new assay (`sketch`), consisting of 50,000 cells, but these cells are now stored in-memory. Users can now easily switch between the in-memory and on-disk representation just by changing the default assay.

```
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- SketchData(
  object = obj,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
obj
```

```
## An object of class Seurat 
## 54564 features across 1306127 samples within 2 assays 
## Active assay: sketch (27282 features, 2000 variable features)
##  2 layers present: counts, data
##  1 other assay present: RNA
```

```
# switch to analyzing the full dataset (on-disk)
DefaultAssay(obj) <- "RNA"
# switch to analyzing the sketched dataset (in-memory)
DefaultAssay(obj) <- "sketch"
```

## Perform clustering on the sketched dataset

Now that we have compressed the dataset, we can perform standard clustering and visualization of a 50,000 cell dataset. After clustering, we can see groups of cells that clearly correspond to precursors of distinct lineages, including endothelial cells (Igfbp7), Excitatory (Neurod6) and Inhibitory (Dlx2) neurons, Intermediate Progenitors (Eomes), Radial Glia (Vim), Cajal-Retzius cells (Reln), Oligodendroytes (Olig1), and extremely rare populations of macrophages (C1qa) that were oversampled in our sketched data.

```
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 2)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 50000
## Number of edges: 1992657
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8766
## Number of communities: 51
## Elapsed time: 19 seconds
```

```
obj <- RunUMAP(obj, dims = 1:50, return.model = T)
DimPlot(obj, label = T, label.size = 3, reduction = "umap") + NoLegend()
```

![](seurat5_sketch_analysis_files/figure-html/unnamed-chunk-3-1.png)

```
FeaturePlot(
  object = obj,
  features = c(
    "Igfbp7", "Neurod6", "Dlx2", "Gad2",
    "Eomes", "Vim", "Reln", "Olig1", "C1qa"
  ),
  ncol = 3
)
```

![](seurat5_sketch_analysis_files/figure-html/unnamed-chunk-4-1.png)

## Extend results to the full datasets

We can now extend the cluster labels and dimensional reductions learned on the sketched cells to the full dataset. The `ProjectData` function projects the on-disk data, onto the `sketch` assay. It returns a Seurat object that includes a

* Dimensional reduction (PCA): The `pca.full` dimensional reduction extends the `pca` reduction on the sketched cells to all cells in the dataset
* Dimensional reduction (UMAP): The `full.umap` dimensional reduction extends the `umap` reduction on the sketched cells to all cells in the dataset
* Cluster labels: The `cluster_full` column in the object metadata now labels all cells in the dataset with one of the cluster labels derived from the sketched cells

```
obj <- ProjectData(
  object = obj,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(cluster_full = "seurat_clusters")
)
# now that we have projected the full dataset, switch back to analyzing all cells
DefaultAssay(obj) <- "RNA"
```

```
DimPlot(obj, label = T, label.size = 3, reduction = "full.umap", group.by = "cluster_full", alpha = 0.1) + NoLegend()
```

![](seurat5_sketch_analysis_files/figure-html/unnamed-chunk-6-1.png)

```
# visualize gene expression on the sketched cells (fast) and the full dataset (slower)
DefaultAssay(obj) <- "sketch"
x1 <- FeaturePlot(obj, "C1qa")
DefaultAssay(obj) <- "RNA"
x2 <- FeaturePlot(obj, "C1qa")
x1 | x2
```

![](seurat5_sketch_analysis_files/figure-html/unnamed-chunk-7-1.png)

## Perform iterative sub-clustering

Now that we have performed an initial analysis of the dataset, we can iteratively ‘zoom-in’ on a cell subtype of interest, extract all cells of this type, and perform iterative sub-clustering. For example, we can see that Dlx2+ interneuron precursors are defined by clusters 2, 15, 18, 28 and 40.

```
DefaultAssay(obj) <- "sketch"
VlnPlot(obj, "Dlx2")
```

![](seurat5_sketch_analysis_files/figure-html/unnamed-chunk-8-1.png)

We therefore extract all cells from the full on-disk dataset that are present in these clusters. There are 200,892 of them. Since this is a manageable number, we can convert these data from on-disk storage into in-memory storage. We can then proceed with standard clustering.

```
# subset cells in these clusters. Note that the data remains on-disk after subsetting
obj.sub <- subset(obj, subset = cluster_full %in% c(2, 15, 18, 28, 40))
DefaultAssay(obj.sub) <- "RNA"

# now convert the RNA assay (previously on-disk) into an in-memory representation (sparse Matrix)
# we only convert the data layer, and keep the counts on-disk
obj.sub[["RNA"]]$data <- as(obj.sub[["RNA"]]$data, Class = "dgCMatrix")

# recluster the cells
obj.sub <- FindVariableFeatures(obj.sub)
obj.sub <- ScaleData(obj.sub)
obj.sub <- RunPCA(obj.sub)
obj.sub <- RunUMAP(obj.sub, dims = 1:30)
obj.sub <- FindNeighbors(obj.sub, dims = 1:30)
obj.sub <- FindClusters(obj.sub)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 236276
## Number of edges: 5470656
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8577
## Number of communities: 28
## Elapsed time: 191 seconds
```

```
DimPlot(obj.sub, label = T, label.size = 3) + NoLegend()
```

![](seurat5_sketch_analysis_files/figure-html/unnamed-chunk-10-1.png)

Note that we can start to see distinct interneuron lineages emerging in this dataset. We can see a clear separation of interneuron precursors that originated from the medial ganglionic eminence (Lhx6) or caudal ganglionic eminence (Nr2f2). We can further see the emergence of Sst (Sst) and Pvalb (Mef2c)-committed interneurons, and a CGE-derived Meis2-expressing progenitor population. These results closely mirror our findings from [Mayer*, Hafemeister*, Bandler\* et al, Nature 2018](https://www.nature.com/articles/nature25999), where we enriched for interneuron precursors using a Dlx6a-cre fate-mapping strategy. Here, we obtain similar results using only computational enrichment, enabled by the large size of the original dataset.

```
FeaturePlot(
  object = obj.sub,
  features = c(
    "Dlx2", "Gad2", "Lhx6", "Nr2f2", "Sst",
    "Mef2c", "Meis2", "Id2", "Dlx6os1"
  ),
  ncol = 3
)
```

![](seurat5_sketch_analysis_files/figure-html/unnamed-chunk-11-1.png)

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
## [1] ggplot2_3.4.4      BPCells_0.1.0      Seurat_5.0.0       SeuratObject_5.0.0
## [5] sp_2.1-1          
## 
## loaded via a namespace (and not attached):
##   [1] spam_2.10-0            systemfonts_1.0.4      plyr_1.8.9            
##   [4] igraph_1.5.1           lazyeval_0.2.2         splines_4.2.2         
##   [7] RcppHNSW_0.5.0         listenv_0.9.0          scattermore_1.2       
##  [10] GenomeInfoDb_1.34.9    digest_0.6.33          htmltools_0.5.6.1     
##  [13] fansi_1.0.5            magrittr_2.0.3         memoise_2.0.1         
##  [16] tensor_1.5             cluster_2.1.4          ROCR_1.0-11           
##  [19] globals_0.16.2         matrixStats_1.0.0      R.utils_2.12.2        
##  [22] pkgdown_2.0.7          spatstat.sparse_3.0-3  colorspace_2.1-0      
##  [25] ggrepel_0.9.4          textshaping_0.3.6      xfun_0.40             
##  [28] dplyr_1.1.3            RCurl_1.98-1.12        jsonlite_1.8.7        
##  [31] progressr_0.14.0       spatstat.data_3.0-3    survival_3.5-7        
##  [34] zoo_1.8-12             glue_1.6.2             polyclip_1.10-6       
##  [37] gtable_0.3.4           zlibbioc_1.44.0        XVector_0.38.0        
##  [40] leiden_0.4.3           R.cache_0.16.0         future.apply_1.11.0   
##  [43] BiocGenerics_0.44.0    abind_1.4-5            scales_1.2.1          
##  [46] spatstat.random_3.2-1  miniUI_0.1.1.1         Rcpp_1.0.11           
##  [49] viridisLite_0.4.2      xtable_1.8-4           reticulate_1.34.0     
##  [52] dotCall64_1.1-0        stats4_4.2.2           htmlwidgets_1.6.2     
##  [55] httr_1.4.7             RColorBrewer_1.1-3     ellipsis_0.3.2        
##  [58] ica_1.0-3              farver_2.1.1           pkgconfig_2.0.3       
##  [61] R.methodsS3_1.8.2      sass_0.4.7             uwot_0.1.16           
##  [64] deldir_1.0-9           utf8_1.2.4             labeling_0.4.3        
##  [67] tidyselect_1.2.0       rlang_1.1.1            reshape2_1.4.4        
##  [70] later_1.3.1            munsell_0.5.0          tools_4.2.2           
##  [73] cachem_1.0.8           cli_3.6.1              generics_0.1.3        
##  [76] ggridges_0.5.4         evaluate_0.22          stringr_1.5.0         
##  [79] fastmap_1.1.1          yaml_2.3.7             ragg_1.2.5            
##  [82] goftest_1.2-3          knitr_1.45             fs_1.6.3              
##  [85] fitdistrplus_1.1-11    purrr_1.0.2            RANN_2.6.1            
##  [88] pbapply_1.7-2          future_1.33.0          nlme_3.1-162          
##  [91] mime_0.12              ggrastr_1.0.1          R.oo_1.25.0           
##  [94] compiler_4.2.2         beeswarm_0.4.0         plotly_4.10.3         
##  [97] png_0.1-8              spatstat.utils_3.0-4   tibble_3.2.1          
## [100] bslib_0.5.1            stringi_1.7.12         highr_0.10            
## [103] desc_1.4.2             RSpectra_0.16-1        lattice_0.21-9        
## [106] Matrix_1.6-1.1         styler_1.10.2          vctrs_0.6.4           
## [109] pillar_1.9.0           lifecycle_1.0.3        spatstat.geom_3.2-7   
## [112] lmtest_0.9-40          jquerylib_0.1.4        RcppAnnoy_0.0.21      
## [115] bitops_1.0-7           data.table_1.14.8      cowplot_1.1.1         
## [118] irlba_2.3.5.1          httpuv_1.6.12          patchwork_1.1.3       
## [121] GenomicRanges_1.50.2   R6_2.5.1               promises_1.2.1        
## [124] KernSmooth_2.23-22     gridExtra_2.3          vipor_0.4.5           
## [127] IRanges_2.32.0         parallelly_1.36.0      codetools_0.2-19      
## [130] fastDummies_1.7.3      MASS_7.3-58.2          rprojroot_2.0.3       
## [133] withr_2.5.2            sctransform_0.4.1      GenomeInfoDbData_1.2.9
## [136] S4Vectors_0.36.2       parallel_4.2.2         grid_4.2.2            
## [139] tidyr_1.3.0            rmarkdown_2.25         Rtsne_0.16            
## [142] spatstat.explore_3.2-5 shiny_1.7.5.1          ggbeeswarm_0.7.1
```
