---
source_url: https://satijalab.org/seurat/articles/parsebio_sketch_integration
download_date: 2025-11-20
---

# Sketch integration using a 1 million cell dataset from Parse Biosciences

#### Compiled: 2023-11-15

Source: [`vignettes/ParseBio_sketch_integration.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/ParseBio_sketch_integration.Rmd)

`ParseBio_sketch_integration.Rmd`

The recent increase in publicly available single-cell datasets poses a significant challenge for integrative analysis. For example, multiple tissues have now been profiled across dozens of studies, representing hundreds of individuals and millions of cells. In [Hao et al, 2023](https://www.nature.com/articles/s41587-023-01767-y) proposed a dictionary learning based method, atomic sketch integration, could also enable efficient and large-scale integrative analysis. Our procedure enables the integration of large compendiums of datasets without ever needing to load the full scale of data into memory. In [our manuscript](https://www.nature.com/articles/s41587-023-01767-y) we use atomic sketch integration to integrate millions of scRNA-seq from human lung and human PBMC.

In this vignette, we demonstrate how to use atomic sketch integration to harmonize scRNA-seq experiments 1M cells, though we have used this procedure to integrate datasets of 10M+ cells as well. We analyze a dataset from Parse Biosciences, in which PBMC from 24 human samples (12 healthy donors, 12 Type-1 diabetes donors), which is available [here](https://cdn.parsebiosciences.com/1M_PBMC_T1D_Parse.zip).

* Sample a representative subset of cells (‘atoms’) from each dataset
* Integrate the atoms from each dataset, and define a set of cell states
* Reconstruct (integrate) the full datasets, based on the atoms
* Annotate all cells in the full datasets
* Identify cell-type specific differences between healthy and diabetic patients

Prior to running this vignette, please [install Seurat v5](/seurat/articles/install), as well as the [BPCells](https://github.com/bnprks/BPCells) package, which we use for on-disk storage. You can read more about using BPCells in Seurat v5 [here](/seurat/articles/seurat5_bpcells_interaction_vignette). We also recommend reading the [Sketch-based analysis in Seurat v5](/seurat/articles/seurat5_sketch_analysis) vignette, which introduces the concept of on-disk and in-memory storage in Seurat v5.

```
library(Seurat)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)
```

## Create a Seurat object containing data from 24 patients

We downloaded the original dataset and donor metadata from [Parse Biosciences](https://cdn.parsebiosciences.com/1M_PBMC_T1D_Parse.zip). While the BPCells package can work directly with h5ad files, for optimal performance, we converted the dataset to the compressed sparse format used by BPCells, as described [here](/seurat/articles/seurat5_bpcells_interaction_vignette). We create a Seurat object for this dataset. Since the input to `CreateSeuratObject` is a BPCells matrix, the data remains on-disk and is not loaded into memory. After creating the object, we split the dataset into 24 [layers](/seurat/articles/seurat5_essential_commands), one for each sample (i.e. patient), to facilitate integration.

```
parse.mat <- open_matrix_dir(dir = "/brahms/hartmana/vignette_data/bpcells/parse_1m_pbmc")
# need to move
metadata <- readRDS("/brahms/haoy/vignette_data/ParseBio_PBMC_meta.rds")
object <- CreateSeuratObject(counts = parse.mat, meta.data = metadata)

object <- NormalizeData(object)
# split assay into 24 layers
object[["RNA"]] <- split(object[["RNA"]], f = object$sample)
object <- FindVariableFeatures(object, verbose = FALSE)
```

## Sample representative cells from each dataset

Inspired by pioneering work aiming to identify [‘sketches’](https://www.sciencedirect.com/science/article/pii/S2405471219301528) of scRNA-seq data, our first step is to sample a representative set of cells from each dataset. We compute a leverage score (estimate of [‘statistical leverage’](https://arxiv.org/abs/1109.3843)) for each cell, which helps to identify cells that are likely to be member of rare subpopulations and ensure that these are included in our representative sample. Importantly, the estimation of leverage scores only requires data normalization, can be computed efficiently for sparse datasets, and does not require any intensive computation or dimensional reduction steps. We load each object separately, perform basic preprocessing (normalization and variable feature selection), and select and store 5,000 representative cells from each dataset. Since there are 24 datasets, the sketched dataset now contains 120,000 cells. These cells are stored in a new `sketch` assay, and are loaded in-memory.

```
object <- SketchData(object = object, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
object
```

```
## An object of class Seurat 
## 121328 features across 966344 samples within 2 assays 
## Active assay: sketch (60664 features, 2000 variable features)
##  48 layers present: counts.H_3128, counts.H_3060, counts.H_409, counts.H_7108, counts.H_120, counts.D_504, counts.H_2928, counts.H_727, counts.H_4579, counts.H_6625, counts.H_630, counts.D_500, counts.D_501, counts.D_610, counts.H_4119, counts.D_609, counts.D_644, counts.D_498, counts.H_1334, counts.D_497, counts.D_505, counts.D_506, counts.D_502, counts.D_503, data.H_3128, data.H_3060, data.H_409, data.H_7108, data.H_120, data.D_504, data.H_2928, data.H_727, data.H_4579, data.H_6625, data.H_630, data.D_500, data.D_501, data.D_610, data.H_4119, data.D_609, data.D_644, data.D_498, data.H_1334, data.D_497, data.D_505, data.D_506, data.D_502, data.D_503
##  1 other assay present: RNA
```

## Perform integration on the sketched cells across samples

Next we perform integrative analysis on the ‘atoms’ from each of the datasets. Here, we perform integration using the streamlined [Seurat v5 integration worfklow](/seurat/articles/seurat5_integration), and utilize the reference-based `RPCAIntegration` method. The function performs all corrections in low-dimensional space (rather than on the expression values themselves) to further improve speed and memory usage, and outputs a merged Seurat object where all cells have been placed in an integrated low-dimensional space (stored as `integrated.rpca`). However, we emphasize that you can perform integration here using any analysis technique that places cells across datasets into a shared space. This includes CCA Integration, Harmony, and scVI. We demonstrate how to use these tools in Seurat v5 [here](/seurat/articles/seurat5_integration).

```
DefaultAssay(object) <- "sketch"
object <- FindVariableFeatures(object, verbose = F)
object <- ScaleData(object, verbose = F)
object <- RunPCA(object, verbose = F)
# integrate the datasets
object <- IntegrateLayers(object, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca",
    dims = 1:30, k.anchor = 20, reference = which(Layers(object, search = "data") %in% c("data.H_3060")),
    verbose = F)
# cluster the integrated data
object <- FindNeighbors(object, reduction = "integrated.rpca", dims = 1:30)
object <- FindClusters(object, resolution = 2)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 120000
## Number of edges: 5137782
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8765
## Number of communities: 47
## Elapsed time: 83 seconds
```

```
object <- RunUMAP(object, reduction = "integrated.rpca", dims = 1:30, return.model = T, verbose = F)
```

```
# you can now rejoin the layers in the sketched assay this is required to perform differential
# expression
object[["sketch"]] <- JoinLayers(object[["sketch"]])
c10_markers <- FindMarkers(object = object, ident.1 = 10, max.cells.per.ident = 500, only.pos = TRUE)
head(c10_markers)
```

```
##                  p_val avg_log2FC pct.1 pct.2     p_val_adj
## BANK1    8.652285e-142   3.183046 0.980 0.183 5.248822e-137
## RALGPS2  1.870876e-134   3.168787 0.949 0.215 1.134948e-129
## OSBPL10  4.406924e-133   3.556864 0.890 0.118 2.673417e-128
## ARHGAP24 1.666152e-121   2.957065 0.916 0.263 1.010754e-116
## AFF3     6.864595e-107   2.248398 0.978 0.306 4.164338e-102
## EBF1     7.805461e-107   2.994191 0.841 0.126 4.735105e-102
```

```
# You can now annotate clusters using marker genes.  We performed this step, and include the
# results in the 'sketch.celltype' metadata column

plot.s1 <- DimPlot(object, group.by = "sample", reduction = "umap")
plot.s2 <- DimPlot(object, group.by = "celltype.manual", reduction = "umap")
```

```
plot.s1 + plot.s2 + plot_layout(ncol = 1)
```

![](ParseBio_sketch_integration_files/figure-html/unnamed-chunk-6-1.png)

## Integrate the full datasets

Now that we have integrated the subset of atoms of each dataset, placing them each in an integrated low-dimensional space, we can now place each cell from each dataset in this space as well. We load the full datasets back in individually, and use the `ProjectIntegration` function to integrate all cells. After this function is run, the `integrated.rpca.full` space now embeds all cells in the dataset.Even though all cells in the dataset have been integrated together, the non-sketched cells are not loaded into memory. Users can still switch between the `sketch` (sketched cells, in-memory) and `RNA` (full dataset, on disk) for analysis. After integration, we can also project cell type labels from the sketched cells onto the full dataset using `ProjectData`.

```
# resplit the sketched cell assay into layers this is required to project the integration onto
# all cells
object[["sketch"]] <- split(object[["sketch"]], f = object$sample)

object <- ProjectIntegration(object = object, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


object <- ProjectData(object = object, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
    full.reduction = "integrated.rpca.full", dims = 1:30, refdata = list(celltype.full = "celltype.manual"))
```

```
object <- RunUMAP(object, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full",
    reduction.key = "UMAP_full_")
```

```
p1 <- DimPlot(object, reduction = "umap.full", group.by = "sample", alpha = 0.1)
p2 <- DimPlot(object, reduction = "umap.full", group.by = "celltype.full", alpha = 0.1)
p1 + p2 + plot_layout(ncol = 1)
```

![](ParseBio_sketch_integration_files/figure-html/unnamed-chunk-9-1.png)

## Compare healthy and diabetic samples

By integrating all samples together, we can now compare healthy and diabetic cells in matched cell states. To maximize statistical power, we want to use all cells - not just the sketched cells - to perform this analysis. As recommended by [Soneson et all.](https://www.nature.com/articles/nmeth.4612) and [Crowell et al.](https://www.nature.com/articles/s41467-020-19894-4), we use an aggregation-based (pseudobulk) workflow. We aggregate all cells within the same cell type and sample using the `AggregateExpression` function. This returns a Seurat object where each ‘cell’ represents the pseudobulk profile of one cell type in one individual.

After we aggregate cells, we can perform celltype-specific differential expression between healthy and diabetic samples using DESeq2. We demonstrate this for CD14 monocytes.

```
bulk <- AggregateExpression(object, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("celltype.full",
    "sample", "disease"))
```

```
# each sample is an individual-specific celltype-specific pseudobulk profile
tail(Cells(bulk))
```

```
## [1] "Treg_H-4119_H" "Treg_H-4579_H" "Treg_H-630_H"  "Treg_H-6625_H"
## [5] "Treg_H-7108_H" "Treg_H-727_H"
```

```
  cd14.bulk <- subset(bulk, celltype.full == "CD14 Mono")
  Idents(cd14.bulk) <- "disease"
  de_markers <- FindMarkers(cd14.bulk, ident.1 = "D", ident.2 = "H", slot = "counts", test.use = "DESeq2",
      verbose = F)
  de_markers$gene <- rownames(de_markers)
  ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
      ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,
      "")), colour = "red", size = 3)
```

![](ParseBio_sketch_integration_files/figure-html/unnamed-chunk-12-1.png)

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
## [1] patchwork_1.1.3    ggrepel_0.9.4      ggplot2_3.4.4      dplyr_1.1.3       
## [5] BPCells_0.1.0      Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-1          
## 
## loaded via a namespace (and not attached):
##   [1] spam_2.10-0            systemfonts_1.0.4      plyr_1.8.9            
##   [4] igraph_1.5.1           lazyeval_0.2.2         splines_4.2.2         
##   [7] RcppHNSW_0.5.0         listenv_0.9.0          scattermore_1.2       
##  [10] GenomeInfoDb_1.34.9    digest_0.6.33          htmltools_0.5.6.1     
##  [13] fansi_1.0.5            magrittr_2.0.3         memoise_2.0.1         
##  [16] tensor_1.5             cluster_2.1.4          ROCR_1.0-11           
##  [19] limma_3.54.1           globals_0.16.2         matrixStats_1.0.0     
##  [22] pkgdown_2.0.7          spatstat.sparse_3.0-3  colorspace_2.1-0      
##  [25] textshaping_0.3.6      xfun_0.40              RCurl_1.98-1.12       
##  [28] jsonlite_1.8.7         progressr_0.14.0       spatstat.data_3.0-3   
##  [31] survival_3.5-7         zoo_1.8-12             glue_1.6.2            
##  [34] polyclip_1.10-6        gtable_0.3.4           zlibbioc_1.44.0       
##  [37] XVector_0.38.0         leiden_0.4.3           future.apply_1.11.0   
##  [40] BiocGenerics_0.44.0    abind_1.4-5            scales_1.2.1          
##  [43] spatstat.random_3.2-1  miniUI_0.1.1.1         Rcpp_1.0.11           
##  [46] viridisLite_0.4.2      xtable_1.8-4           reticulate_1.34.0     
##  [49] dotCall64_1.1-0        stats4_4.2.2           htmlwidgets_1.6.2     
##  [52] httr_1.4.7             RColorBrewer_1.1-3     ellipsis_0.3.2        
##  [55] ica_1.0-3              farver_2.1.1           pkgconfig_2.0.3       
##  [58] sass_0.4.7             uwot_0.1.16            deldir_1.0-9          
##  [61] utf8_1.2.4             labeling_0.4.3         tidyselect_1.2.0      
##  [64] rlang_1.1.1            reshape2_1.4.4         later_1.3.1           
##  [67] munsell_0.5.0          tools_4.2.2            cachem_1.0.8          
##  [70] cli_3.6.1              generics_0.1.3         ggridges_0.5.4        
##  [73] evaluate_0.22          stringr_1.5.0          fastmap_1.1.1         
##  [76] yaml_2.3.7             ragg_1.2.5             goftest_1.2-3         
##  [79] knitr_1.45             fs_1.6.3               fitdistrplus_1.1-11   
##  [82] purrr_1.0.2            RANN_2.6.1             pbapply_1.7-2         
##  [85] future_1.33.0          nlme_3.1-162           mime_0.12             
##  [88] formatR_1.14           compiler_4.2.2         plotly_4.10.3         
##  [91] png_0.1-8              spatstat.utils_3.0-4   tibble_3.2.1          
##  [94] bslib_0.5.1            stringi_1.7.12         highr_0.10            
##  [97] desc_1.4.2             RSpectra_0.16-1        lattice_0.21-9        
## [100] Matrix_1.6-1.1         vctrs_0.6.4            pillar_1.9.0          
## [103] lifecycle_1.0.4        spatstat.geom_3.2-7    lmtest_0.9-40         
## [106] jquerylib_0.1.4        RcppAnnoy_0.0.21       data.table_1.14.8     
## [109] cowplot_1.1.1          bitops_1.0-7           irlba_2.3.5.1         
## [112] httpuv_1.6.12          GenomicRanges_1.50.2   R6_2.5.1              
## [115] promises_1.2.1         KernSmooth_2.23-22     gridExtra_2.3         
## [118] IRanges_2.32.0         parallelly_1.36.0      codetools_0.2-19      
## [121] fastDummies_1.7.3      MASS_7.3-58.2          rprojroot_2.0.4       
## [124] withr_2.5.2            presto_1.0.0           sctransform_0.4.1     
## [127] S4Vectors_0.36.2       GenomeInfoDbData_1.2.9 parallel_4.2.2        
## [130] grid_4.2.2             tidyr_1.3.0            rmarkdown_2.25        
## [133] Rtsne_0.16             spatstat.explore_3.2-5 shiny_1.7.5.1
```
