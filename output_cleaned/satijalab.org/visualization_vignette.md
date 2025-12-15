---
source_url: https://satijalab.org/seurat/articles/visualization_vignette
download_date: 2025-11-20
---

# Data visualization methods in Seurat

#### Compiled: 2023-10-31

Source: [`vignettes/visualization_vignette.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/visualization_vignette.Rmd)

`visualization_vignette.Rmd`

We’ll demonstrate visualization techniques in Seurat using our previously computed Seurat object from the 2,700 PBMC tutorial. You can download this dataset from [SeuratData](https://github.com/satijalab/seurat-data)

```
SeuratData::InstallData("pbmc3k")
```

```
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
pbmc3k.final <- LoadData("pbmc3k", type = "pbmc3k.final")
pbmc3k.final$groups <- sample(c("group1", "group2"), size = ncol(pbmc3k.final), replace = TRUE)
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
pbmc3k.final
```

```
## An object of class Seurat 
## 13714 features across 2638 samples within 1 assay 
## Active assay: RNA (13714 features, 2000 variable features)
##  3 layers present: data, counts, scale.data
##  2 dimensional reductions calculated: pca, umap
```

## Five visualizations of marker feature expression

```
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(pbmc3k.final, features = features, ncol = 2)
```

![](visualization_vignette_files/figure-html/visualization_smorgasbord-1.png)

```
# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(pbmc3k.final, features = features)
```

![](visualization_vignette_files/figure-html/visualization_smorgasbord-2.png)

```
# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(pbmc3k.final, features = features)
```

![](visualization_vignette_files/figure-html/visualization_smorgasbord-3.png)

```
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(pbmc3k.final, features = features) + RotatedAxis()
```

![](visualization_vignette_files/figure-html/visualization_smorgasbord2-1.png)

```
# Single cell heatmap of feature expression
DoHeatmap(subset(pbmc3k.final, downsample = 100), features = features, size = 3)
```

![](visualization_vignette_files/figure-html/visualization_smorgasbord2-2.png)

## New additions to `FeaturePlot`

```
# Plot a legend to map colors to expression levels
FeaturePlot(pbmc3k.final, features = "MS4A1")
```

![](visualization_vignette_files/figure-html/featureplot-1.png)

```
# Adjust the contrast in the plot
FeaturePlot(pbmc3k.final, features = "MS4A1", min.cutoff = 1, max.cutoff = 3)
```

![](visualization_vignette_files/figure-html/featureplot-2.png)

```
# Calculate feature-specific contrast levels based on quantiles of non-zero expression.
# Particularly useful when plotting multiple markers
FeaturePlot(pbmc3k.final, features = c("MS4A1", "PTPRCAP"), min.cutoff = "q10", max.cutoff = "q90")
```

![](visualization_vignette_files/figure-html/featureplot2-1.png)

```
# Visualize co-expression of two features simultaneously
FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), blend = TRUE)
```

![](visualization_vignette_files/figure-html/featureplot2-2.png)

```
# Split visualization to view expression by groups (replaces FeatureHeatmap)
FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), split.by = "groups")
```

![](visualization_vignette_files/figure-html/featureplot.split-1.png)

## Updated and expanded visualization functions

In addition to changes to `FeaturePlot()`, several other plotting functions have been updated and expanded with new features and taking over the role of now-deprecated functions

```
# Violin plots can also be split on some variable. Simply add the splitting variable to object
# metadata and pass it to the split.by argument
VlnPlot(pbmc3k.final, features = "percent.mt", split.by = "groups")
```

![](visualization_vignette_files/figure-html/new_functions-1.png)

```
# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(pbmc3k.final, features = features, split.by = "groups") + RotatedAxis()
```

![](visualization_vignette_files/figure-html/new_functions-2.png)

```
# DimPlot replaces TSNEPlot, PCAPlot, etc. In addition, it will plot either 'umap', 'tsne', or
# 'pca' by default, in that order
DimPlot(pbmc3k.final)
```

![](visualization_vignette_files/figure-html/new_functions-3.png)

```
pbmc3k.final.no.umap <- pbmc3k.final
pbmc3k.final.no.umap[["umap"]] <- NULL
DimPlot(pbmc3k.final.no.umap) + RotatedAxis()
```

![](visualization_vignette_files/figure-html/new_functions-4.png)

```
# DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. This can
# be changed with the `group.by` parameter
DoHeatmap(pbmc3k.final, features = VariableFeatures(pbmc3k.final)[1:100], cells = 1:500, size = 4,
    angle = 90) + NoLegend()
```

![](visualization_vignette_files/figure-html/new2-1.png)

## Applying themes to plots

With Seurat, all plotting functions return ggplot2-based plots by default, allowing one to easily capture and manipulate plots just like any other ggplot2-based plot.

```
baseplot <- DimPlot(pbmc3k.final, reduction = "umap")
# Add custom labels and titles
baseplot + labs(title = "Clustering of 2,700 PBMCs")
```

![](visualization_vignette_files/figure-html/themeing-1.png)

```
# Use community-created themes, overwriting the default Seurat-applied theme Install ggmin
# with remotes::install_github('sjessa/ggmin')
baseplot + ggmin::theme_powerpoint()
```

![](visualization_vignette_files/figure-html/themeing-2.png)

```
# Seurat also provides several built-in themes, such as DarkTheme; for more details see
# ?SeuratTheme
baseplot + DarkTheme()
```

![](visualization_vignette_files/figure-html/themeing-3.png)

```
# Chain themes together
baseplot + FontSize(x.title = 20, y.title = 20) + NoLegend()
```

![](visualization_vignette_files/figure-html/themeing-4.png)

## Interactive plotting features

Seurat utilizes R’s plotly graphing library to create interactive plots. This interactive plotting feature works with any ggplot2-based scatter plots (requires a `geom_point` layer). To use, simply make a ggplot2-based scatter plot (such as `DimPlot()` or `FeaturePlot()`) and pass the resulting plot to `HoverLocator()`

```
# Include additional data to display alongside cell names by passing in a data frame of
# information.  Works well when using FetchData
plot <- FeaturePlot(pbmc3k.final, features = "MS4A1")
HoverLocator(plot = plot, information = FetchData(pbmc3k.final, vars = c("ident", "PC_1", "nFeature_RNA")))
```

Another interactive feature provided by Seurat is being able to manually select cells for further investigation. We have found this particularly useful for small clusters that do not always separate using unbiased clustering, but which look tantalizingly distinct. You can now select these cells by creating a ggplot2-based scatter plot (such as with `DimPlot()` or `FeaturePlot()`, and passing the returned plot to `CellSelector()`. `CellSelector()` will return a vector with the names of the points selected, so that you can then set them to a new identity class and perform differential expression.

For example, let’s pretend that DCs had merged with monocytes in the clustering, but we wanted to see what was unique about them based on their position in the tSNE plot.

```
pbmc3k.final <- RenameIdents(pbmc3k.final, DC = "CD14+ Mono")
plot <- DimPlot(pbmc3k.final, reduction = "umap")
select.cells <- CellSelector(plot = plot)
```

![](assets/pbmc_select.gif)

We can then change the identity of these cells to turn them into their own mini-cluster.

```
head(select.cells)
```

```
## [1] "AAGATTACCGCCTT" "AAGCCATGAACTGC" "AATTACGAATTCCT" "ACCCGTTGCTTCTA"
## [5] "ACGAGGGACAGGAG" "ACGTGATGCCATGA"
```

```
Idents(pbmc3k.final, cells = select.cells) <- "NewCells"

# Now, we find markers that are specific to the new cells, and find clear DC markers
newcells.markers <- FindMarkers(pbmc3k.final, ident.1 = "NewCells", ident.2 = "CD14+ Mono", min.diff.pct = 0.3,
    only.pos = TRUE)
head(newcells.markers)
```

```
##                 p_val avg_log2FC pct.1 pct.2    p_val_adj
## FCER1A   3.239004e-69   6.504163 0.800 0.017 4.441970e-65
## SERPINF1 7.761413e-36   5.456560 0.457 0.013 1.064400e-31
## HLA-DQB2 1.721094e-34   4.752397 0.429 0.010 2.360309e-30
## CD1C     2.304106e-33   4.929607 0.514 0.025 3.159851e-29
## ENHO     5.099765e-32   5.076634 0.400 0.010 6.993818e-28
## ITM2C    4.299994e-29   5.660234 0.371 0.010 5.897012e-25
```

Using `CellSelector` to Automatically Assign Cell Identities

In addition to returning a vector of cell names, `CellSelector()` can also take the selected cells and assign a new identity to them, returning a Seurat object with the identity classes already set. This is done by passing the Seurat object used to make the plot into `CellSelector()`, as well as an identity class. As an example, we’re going to select the same set of cells as before, and set their identity class to “selected”

```
pbmc3k.final <- CellSelector(plot = plot, object = pbmc3k.final, ident = "selected")
```

![](assets/pbmc_select.gif)

```
levels(pbmc3k.final)
```

```
## [1] "selected"     "Naive CD4 T"  "Memory CD4 T" "CD14+ Mono"   "B"           
## [6] "CD8 T"        "FCGR3A+ Mono" "NK"           "Platelet"
```

## Plotting Accessories

Along with new functions add interactive functionality to plots, Seurat provides new accessory functions for manipulating and combining plots.

```
# LabelClusters and LabelPoints will label clusters (a coloring variable) or individual points
# on a ggplot2-based scatter plot
plot <- DimPlot(pbmc3k.final, reduction = "pca") + NoLegend()
LabelClusters(plot = plot, id = "ident")
```

![](visualization_vignette_files/figure-html/labelling-1.png)

```
# Both functions support `repel`, which will intelligently stagger labels and draw connecting
# lines from the labels to the points or clusters
LabelPoints(plot = plot, points = TopCells(object = pbmc3k.final[["pca"]]), repel = TRUE)
```

![](visualization_vignette_files/figure-html/labelling-2.png)

Plotting multiple plots was previously achieved with the `CombinePlot()` function. We are deprecating this functionality in favor of the [patchwork](https://patchwork.data-imaginist.com/) system. Below is a brief demonstration but please see the patchwork package website [here](https://patchwork.data-imaginist.com/) for more details and examples.

```
plot1 <- DimPlot(pbmc3k.final)
# Create scatter plot with the Pearson correlation value as the title
plot2 <- FeatureScatter(pbmc3k.final, feature1 = "LYZ", feature2 = "CCL5")
# Combine two plots
plot1 + plot2
```

![](visualization_vignette_files/figure-html/combining_plots-1.png)

```
# Remove the legend from all plots
(plot1 + plot2) & NoLegend()
```

![](visualization_vignette_files/figure-html/combining_plots-2.png)

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
##  [3] thp1.eccite.SeuratData_3.1.5  stxBrain.SeuratData_0.1.1    
##  [5] ssHippo.SeuratData_3.1.4      pbmcsca.SeuratData_3.0.0     
##  [7] pbmcref.SeuratData_1.0.0      pbmcMultiome.SeuratData_0.1.4
##  [9] pbmc3k.SeuratData_3.1.4       panc8.SeuratData_3.0.2       
## [11] ifnb.SeuratData_3.0.0         hcabm40k.SeuratData_3.0.0    
## [13] cbmc.SeuratData_3.1.4         bmcite.SeuratData_0.3.0      
## [15] SeuratData_0.2.2.9001         Seurat_5.0.0                 
## [17] SeuratObject_5.0.0            sp_2.1-1                     
## 
## loaded via a namespace (and not attached):
##   [1] spam_2.10-0            systemfonts_1.0.4      plyr_1.8.9            
##   [4] igraph_1.5.1           lazyeval_0.2.2         splines_4.2.2         
##   [7] RcppHNSW_0.5.0         crosstalk_1.2.0        listenv_0.9.0         
##  [10] scattermore_1.2        digest_0.6.33          htmltools_0.5.6.1     
##  [13] fansi_1.0.5            magrittr_2.0.3         memoise_2.0.1         
##  [16] tensor_1.5             cluster_2.1.4          ROCR_1.0-11           
##  [19] limma_3.54.1           globals_0.16.2         matrixStats_1.0.0     
##  [22] pkgdown_2.0.7          spatstat.sparse_3.0-3  colorspace_2.1-0      
##  [25] rappdirs_0.3.3         ggrepel_0.9.4          textshaping_0.3.6     
##  [28] xfun_0.40              dplyr_1.1.3            crayon_1.5.2          
##  [31] jsonlite_1.8.7         progressr_0.14.0       spatstat.data_3.0-3   
##  [34] survival_3.5-7         zoo_1.8-12             glue_1.6.2            
##  [37] polyclip_1.10-6        gtable_0.3.4           leiden_0.4.3          
##  [40] future.apply_1.11.0    abind_1.4-5            scales_1.2.1          
##  [43] spatstat.random_3.2-1  miniUI_0.1.1.1         Rcpp_1.0.11           
##  [46] viridisLite_0.4.2      xtable_1.8-4           reticulate_1.34.0     
##  [49] dotCall64_1.1-0        ggmin_0.0.0.9000       htmlwidgets_1.6.2     
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
##  [88] formatR_1.14           ggrastr_1.0.1          compiler_4.2.2        
##  [91] beeswarm_0.4.0         plotly_4.10.3          png_0.1-8             
##  [94] spatstat.utils_3.0-4   tibble_3.2.1           bslib_0.5.1           
##  [97] stringi_1.7.12         highr_0.10             desc_1.4.2            
## [100] RSpectra_0.16-1        lattice_0.21-9         Matrix_1.6-1.1        
## [103] vctrs_0.6.4            pillar_1.9.0           lifecycle_1.0.3       
## [106] spatstat.geom_3.2-7    lmtest_0.9-40          jquerylib_0.1.4       
## [109] RcppAnnoy_0.0.21       data.table_1.14.8      cowplot_1.1.1         
## [112] irlba_2.3.5.1          httpuv_1.6.12          R6_2.5.1              
## [115] promises_1.2.1         KernSmooth_2.23-22     gridExtra_2.3         
## [118] vipor_0.4.5            parallelly_1.36.0      codetools_0.2-19      
## [121] fastDummies_1.7.3      MASS_7.3-58.2          rprojroot_2.0.3       
## [124] withr_2.5.2            presto_1.0.0           sctransform_0.4.1     
## [127] parallel_4.2.2         grid_4.2.2             tidyr_1.3.0           
## [130] rmarkdown_2.25         Rtsne_0.16             spatstat.explore_3.2-5
## [133] shiny_1.7.5.1          ggbeeswarm_0.7.1
```
