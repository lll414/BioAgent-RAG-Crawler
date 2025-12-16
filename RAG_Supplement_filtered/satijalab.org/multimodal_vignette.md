---
source_url: https://satijalab.org/seurat/articles/multimodal_vignette
download_date: 2025-11-20
---

# Using Seurat with multimodal data

#### Compiled: 2023-11-16

Source: [`vignettes/multimodal_vignette.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/multimodal_vignette.Rmd)

`multimodal_vignette.Rmd`

## Load in the data

The ability to make simultaneous measurements of multiple data types from the same cell, known as multimodal analysis, represents a new and exciting frontier for single-cell genomics. For example, [CITE-seq](http://www.nature.com/nmeth/journal/v14/n9/full/nmeth.4380.html) enables the simultaneous measurements of transcriptomes and cell-surface proteins from the same cell. Other exciting multimodal technologies, such as the [10x multiome kit](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression) allow for the paired measurements of cellular transcriptome and chromatin accessibility (i.e scRNA-seq+scATAC-seq). Other modalities that can be measured alongside cellular transcriptomes include genetic perturbations, cellular methylomes, and hashtag oligos from [Cell Hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1). We have designed Seurat to enable for the seamless storage, analysis, and exploration of diverse multimodal single-cell datasets.

In this vignette, we present an introductory workflow for creating a multimodal Seurat object and performing an initial analysis. For example, we demonstrate how to cluster a CITE-seq dataset on the basis of the measured cellular transcriptomes, and subsequently discover cell surface proteins that are enriched in each cluster. We note that Seurat also enables more advanced techniques for the analysis of multimodal data, in particular the application of our [Weighted Nearest Neighbors (WNN) approach](https://doi.org/10.1016/j.cell.2021.04.048) that enables simultaneous clustering of cells based on a weighted combination of both modalities, and you can explore this functionality [here](/seurat/articles/weighted_nearest_neighbor_analysis).

Here, we analyze a dataset of 8,617 cord blood mononuclear cells (CBMCs), where transcriptomic measurements are paired with abundance estimates for 11 surface proteins, whose levels are quantified with DNA-barcoded antibodies. First, we load in two count matrices : one for the RNA measurements, and one for the antibody-derived tags (ADT). You can download the ADT file [here](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz) and the RNA file [here](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz)

```
library(Seurat)
library(ggplot2)
library(patchwork)
```

```
# Load in the RNA UMI matrix

# Note that this dataset also contains ~5% of mouse cells, which we can use as negative
# controls for the protein measurements. For this reason, the gene expression matrix has
# HUMAN_ or MOUSE_ appended to the beginning of each gene.
cbmc.rna <- as.sparse(read.csv(file = "/brahms/shared/vignette-data/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz",
    sep = ",", header = TRUE, row.names = 1))

# To make life a bit easier going forward, we're going to discard all but the top 100 most
# highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Load in the ADT UMI matrix
cbmc.adt <- as.sparse(read.csv(file = "/brahms/shared/vignette-data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz",
    sep = ",", header = TRUE, row.names = 1))

# Note that since measurements were made in the same cells, the two matrices have identical
# column names
all.equal(colnames(cbmc.rna), colnames(cbmc.adt))
```

```
## [1] TRUE
```

## Setup a Seurat object, add the RNA and protein data

Now we create a Seurat object, and add the ADT data as a second assay

```
# creates a Seurat object based on the scRNA-seq data
cbmc <- CreateSeuratObject(counts = cbmc.rna)

# We can see that by default, the cbmc object contains an assay storing RNA measurement
Assays(cbmc)
```

```
## [1] "RNA"
```

```
# create a new assay to store ADT information
adt_assay <- CreateAssay5Object(counts = cbmc.adt)

# add this assay to the previously created Seurat object
cbmc[["ADT"]] <- adt_assay

# Validate that the object now contains multiple assays
Assays(cbmc)
```

```
## [1] "RNA" "ADT"
```

```
# Extract a list of features measured in the ADT assay
rownames(cbmc[["ADT"]])
```

```
##  [1] "CD3"    "CD4"    "CD8"    "CD45RA" "CD56"   "CD16"   "CD10"   "CD11c" 
##  [9] "CD14"   "CD19"   "CD34"   "CCR5"   "CCR7"
```

```
# Note that we can easily switch back and forth between the two assays to specify the default
# for visualization and analysis

# List the current default assay
DefaultAssay(cbmc)
```

```
## [1] "RNA"
```

```
# Switch the default to ADT
DefaultAssay(cbmc) <- "ADT"
DefaultAssay(cbmc)
```

```
## [1] "ADT"
```

## Cluster cells on the basis of their scRNA-seq profiles

The steps below represent a quick clustering of the PBMCs based on the scRNA-seq data. For more detail on individual steps or more advanced options, see our PBMC clustering guided tutorial [here](/seurat/articles/pbmc3k_tutorial)

```
# Note that all operations below are performed on the RNA assay Set and verify that the
# default assay is RNA
DefaultAssay(cbmc) <- "RNA"
DefaultAssay(cbmc)
```

```
## [1] "RNA"
```

```
# perform visualization and clustering steps
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
cbmc <- FindNeighbors(cbmc, dims = 1:30)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:30)
DimPlot(cbmc, label = TRUE)
```

![](multimodal_vignette_files/figure-html/cluster1-1.png)

## Visualize multiple modalities side-by-side

Now that we have obtained clusters from scRNA-seq profiles, we can visualize the expression of either protein or RNA molecules in our dataset. Importantly, Seurat provides a couple ways to switch between modalities, and specify which modality you are interested in analyzing or visualizing. This is particularly important as, in some cases, the same feature can be present in multiple modalities - for example this dataset contains independent measurements of the B cell marker CD19 (both protein and RNA levels).

```
# Normalize ADT data,
DefaultAssay(cbmc) <- "ADT"
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2)
DefaultAssay(cbmc) <- "RNA"

# Note that the following command is an alternative but returns the same result
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2, assay = "ADT")

# Now, we will visualize CD14 levels for RNA and protein By setting the default assay, we can
# visualize one or the other
DefaultAssay(cbmc) <- "ADT"
p1 <- FeaturePlot(cbmc, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")

# place plots side-by-side
p1 | p2
```

![](multimodal_vignette_files/figure-html/vis-1.png)

```
# Alternately, we can use specific assay keys to specify a specific modality Identify the key
# for the RNA and protein assays
Key(cbmc[["RNA"]])
```

```
## [1] "rna_"
```

```
Key(cbmc[["ADT"]])
```

```
## [1] "adt_"
```

```
# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(cbmc, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(cbmc, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2
```

![](multimodal_vignette_files/figure-html/vis-2.png)

## Identify cell surface markers for scRNA-seq clusters

We can leverage our paired CITE-seq measurements to help annotate clusters derived from scRNA-seq, and to identify both protein and RNA markers.

```
# as we know that CD19 is a B cell marker, we can identify cluster 6 as expressing CD19 on the
# surface
VlnPlot(cbmc, "adt_CD19")
```

![](multimodal_vignette_files/figure-html/markers-1.png)

```
# we can also identify alternative protein and RNA markers for this cluster through
# differential expression
adt_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "ADT")
rna_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "RNA")

head(adt_markers)
```

```
##                p_val avg_log2FC pct.1 pct.2     p_val_adj
## CD19   2.067533e-215  2.5741873     1     1 2.687793e-214
## CD45RA 8.108073e-109  0.5300346     1     1 1.054049e-107
## CD4    1.123162e-107 -1.6707420     1     1 1.460110e-106
## CD14   7.212876e-106 -1.0332070     1     1 9.376739e-105
## CD3     1.639633e-87 -1.5823056     1     1  2.131523e-86
## CCR5    2.552859e-63  0.3753989     1     1  3.318716e-62
```

```
head(rna_markers)
```

```
##       p_val avg_log2FC pct.1 pct.2 p_val_adj
## IGHM      0   6.660187 0.977 0.044         0
## CD79A     0   6.748356 0.965 0.045         0
## TCL1A     0   7.428099 0.904 0.028         0
## CD79B     0   5.525568 0.944 0.089         0
## IGHD      0   7.811884 0.857 0.015         0
## MS4A1     0   7.523215 0.851 0.016         0
```

## Additional visualizations of multimodal data

```
# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if
# desired by using HoverLocator and FeatureLocator
FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")
```

![](multimodal_vignette_files/figure-html/viz.cite.two-1.png)

```
# view relationship between protein and RNA
FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")
```

![](multimodal_vignette_files/figure-html/viz.cite.two-2.png)

```
FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8")
```

![](multimodal_vignette_files/figure-html/viz.cite.two-3.png)

```
# Let's look at the raw (non-normalized) ADT counts. You can see the values are quite high,
# particularly in comparison to RNA values. This is due to the significantly higher protein
# copy number in cells, which significantly reduces 'drop-out' in ADT data
FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")
```

![](multimodal_vignette_files/figure-html/viz.cite.two-4.png)

## Loading data from 10X multi-modal experiments

Seurat is also able to analyze data from multimodal 10X experiments processed using CellRanger v3; as an example, we recreate the plots above using a dataset of 7,900 peripheral blood mononuclear cells (PBMC), freely available from 10X Genomics [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3).

```
pbmc10k.data <- Read10X(data.dir = "/brahms/shared/vignette-data/pbmc10k/filtered_feature_bc_matrix/")
rownames(x = pbmc10k.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
    x = rownames(x = pbmc10k.data[["Antibody Capture"]]))

pbmc10k <- CreateSeuratObject(counts = pbmc10k.data[["Gene Expression"]], min.cells = 3, min.features = 200)
pbmc10k <- NormalizeData(pbmc10k)
pbmc10k[["ADT"]] <- CreateAssayObject(pbmc10k.data[["Antibody Capture"]][, colnames(x = pbmc10k)])
pbmc10k <- NormalizeData(pbmc10k, assay = "ADT", normalization.method = "CLR")

plot1 <- FeatureScatter(pbmc10k, feature1 = "adt_CD19", feature2 = "adt_CD3", pt.size = 1)
plot2 <- FeatureScatter(pbmc10k, feature1 = "adt_CD4", feature2 = "adt_CD8a", pt.size = 1)
plot3 <- FeatureScatter(pbmc10k, feature1 = "adt_CD3", feature2 = "CD3E", pt.size = 1)
(plot1 + plot2 + plot3) & NoLegend()
```

![](multimodal_vignette_files/figure-html/pbmc10x-1.png)

```
plot <- FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3") + NoLegend() + theme(axis.title = element_text(size = 18),
    legend.text = element_text(size = 18))
ggsave(filename = "../output/images/citeseq_plot.jpg", height = 7, width = 12, plot = plot, quality = 50)
```

## Additional functionality for multimodal data in Seurat

Seurat v4 also includes additional functionality for the analysis, visualization, and integration of multimodal datasets. For more information, please explore the resources below:

* Defining cellular identity from multimodal data using WNN analysis in Seurat v4 [vignette](/seurat/articles/weighted_nearest_neighbor_analysis)
* Mapping scRNA-seq data onto CITE-seq references [vignette](/seurat/articles/multimodal_reference_mapping)
* Introduction to the analysis of spatial transcriptomics analysis [vignette](/seurat/articles/spatial_vignette)
* Analysis of 10x multiome (paired scRNA-seq + ATAC) using WNN analysis [vignette](/seurat/articles/weighted_nearest_neighbor_analysis)
* Signac: Analysis, interpretation, and exploration of single-cell chromatin datasets [package](https://satijalab.org/signac/)
* Mixscape: an analytical toolkit for pooled single-cell genetic screens [vignette](/seurat/articles/mixscape_vignette)
* Bridge integration: mapping multi-omic datasets across molecular modalities [vignette](https://satijalab.org/seurat/articles/seurat5_integration_bridge)

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
## [1] patchwork_1.1.3    ggplot2_3.4.4      Seurat_5.0.1  testthat_3.2.0    
## [5] SeuratObject_5.0.1 sp_2.1-1          
## 
## loaded via a namespace (and not attached):
##   [1] spam_2.10-0            systemfonts_1.0.4      plyr_1.8.9            
##   [4] igraph_1.5.1           lazyeval_0.2.2         splines_4.2.2         
##   [7] RcppHNSW_0.5.0         listenv_0.9.0          scattermore_1.2       
##  [10] usethis_2.1.6          digest_0.6.33          htmltools_0.5.6.1     
##  [13] fansi_1.0.5            magrittr_2.0.3         memoise_2.0.1         
##  [16] tensor_1.5             cluster_2.1.4          ROCR_1.0-11           
##  [19] limma_3.54.1           remotes_2.4.2.1        globals_0.16.2        
##  [22] matrixStats_1.0.0      R.utils_2.12.2         pkgdown_2.0.7         
##  [25] spatstat.sparse_3.0-3  prettyunits_1.2.0      colorspace_2.1-0      
##  [28] ggrepel_0.9.4          textshaping_0.3.6      xfun_0.40             
##  [31] dplyr_1.1.3            callr_3.7.3            crayon_1.5.2          
##  [34] jsonlite_1.8.7         progressr_0.14.0       spatstat.data_3.0-3   
##  [37] survival_3.5-7         zoo_1.8-12             glue_1.6.2            
##  [40] polyclip_1.10-6        gtable_0.3.4           leiden_0.4.3          
##  [43] pkgbuild_1.4.2         future.apply_1.11.0    abind_1.4-5           
##  [46] scales_1.2.1           spatstat.random_3.2-1  miniUI_0.1.1.1        
##  [49] Rcpp_1.0.11            viridisLite_0.4.2      xtable_1.8-4          
##  [52] reticulate_1.34.0      dotCall64_1.1-0        profvis_0.3.7         
##  [55] htmlwidgets_1.6.2      httr_1.4.7             RColorBrewer_1.1-3    
##  [58] ellipsis_0.3.2         ica_1.0-3              R.methodsS3_1.8.2     
##  [61] farver_2.1.1           urlchecker_1.0.1       pkgconfig_2.0.3       
##  [64] uwot_0.1.16            sass_0.4.7             deldir_1.0-9          
##  [67] utf8_1.2.4             labeling_0.4.3         tidyselect_1.2.0      
##  [70] rlang_1.1.1            reshape2_1.4.4         later_1.3.1           
##  [73] munsell_0.5.0          tools_4.2.2            cachem_1.0.8          
##  [76] cli_3.6.1              generics_0.1.3         devtools_2.4.5        
##  [79] ggridges_0.5.4         evaluate_0.22          stringr_1.5.0         
##  [82] fastmap_1.1.1          goftest_1.2-3          yaml_2.3.7            
##  [85] ragg_1.2.5             processx_3.8.2         knitr_1.45            
##  [88] fs_1.6.3               fitdistrplus_1.1-11    purrr_1.0.2           
##  [91] RANN_2.6.1             nlme_3.1-162           pbapply_1.7-2         
##  [94] future_1.33.0          mime_0.12              formatR_1.14          
##  [97] R.oo_1.25.0            ggrastr_1.0.1          brio_1.1.3            
## [100] compiler_4.2.2         rstudioapi_0.14        beeswarm_0.4.0        
## [103] plotly_4.10.3          png_0.1-8              spatstat.utils_3.0-4  
## [106] tibble_3.2.1           bslib_0.5.1            stringi_1.7.12        
## [109] highr_0.10             ps_1.7.5               desc_1.4.2            
## [112] RSpectra_0.16-1        lattice_0.21-9         Matrix_1.6-1.1        
## [115] vctrs_0.6.4            pillar_1.9.0           lifecycle_1.0.4       
## [118] spatstat.geom_3.2-7    lmtest_0.9-40          jquerylib_0.1.4       
## [121] RcppAnnoy_0.0.21       data.table_1.14.8      cowplot_1.1.1         
## [124] irlba_2.3.5.1          httpuv_1.6.12          R6_2.5.1              
## [127] promises_1.2.1         KernSmooth_2.23-22     gridExtra_2.3         
## [130] vipor_0.4.5            parallelly_1.36.0      sessioninfo_1.2.2     
## [133] codetools_0.2-19       fastDummies_1.7.3      MASS_7.3-58.2         
## [136] pkgload_1.3.3          rprojroot_2.0.4        withr_2.5.2           
## [139] presto_1.0.0           sctransform_0.4.1      parallel_4.2.2        
## [142] grid_4.2.2             tidyr_1.3.0            rmarkdown_2.25        
## [145] Rtsne_0.16             spatstat.explore_3.2-5 shiny_1.7.5.1         
## [148] ggbeeswarm_0.7.1
```
