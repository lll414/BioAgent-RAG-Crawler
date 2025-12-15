---
source_url: https://satijalab.org/seurat/articles/covid_sctmapping
download_date: 2025-11-20
---

# Map COVID PBMC datasets to a healthy reference

#### Compiled: 2023-11-15

Source: [`vignettes/COVID_SCTMapping.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/COVID_SCTMapping.Rmd)

`COVID_SCTMapping.Rmd`

```
library(Seurat)
library(BPCells)
library(dplyr)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 1e9)
```

## Introduction: Reference mapping analysis in Seurat v5

In Seurat v5, we introduce a scalable approach for reference mapping datasets from separate studies or individuals. Reference mapping is a powerful approach to identify consistent labels across studies and perform cross-dataset analysis. We emphasize that while individual datasets are manageable in size, the aggregate of many datasets often amounts to millions of cell which do not fit in-memory. Furthermore, cross-dataset analysis is often challenged by disparate or unique cell type labels. Through reference mapping, we annotate all cells with a common reference for consistent cell type labels. Importantly, we never simultaneously load all of the cells in-memory to maintain low memory usage.

In this vignette, we reference map three publicly available datasets totaling 1,498,064 cells and 277 donors which are available through [CZI cellxgene collections](https://cellxgene.cziscience.com/collections): [Ahern, et al., Nature 2022](https://cellxgene.cziscience.com/collections/8f126edf-5405-4731-8374-b5ce11f53e82), [Jin, et al., Science 2021](https://cellxgene.cziscience.com/collections/b9fc3d70-5a72-4479-a046-c2cc1ab19efc), and [Yoshida, et al., Nature 2022](https://cellxgene.cziscience.com/collections/03f821b4-87be-4ff4-b65a-b5fc00061da7). Each dataset consists of PBMCs from both healthy donors and donors diagnosed with COVID-19. Using the harmonized annotations, we demonstrate how to prepare a pseudobulk object to perform differential expression analysis across disease within cell types.

Prior to running this vignette, please [install Seurat v5](/seurat/articles/install) and see the [BPCells vignette](/seurat/articles/seurat5_bpcells_interaction_vignette) to construct the on-disk object used in this vignette. Additionally, we map to our annotated CITE-seq reference containing 162,000 cells and 228 antibodies ([Hao*, Hao*, et al., Cell 2021](https://doi.org/10.1016/j.cell.2021.04.048)) which is available for download [here](https://zenodo.org/record/7779017#.ZCMojezMJqs).

## Load the PBMC Reference Dataset and Query Datasets

We first load the reference (available [here](https://zenodo.org/record/7779017#.ZCMojezMJqs)) and normalize the query Seurat object prepared in the [BPCells interaction vignette](/seurat/articles/seurat5_bpcells_interaction_vignette). The query object consists of datasets from three different studies constructed using the `CreateSeuratObject` function, which accepts a list of BPCells matrices as input. Within the Seurat object, the three datasets reside in the `RNA` assay in three separate `layers` on-disk.

```
reference <- readRDS("/brahms/hartmana/vignette_data/pbmc_multimodal_2023.rds")
object <- readRDS("/brahms/mollag/seurat_v5/vignette_data/merged_covid_object.rds")
object <- NormalizeData(object, verbose = FALSE)
```

## Mapping

Using the same code from the [v4 reference mapping vignette](articles/multimodal_reference_mapping.html), we find anchors between the reference and query in the precomputed supervised PCA. We recommend the use of supervised PCA for CITE-seq reference datasets, and demonstrate how to compute this transformation in [v4 reference mapping vignette](articles/multimodal_reference_mapping.html). In Seurat v5, we only need to call `FindTransferAnchors` and `MapQuery` once to map all three datasets as they are all contained within the query object. Furthermore, utilizing the on-disk capabilities of [BPCells](https://github.com/bnprks/BPCells), we map 1.5 million cells without ever loading them all into memory.

```
anchor <- FindTransferAnchors(
  reference = reference,
  query = object,
  reference.reduction = "spca",
  normalization.method = "SCT",
  dims = 1:50
)
object <- MapQuery(
  anchorset = anchor,
  query = object,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2"
  ),
  reduction.model = "wnn.umap"
)
```

## Explore the mapping results

Next, we visualize all cells from the three studies which have been projected into a UMAP-space defined by the reference. Each cell is annotated at two levels of granularity (`predicted.celltype.l1` and `predicted.celltype.l2`). We can compare the differing ontologies used in the original annotations (`cell_type`) to the now harmonized annotations (`predicted.celltype.l2`, for example) that were predicted from reference-mapping. Previously, the lack of standardization prevented us from directly performing integrative analysis across studies, but now we can easily compare.

```
DimPlot(object, reduction = "ref.umap", group.by = "cell_type", alpha = 0.1, label = TRUE, split.by = "publication", ncol = 3, label.size = 3) + NoLegend()
```

![](COVID_SCTMapping_files/figure-html/unnamed-chunk-3-1.png)

```
DimPlot(object, reduction = "ref.umap", group.by = "predicted.celltype.l2", alpha = 0.1, label = TRUE, split.by = "publication", ncol = 3, label.size = 3) + NoLegend()
```

![](COVID_SCTMapping_files/figure-html/unnamed-chunk-4-1.png)

## Differential composition analysis

We utilize our harmonized annotations to identify differences in the proportion of different cell types between healthy individuals and COVID-19 patients. For example, we noticed a reduction in MAIT cells as well as an increase in plasmablasts among COVID-19 patients.

```
df_comp <- as.data.frame.matrix(table(object$donor_id, object$predicted.celltype.l2))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_disease <- as.data.frame.matrix(table(object$donor_id, object$disease))[select.donors, ]

df_comp_relative$disease <- "other"
df_comp_relative$disease[df_disease$normal != 0] <- "normal"
df_comp_relative$disease[df_disease$`COVID-19` != 0] <- "COVID-19"
df_comp_relative$disease <- factor(df_comp_relative$disease, levels = c("normal", "COVID-19", "other"))
df_comp_relative <- df_comp_relative[df_comp_relative$disease %in% c("normal", "COVID-19"), ]
```

```
p1 <- ggplot(data = df_comp_relative, mapping = aes(x = disease, y = MAIT, fill = disease)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
  xlab("") +
  ylab("relative abundance") +
  ggtitle("MAIT") +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold")
  )

p2 <- ggplot(data = df_comp_relative, mapping = aes(x = disease, y = Plasmablast, fill = disease)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
  xlab("") +
  ylab("relative abundance") +
  ggtitle("Plasmablast") +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 15, hjust = 0.5, face = "bold")
  )

p1 + p2 + plot_layout(ncol = 2)
```

![](COVID_SCTMapping_files/figure-html/unnamed-chunk-6-1.png)

## Differential expression analysis

In addition to composition analysis, we use an aggregation-based (pseudobulk) workflow to explore differential genes between healthy individuals and COVID-19 donors. We aggregate all cells within the same cell type and donor using the `AggregateExpression` function. This returns a Seurat object where each ‘cell’ represents the pseudobulk profile of one cell type in one individual.

```
bulk <- AggregateExpression(object,
  return.seurat = TRUE,
  assays = "RNA",
  group.by = c("predicted.celltype.l2", "donor_id", "disease")
)
```

```
bulk <- subset(bulk, subset = disease %in% c("normal", "COVID-19"))
bulk <- subset(bulk, subset = predicted.celltype.l2 != "Doublet")
bulk$disease <- factor(bulk$disease, levels = c("normal", "COVID-19"))
```

Once a pseudobulk object is created, we can perform cell type-specific differential expression analysis between healthy individuals and COVID-19 donors. Here, we only visualize certain interferon-stimulated genes which are often upregulated during viral infection.

```
p1 <- VlnPlot(
  object = bulk, features = "IFI6", group.by = "predicted.celltype.l2",
  split.by = "disease", cols = c("#377eb8", "#e41a1c")
)
p2 <- VlnPlot(
  object = bulk, features = c("ISG15"), group.by = "predicted.celltype.l2",
  split.by = "disease", cols = c("#377eb8", "#e41a1c")
)
p3 <- VlnPlot(
  object = bulk, features = c("IFIT5"), group.by = "predicted.celltype.l2",
  split.by = "disease", cols = c("#377eb8", "#e41a1c")
)
p1 + p2 + p3 + plot_layout(ncol = 1)
```

![](COVID_SCTMapping_files/figure-html/unnamed-chunk-10-1.png)

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
## [1] ggplot2_3.4.4      patchwork_1.1.3    dplyr_1.1.3        BPCells_0.1.0     
## [5] Seurat_5.0.1  testthat_3.2.0     SeuratObject_5.0.1 sp_2.1-1          
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.4             spatstat.explore_3.2-5 reticulate_1.34.0     
##   [4] R.utils_2.12.2         tidyselect_1.2.0       htmlwidgets_1.6.2     
##   [7] grid_4.2.2             Rtsne_0.16             devtools_2.4.5        
##  [10] munsell_0.5.0          codetools_0.2-19       ragg_1.2.5            
##  [13] ica_1.0-3              future_1.33.0          miniUI_0.1.1.1        
##  [16] withr_2.5.2            spatstat.random_3.2-1  colorspace_2.1-0      
##  [19] progressr_0.14.0       highr_0.10             knitr_1.45            
##  [22] rstudioapi_0.14        stats4_4.2.2           ROCR_1.0-11           
##  [25] tensor_1.5             listenv_0.9.0          labeling_0.4.3        
##  [28] GenomeInfoDbData_1.2.9 polyclip_1.10-6        farver_2.1.1          
##  [31] rprojroot_2.0.4        parallelly_1.36.0      vctrs_0.6.4           
##  [34] generics_0.1.3         xfun_0.40              R6_2.5.1              
##  [37] GenomeInfoDb_1.34.9    ggbeeswarm_0.7.1       bitops_1.0-7          
##  [40] spatstat.utils_3.0-4   cachem_1.0.8           promises_1.2.1        
##  [43] scales_1.2.1           beeswarm_0.4.0         gtable_0.3.4          
##  [46] globals_0.16.2         processx_3.8.2         goftest_1.2-3         
##  [49] spam_2.10-0            rlang_1.1.1            systemfonts_1.0.4     
##  [52] splines_4.2.2          lazyeval_0.2.2         spatstat.geom_3.2-7   
##  [55] yaml_2.3.7             reshape2_1.4.4         abind_1.4-5           
##  [58] httpuv_1.6.12          tools_4.2.2            usethis_2.1.6         
##  [61] ellipsis_0.3.2         jquerylib_0.1.4        RColorBrewer_1.1-3    
##  [64] BiocGenerics_0.44.0    sessioninfo_1.2.2      ggridges_0.5.4        
##  [67] Rcpp_1.0.11            plyr_1.8.9             zlibbioc_1.44.0       
##  [70] purrr_1.0.2            RCurl_1.98-1.12        ps_1.7.5              
##  [73] prettyunits_1.2.0      deldir_1.0-9           pbapply_1.7-2         
##  [76] cowplot_1.1.1          urlchecker_1.0.1       S4Vectors_0.36.2      
##  [79] zoo_1.8-12             ggrepel_0.9.4          cluster_2.1.4         
##  [82] fs_1.6.3               magrittr_2.0.3         data.table_1.14.8     
##  [85] RSpectra_0.16-1        scattermore_1.2        lmtest_0.9-40         
##  [88] RANN_2.6.1             fitdistrplus_1.1-11    R.cache_0.16.0        
##  [91] matrixStats_1.0.0      pkgload_1.3.3          mime_0.12             
##  [94] evaluate_0.22          xtable_1.8-4           fastDummies_1.7.3     
##  [97] IRanges_2.32.0         gridExtra_2.3          compiler_4.2.2        
## [100] tibble_3.2.1           KernSmooth_2.23-22     crayon_1.5.2          
## [103] R.oo_1.25.0            htmltools_0.5.6.1      later_1.3.1           
## [106] tidyr_1.3.0            MASS_7.3-58.2          Matrix_1.6-1.1        
## [109] brio_1.1.3             cli_3.6.1              R.methodsS3_1.8.2     
## [112] parallel_4.2.2         dotCall64_1.1-0        igraph_1.5.1          
## [115] GenomicRanges_1.50.2   pkgconfig_2.0.3        pkgdown_2.0.7         
## [118] plotly_4.10.3          spatstat.sparse_3.0-3  vipor_0.4.5           
## [121] bslib_0.5.1            XVector_0.38.0         stringr_1.5.0         
## [124] callr_3.7.3            digest_0.6.33          sctransform_0.4.1     
## [127] RcppAnnoy_0.0.21       spatstat.data_3.0-3    rmarkdown_2.25        
## [130] leiden_0.4.3           uwot_0.1.16            shiny_1.7.5.1         
## [133] lifecycle_1.0.4        nlme_3.1-162           jsonlite_1.8.7        
## [136] desc_1.4.2             viridisLite_0.4.2      fansi_1.0.5           
## [139] pillar_1.9.0           lattice_0.21-9         ggrastr_1.0.1         
## [142] fastmap_1.1.1          httr_1.4.7             pkgbuild_1.4.2        
## [145] survival_3.5-7         glue_1.6.2             remotes_2.4.2.1       
## [148] png_0.1-8              stringi_1.7.12         sass_0.4.7            
## [151] profvis_0.3.7          textshaping_0.3.6      RcppHNSW_0.5.0        
## [154] memoise_2.0.1          styler_1.10.2          irlba_2.3.5.1         
## [157] future.apply_1.11.0
```