---
source_url: https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette
download_date: 2025-11-20
---

# Integrating scRNA-seq and scATAC-seq data

#### Compiled: October 31, 2023

Source: [`vignettes/seurat5_atacseq_integration_vignette.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/seurat5_atacseq_integration_vignette.Rmd)

`seurat5_atacseq_integration_vignette.Rmd`

Single-cell transcriptomics has transformed our ability to characterize cell states, but deep biological understanding requires more than a taxonomic listing of clusters. As new methods arise to measure distinct cellular modalities, a key analytical challenge is to integrate these datasets to better understand cellular identity and function. For example, users may perform scRNA-seq and scATAC-seq experiments on the same biological system and to consistently annotate both datasets with the same set of cell type labels. This analysis is particularly challenging as scATAC-seq datasets are difficult to annotate, due to both the sparsity of genomic data collected at single-cell resolution, and the lack of interpretable gene markers in scRNA-seq data.

In [Stuart\*, Butler\* et al, 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8), we introduce methods to integrate scRNA-seq and scATAC-seq datasets collected from the same biological system, and demonstrate these methods in this vignette. In particular, we demonstrate the following analyses:

* How to use an annotated scRNA-seq dataset to label cells from an scATAC-seq experiment
* How to co-visualize (co-embed) cells from scRNA-seq and scATAC-seq
* How to project scATAC-seq cells onto a UMAP derived from an scRNA-seq experiment

This vignette makes extensive use of the [Signac package](https://stuartlab.org/signac/), recently developed for the analysis of chromatin datasets collected at single-cell resolution, including scATAC-seq. Please see the Signac website for additional [vignettes](https://stuartlab.org/signac/articles/pbmc_vignette.html) and documentation for analyzing scATAC-seq data.

We demonstrate these methods using a publicly available ~12,000 human PBMC ‘multiome’ dataset from 10x Genomics. In this dataset, scRNA-seq and scATAC-seq profiles were simultaneously collected in the same cells. For the purposes of this vignette, we treat the datasets as originating from two different experiments and integrate them together. Since they were originally measured in the same cells, this provides a ground truth that we can use to assess the accuracy of integration. We emphasize that our use of the multiome dataset here is for demonstration and evaluation purposes, and that users should apply these methods to scRNA-seq and scATAC-seq datasets that are collected separately. We provide a separate [weighted nearest neighbors vignette (WNN)](/seurat/articles/weighted_nearest_neighbor_analysis) that describes analysis strategies for multi-omic single-cell data.

## Load in data and process each modality individually

The PBMC multiome dataset is available from [10x genomics](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k). To facilitate easy loading and exploration, it is also available as part of our SeuratData package. We load the RNA and ATAC data in separately, and pretend that these profiles were measured in separate experiments. We annotated these cells in our [WNN](/seurat/articles/weighted_nearest_neighbor_analysis) vignette, and the annotations are also included in SeuratData.

```
library(SeuratData)
# install the dataset and load requirements
InstallData("pbmcMultiome")
```

```
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
```

```
# load both modalities
pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")

pbmc.rna[["RNA"]] <- as(pbmc.rna[["RNA"]], Class = "Assay5")
# repeat QC steps performed in the WNN vignette
pbmc.rna <- subset(pbmc.rna, seurat_annotations != "filtered")
pbmc.atac <- subset(pbmc.atac, seurat_annotations != "filtered")

# Perform standard analysis of each modality independently RNA analysis
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(pbmc.atac) <- annotations

# We exclude the first dimension as this is typically correlated with sequencing depth
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
pbmc.atac <- RunSVD(pbmc.atac)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
```

Now we plot the results from both modalities. Cells have been previously annotated based on transcriptomic state. We will predict annotations for the scATAC-seq cells.

```
p1 <- DimPlot(pbmc.rna, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2
```

![](seurat5_atacseq_integration_vignette_files/figure-html/viz1-1.png)

```
plot <- (p1 + p2) & xlab("UMAP 1") & ylab("UMAP 2") & theme(axis.title = element_text(size = 18))
ggsave(filename = "../output/images/atacseq_integration_vignette.jpg", height = 7, width = 12, plot = plot,
    quality = 50)
```

## Identifying anchors between scRNA-seq and scATAC-seq datasets

In order to identify ‘anchors’ between scRNA-seq and scATAC-seq experiments, we first generate a rough estimate of the transcriptional activity of each gene by quantifying ATAC-seq counts in the 2 kb-upstream region and gene body, using the `GeneActivity()` function in the Signac package. The ensuing gene activity scores from the scATAC-seq data are then used as input for canonical correlation analysis, along with the gene expression quantifications from scRNA-seq. We perform this quantification for all genes identified as being highly variable from the scRNA-seq dataset.

```
# quantify gene activity
gene.activities <- GeneActivity(pbmc.atac, features = VariableFeatures(pbmc.rna))
```

```
# add gene activities as a new assay
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))
```

```
# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna),
    reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
```

## Annotate scATAC-seq cells via label transfer

After identifying anchors, we can transfer annotations from the scRNA-seq dataset onto the scATAC-seq cells. The annotations are stored in the `seurat_annotations` field, and are provided as input to the `refdata` parameter. The output will contain a matrix with predictions and confidence scores for each ATAC-seq cell.

```
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$seurat_annotations,
    weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)

pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
```

**Why do you choose different (non-default) values for reduction and weight.reduction?**

In `FindTransferAnchors()`, we typically project the PCA structure from the reference onto the query when transferring between scRNA-seq datasets. However, when transferring across modalities we find that CCA better captures the shared feature correlation structure and therefore set `reduction = 'cca'` here. Additionally, by default in `TransferData()` we use the same projected PCA structure to compute the weights of the local neighborhood of anchors that influence each cell’s prediction. In the case of scRNA-seq to scATAC-seq transfer, we use the low dimensional space learned by computing an LSI on the ATAC-seq data to compute these weights as this better captures the internal structure of the ATAC-seq data.

After performing transfer, the ATAC-seq cells have predicted annotations (transferred from the scRNA-seq dataset) stored in the `predicted.id` field. Since these cells were measured with the multiome kit, we also have a ground-truth annotation that can be used for evaluation. You can see that the predicted and actual annotations are extremely similar.

```
pbmc.atac$annotation_correct <- pbmc.atac$predicted.id == pbmc.atac$seurat_annotations
p1 <- DimPlot(pbmc.atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(pbmc.atac, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2
```

![](seurat5_atacseq_integration_vignette_files/figure-html/viz.label.accuracy-1.png)

In this example, the annotation for an scATAC-seq profile is correctly predicted via scRNA-seq integration ~90% of the time. In addition, the `prediction.score.max` field quantifies the uncertainty associated with our predicted annotations. We can see that cells that are correctly annotated are typically associated with high prediction scores (>90%), while cells that are incorrectly annotated are associated with sharply lower prediction scores (<50%). Incorrect assignments also tend to reflect closely related cell types (i.e. Intermediate vs. Naive B cells).

```
predictions <- table(pbmc.atac$seurat_annotations, pbmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(pbmc.atac$seurat_annotations == pbmc.atac$predicted.id))
incorrect <- length(which(pbmc.atac$seurat_annotations != pbmc.atac$predicted.id))
data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
    geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2
```

![](seurat5_atacseq_integration_vignette_files/figure-html/score.viz-1.png)

## Co-embedding scRNA-seq and scATAC-seq datasets

In addition to transferring labels across modalities, it is also possible to visualize scRNA-seq and scATAC-seq cells on the same plot. We emphasize that this step is primarily for visualization, and is an optional step. Typically, when we perform integrative analysis between scRNA-seq and scATAC-seq datasets, we focus primarily on label transfer as described above. We demonstrate our workflows for co-embedding below, and again highlight that this is for demonstration purposes, especially as in this particular case both the scRNA-seq profiles and scATAC-seq profiles were actually measured in the same cells.

In order to perform co-embedding, we first ‘impute’ RNA expression into the scATAC-seq cells based on the previously computed anchors, and then merge the datasets.

```
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(pbmc.rna)
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]],
    dims = 2:30)
pbmc.atac[["RNA"]] <- imputation

coembed <- merge(x = pbmc.rna, y = pbmc.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

DimPlot(coembed, group.by = c("orig.ident", "seurat_annotations"))
```

![](seurat5_atacseq_integration_vignette_files/figure-html/coembed-1.png)

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
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] cowplot_1.1.1                 ggplot2_3.4.4                
##  [3] EnsDb.Hsapiens.v86_2.99.0     ensembldb_2.22.0             
##  [5] AnnotationFilter_1.22.0       GenomicFeatures_1.50.4       
##  [7] AnnotationDbi_1.60.2          Biobase_2.58.0               
##  [9] GenomicRanges_1.50.2          GenomeInfoDb_1.34.9          
## [11] IRanges_2.32.0                S4Vectors_0.36.2             
## [13] BiocGenerics_0.44.0           Signac_1.12.9000             
## [15] Seurat_5.0.1                  SeuratObject_5.0.1           
## [17] sp_2.1-1                      thp1.eccite.SeuratData_3.1.5 
## [19] stxBrain.SeuratData_0.1.1     ssHippo.SeuratData_3.1.4     
## [21] pbmcsca.SeuratData_3.0.0      pbmcref.SeuratData_1.0.0     
## [23] pbmcMultiome.SeuratData_0.1.4 pbmc3k.SeuratData_3.1.4      
## [25] panc8.SeuratData_3.0.2        ifnb.SeuratData_3.0.0        
## [27] hcabm40k.SeuratData_3.0.0     cbmc.SeuratData_3.1.4        
## [29] bmcite.SeuratData_0.3.0       SeuratData_0.2.2.9001        
## 
## loaded via a namespace (and not attached):
##   [1] rappdirs_0.3.3              rtracklayer_1.58.0         
##   [3] scattermore_1.2             ragg_1.2.5                 
##   [5] tidyr_1.3.0                 bit64_4.0.5                
##   [7] knitr_1.45                  irlba_2.3.5.1              
##   [9] DelayedArray_0.24.0         data.table_1.14.8          
##  [11] rpart_4.1.19                KEGGREST_1.38.0            
##  [13] RCurl_1.98-1.12             generics_0.1.3             
##  [15] RSQLite_2.3.1               RANN_2.6.1                 
##  [17] future_1.33.0               bit_4.0.5                  
##  [19] spatstat.data_3.0-3         xml2_1.3.5                 
##  [21] httpuv_1.6.12               SummarizedExperiment_1.28.0
##  [23] xfun_0.40                   hms_1.1.3                  
##  [25] jquerylib_0.1.4             evaluate_0.22              
##  [27] promises_1.2.1              fansi_1.0.5                
##  [29] restfulr_0.0.15             progress_1.2.2             
##  [31] dbplyr_2.3.4                igraph_1.5.1               
##  [33] DBI_1.1.3                   htmlwidgets_1.6.2          
##  [35] spatstat.geom_3.2-7         purrr_1.0.2                
##  [37] ellipsis_0.3.2              RSpectra_0.16-1            
##  [39] dplyr_1.1.3                 backports_1.4.1            
##  [41] biomaRt_2.54.1              deldir_1.0-9               
##  [43] MatrixGenerics_1.10.0       vctrs_0.6.4                
##  [45] ROCR_1.0-11                 abind_1.4-5                
##  [47] cachem_1.0.8                withr_2.5.2                
##  [49] BSgenome_1.66.3             progressr_0.14.0           
##  [51] checkmate_2.2.0             sctransform_0.4.1          
##  [53] GenomicAlignments_1.34.1    prettyunits_1.1.1          
##  [55] goftest_1.2-3               cluster_2.1.4              
##  [57] dotCall64_1.1-0             lazyeval_0.2.2             
##  [59] crayon_1.5.2                spatstat.explore_3.2-5     
##  [61] labeling_0.4.3              pkgconfig_2.0.3            
##  [63] nlme_3.1-162                ProtGenerics_1.30.0        
##  [65] nnet_7.3-18                 rlang_1.1.1                
##  [67] globals_0.16.2              lifecycle_1.0.3            
##  [69] miniUI_0.1.1.1              filelock_1.0.2             
##  [71] fastDummies_1.7.3           BiocFileCache_2.6.1        
##  [73] dichromat_2.0-0.1           rprojroot_2.0.3            
##  [75] polyclip_1.10-6             RcppHNSW_0.5.0             
##  [77] matrixStats_1.0.0           lmtest_0.9-40              
##  [79] Matrix_1.6-1.1              zoo_1.8-12                 
##  [81] base64enc_0.1-3             ggridges_0.5.4             
##  [83] png_0.1-8                   viridisLite_0.4.2          
##  [85] rjson_0.2.21                bitops_1.0-7               
##  [87] KernSmooth_2.23-22          spam_2.10-0                
##  [89] Biostrings_2.66.0           blob_1.2.4                 
##  [91] stringr_1.5.0               parallelly_1.36.0          
##  [93] spatstat.random_3.2-1       scales_1.2.1               
##  [95] memoise_2.0.1               magrittr_2.0.3             
##  [97] plyr_1.8.9                  ica_1.0-3                  
##  [99] zlibbioc_1.44.0             compiler_4.2.2             
## [101] BiocIO_1.8.0                RColorBrewer_1.1-3         
## [103] fitdistrplus_1.1-11         Rsamtools_2.14.0           
## [105] cli_3.6.1                   XVector_0.38.0             
## [107] listenv_0.9.0               patchwork_1.1.3            
## [109] pbapply_1.7-2               htmlTable_2.4.1            
## [111] formatR_1.14                Formula_1.2-5              
## [113] MASS_7.3-58.2               tidyselect_1.2.0           
## [115] stringi_1.7.12              textshaping_0.3.6          
## [117] highr_0.10                  yaml_2.3.7                 
## [119] ggrepel_0.9.4               grid_4.2.2                 
## [121] VariantAnnotation_1.44.1    sass_0.4.7                 
## [123] fastmatch_1.1-4             tools_4.2.2                
## [125] future.apply_1.11.0         parallel_4.2.2             
## [127] rstudioapi_0.14             foreign_0.8-84             
## [129] gridExtra_2.3               farver_2.1.1               
## [131] Rtsne_0.16                  digest_0.6.33              
## [133] shiny_1.7.5.1               Rcpp_1.0.11                
## [135] later_1.3.1                 RcppAnnoy_0.0.21           
## [137] httr_1.4.7                  biovizBase_1.46.0          
## [139] colorspace_2.1-0            XML_3.99-0.14              
## [141] fs_1.6.3                    tensor_1.5                 
## [143] reticulate_1.34.0           splines_4.2.2              
## [145] uwot_0.1.16                 RcppRoll_0.3.0             
## [147] spatstat.utils_3.0-4        pkgdown_2.0.7              
## [149] plotly_4.10.3               systemfonts_1.0.4          
## [151] xtable_1.8-4                jsonlite_1.8.7             
## [153] R6_2.5.1                    Hmisc_5.1-1                
## [155] pillar_1.9.0                htmltools_0.5.6.1          
## [157] mime_0.12                   glue_1.6.2                 
## [159] fastmap_1.1.1               BiocParallel_1.32.6        
## [161] codetools_0.2-19            utf8_1.2.4                 
## [163] lattice_0.21-9              bslib_0.5.1                
## [165] spatstat.sparse_3.0-3       tibble_3.2.1               
## [167] curl_5.1.0                  leiden_0.4.3               
## [169] survival_3.5-7              rmarkdown_2.25             
## [171] desc_1.4.2                  munsell_0.5.0              
## [173] GenomeInfoDbData_1.2.9      reshape2_1.4.4             
## [175] gtable_0.3.4
```
