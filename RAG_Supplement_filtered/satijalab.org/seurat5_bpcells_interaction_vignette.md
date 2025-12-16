---
source_url: https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette
download_date: 2025-11-20
---

# Using BPCells with Seurat Objects

#### Compiled: October 20, 2023

Source: [`vignettes/v5/seurat5_bpcells_interaction_vignette.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/v5/seurat5_bpcells_interaction_vignette.Rmd)

`seurat5_bpcells_interaction_vignette.Rmd`

---

BPCells is an [R package](https://github.com/bnprks/BPCells) that allows for computationally efficient single-cell analysis. It utilizes bit-packing compression to store counts matrices on disk and C++ code to cache operations.

We leverage the high performance capabilities of BPCells to work with Seurat objects in memory while accessing the counts on disk. In this vignette, we show how to use BPCells to load data, work with a Seurat objects in a more memory-efficient way, and write out Seurat objects with BPCells matrices.

We will show the methods for interacting with both a single dataset in one file or multiple datasets across multiple files using BPCells. BPCells allows us to easily analyze these large datasets in memory, and we encourage users to check out some of our other vignettes [here](/seurat/articles/seurat5_sketch_analysis) and [here](/seurat/articles/parsebio_sketch_integration) to see further applications.

```
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)
```

We use BPCells functionality to both load in our data and write the counts layers to bitpacked compressed binary files on disk to improve computation speeds. BPCells has multiple functions for reading in files.

## Load Data

### Load Data from one h5 file

In this section, we will load a dataset of mouse brain cells freely available from 10x Genomics. This includes 1.3 Million single cells that were sequenced on the Illumina NovaSeq 6000. The raw data can be found [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.3.0/1M_neurons?).

To read in the file, we will use open\_matrix\_10x\_hdf5, a BPCells function written to read in feature matrices from 10x. We then write a matrix directory, load the matrix, and create a Seurat object.

```
brain.data <- open_matrix_10x_hdf5(
  path = "/brahms/hartmana/vignette_data/1M_neurons_filtered_gene_bc_matrices_h5.h5"
)
# Write the matrix to a directory
 write_matrix_dir(
   mat = brain.data,
   dir = '/brahms/hartmana/vignette_data/bpcells/brain_counts')
# Now that we have the matrix on disk, we can load it
brain.mat <- open_matrix_dir(dir = "/brahms/hartmana/vignette_data/bpcells/brain_counts")
brain.mat <- Azimuth:::ConvertEnsembleToSymbol(mat = brain.mat, species = "mouse")

# Create Seurat Object
brain <- CreateSeuratObject(counts = brain.mat)
```

**What if I already have a Seurat Object?**

You can use BPCells to convert the matrices in your already created Seurat objects to on-disk matrices. Note, that this is only possible for V5 assays. As an example, if you’d like to convert the counts matrix of your RNA assay to a BPCells matrix, you can use the following:

```
# obj <- readRDS("/path/to/reference.rds")

# # Write the counts layer to a directory
 write_matrix_dir(mat = obj[["RNA"]]$counts, dir = '/brahms/hartmana/vignette_data/bpcells/brain_counts')
 counts.mat <- open_matrix_dir(dir = "/brahms/hartmana/vignette_data/bpcells/brain_counts")

# obj[["RNA"]]$counts <- counts.mat
```

#### Example Analsyis

Once this conversion is done, you can perform typical Seurat functions on the object. For example, we can normalize data and visualize features by automatically accessing the on-disk counts.

```
VlnPlot(brain, features = c("Sox10", "Slc17a7", "Aif1"), ncol = 3, layer = "counts", alpha = 0.1)
```

![](seurat5_bpcells_interaction_vignette_files/figure-html/unnamed-chunk-3-1.png)

```
# We then normalize and visualize again
brain <- NormalizeData(brain, normalization.method = "LogNormalize")
VlnPlot(brain, features = c("Sox10", "Slc17a7", "Aif1"), ncol = 3, layer = "data", alpha = 0.1)
```

![](seurat5_bpcells_interaction_vignette_files/figure-html/unnamed-chunk-3-2.png)

#### Saving Seurat objects with on-disk layers

If you save your object and load it in in the future, Seurat will access the on-disk matrices by their path, which is stored in the assay level data. `saveRDS()` can still be used to save your Seurat objects with on-disk matrices as shown below. Note, if you move the object across computers or to a place that no longer has access to the paths where the on-disk matrices are stored, you will have to move the folders that contain the on-disk matrices as well. We are actively working on functionality to make this process easier for users.

```
saveRDS(
   object = brain,
   file = "obj.Rds"
)
```

If needed, a layer with an on-disk matrix can be converted to an in-memory matrix using the `as()` function. For the purposes of this demo, we’ll subset the object so that it takes up less space in memory.

```
brain <- subset(brain, downsample = 1000)
brain[["RNA"]]$counts <- as(object = brain[["RNA"]]$counts, Class = "dgCMatrix")
```

### Load data from multiple h5ad files

You can also download data from multiple matrices. In this section, we create a Seurat object using multiple peripheral blood mononuclear cell (PBMC) samples that are freely available for download from CZI [here](https://cellxgene.cziscience.com/collections). We download data from [Ahern et al. (2022) Nature](https://cellxgene.cziscience.com/collections/8f126edf-5405-4731-8374-b5ce11f53e82), [Jin et al. (2021) Science](https://cellxgene.cziscience.com/collections/b9fc3d70-5a72-4479-a046-c2cc1ab19efc), and [Yoshida et al. (2022) Nature](https://cellxgene.cziscience.com/collections/03f821b4-87be-4ff4-b65a-b5fc00061da7). We use the BPCells function to read h5ad files.

```
file.dir <- "/brahms/hartmana/vignette_data/h5ad_files/"
files.set <- c("ahern_pbmc.h5ad", "jin_pbmc.h5ad", "yoshida_pbmc.h5ad")

# Loop through h5ad files and output BPCells matrices on-disk
data.list <- c()
metadata.list <- c()

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i])
  data <- open_matrix_anndata_hdf5(path)
   write_matrix_dir(
     mat = data,
     dir = paste0(gsub(".h5ad", "", path), "_BP")
   )
  # Load in BP matrices
  mat <- open_matrix_dir(dir = paste0(gsub(".h5ad", "", path), "_BP"))
  mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
  # Get metadata
  metadata.list[[i]] <- LoadH5ADobs(path = path)
  data.list[[i]] <- mat
}
# Name layers
names(data.list) <- c("ahern", "jin", "yoshida")

# Add Metadata
for (i in 1:length(metadata.list)) {
  metadata.list[[i]]$publication <- names(data.list)[i]
}
metadata.list <- lapply(metadata.list, function(x) {
  x <- x[, c("publication", "sex", "cell_type", "donor_id", "disease")]
  return(x)
})
metadata <- Reduce(rbind, metadata.list)
```

When we create the Seurat object with the list of matrices from each publication, we can then see that multiple counts layers exist that represent each dataset. This object contains over a million cells, yet only takes up minimal space in memory!

```
merged.object <- CreateSeuratObject(counts = data.list, meta.data = metadata)
merged.object
```

```
## An object of class Seurat 
## 26143 features across 1498064 samples within 1 assay 
## Active assay: RNA (26143 features, 0 variable features)
##  3 layers present: counts.ahern, counts.jin, counts.yoshida
```

```
saveRDS(
  object = merged.object,
  file = "obj.Rds"
)
```

### Parse Biosciences

Here, we show how to load a 1 million cell data set from Parse Biosciences and create a Seurat Object. The data is available for download [here](https://support.parsebiosciences.com/hc/en-us/articles/7704577188500-How-to-analyze-a-1-million-cell-data-set-using-Scanpy-and-Harmony)

```
parse.data <- open_matrix_anndata_hdf5(
  "/brahms/hartmana/vignette_data/h5ad_files/ParseBio_PBMC.h5ad"
)
```

```
write_matrix_dir(mat = parse.data, dir = "/brahms/hartmana/vignette_data/bpcells/parse_1m_pbmc")
```

```
parse.mat <- open_matrix_dir(dir = "/brahms/hartmana/vignette_data/bpcells/parse_1m_pbmc")
metadata <- readRDS("/brahms/hartmana/vignette_data/ParseBio_PBMC_meta.rds")
metadata$disease <- sapply(strsplit(x = metadata$sample, split = "_"), "[", 1)
parse.object <- CreateSeuratObject(counts = parse.mat, meta.data = metadata)
```

```
saveRDS(
  object = parse.object,
  file = "obj.Rds"
)
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
## [1] dplyr_1.1.3             biomaRt_2.54.1          Azimuth_0.4.6.9004     
## [4] shinyBS_0.61.1          SeuratDisk_0.0.0.9020   Seurat_4.9.9.9085      
## [7] SeuratObject_4.9.9.9092 sp_2.1-1                BPCells_0.1.0          
## 
## loaded via a namespace (and not attached):
##   [1] rappdirs_0.3.3                    rtracklayer_1.58.0               
##   [3] scattermore_1.2                   R.methodsS3_1.8.2                
##   [5] ragg_1.2.5                        tidyr_1.3.0                      
##   [7] JASPAR2020_0.99.10                ggplot2_3.4.4                    
##   [9] bit64_4.0.5                       knitr_1.44                       
##  [11] irlba_2.3.5.1                     DelayedArray_0.24.0              
##  [13] R.utils_2.12.2                    styler_1.10.2                    
##  [15] data.table_1.14.8                 KEGGREST_1.38.0                  
##  [17] TFBSTools_1.36.0                  RCurl_1.98-1.12                  
##  [19] AnnotationFilter_1.22.0           generics_0.1.3                   
##  [21] BiocGenerics_0.44.0               GenomicFeatures_1.50.4           
##  [23] cowplot_1.1.1                     RSQLite_2.3.1                    
##  [25] RANN_2.6.1                        future_1.33.0                    
##  [27] bit_4.0.5                         tzdb_0.4.0                       
##  [29] spatstat.data_3.0-1               xml2_1.3.5                       
##  [31] httpuv_1.6.11                     SummarizedExperiment_1.28.0      
##  [33] DirichletMultinomial_1.40.0       gargle_1.5.2                     
##  [35] xfun_0.40                         hms_1.1.3                        
##  [37] jquerylib_0.1.4                   evaluate_0.22                    
##  [39] promises_1.2.1                    fansi_1.0.5                      
##  [41] restfulr_0.0.15                   progress_1.2.2                   
##  [43] caTools_1.18.2                    dbplyr_2.3.4                     
##  [45] igraph_1.5.1                      DBI_1.1.3                        
##  [47] htmlwidgets_1.6.2                 spatstat.geom_3.2-5              
##  [49] googledrive_2.1.1                 stats4_4.2.2                     
##  [51] purrr_1.0.2                       ellipsis_0.3.2                   
##  [53] RSpectra_0.16-1                   annotate_1.76.0                  
##  [55] deldir_1.0-9                      MatrixGenerics_1.10.0            
##  [57] vctrs_0.6.4                       Biobase_2.58.0                   
##  [59] Cairo_1.6-0                       ensembldb_2.22.0                 
##  [61] ROCR_1.0-11                       abind_1.4-5                      
##  [63] cachem_1.0.8                      withr_2.5.1                      
##  [65] BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.66.3                  
##  [67] progressr_0.14.0                  presto_1.0.0                     
##  [69] sctransform_0.4.1                 GenomicAlignments_1.34.1         
##  [71] prettyunits_1.1.1                 goftest_1.2-3                    
##  [73] cluster_2.1.4                     dotCall64_1.1-0                  
##  [75] lazyeval_0.2.2                    seqLogo_1.64.0                   
##  [77] crayon_1.5.2                      hdf5r_1.3.8                      
##  [79] spatstat.explore_3.2-3            labeling_0.4.3                   
##  [81] pkgconfig_2.0.3                   GenomeInfoDb_1.34.9              
##  [83] vipor_0.4.5                       nlme_3.1-162                     
##  [85] ProtGenerics_1.30.0               rlang_1.1.1                      
##  [87] globals_0.16.2                    lifecycle_1.0.3                  
##  [89] miniUI_0.1.1.1                    filelock_1.0.2                   
##  [91] fastDummies_1.7.3                 BiocFileCache_2.6.1              
##  [93] SeuratData_0.2.2.9001             ggrastr_1.0.1                    
##  [95] cellranger_1.1.0                  rprojroot_2.0.3                  
##  [97] polyclip_1.10-6                   RcppHNSW_0.5.0                   
##  [99] matrixStats_1.0.0                 lmtest_0.9-40                    
## [101] Matrix_1.6-1.1                    Rhdf5lib_1.20.0                  
## [103] zoo_1.8-12                        beeswarm_0.4.0                   
## [105] googlesheets4_1.1.1               ggridges_0.5.4                   
## [107] png_0.1-8                         viridisLite_0.4.2                
## [109] rjson_0.2.21                      shinydashboard_0.7.2             
## [111] bitops_1.0-7                      R.oo_1.25.0                      
## [113] rhdf5filters_1.10.1               KernSmooth_2.23-22               
## [115] spam_2.9-1                        Biostrings_2.66.0                
## [117] blob_1.2.4                        stringr_1.5.0                    
## [119] parallelly_1.36.0                 spatstat.random_3.1-6            
## [121] R.cache_0.16.0                    readr_2.1.4                      
## [123] S4Vectors_0.36.2                  CNEr_1.34.0                      
## [125] scales_1.2.1                      memoise_2.0.1                    
## [127] magrittr_2.0.3                    plyr_1.8.9                       
## [129] ica_1.0-3                         zlibbioc_1.44.0                  
## [131] compiler_4.2.2                    BiocIO_1.8.0                     
## [133] RColorBrewer_1.1-3                fitdistrplus_1.1-11              
## [135] Rsamtools_2.14.0                  cli_3.6.1                        
## [137] XVector_0.38.0                    listenv_0.9.0                    
## [139] patchwork_1.1.3                   pbapply_1.7-2                    
## [141] MASS_7.3-58.2                     tidyselect_1.2.0                 
## [143] stringi_1.7.12                    textshaping_0.3.6                
## [145] yaml_2.3.7                        ggrepel_0.9.4                    
## [147] grid_4.2.2                        sass_0.4.7                       
## [149] fastmatch_1.1-4                   EnsDb.Hsapiens.v86_2.99.0        
## [151] tools_4.2.2                       future.apply_1.11.0              
## [153] parallel_4.2.2                    rstudioapi_0.14                  
## [155] TFMPvalue_0.0.9                   gridExtra_2.3                    
## [157] farver_2.1.1                      Rtsne_0.16                       
## [159] digest_0.6.33                     pracma_2.4.2                     
## [161] shiny_1.7.5.1                     Rcpp_1.0.11                      
## [163] GenomicRanges_1.50.2              later_1.3.1                      
## [165] RcppAnnoy_0.0.21                  httr_1.4.7                       
## [167] AnnotationDbi_1.60.2              colorspace_2.1-0                 
## [169] XML_3.99-0.14                     fs_1.6.3                         
## [171] tensor_1.5                        reticulate_1.34.0                
## [173] IRanges_2.32.0                    splines_4.2.2                    
## [175] uwot_0.1.16                       RcppRoll_0.3.0                   
## [177] spatstat.utils_3.0-3              pkgdown_2.0.7                    
## [179] plotly_4.10.2                     systemfonts_1.0.4                
## [181] xtable_1.8-4                      jsonlite_1.8.7                   
## [183] poweRlaw_0.70.6                   R6_2.5.1                         
## [185] pillar_1.9.0                      htmltools_0.5.6.1                
## [187] mime_0.12                         glue_1.6.2                       
## [189] fastmap_1.1.1                     DT_0.30                          
## [191] BiocParallel_1.32.6               codetools_0.2-19                 
## [193] Signac_1.12.9000                  utf8_1.2.3                       
## [195] lattice_0.21-9                    bslib_0.5.1                      
## [197] spatstat.sparse_3.0-2             tibble_3.2.1                     
## [199] ggbeeswarm_0.7.1                  curl_5.1.0                       
## [201] leiden_0.4.3                      gtools_3.9.4                     
## [203] shinyjs_2.1.0                     GO.db_3.16.0                     
## [205] survival_3.5-7                    rmarkdown_2.25                   
## [207] desc_1.4.2                        munsell_0.5.0                    
## [209] rhdf5_2.42.1                      GenomeInfoDbData_1.2.9           
## [211] reshape2_1.4.4                    gtable_0.3.4
```
