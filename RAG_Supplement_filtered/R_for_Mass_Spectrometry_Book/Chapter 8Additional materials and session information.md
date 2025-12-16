Source URL: https://rformassspectrometry.github.io/book/sec-si.html
Date Scraped: 2025-12-15

---

# Chapter 8 Additional materials and session information

## 8.1 Additional materials

* The [Single-cell proteomics data analysis using `QFeatures` and
  `scp`](https://lgatto.github.io/QFeaturesScpWorkshop2021/) workshop
  is provided as two vignettes. The first one provides a general
  introduction to the `QFeatures` class in the general context of mass
  spectrometry-based proteomics data manipulation. The second vignette
  focuses on single-cell application and introduces the `scp` package
  (Vanderaa and Gatto 2021) as an extension of `QFeatures`. This second
  vignette also provides exercises that give the attendee the
  opportunity to apply the learned concepts to reproduce a published
  analysis on a subset of a real data set. A [recent
  workshop](https://github.com/lgatto/2024_scpworkshop_EUBIC), offered
  at the 2024 EuBIC winter school, provides teaching material for the
  new *scplainer* analysis workflow.
* A tutorial presenting [Use Cases and Examples for Annotation of
  Untargeted Metabolomics
  Data](https://jorainer.github.io/MetaboAnnotationTutorials/) using
  the `MetaboAnnotation` and `MetaboCoreUtils` packages
  (Rainer et al. 2022).
* [Exploring and analyzing LC-MS data with Spectra and
  xcms](https://jorainer.github.io/xcmsTutorials/) provides an
  overview of recent developments in Bioconductor to work with mass
  spectrometry
  ([MsExperiment](https://github.com/RforMassSpectrometry/MsExperiment),
  [Spectra](https://github.com/RforMassSpectrometry/Spectra)) and
  specifically LC-MS data ([xcms](https://github.com/sneumann/xcms))
  and walks through the preprocessing of a small data set emphasizing
  on selection of data-dependent settings for the individual
  pre-processing steps.

The [SpectraTutorials](https://jorainer.github.io/SpectraTutorials) package provides three different vignettes:

* [Seamless Integration of Mass Spectrometry Data from Different
  Sources](https://jorainer.github.io/SpectraTutorials/articles/analyzing-MS-data-from-different-sources-with-Spectra.html):
  describes import/export of MS data from/to files in different format as well
  as processing and handling of MS data with the *Spectra* package.
* [Spectra: an Expandable Infrastructure to Handle Mass Spectrometry
  Data](https://jorainer.github.io/SpectraTutorials/articles/Spectra-backends.html):
  explains the concept of backends in *Spectra*, their properties, use cases
  along with performance considerations.
* [MS/MS Spectra Matching with the MetaboAnnotation
  Package](https://jorainer.github.io/SpectraTutorials/articles/Spectra-matching-with-MetaboAnnotation.html):
  explains how the *Spectra* package can be used together with the
  *[MetaboAnnotation](https://bioconductor.org/packages/3.22/MetaboAnnotation)* package in LC-MS/MS annotation
  workflows for untargeted metabolomics data.

## 8.2 Compiling the book locally

To compile and render the teaching material, you will also need
the *[BiocStyle](https://bioconductor.org/packages/3.22/BiocStyle)* package and the (slighly
modified) [Modern Statistics for Model Biology (msmb) HTML Book
Style](https://www-huber.embl.de/users/msmith/msmbstyle/) by Mike
Smith:

```
BiocManager::install(c("bookdown", "BiocStyle", "lgatto/msmbstyle"))
```

Clone the [book
repository](https://github.com/Rformassspectrometry/book) and render
the book with

```
bookdown::render_book(".")
```

## 8.3 Session information

The following packages have been used to generate this document.

```
sessionInfo()
```

```
## R version 4.5.0 (2025-04-11)
## Platform: x86_64-pc-linux-gnu
## Running under: Ubuntu 24.04.3 LTS
## 
## Matrix products: default
## BLAS:   /opt/R-4.5/lib/R/lib/libRblas.so 
## LAPACK: /opt/R-4.5/lib/R/lib/libRlapack.so;  LAPACK version 3.12.1
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Europe/Brussels
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] mzID_1.48.0                 patchwork_1.3.2            
##  [3] factoextra_1.0.7            gplots_3.2.0               
##  [5] limma_3.66.0                lubridate_1.9.4            
##  [7] forcats_1.0.1               stringr_1.6.0              
##  [9] purrr_1.2.0                 readr_2.1.6                
## [11] tidyr_1.3.1                 tibble_3.3.0               
## [13] tidyverse_2.0.0             MSnID_1.44.0               
## [15] cleaver_1.48.0              Biostrings_2.78.0          
## [17] XVector_0.50.0              ggplot2_4.0.1              
## [19] dplyr_1.1.4                 msdata_0.50.0              
## [21] MsDataHub_1.10.0            rpx_2.18.0                 
## [23] MsCoreUtils_1.21.0          QFeatures_1.20.0           
## [25] MultiAssayExperiment_1.36.0 SummarizedExperiment_1.40.0
## [27] Biobase_2.70.0              GenomicRanges_1.62.0       
## [29] Seqinfo_1.0.0               IRanges_2.44.0             
## [31] MatrixGenerics_1.22.0       matrixStats_1.5.0          
## [33] Spectra_1.20.0              BiocParallel_1.44.0        
## [35] S4Vectors_0.48.0            BiocGenerics_0.56.0        
## [37] generics_0.1.4              mzR_2.44.0                 
## [39] Rcpp_1.1.0                  BiocStyle_2.38.0           
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3     rstudioapi_0.17.1      jsonlite_2.0.0        
##   [4] magrittr_2.0.4         MALDIquant_1.22.3      farver_2.1.2          
##   [7] rmarkdown_2.30         fs_1.6.6               vctrs_0.6.5           
##  [10] memoise_2.0.1          RCurl_1.98-1.17        rstatix_0.7.3         
##  [13] htmltools_0.5.8.1      S4Arrays_1.10.0        BiocBaseUtils_1.12.0  
##  [16] AnnotationHub_4.0.0    curl_7.0.0             broom_1.0.10          
##  [19] Rhdf5lib_1.32.0        Formula_1.2-5          SparseArray_1.10.1    
##  [22] rhdf5_2.54.0           sass_0.4.10            KernSmooth_2.23-26    
##  [25] bslib_0.9.0            plyr_1.8.9             httr2_1.2.1           
##  [28] impute_1.84.0          cachem_1.1.0           igraph_2.2.1          
##  [31] lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3       
##  [34] Matrix_1.7-4           R6_2.6.1               fastmap_1.2.0         
##  [37] clue_0.3-66            digest_0.6.38          pcaMethods_2.2.0      
##  [40] AnnotationDbi_1.72.0   ExperimentHub_3.0.0    RSQLite_2.4.4         
##  [43] ggpubr_0.6.2           filelock_1.0.3         labeling_0.4.3        
##  [46] timechange_0.3.0       httr_1.4.7             abind_1.4-8           
##  [49] compiler_4.5.0         bit64_4.6.0-1          withr_3.0.2           
##  [52] doParallel_1.0.17      backports_1.5.0        S7_0.2.1              
##  [55] carData_3.0-5          DBI_1.2.3              R.utils_2.13.0        
##  [58] ggsignif_0.6.4         MASS_7.3-65            rappdirs_0.3.3        
##  [61] DelayedArray_0.36.0    caTools_1.18.3         gtools_3.9.5          
##  [64] tools_4.5.0            PSMatch_1.14.0         R.oo_1.27.1           
##  [67] glue_1.8.0             R.cache_0.17.0         rhdf5filters_1.22.0   
##  [70] grid_4.5.0             cluster_2.1.8.1        reshape2_1.4.5        
##  [73] gtable_0.3.6           msmbstyle_0.0.22       tzdb_0.5.0            
##  [76] preprocessCore_1.72.0  R.methodsS3_1.8.2      hms_1.1.4             
##  [79] data.table_1.17.8      MetaboCoreUtils_1.18.0 car_3.1-3             
##  [82] xml2_1.5.0             utf8_1.2.6             ggrepel_0.9.6         
##  [85] BiocVersion_3.22.0     foreach_1.5.2          pillar_1.11.1         
##  [88] BiocFileCache_3.0.0    lattice_0.22-7         bit_4.6.0             
##  [91] tidyselect_1.2.1       knitr_1.50             bookdown_0.34.2       
##  [94] ProtGenerics_1.42.0    xfun_0.54              statmod_1.5.1         
##  [97] MSnbase_2.36.0         pheatmap_1.0.13        stringi_1.8.7         
## [100] lazyeval_0.2.2        
##  [ reached 'max' / getOption("max.print") -- omitted 21 entries ]
```