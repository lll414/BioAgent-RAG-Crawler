Source URL: https://bioconductor.org/books/release/SingleRBook/introduction.html
Date Scraped: 2025-12-15

---

# Chapter 1 Introduction



## 1.1 Motivation

The Bioconductor package *[SingleR](https://bioconductor.org/packages/3.22/SingleR)* implements an automatic annotation method
for single-cell RNA sequencing (scRNA-seq) data (Aran et al. [2019](#ref-aran2019reference)).
Given a reference dataset of samples (single-cell or bulk) with known labels,
it assigns those labels to new cells from a test dataset based on similarities in their expression profiles.
This provides a convenient way of transferring biological knowledge across datasets,
allowing users to leverage the domain expertise implicit in the creation of each reference.
The most common application of *[SingleR](https://bioconductor.org/packages/3.22/SingleR)* involves predicting cell type (or “state”, or “kind”) in a new dataset,
a process that is facilitated by the availability of curated references and compatibility with user-supplied datasets.
In this manner, the burden of manually interpreting clusters and defining marker genes only has to be done once, for the reference dataset, and this knowledge can be propagated to new datasets in an automated manner.

## 1.2 Method description

*[SingleR](https://bioconductor.org/packages/3.22/SingleR)* can be considered a robust variant of nearest-neighbors classification,
with some tweaks to improve resolution for closely related labels.
For each test cell:

1. We compute the Spearman correlation between its expression profile and that of each reference sample.
   The use of Spearman’s correlation provides a measure of robustness to batch effects across datasets.
   The calculation only uses the union of marker genes identified by pairwise comparisons between labels in the reference data,
   so as to improve resolution of separation between labels.
2. We define the per-label score as a fixed quantile (by default, 0.8) of the correlations across all samples with that label.
   This accounts for differences in the number of reference samples for each label,
   which interferes with simpler flavors of nearest neighbor classification;
   it also avoids penalizing classifications to heterogeneous labels by only requiring a good match to a minority of samples.
3. We repeat the score calculation for all labels in the reference dataset.
   The label with the highest score is used as *[SingleR](https://bioconductor.org/packages/3.22/SingleR)*’s prediction for this cell.
4. We optionally perform a fine-tuning step to improve resolution between closely related labels.
   The reference dataset is subsetted to only include labels with scores close to the maximum;
   scores are recomputed using only marker genes for the subset of labels, thus focusing on the most relevant features;
   and this process is iterated until only one label remains.

## 1.3 Quick start

We will demonstrate the use of `SingleR()` on a well-known 10X Genomics dataset (Zheng et al. [2017](#ref-zheng2017massively))
with the Human Primary Cell Atlas dataset (Mabbott et al. [2013](#ref-hpcaRef)) as the reference.

```
# Loading test data.
library(TENxPBMCData)
new.data <- TENxPBMCData("pbmc4k")

# Loading reference data with Ensembl annotations.
library(celldex)
ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE)

# Performing predictions.
library(SingleR)
predictions <- SingleR(test=new.data, assay.type.test=1, 
    ref=ref.data, labels=ref.data$label.main)

table(predictions$labels)
```

```
## 
##           B_cell              CMP               DC              GMP 
##              606                8                1                2 
##         Monocyte          NK_cell        Platelets Pre-B_cell_CD34- 
##             1164              217                3               46 
##          T_cells 
##             2293
```

And that’s it, really.

## 1.4 Where to get help

Questions on the general use of *[SingleR](https://bioconductor.org/packages/3.22/SingleR)* should be posted to
the [Bioconductor support site](https://support.bioconductor.org).
Please send requests for general assistance and advice to the
support site rather than to the individual authors.
Bug reports or feature requests should be made to the [GitHub repository](https://github.com/LTLA/SingleR);
well-considered suggestions for improvements are always welcome.

## Session information

View session info

```
R version 4.5.1 Patched (2025-08-23 r88802)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.3 LTS

Matrix products: default
BLAS:   /home/biocbuild/bbs-3.22-bioc/R/lib/libRblas.so 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB              LC_COLLATE=C              
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/New_York
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] SingleR_2.12.0              ensembldb_2.34.0           
 [3] AnnotationFilter_1.34.0     GenomicFeatures_1.62.0     
 [5] AnnotationDbi_1.72.0        celldex_1.20.0             
 [7] TENxPBMCData_1.28.0         HDF5Array_1.38.0           
 [9] h5mread_1.2.0               rhdf5_2.54.0               
[11] DelayedArray_0.36.0         SparseArray_1.10.1         
[13] S4Arrays_1.10.0             abind_1.4-8                
[15] Matrix_1.7-4                SingleCellExperiment_1.32.0
[17] SummarizedExperiment_1.40.0 Biobase_2.70.0             
[19] GenomicRanges_1.62.0        Seqinfo_1.0.0              
[21] IRanges_2.44.0              S4Vectors_0.48.0           
[23] BiocGenerics_0.56.0         generics_0.1.4             
[25] MatrixGenerics_1.22.0       matrixStats_1.5.0          
[27] BiocStyle_2.38.0            rebook_1.20.0              

loaded via a namespace (and not attached):
 [1] DBI_1.2.3                 bitops_1.0-9             
 [3] httr2_1.2.1               CodeDepends_0.6.6        
 [5] rlang_1.1.6               magrittr_2.0.4           
 [7] gypsum_1.6.0              compiler_4.5.1           
 [9] RSQLite_2.4.3             dir.expiry_1.18.0        
[11] DelayedMatrixStats_1.32.0 png_0.1-8                
[13] vctrs_0.6.5               ProtGenerics_1.42.0      
[15] pkgconfig_2.0.3           crayon_1.5.3             
[17] fastmap_1.2.0             dbplyr_2.5.1             
[19] XVector_0.50.0            Rsamtools_2.26.0         
[21] rmarkdown_2.30            UCSC.utils_1.6.0         
[23] graph_1.88.0              purrr_1.1.0              
[25] bit_4.6.0                 xfun_0.54                
[27] beachmat_2.26.0           cachem_1.1.0             
[29] cigarillo_1.0.0           GenomeInfoDb_1.46.0      
[31] jsonlite_2.0.0            blob_1.2.4               
[33] rhdf5filters_1.22.0       BiocParallel_1.44.0      
[35] Rhdf5lib_1.32.0           parallel_4.5.1           
[37] R6_2.6.1                  bslib_0.9.0              
[39] rtracklayer_1.70.0        jquerylib_0.1.4          
[41] Rcpp_1.1.0                bookdown_0.45            
[43] knitr_1.50                tidyselect_1.2.1         
[45] yaml_2.3.10               codetools_0.2-20         
[47] curl_7.0.0                lattice_0.22-7           
[49] tibble_3.3.0              withr_3.0.2              
[51] KEGGREST_1.50.0           evaluate_1.0.5           
[53] BiocFileCache_3.0.0       alabaster.schemas_1.10.0 
[55] ExperimentHub_3.0.0       Biostrings_2.78.0        
[57] pillar_1.11.1             BiocManager_1.30.26      
[59] filelock_1.0.3            RCurl_1.98-1.17          
[61] BiocVersion_3.22.0        sparseMatrixStats_1.22.0 
[63] alabaster.base_1.10.0     glue_1.8.0               
[65] alabaster.ranges_1.10.0   lazyeval_0.2.2           
[67] alabaster.matrix_1.10.0   tools_4.5.1              
[69] beachmat.hdf5_1.8.0       AnnotationHub_4.0.0      
[71] BiocIO_1.20.0             BiocNeighbors_2.4.0      
[73] GenomicAlignments_1.46.0  XML_3.99-0.19            
[75] grid_4.5.1                restfulr_0.0.16          
[77] cli_3.6.5                 rappdirs_0.3.3           
[79] dplyr_1.1.4               alabaster.se_1.10.0      
[81] sass_0.4.10               digest_0.6.37            
[83] rjson_0.2.23              memoise_2.0.1            
[85] htmltools_0.5.8.1         lifecycle_1.0.4          
[87] httr_1.4.7                bit64_4.6.0-1
```

### Bibliography

Aran, D., A. P. Looney, L. Liu, E. Wu, V. Fong, A. Hsu, S. Chak, et al. 2019. “Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.” *Nat. Immunol.* 20 (2): 163–72.

Mabbott, Neil A., J. K. Baillie, Helen Brown, Tom C. Freeman, and David A. Hume. 2013. “An expression atlas of human primary cells: Inference of gene function from coexpression networks.” *BMC Genomics* 14. <https://doi.org/10.1186/1471-2164-14-632>.

Zheng, G. X., J. M. Terry, P. Belgrader, P. Ryvkin, Z. W. Bent, R. Wilson, S. B. Ziraldo, et al. 2017. “Massively parallel digital transcriptional profiling of single cells.” *Nat Commun* 8 (January): 14049.