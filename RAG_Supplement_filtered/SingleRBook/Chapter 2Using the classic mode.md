Source URL: https://bioconductor.org/books/release/SingleRBook/classic-mode.html
Date Scraped: 2025-12-15

---

# Chapter 2 Using the classic mode



## 2.1 Overview

*[SingleR](https://bioconductor.org/packages/3.22/SingleR)* detects markers in a pairwise manner between labels in the reference dataset.
Specifically, for each label of interest, it performs pairwise comparisons to every other label in the reference
and identifies the genes that are upregulated in the label of interest for each comparison.
The initial score calculation is then performed on the union of marker genes across all comparisons for all label.
This approach ensures that the selected subset of features will contain genes that distinguish each label from any other label.
(In contrast, other approaches that treat the “other” labels as a single group do not offer this guarantee;
see [here](https://osca.bioconductor.org/marker-detection.html#standard-application) for a discussion.)
It also allows the fine-tuning step to aggressively improve resolution by only using marker genes
from comparisons where both labels have scores close to the maximum.

The original (“classic”) marker detection algorithm used in Aran et al. ([2019](#ref-aran2019reference)) identified marker genes
based on their log-fold changes in each pairwise comparison.
Specifically, it used the genes with the largest positive differences in the per-label median log-expression values between labels.
The number of genes taken from each pairwise comparison was defined as 500(23)log2(n),
where n is the number of unique labels in the reference;
this scheme aimed to reduce the number of genes (and thus the computational time) as the number of labels and pairwise comparisons increased.
Classic mode is primarily intended for reference datasets that have little or no replication,
a description that covers many of the bulk-derived references
and precludes more complicated marker detection procedures (Chapter [3](sc-mode.html#sc-mode)).

## 2.2 Annotating the test dataset

For demonstration purposes, we will use the Grun et al. ([2016](#ref-grun2016denovo)) haematopoietic stem cell (HSC)
dataset from the *[scRNAseq](https://bioconductor.org/packages/3.22/scRNAseq)* package.
The `GrunHSCData()` function conveniently returns a `SingleCellExperiment`
object containing the count matrix for this dataset.

```
library(scRNAseq)
sce <- GrunHSCData(ensembl=TRUE)
sce
```

```
## class: SingleCellExperiment 
## dim: 21817 1915 
## metadata(0):
## assays(1): counts
## rownames(21817): ENSMUSG00000109644 ENSMUSG00000007777 ...
##   ENSMUSG00000055670 ENSMUSG00000039068
## rowData names(3): symbol chr originalName
## colnames(1915): JC4_349_HSC_FE_S13_ JC4_350_HSC_FE_S13_ ...
##   JC48P6_1203_HSC_FE_S8_ JC48P6_1204_HSC_FE_S8_
## colData names(2): sample protocol
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```

Our aim is to annotate each cell with the ImmGen reference dataset (Heng et al. [2008](#ref-ImmGenRef)) from the *[celldex](https://bioconductor.org/packages/3.22/celldex)* package.
(Some further comments on the choice of reference are provided below in Section [2.5](classic-mode.html#reference-choice).)
Calling the `ImmGenData()` function returns a `SummarizedExperiment` object
containing a matrix of log-expression values with sample-level labels.
We also set `ensembl=TRUE` to match the reference’s gene annotation with that in the `sce` object -
the default behavior is to use the gene symbol.

```
library(celldex)
immgen <- ImmGenData(ensembl=TRUE)
immgen
```

```
## class: SummarizedExperiment 
## dim: 21352 830 
## metadata(0):
## assays(1): logcounts
## rownames(21352): ENSMUSG00000079681 ENSMUSG00000066372 ...
##   ENSMUSG00000034640 ENSMUSG00000036940
## rowData names(0):
## colnames(830):
##   GSM1136119_EA07068_260297_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_1.CEL
##   GSM1136120_EA07068_260298_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_2.CEL ...
##   GSM920654_EA07068_201214_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_1.CEL
##   GSM920655_EA07068_201215_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_2.CEL
## colData names(3): label.main label.fine label.ont
```

Each *[celldex](https://bioconductor.org/packages/3.22/celldex)* dataset actually has three sets of labels that primarily differ in their resolution.
For the purposes of this demonstration, we will use the “fine” labels in the `label.fine` metadata field,
which represents the highest resolution of annotation available for this dataset.

```
head(immgen$label.fine)
```

```
## [1] "Macrophages (MF.11C-11B+)" "Macrophages (MF.11C-11B+)"
## [3] "Macrophages (MF.11C-11B+)" "Macrophages (MF.ALV)"     
## [5] "Macrophages (MF.ALV)"      "Macrophages (MF.ALV)"
```

We perform annotation by calling `SingleR()` on our test (Grun) dataset and the reference (ImmGen) dataset,
leaving the default of `de.method="classic"` to use the original marker detection scheme.
This applies the algorithm described in Section [1.2](introduction.html#method-description),
returning a `DataFrame` where each row contains prediction results for a single cell in the `sce` object.
Labels are provided in the `labels` column; some of the other fields are discussed in more detail in Chapter [4](annotation-diagnostics.html#annotation-diagnostics).

```
library(SingleR)

# See 'Choices of assay data' for 'assay.type.test=' explanation.
pred <- SingleR(test = sce, ref = immgen, 
    labels = immgen$label.fine, assay.type.test=1)
colnames(pred)
```

```
## [1] "scores"        "labels"        "delta.next"    "pruned.labels"
```

## 2.3 Interaction with quality control

Upon examining the distribution of assigned labels, we see that many of them are related to stem cells.
However, there are quite a large number of more differentiated labels mixed in,
which is not what we expect from a sorted population of HSCs.

```
head(sort(table(pred$labels), decreasing=TRUE))
```

```
## 
##   Stem cells (SC.MEP) Neutrophils (GN.ARTH)      Macrophages (MF) 
##                   362                   300                   167 
##  Stem cells (SC.STSL)    B cells (proB.FrA) Stem cells (SC.LT34F) 
##                   132                   128                   109
```

This is probably because - despite what its name might suggest -
the dataset obtained by `GrunHSCData()` actually contains more than HSCs.
If we restrict our analysis to the sorted HSCs (obviously) and remove one low-quality batch
(see [the analysis here](https://osca.bioconductor.org/merged-hcsc.html#quality-control-12) for the rationale)
we can see that the distribution of cell type labels is more similar to what we might expect.
Low-quality cells lack information for accurate label assignment and need to be removed to enable interpretation of the results.

```
actual.hsc <- pred$labels[sce$protocol=="sorted hematopoietic stem cells" & sce$sample!="JC4"]
head(sort(table(actual.hsc), decreasing=TRUE))
```

```
## actual.hsc
##       Stem cells (SC.LT34F)        Stem cells (SC.STSL) 
##                         104                         100 
##       Stem cells (SC.ST34F) Stem cells (SC.CD150-CD48-) 
##                          40                          15 
##          Stem cells (LTHSC)            Stem cells (MLP) 
##                          12                           8
```

Filtering the annotation results in the above manner is valid because `SingleR()` operates independently on each test cell.
The annotation is orthogonal to any decisions about the relative quality of the cells in the test dataset;
the same results will be obtained regardless of whether *[SingleR](https://bioconductor.org/packages/3.22/SingleR)* is run before or after quality control.
This is logistically convenient as it means that the annotation does not have to be repeated
if the quality control scheme (or any other downstream step, like clustering) changes throughout the lifetime of the analysis.

## 2.4 Choices of assay data

For the reference dataset, the assay matrix *must* contain log-transformed normalized expression values.
This is because the default marker detection scheme computes log-fold changes by subtracting the medians,
which makes little sense unless the input expression values are already log-transformed.
For alternative schemes, this requirement may be relaxed (e.g., Wilcoxon rank sum tests do not require transformation);
similarly, if pre-defined markers are supplied, no transformation or normalization is necessary.

For the test data, the assay data need not be log-transformed or even (scale) normalized.
This is because `SingleR()` computes Spearman correlations within each cell,
which is unaffected by monotonic transformations like cell-specific scaling or log-transformation.
It is perfectly satisfactory to provide the raw counts for the test dataset to `SingleR()`,
which is the reason for setting `assay.type.test=1` in our previous `SingleR()` call for the Grun dataset.

The exception to this rule occurs when comparing data from full-length technologies to the *[celldex](https://bioconductor.org/packages/3.22/celldex)* references.
These references are intended to be comparable to data from unique molecular identifier (UMI) protocols
where the expression values are less sensitive to differences in gene length.
Thus, when annotating Smart-seq2 test datasets against the *[celldex](https://bioconductor.org/packages/3.22/celldex)* references,
better performance can often be achieved by processing the test counts to transcripts-per-million values.

We demonstrate below using another HSC dataset that was generated using the Smart-seq2 protocol (Nestorowa et al. [2016](#ref-nestorowa2016singlecell)).
Again, we see that most of the predicted labels are related to stem cells, which is comforting.

```
sce.nest <- NestorowaHSCData()

# Getting the exonic gene lengths.
library(AnnotationHub)
mm.db <- AnnotationHub()[["AH73905"]]
mm.exons <- exonsBy(mm.db, by="gene")
mm.exons <- reduce(mm.exons)
mm.len <- sum(width(mm.exons))

# Computing the TPMs with a simple scaling by gene length.
library(scater)
keep <- intersect(names(mm.len), rownames(sce.nest))
tpm.nest <- calculateTPM(sce.nest[keep,], lengths=mm.len[keep])

# Performing the assignment.
pred <- SingleR(test = tpm.nest, ref = immgen, labels = immgen$label.fine)
head(sort(table(pred$labels), decreasing=TRUE), 10)
```

```
## 
##         Stem cells (SC.MEP)       Stem cells (SC.ST34F) 
##                         409                         357 
##      Stem cells (SC.MPP34F)      Stem cells (SC.CMP.DR) 
##                         329                         298 
##            Stem cells (MLP)            Stem cells (GMP) 
##                         167                         102 
##        Stem cells (SC.STSL)         Stem cells (SC.MDP) 
##                          71                          66 
## Stem cells (SC.CD150-CD48-)       Stem cells (SC.LT34F) 
##                          55                          37
```

## 2.5 Comments on choice of references

Unsurprisingly, the choice of reference has a major impact on the annotation results.
We need to pick a reference that contains a superset of the labels that we expect to be present in our test dataset.
Whether the original authors assigned appropriate labels to the reference samples is largely taken as a matter of faith;
it is not entirely unexpected that some references are “better” than others depending on the quality of sample preparation.
We would also prefer a reference that is generated from a similar technology or protocol as our test dataset,
though this is usually not an issue when using `SingleR()` to annotate well-defined cell types.

Users are advised to read the [relevant vignette](https://bioconductor.org/packages/3.22/celldex/vignettes/userguide.html)
for more details about the available references as well as some recommendations on which to use.
(As an aside, the ImmGen dataset and other references were originally supplied along with *[SingleR](https://bioconductor.org/packages/3.22/SingleR)* itself
but have since been migrated to the separate *[celldex](https://bioconductor.org/packages/3.22/celldex)* package for more general use throughout Bioconductor.)
Of course, as we shall see in the next Chapter, it is entirely possible to supply your own reference datasets instead;
all we need are log-expression values and a set of labels for the cells or samples.

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
 [1] scater_1.38.0               ggplot2_4.0.0              
 [3] scuttle_1.20.0              AnnotationHub_4.0.0        
 [5] BiocFileCache_3.0.0         dbplyr_2.5.1               
 [7] SingleR_2.12.0              celldex_1.20.0             
 [9] ensembldb_2.34.0            AnnotationFilter_1.34.0    
[11] GenomicFeatures_1.62.0      AnnotationDbi_1.72.0       
[13] scRNAseq_2.24.0             SingleCellExperiment_1.32.0
[15] SummarizedExperiment_1.40.0 Biobase_2.70.0             
[17] GenomicRanges_1.62.0        Seqinfo_1.0.0              
[19] IRanges_2.44.0              S4Vectors_0.48.0           
[21] BiocGenerics_0.56.0         generics_0.1.4             
[23] MatrixGenerics_1.22.0       matrixStats_1.5.0          
[25] BiocStyle_2.38.0            rebook_1.20.0              

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3        jsonlite_2.0.0           
  [3] CodeDepends_0.6.6         magrittr_2.0.4           
  [5] ggbeeswarm_0.7.2          gypsum_1.6.0             
  [7] farver_2.1.2              rmarkdown_2.30           
  [9] BiocIO_1.20.0             vctrs_0.6.5              
 [11] memoise_2.0.1             Rsamtools_2.26.0         
 [13] DelayedMatrixStats_1.32.0 RCurl_1.98-1.17          
 [15] htmltools_0.5.8.1         S4Arrays_1.10.0          
 [17] curl_7.0.0                BiocNeighbors_2.4.0      
 [19] Rhdf5lib_1.32.0           SparseArray_1.10.1       
 [21] rhdf5_2.54.0              sass_0.4.10              
 [23] alabaster.base_1.10.0     bslib_0.9.0              
 [25] alabaster.sce_1.10.0      httr2_1.2.1              
 [27] cachem_1.1.0              GenomicAlignments_1.46.0 
 [29] lifecycle_1.0.4           pkgconfig_2.0.3          
 [31] rsvd_1.0.5                Matrix_1.7-4             
 [33] R6_2.6.1                  fastmap_1.2.0            
 [35] digest_0.6.37             irlba_2.3.5.1            
 [37] ExperimentHub_3.0.0       RSQLite_2.4.3            
 [39] beachmat_2.26.0           filelock_1.0.3           
 [41] httr_1.4.7                abind_1.4-8              
 [43] compiler_4.5.1            bit64_4.6.0-1            
 [45] withr_3.0.2               S7_0.2.0                 
 [47] BiocParallel_1.44.0       viridis_0.6.5            
 [49] DBI_1.2.3                 HDF5Array_1.38.0         
 [51] alabaster.ranges_1.10.0   alabaster.schemas_1.10.0 
 [53] rappdirs_0.3.3            DelayedArray_0.36.0      
 [55] rjson_0.2.23              tools_4.5.1              
 [57] vipor_0.4.7               beeswarm_0.4.0           
 [59] glue_1.8.0                h5mread_1.2.0            
 [61] restfulr_0.0.16           rhdf5filters_1.22.0      
 [63] grid_4.5.1                gtable_0.3.6             
 [65] BiocSingular_1.26.0       ScaledMatrix_1.18.0      
 [67] XVector_0.50.0            ggrepel_0.9.6            
 [69] BiocVersion_3.22.0        pillar_1.11.1            
 [71] dplyr_1.1.4               lattice_0.22-7           
 [73] rtracklayer_1.70.0        bit_4.6.0                
 [75] tidyselect_1.2.1          Biostrings_2.78.0        
 [77] knitr_1.50                gridExtra_2.3            
 [79] bookdown_0.45             ProtGenerics_1.42.0      
 [81] xfun_0.54                 UCSC.utils_1.6.0         
 [83] lazyeval_0.2.2            yaml_2.3.10              
 [85] evaluate_1.0.5            codetools_0.2-20         
 [87] cigarillo_1.0.0           tibble_3.3.0             
 [89] alabaster.matrix_1.10.0   BiocManager_1.30.26      
 [91] graph_1.88.0              cli_3.6.5                
 [93] jquerylib_0.1.4           dichromat_2.0-0.1        
 [95] Rcpp_1.1.0                GenomeInfoDb_1.46.0      
 [97] dir.expiry_1.18.0         png_0.1-8                
 [99] XML_3.99-0.19             parallel_4.5.1           
[101] blob_1.2.4                sparseMatrixStats_1.22.0 
[103] bitops_1.0-9              viridisLite_0.4.2        
[105] alabaster.se_1.10.0       scales_1.4.0             
[107] purrr_1.1.0               crayon_1.5.3             
[109] rlang_1.1.6               KEGGREST_1.50.0
```

### Bibliography

Aran, D., A. P. Looney, L. Liu, E. Wu, V. Fong, A. Hsu, S. Chak, et al. 2019. “Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.” *Nat. Immunol.* 20 (2): 163–72.

Grun, D., M. J. Muraro, J. C. Boisset, K. Wiebrands, A. Lyubimova, G. Dharmadhikari, M. van den Born, et al. 2016. “De Novo Prediction of Stem Cell Identity using Single-Cell Transcriptome Data.” *Cell Stem Cell* 19 (2): 266–77.

Heng, Tracy S. P., Michio W. Painter, Kutlu Elpek, Veronika Lukacs-Kornek, Nora Mauermann, Shannon J. Turley, Daphne Koller, et al. 2008. “The immunological genome project: Networks of gene expression in immune cells.” *Nature Immunology* 9 (10): 1091–4. <https://doi.org/10.1038/ni1008-1091>.

Nestorowa, S., F. K. Hamey, B. Pijuan Sala, E. Diamanti, M. Shepherd, E. Laurenti, N. K. Wilson, D. G. Kent, and B. Gottgens. 2016. “A single-cell resolution map of mouse hematopoietic stem and progenitor cell differentiation.” *Blood* 128 (8): 20–31.