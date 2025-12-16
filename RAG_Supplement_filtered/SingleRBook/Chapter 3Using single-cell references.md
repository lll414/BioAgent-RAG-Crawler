Source URL: https://bioconductor.org/books/release/SingleRBook/sc-mode.html
Date Scraped: 2025-12-15

---

# Chapter 3 Using single-cell references



## 3.1 Overview

*[SingleR](https://bioconductor.org/packages/3.22/SingleR)* can also be applied to reference datasets derived from single-cell RNA-seq experiments.
This involves the same algorithm as that used in the classic mode (Chapter [2](classic-mode.html#classic-mode))
but performs marker detection with conventional statistical tests instead of the log-fold change.
In particular, we identify top-ranked markers based on pairwise Wilcoxon rank sum tests or t-tests between labels;
this allows us to account for the variability across cells to choose genes that are robustly upregulated in each label.
Users can also supply their own custom marker lists to `SingleR()`,
facilitating incorporation of prior biological knowledge into the annotation process.
We will demonstrate these capabilities below in this chapter.

## 3.2 Annotation with test-based marker detection

To demonstrate, we will use two human pancreas scRNA-seq datasets from the *[scRNAseq](https://bioconductor.org/packages/3.22/scRNAseq)* package.
The aim is to use one pre-labelled dataset to annotate the other unlabelled dataset.
First, we set up the Muraro et al. ([2016](#ref-muraro2016singlecell)) dataset to be our reference,
computing log-normalized expression values as discussed in Section [2.4](classic-mode.html#choices-of-assay-data).

```
library(scRNAseq)
sceM <- MuraroPancreasData()

# Removing unlabelled cells or cells without a clear label.
sceM <- sceM[,!is.na(sceM$label) & sceM$label!="unclear"] 

library(scater)
sceM <- logNormCounts(sceM)
sceM
```

```
## class: SingleCellExperiment 
## dim: 19059 2122 
## metadata(0):
## assays(2): counts logcounts
## rownames(19059): A1BG-AS1__chr19 A1BG__chr19 ... ZZEF1__chr17
##   ZZZ3__chr1
## rowData names(2): symbol chr
## colnames(2122): D28-1_1 D28-1_2 ... D30-8_93 D30-8_94
## colData names(4): label donor plate sizeFactor
## reducedDimNames(0):
## mainExpName: endogenous
## altExpNames(1): ERCC
```

```
# Seeing the available labels in this dataset.
table(sceM$label)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         219         812         448         193         245          21 
##     epsilon mesenchymal          pp 
##           3          80         101
```

We then set up our test dataset from Grun et al. ([2016](#ref-grun2016denovo)), applying some basic quality control as discusssed
[here](http://bioconductor.org/books/release/OSCA.workflows/grun-human-pancreas-cel-seq2.html#quality-control-4)
and in Section [2.3](classic-mode.html#interaction-with-quality-control).
We also compute the log-transformed values here, not because it is strictly necessary
but so that we don’t have to keep on typing `assay.type.test=1` in later calls to `SingleR()`.

```
sceG <- GrunPancreasData()

sceG <- addPerCellQC(sceG)
qc <- quickPerCellQC(colData(sceG), 
    percent_subsets="altexps_ERCC_percent",
    batch=sceG$donor,
    subset=sceG$donor %in% c("D17", "D7", "D2"))
sceG <- sceG[,!qc$discard]

sceG <- logNormCounts(sceG)
sceG
```

```
## class: SingleCellExperiment 
## dim: 20064 1064 
## metadata(0):
## assays(2): counts logcounts
## rownames(20064): A1BG-AS1__chr19 A1BG__chr19 ... ZZEF1__chr17
##   ZZZ3__chr1
## rowData names(2): symbol chr
## colnames(1064): D2ex_1 D2ex_2 ... D17TGFB_94 D17TGFB_95
## colData names(9): donor sample ... total sizeFactor
## reducedDimNames(0):
## mainExpName: endogenous
## altExpNames(1): ERCC
```

We run `SingleR()` as described previously but with a marker detection mode that considers the variance of expression across cells.
Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels.
This is slower but more appropriate for single-cell data compared to the default marker detection algorithm,
as the latter may fail for low-coverage data where the median for each label is often zero.

```
library(SingleR)
pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
table(pred.grun$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         289         201         178          54         295           5 
##     epsilon mesenchymal          pp 
##           1          23          18
```

By default, the function will take the top `de.n` (default: 10) genes from each pairwise comparison between labels.
A larger number of markers increases the robustness of the annotation by ensuring that relevant genes are not omitted,
especially if the reference dataset has study-specific effects that cause uninteresting genes to dominate the top set.
However, this comes at the cost of increasing noise and computational time.

```
library(SingleR)
pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, 
    de.method="wilcox", de.n=50)
table(pred.grun$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         274         201         176          56         309           5 
##     epsilon mesenchymal          pp 
##           1          23          19
```

## 3.3 Defining custom markers

The marker detection in `SingleR()` is based on *[scrapper](https://bioconductor.org/packages/3.22/scrapper)*,
so most options in `?scoreMarkers` and friends can be applied via the `de.args=` option.
For example, we could apply a log-fold change threshold with `de.args=list(threshold=1)`.

```
library(SingleR)
pred.grun2 <- SingleR(test=sceG, ref=sceM, labels=sceM$label, 
    de.method="t", de.args=list(threshold=1))
table(pred.grun2$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         286         201         179          53         296           5 
##     epsilon mesenchymal          pp 
##           1          24          19
```

However, users can also construct their own marker lists with any DE testing machinery.
For example, we can perform pairwise binomial tests to identify genes that are differentially detected
(i.e., have differences in the proportion of cells with non-zero counts) between labels in the reference Muraro dataset.
We then take the top 10 marker genes from each pairwise comparison,
obtaining a list of lists of character vectors containing the identities of the markers for that comparison.

```
library(scran)
out <- pairwiseBinom(counts(sceM), sceM$label, direction="up")
markers <- getTopMarkers(out$statistics, out$pairs, n=10)

# Upregulated in acinar compared to alpha:
markers$acinar$alpha
```

```
##  [1] "KCNQ1__chr11"  "FAM129A__chr1" "KLK1__chr19"   "NTN4__chr12"  
##  [5] "RASEF__chr9"   "CTRL__chr16"   "LGALS2__chr22" "NUPR1__chr16" 
##  [9] "LGALS3__chr14" "NR5A2__chr1"
```

```
# Upregulated in alpha compared to acinar:
markers$alpha$acinar
```

```
##  [1] "SLC38A4__chr12" "ARX__chrX"      "CRYBA2__chr2"   "FSTL5__chr4"   
##  [5] "GNG2__chr14"    "NOL4__chr18"    "IRX2__chr5"     "KCNMB2__chr3"  
##  [9] "CFC1__chr2"     "KCNJ6__chr21"
```

Once we have this list of lists, we supply it to `SingleR()` via the `genes=` argument,
which causes the function to bypass the internal marker detection to use the supplied gene sets instead.
The most obvious benefit of this approach is that the user can achieve greater control of the markers,
allowing integration of prior biological knowledge to obtain more relevant genes and a more robust annotation.

```
pred.grun2b <- SingleR(test=sceG, ref=sceM, labels=sceM$label, genes=markers)
table(pred.grun2b$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         276         202         175          54         302           6 
##     epsilon mesenchymal          pp 
##           2          24          23
```

In some cases, markers may only be available for specific labels rather than for pairwise comparisons between labels.
This is accommodated by supplying a named list of character vectors to `genes`.
Note that this is likely to be less powerful than the list-of-lists approach as information about pairwise differences is discarded.

```
# Creating label-specific markers.
label.markers <- lapply(markers, unlist)
label.markers <- lapply(label.markers, unique)
str(label.markers)
```

```
## List of 9
##  $ acinar     : chr [1:40] "KCNQ1__chr11" "FAM129A__chr1" "KLK1__chr19" "NTN4__chr12" ...
##  $ alpha      : chr [1:41] "SLC38A4__chr12" "ARX__chrX" "CRYBA2__chr2" "FSTL5__chr4" ...
##  $ beta       : chr [1:47] "ELAVL4__chr1" "PRUNE2__chr9" "NMNAT2__chr1" "PLCB4__chr20" ...
##  $ delta      : chr [1:44] "NOL4__chr18" "CABP7__chr22" "UNC80__chr2" "HEPACAM2__chr7" ...
##  $ duct       : chr [1:50] "ADCY5__chr3" "PDE3A__chr12" "SLC3A1__chr2" "BICC1__chr10" ...
##  $ endothelial: chr [1:26] "GPR4__chr19" "TMEM204__chr16" "GPR116__chr6" "CYYR1__chr21" ...
##  $ epsilon    : chr [1:14] "BHMT__chr5" "JPH3__chr16" "SERPINA10__chr14" "UGT2B4__chr4" ...
##  $ mesenchymal: chr [1:34] "TNFAIP6__chr2" "THBS2__chr6" "CDH11__chr16" "SRPX2__chrX" ...
##  $ pp         : chr [1:44] "SERTM1__chr13" "ETV1__chr7" "ARX__chrX" "ELAVL4__chr1" ...
```

```
pred.grun2c <- SingleR(test=sceG, ref=sceM, labels=sceM$label, genes=label.markers)
table(pred.grun2c$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         262         204         169          59         317           6 
##     epsilon mesenchymal          pp 
##           2          24          21
```

## 3.4 Pseudo-bulk aggregation

Single-cell reference datasets provide a like-for-like comparison to our test single-cell datasets,
yielding a more accurate classification of the cells in the latter (hopefully).
However, there are frequently many more samples in single-cell references compared to bulk references,
increasing the computational work involved in classification.
We overcome this by aggregating cells into one “pseudo-bulk” sample per label (e.g., by averaging across log-expression values)
and using that as the reference profile, which allows us to achieve the same efficiency as the use of bulk references.

The obvious cost of this approach is that we discard potentially useful information
about the distribution of cells within each label.
Cells that belong to a heterogeneous population may not be correctly assigned if they are far from the population center.
To preserve some of this information, we perform k-means clustering within each label
to create pseudo-bulk samples that are representative of a particular region of the expression space (i.e., vector quantization).
We create √N clusters given a label with N cells,
which provides a reasonable compromise between reducing computational work and preserving the label’s internal distribution.

To enable this aggregation, we simply set `aggr.ref=TRUE` in the `SingleR()` call.
This uses the `aggregateReference()` function to perform k-means clustering *within* each label
(typically after principal components analysis on the log-expression matrix, for greater speed)
and average expression values for each within-label cluster.
Note that marker detection is still performed on the unaggregated data
so as to make full use of the distribution of expression values across cells.

```
set.seed(100) # for the k-means step.
pred.grun3 <- SingleR(test=sceG, ref=sceM, labels=sceM$label, 
    de.method="wilcox", aggr.ref=TRUE)
table(pred.grun3$labels)
```

```
## 
##      acinar       alpha        beta       delta        duct endothelial 
##         285         201         179          53         298           5 
##     epsilon mesenchymal          pp 
##           1          23          19
```

Advanced users can also perform the aggregation manually by calling the `aggregateReference()` function.
This returns a `SummarizedExperiment` object for use as `ref=` in the `SingleR()` function.

```
set.seed(100)
aggregated <- aggregateReference(sceM, sceM$label)
aggregated
```

```
## class: SummarizedExperiment 
## dim: 19059 114 
## metadata(0):
## assays(1): logcounts
## rownames(19059): A1BG-AS1__chr19 A1BG__chr19 ... ZZEF1__chr17
##   ZZZ3__chr1
## rowData names(0):
## colnames(114): acinar.1 acinar.2 ... pp.9 pp.10
## colData names(1): label
```

Obviously, the aggregation itself requires computational work so setting `aggr.ref=TRUE` in `SingleR()` itself may not improve speed.
Rather, the real power of this approach lies in pre-aggregating the reference dataset
so that it can be repeatedly applied to quickly annotate multiple test datasets.
This approach is discussed in more detail in Chapter [7](advanced-options.html#advanced-options).

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
 [1] scran_1.38.0                SingleR_2.12.0             
 [3] scater_1.38.0               ggplot2_4.0.0              
 [5] scuttle_1.20.0              scRNAseq_2.24.0            
 [7] SingleCellExperiment_1.32.0 SummarizedExperiment_1.40.0
 [9] Biobase_2.70.0              GenomicRanges_1.62.0       
[11] Seqinfo_1.0.0               IRanges_2.44.0             
[13] S4Vectors_0.48.0            BiocGenerics_0.56.0        
[15] generics_0.1.4              MatrixGenerics_1.22.0      
[17] matrixStats_1.5.0           BiocStyle_2.38.0           
[19] rebook_1.20.0              

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3        jsonlite_2.0.0           
  [3] CodeDepends_0.6.6         magrittr_2.0.4           
  [5] ggbeeswarm_0.7.2          GenomicFeatures_1.62.0   
  [7] gypsum_1.6.0              farver_2.1.2             
  [9] rmarkdown_2.30            BiocIO_1.20.0            
 [11] vctrs_0.6.5               DelayedMatrixStats_1.32.0
 [13] memoise_2.0.1             Rsamtools_2.26.0         
 [15] RCurl_1.98-1.17           htmltools_0.5.8.1        
 [17] S4Arrays_1.10.0           AnnotationHub_4.0.0      
 [19] curl_7.0.0                BiocNeighbors_2.4.0      
 [21] Rhdf5lib_1.32.0           SparseArray_1.10.1       
 [23] rhdf5_2.54.0              sass_0.4.10              
 [25] alabaster.base_1.10.0     bslib_0.9.0              
 [27] alabaster.sce_1.10.0      httr2_1.2.1              
 [29] cachem_1.1.0              GenomicAlignments_1.46.0 
 [31] igraph_2.2.1              lifecycle_1.0.4          
 [33] pkgconfig_2.0.3           rsvd_1.0.5               
 [35] Matrix_1.7-4              R6_2.6.1                 
 [37] fastmap_1.2.0             digest_0.6.37            
 [39] AnnotationDbi_1.72.0      dqrng_0.4.1              
 [41] irlba_2.3.5.1             ExperimentHub_3.0.0      
 [43] RSQLite_2.4.3             beachmat_2.26.0          
 [45] filelock_1.0.3            httr_1.4.7               
 [47] abind_1.4-8               compiler_4.5.1           
 [49] bit64_4.6.0-1             withr_3.0.2              
 [51] S7_0.2.0                  BiocParallel_1.44.0      
 [53] viridis_0.6.5             DBI_1.2.3                
 [55] HDF5Array_1.38.0          alabaster.ranges_1.10.0  
 [57] alabaster.schemas_1.10.0  rappdirs_0.3.3           
 [59] DelayedArray_0.36.0       bluster_1.20.0           
 [61] rjson_0.2.23              tools_4.5.1              
 [63] vipor_0.4.7               beeswarm_0.4.0           
 [65] glue_1.8.0                h5mread_1.2.0            
 [67] restfulr_0.0.16           rhdf5filters_1.22.0      
 [69] grid_4.5.1                cluster_2.1.8.1          
 [71] gtable_0.3.6              ensembldb_2.34.0         
 [73] metapod_1.18.0            BiocSingular_1.26.0      
 [75] ScaledMatrix_1.18.0       XVector_0.50.0           
 [77] ggrepel_0.9.6             BiocVersion_3.22.0       
 [79] pillar_1.11.1             limma_3.66.0             
 [81] dplyr_1.1.4               BiocFileCache_3.0.0      
 [83] lattice_0.22-7            rtracklayer_1.70.0       
 [85] bit_4.6.0                 tidyselect_1.2.1         
 [87] locfit_1.5-9.12           Biostrings_2.78.0        
 [89] knitr_1.50                gridExtra_2.3            
 [91] scrapper_1.4.0            bookdown_0.45            
 [93] ProtGenerics_1.42.0       edgeR_4.8.0              
 [95] xfun_0.54                 statmod_1.5.1            
 [97] UCSC.utils_1.6.0          lazyeval_0.2.2           
 [99] yaml_2.3.10               evaluate_1.0.5           
[101] codetools_0.2-20          cigarillo_1.0.0          
[103] tibble_3.3.0              alabaster.matrix_1.10.0  
[105] BiocManager_1.30.26       graph_1.88.0             
[107] cli_3.6.5                 jquerylib_0.1.4          
[109] dichromat_2.0-0.1         Rcpp_1.1.0               
[111] GenomeInfoDb_1.46.0       dir.expiry_1.18.0        
[113] dbplyr_2.5.1              png_0.1-8                
[115] XML_3.99-0.19             parallel_4.5.1           
[117] blob_1.2.4                AnnotationFilter_1.34.0  
[119] sparseMatrixStats_1.22.0  bitops_1.0-9             
[121] viridisLite_0.4.2         alabaster.se_1.10.0      
[123] scales_1.4.0              crayon_1.5.3             
[125] rlang_1.1.6               KEGGREST_1.50.0
```

### Bibliography

Grun, D., M. J. Muraro, J. C. Boisset, K. Wiebrands, A. Lyubimova, G. Dharmadhikari, M. van den Born, et al. 2016. “De Novo Prediction of Stem Cell Identity using Single-Cell Transcriptome Data.” *Cell Stem Cell* 19 (2): 266–77.

Muraro, M. J., G. Dharmadhikari, D. Grun, N. Groen, T. Dielen, E. Jansen, L. van Gurp, et al. 2016. “A Single-Cell Transcriptome Atlas of the Human Pancreas.” *Cell Syst* 3 (4): 385–94.