Source URL: https://bioconductor.org/books/release/SingleRBook/annotation-diagnostics.html
Date Scraped: 2025-12-15

---

# Chapter 4 Annotation diagnostics



## 4.1 Overview

In addition to the labels, `SingleR()` returns a number of helpful diagnostics about the annotation process
that can be used to determine whether the assignments are appropriate.
Unambiguous assignments corroborated by expression of canonical markers add confidence to the results;
conversely, low-confidence assignments can be pruned out to avoid adding noise to downstream analyses.
This chapter will demonstrate some of these common sanity checks on the pancreas datasets
from Chapter [3](sc-mode.html#sc-mode) (Muraro et al. [2016](#ref-muraro2016singlecell); Grun et al. [2016](#ref-grun2016denovo)).

View set-up code (Chapter [8](pancreas-case-study.html#pancreas-case-study))

```
#--- loading-muraro ---#
library(scRNAseq)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label) & sceM$label!="unclear"] 

#--- normalize-muraro ---#
library(scater)
sceM <- logNormCounts(sceM)

#--- loading-grun ---#
sceG <- GrunPancreasData()

sceG <- addPerCellQC(sceG)
qc <- quickPerCellQC(colData(sceG), 
    percent_subsets="altexps_ERCC_percent",
    batch=sceG$donor,
    subset=sceG$donor %in% c("D17", "D7", "D2"))
sceG <- sceG[,!qc$discard]

#--- normalize-grun ---#
sceG <- logNormCounts(sceG)

#--- annotation ---#
library(SingleR)
pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
```

## 4.2 Based on the scores within cells

The most obvious diagnostic reported by `SingleR()` is the nested matrix of per-cell scores in the `scores` field.
This contains the correlation-based scores prior to any fine-tuning for each cell (row) and reference label (column).
Ideally, we would see unambiguous assignments where, for any given cell, one label’s score is clearly larger than the others.

```
pred.grun$scores[1:10,]
```

```
##       acinar  alpha   beta  delta   duct endothelial epsilon mesenchymal     pp
##  [1,] 0.7312 0.1971 0.2247 0.1912 0.5527      0.3596 0.08934      0.2802 0.1896
##  [2,] 0.7022 0.1847 0.1889 0.1613 0.5438      0.3536 0.05597      0.2700 0.1461
##  [3,] 0.6987 0.2091 0.2432 0.2290 0.4947      0.3132 0.10350      0.1968 0.2049
##  [4,] 0.6692 0.2427 0.2727 0.2424 0.7407      0.4023 0.16615      0.3553 0.2206
##  [5,] 0.6364 0.2482 0.2816 0.2463 0.7067      0.4539 0.18858      0.3847 0.2436
##  [6,] 0.6258 0.2902 0.2884 0.2798 0.7068      0.4504 0.19564      0.4021 0.2450
##  [7,] 0.5946 0.3154 0.3256 0.2842 0.6820      0.4620 0.21596      0.4244 0.2879
##  [8,] 0.5901 0.3300 0.3179 0.3052 0.6471      0.4010 0.30567      0.3331 0.2840
##  [9,] 0.6740 0.2827 0.2801 0.2668 0.7401      0.4705 0.18671      0.4119 0.2668
## [10,] 0.7331 0.2663 0.2432 0.2487 0.5849      0.3212 0.18233      0.2288 0.2521
```

To check whether this is indeed the case,
we use the `plotScoreHeatmap()` function to visualize the score matrix (Figure [4.1](annotation-diagnostics.html#fig:score-heatmap-grun)).
Here, the key is to examine the spread of scores within each cell, i.e., down the columns of the heatmap.
Similar scores for a group of labels indicates that the assignment is uncertain for those columns,
though this may be acceptable if the uncertainty is distributed across closely related cell types.
(Note that the assigned label for a cell may not be the visually top-scoring label if fine-tuning is applied,
as the only the pre-tuned scores are directly comparable across all labels.)

```
library(SingleR)
plotScoreHeatmap(pred.grun)
```

![Heatmap of normalized scores for the Grun dataset. Each cell is a column while each row is a label in the reference Muraro dataset. The final label (after fine-tuning) for each cell is shown in the top color bar.](diagnostics_files/figure-html/score-heatmap-grun-1.png)

Figure 4.1: Heatmap of normalized scores for the Grun dataset. Each cell is a column while each row is a label in the reference Muraro dataset. The final label (after fine-tuning) for each cell is shown in the top color bar.

We can also display other metadata information for each cell by setting `clusters=` or `annotation_col=`.
This is occasionally useful for examining potential batch effects,
differences in cell type composition between conditions,
relationship to clusters from an unsupervised analysis and so on,.
For example, Figure [4.2](annotation-diagnostics.html#fig:score-heatmap-grun-donor) displays the donor of origin for each cell;
we can see that each cell type has contributions from multiple donors,
which is reassuring as it indicates that our assignments are not (purely) driven by donor effects.

```
plotScoreHeatmap(pred.grun, 
    annotation_col=as.data.frame(colData(sceG)[,"donor",drop=FALSE]))
```

![Heatmap of normalized scores for the Grun dataset, including the donor of origin for each cell.](diagnostics_files/figure-html/score-heatmap-grun-donor-1.png)

Figure 4.2: Heatmap of normalized scores for the Grun dataset, including the donor of origin for each cell.

The `scores` matrix has several caveats associated with its interpretation.
Only the pre-tuned scores are stored in this matrix, as scores after fine-tuning are not comparable across all labels.
This means that the label with the highest score for a cell may not be the cell’s final label if fine-tuning is applied.
Moreover, the magnitude of differences in the scores has no clear interpretation;
indeed, `plotScoreHeatmap()` dispenses with any faithful representation of the scores
and instead adjusts the values to highlight any differences between labels within each cell.

## 4.3 Based on the deltas across cells

We identify poor-quality or ambiguous assignments based on the per-cell “delta”,
i.e., the difference between the score for the assigned label and the median across all labels for each cell.
Our assumption is that most of the labels in the reference are not relevant to any given cell.
Thus, the median across all labels can be used as a measure of the baseline correlation,
while the gap from the assigned label to this baseline can be used as a measure of the assignment confidence.

Low deltas indicate that the assignment is uncertain, possibly because the cell’s true label does not exist in the reference.
An obvious next step is to apply a threshold on the delta to filter out these low-confidence assignments.
We use the delta rather than the assignment score as the latter is more sensitive to technical effects.
For example, changes in library size affect the technical noise and can increase/decrease all scores for a given cell,
while the delta is somewhat more robust as it focuses on the differences between scores within each cell.

`SingleR()` will set a threshold on the delta for each label using an outlier-based strategy.
Specifically, we identify cells with deltas that are small outliers relative to the deltas of other cells with the same label.
This assumes that, for any given label, most cells assigned to that label are correct.
We focus on outliers to avoid difficulties with setting a fixed threshold,
especially given that the magnitudes of the deltas are about as uninterpretable as the scores themselves.
Pruned labels are reported in the `pruned.labels` field where low-quality assignments are replaced with `NA`.

```
to.remove <- is.na(pred.grun$pruned.labels)
table(Label=pred.grun$labels, Removed=to.remove)
```

```
##              Removed
## Label         FALSE TRUE
##   acinar        260   29
##   alpha         200    1
##   beta          177    1
##   delta          52    2
##   duct          291    4
##   endothelial     5    0
##   epsilon         1    0
##   mesenchymal    22    1
##   pp             18    0
```

However, the default pruning parameters may not be appropriate for every dataset.
For example, if one label is consistently misassigned, the assumption that most cells are correctly assigned will not be appropriate.
In such cases, we can revert to a fixed threshold by manually calling the underlying `pruneScores()` function with `min.diff.med=`.
The example below discards cells with deltas below an arbitrary threshold of 0.2,
where higher thresholds correspond to greater assignment certainty.

```
to.remove <- pruneScores(pred.grun, min.diff.med=0.2)
table(Label=pred.grun$labels, Removed=to.remove)
```

```
##              Removed
## Label         FALSE TRUE
##   acinar        259   30
##   alpha         168   33
##   beta          149   29
##   delta          37   17
##   duct          291    4
##   endothelial     5    0
##   epsilon         1    0
##   mesenchymal    22    1
##   pp              5   13
```

This entire process can be visualized using the `plotScoreDistribution()` function,
which displays the per-label distribution of the deltas across cells (Figure [4.3](annotation-diagnostics.html#fig:score-dist-grun)).
We can use this plot to check that outlier detection in `pruneScores()` behaved sensibly.
Labels with especially low deltas may warrant some additional caution in their interpretation.

```
plotDeltaDistribution(pred.grun)
```

![Distribution of deltas for the Grun dataset. Each facet represents a label in the Muraro dataset, and each point represents a cell assigned to that label (colored by whether it was pruned).](diagnostics_files/figure-html/score-dist-grun-1.png)

Figure 4.3: Distribution of deltas for the Grun dataset. Each facet represents a label in the Muraro dataset, and each point represents a cell assigned to that label (colored by whether it was pruned).

If fine-tuning was performed, we can apply an even more stringent filter
based on the difference between the highest and second-highest scores after fine-tuning.
Cells will only pass the filter if they are assigned to a label that is clearly distinguishable from any other label.
In practice, this approach tends to be too conservative as assignments involving closely related labels are heavily penalized.

```
to.remove2 <- pruneScores(pred.grun, min.diff.next=0.1)
table(Label=pred.grun$labels, Removed=to.remove2)
```

```
##              Removed
## Label         FALSE TRUE
##   acinar        242   47
##   alpha         167   34
##   beta          107   71
##   delta          31   23
##   duct          173  122
##   endothelial     5    0
##   epsilon         1    0
##   mesenchymal    22    1
##   pp              6   12
```

## 4.4 Based on marker gene expression

Another simple yet effective diagnostic is to examine the expression of the marker genes for each label in the test dataset.
The marker genes used for each label are reported in the `metadata()` of the `SingleR()` output, so we can simply retrieve them to visualize their (usually log-transformed) expression values across the test dataset.
This is done automatically by the `plotMarkerHeatmap()` function, which visualizes marker expression for a particular label (Figure [4.4](annotation-diagnostics.html#fig:grun-beta-heat)).
To avoid showing too many genes, this function only focuses on the most relevant markers,
i.e., those that are upregulated in the test dataset for the label of interest and thus are responsible for driving the classification of cells to that label.

```
plotMarkerHeatmap(pred.grun, sceG, "beta")
```

![Heatmap of log-expression values in the Grun dataset for the top marker genes upregulated in beta cells in the Muraro reference dataset. Assigned labels for each cell are shown at the top of the plot.](diagnostics_files/figure-html/grun-beta-heat-1.png)

Figure 4.4: Heatmap of log-expression values in the Grun dataset for the top marker genes upregulated in beta cells in the Muraro reference dataset. Assigned labels for each cell are shown at the top of the plot.

If a cell in the test dataset is confidently assigned to a particular label,
we would expect it to have strong expression of that label’s markers.
We would also hope that those label’s markers are biologically meaningful;
in this case, we do observe strong upregulation of insulin (*INS*) in the beta cells,
which is reassuring and gives greater confidence to the correctness of the assignment.
If the identified markers are not meaningful or not consistently upregulated,
some skepticism towards the quality of the assignments is warranted.

We can easily create a diagnostic plot for each label by wrapping the above code in a loop (Figure [4.5](annotation-diagnostics.html#fig:grun-beta-heat-all)).
This allows us to quickly visualize the quality of each label’s assignments.

```
collected <- list()
for (lab in unique(pred.grun$labels)) {
    collected[[lab]] <- plotMarkerHeatmap(pred.grun, sceG, lab,
        main=lab, silent=TRUE)[[4]]
}
do.call(gridExtra::grid.arrange, collected)
```

![Heatmaps of log-expression values in the Grun dataset for the top marker genes upregulated in each label in the Muraro reference dataset. Assigned labels for each cell are shown at the top of each plot.](diagnostics_files/figure-html/grun-beta-heat-all-1.png)

Figure 4.5: Heatmaps of log-expression values in the Grun dataset for the top marker genes upregulated in each label in the Muraro reference dataset. Assigned labels for each cell are shown at the top of each plot.

Users can also customize the visualization by re-using the heatmap configuration from `configureMarkerHeatmap()` with other plotting functions,
e.g., `plotDots()` from *[scater](https://bioconductor.org/packages/3.22/scater)* (to create a *[Seurat](https://CRAN.R-project.org/package=Seurat)*-style dot plot) or `dittoHeatmap()` from *[dittoSeq](https://bioconductor.org/packages/3.22/dittoSeq)*.
For example, Figure [4.6](annotation-diagnostics.html#fig:grun-epsilon-heat) creates a per-label heatmap to ensure that the visualization is not dominated by the most abundant labels.

```
config <- configureMarkerHeatmap(pred.grun, sceG, "epsilon")
mat <- assay(sceG, "logcounts")[head(config$rows, 20), config$columns]
aggregated <- scuttle::summarizeAssayByGroup(mat, config$predictions)
pheatmap::pheatmap(assay(aggregated), cluster_col=FALSE)
```

![Heatmap of log-expression values in the Grun dataset for the top marker genes upregulated in epsilon cells in the Muraro reference dataset. Assigned labels for each cell are shown at the top of the plot.](diagnostics_files/figure-html/grun-epsilon-heat-1.png)

Figure 4.6: Heatmap of log-expression values in the Grun dataset for the top marker genes upregulated in epsilon cells in the Muraro reference dataset. Assigned labels for each cell are shown at the top of the plot.

In general, the marker expression heatmap provides a more interpretable diagnostic visualization than the plots of scores and deltas.
However, it does require more effort to inspect and may not be feasible for large numbers of labels.
It is also difficult to use a heatmap to determine the correctness of assignment for closely related labels.

## 4.5 Comparing to unsupervised clustering

It can also be instructive to compare the assigned labels to the groupings generated from unsupervised clustering algorithms.
The assumption is that the differences between reference labels are also the dominant factor of variation in the test dataset;
this implies that we should expect strong agreement between the clusters and the assigned labels.
To demonstrate, we’ll use the `sceG` from Chapter [8](pancreas-case-study.html#pancreas-case-study)
where clusters have generated using a graph-based method (Xu and Su [2015](#ref-xu2015identification)).

View set-up code (Chapter [8](pancreas-case-study.html#pancreas-case-study))

```
#--- loading-muraro ---#
library(scRNAseq)
sceM <- MuraroPancreasData()
sceM <- sceM[,!is.na(sceM$label) & sceM$label!="unclear"] 

#--- normalize-muraro ---#
library(scater)
sceM <- logNormCounts(sceM)

#--- loading-grun ---#
sceG <- GrunPancreasData()

sceG <- addPerCellQC(sceG)
qc <- quickPerCellQC(colData(sceG), 
    percent_subsets="altexps_ERCC_percent",
    batch=sceG$donor,
    subset=sceG$donor %in% c("D17", "D7", "D2"))
sceG <- sceG[,!qc$discard]

#--- normalize-grun ---#
sceG <- logNormCounts(sceG)

#--- annotation ---#
library(SingleR)
pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")

#--- clustering ---#
library(scran)
decG <- modelGeneVarWithSpikes(sceG, "ERCC")

set.seed(1000100)
sceG <- denoisePCA(sceG, decG)

library(bluster)
sceG$cluster <- clusterRows(reducedDim(sceG), NNGraphParam(k=5))
```

We compare these clusters to the labels generated by *[SingleR](https://bioconductor.org/packages/3.22/SingleR)*.
Any similarity can be quantified with the adjusted rand index (ARI) with `pairwiseRand()` from the *[bluster](https://bioconductor.org/packages/3.22/bluster)* package.
Large ARIs indicate that the two partitionings are in agreement, though an acceptable definition of “large” is difficult to gauge;
experience suggests that a reasonable level of consistency is achieved at ARIs above 0.5.

```
library(bluster)
pairwiseRand(sceG$cluster, pred.grun$labels, mode="index")
```

```
## [1] 0.3881
```

In practice, it is more informative to examine the distribution of cells across each cluster/label combination.
Figure [4.7](annotation-diagnostics.html#fig:grun-label-clusters) shows that most clusters are nested within labels,
a difference in resolution that is likely responsible for reducing the ARI.
Clusters containing multiple labels are particularly interesting for diagnostic purposes,
as this suggests that the differences between labels are not strong enough to drive formation of distinct clusters in the test.

```
tab <- table(cluster=sceG$cluster, label=pred.grun$labels) 
pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing.
```

![Heatmap of the log-transformed number of cells in each combination of label (column) and cluster (row) in the Grun dataset.](diagnostics_files/figure-html/grun-label-clusters-1.png)

Figure 4.7: Heatmap of the log-transformed number of cells in each combination of label (column) and cluster (row) in the Grun dataset.

The underlying assumption is somewhat reasonable in most scenarios where the labels relate to cell type identity.
However, disagreements between the clusters and labels should not be cause for much concern.
The whole point of unsupervised clustering is to identify novel variation that, by definition, is not in the reference.
It is entirely possible for the clustering and labels to be different without compromising the validity or utility of either;
the former captures new heterogeneity while the latter facilitates interpretation in the context of existing knowledge.

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
 [1] bluster_1.20.0              SingleR_2.12.0             
 [3] SummarizedExperiment_1.40.0 Biobase_2.70.0             
 [5] GenomicRanges_1.62.0        Seqinfo_1.0.0              
 [7] IRanges_2.44.0              S4Vectors_0.48.0           
 [9] BiocGenerics_0.56.0         generics_0.1.4             
[11] MatrixGenerics_1.22.0       matrixStats_1.5.0          
[13] BiocStyle_2.38.0            rebook_1.20.0              

loaded via a namespace (and not attached):
 [1] gtable_0.3.6                dir.expiry_1.18.0          
 [3] xfun_0.54                   bslib_0.9.0                
 [5] ggplot2_4.0.0               lattice_0.22-7             
 [7] scrapper_1.4.0              vctrs_0.6.5                
 [9] tools_4.5.1                 parallel_4.5.1             
[11] tibble_3.3.0                cluster_2.1.8.1            
[13] pkgconfig_2.0.3             BiocNeighbors_2.4.0        
[15] pheatmap_1.0.13             Matrix_1.7-4               
[17] RColorBrewer_1.1-3          S7_0.2.0                   
[19] sparseMatrixStats_1.22.0    graph_1.88.0               
[21] lifecycle_1.0.4             compiler_4.5.1             
[23] farver_2.1.2                codetools_0.2-20           
[25] htmltools_0.5.8.1           sass_0.4.10                
[27] yaml_2.3.10                 pillar_1.11.1              
[29] jquerylib_0.1.4             SingleCellExperiment_1.32.0
[31] BiocParallel_1.44.0         DelayedArray_0.36.0        
[33] cachem_1.1.0                viridis_0.6.5              
[35] abind_1.4-8                 tidyselect_1.2.1           
[37] digest_0.6.37               dplyr_1.1.4                
[39] bookdown_0.45               labeling_0.4.3             
[41] fastmap_1.2.0               grid_4.5.1                 
[43] cli_3.6.5                   SparseArray_1.10.1         
[45] magrittr_2.0.4              S4Arrays_1.10.0            
[47] dichromat_2.0-0.1           XML_3.99-0.19              
[49] withr_3.0.2                 DelayedMatrixStats_1.32.0  
[51] filelock_1.0.3              CodeDepends_0.6.6          
[53] scales_1.4.0                rmarkdown_2.30             
[55] XVector_0.50.0              igraph_2.2.1               
[57] gridExtra_2.3               beachmat_2.26.0            
[59] evaluate_1.0.5              knitr_1.50                 
[61] viridisLite_0.4.2           rlang_1.1.6                
[63] Rcpp_1.1.0                  scuttle_1.20.0             
[65] glue_1.8.0                  BiocManager_1.30.26        
[67] jsonlite_2.0.0              R6_2.6.1
```

### Bibliography

Grun, D., M. J. Muraro, J. C. Boisset, K. Wiebrands, A. Lyubimova, G. Dharmadhikari, M. van den Born, et al. 2016. “De Novo Prediction of Stem Cell Identity using Single-Cell Transcriptome Data.” *Cell Stem Cell* 19 (2): 266–77.

Muraro, M. J., G. Dharmadhikari, D. Grun, N. Groen, T. Dielen, E. Jansen, L. van Gurp, et al. 2016. “A Single-Cell Transcriptome Atlas of the Human Pancreas.” *Cell Syst* 3 (4): 385–94.

Xu, C., and Z. Su. 2015. “Identification of cell types from single-cell transcriptomes using a novel clustering method.” *Bioinformatics* 31 (12): 1974–80.