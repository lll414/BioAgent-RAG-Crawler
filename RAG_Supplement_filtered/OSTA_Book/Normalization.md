Source URL: https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html
Date Scraped: 2025-12-15

---

1. [Platform-independent analyses](https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html)
2. [26  Normalization](https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html)

# 26  Normalization

## 26.1 Preamble

### 26.1.1 Introduction

Normalization is a critical step in the analysis of transcriptomics and other high-throughput biological data, aiming to remove technical variability while preserving true biological signals.

In single-cell RNA sequencing, normalization aims at removing differences in sequencing coverage between libraries (cells). Such differences in *library size* are typically accounted for by scaling the counts in each cell by a size factor, which is often computed as the total number of counts for that cell (or a robust variant thereof); see [OSCA](https://bioconductor.org/books/release/OSCA.basic/normalization.html).

Single-cell normalization methods have been widely used for the analysis of spatial transcriptomics (ST) data, given the similarity between the two technologies, which both yield count data. While this is justifiable for sequencing-based ST technologies, it is less clear that this is appropriate for imaging-based ones.

In fact, in imaging-based ST, the total number of counts is driven by the number of probes that successfully hybridize to their target transcripts within each cell and not by the sequencing coverage. It is hence less clear that this leads to the same compositional biases as in sequencing platforms. Moreover, these technologies often employ targeted gene panels rather than whole-transcriptome profiling, which may lead to confounding between total counts and biology, as shown by Atta et al. ([2024](https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html#ref-Atta2024-gene-count)) and Bhuva et al. ([2024](https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html#ref-Bhuva2024-library-size)).

Here, we will review both scaling normalization methods and spatially-aware normalization methods that have been specifically developed for ST data.

### 26.1.2 Dependencies

Code

```
library(scater)
library(ggplot2)
library(patchwork)
library(OSTA.data)
library(SpatialExperiment)
# load data from preceding 
# chapter (post quality control)
(spe <- readRDS("img-spe_qc.rds"))
```

```
##  class: SpatialExperiment 
##  dim: 313 140268 
##  metadata(0):
##  assays(1): counts
##  rownames(313): ABCC11 ACTA2 ... ZEB2 ZNF562
##  rowData names(3): ID Symbol Type
##  colnames(140268): 2 3 ... 167779 167780
##  colData names(12): cell_id transcript_counts ... detected total
##  reducedDimNames(0):
##  mainExpName: NULL
##  altExpNames(0):
##  spatialCoords names(2) : x_centroid y_centroid
##  imgData names(1): sample_id
```

Code

```
# get annotations from 'BiocFileCache'
# (data has been retrieved already)
id <- "Xenium_HumanBreast1_Janesick"
pa <- OSTA.data_load(id, mol=FALSE)
dir.create(td <- tempfile())
unzip(pa, "annotation.csv", exdir=td)
df <- read.csv(list.files(td, full.names=TRUE))
# add annotations as cell metadata
cs <- match(spe$cell_id, df$Barcode)
spe$Label <- df$Annotation[cs]
```

## 26.2 Scaling normalization

As a baseline, we will apply log-library size normalization, as implemented in the *[scater](https://bioconductor.org/packages/3.22/scater)* package ([McCarthy et al. 2017](https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html#ref-McCarthy2017-scater)). Note that the `logNormCounts` function adds a new assay named `logcounts` to the `spe` object, which contains the normalized expression values on the log scale.

Code

```
spe$library_size <- librarySizeFactors(spe)
spe <- logNormCounts(spe, size.factors=spe$library_size)
assayNames(spe)
```

```
##  [1] "counts"    "logcounts"
```

An alternative that has been suggested is normalization by cell area (or volume). In principle, this is equivalent to assuming a constant transcription rate across cells and focusing on transcript density, rather than transcript counts, as a measure of gene expression.

Nevertheless, it is important to highlight that cell area itself is not known a priori, but rather estimated through segmentation algorithms that may introduce errors, which can impact the final results.

To compute area-derived size factors, we can use the `cell_area` column in the `colData` of the `spe` object, which contains the area of each cell as estimated by the segmentation.

Code

```
sizeFactors(spe) <- spe$cell_area / median(spe$cell_area)
spe <- logNormCounts(spe, name="normalized_by_area")
```

We can check the distribution of library and area-derived size factors across cell types, to highlight any potential confounding:

Code

```
df <- data.frame(colData(spe), spatialCoords(spe))
p1 <- ggplot(df, aes(y=library_size)) + ggtitle("Library size factors") 
p2 <- ggplot(df, aes(y=sizeFactor)) + ggtitle("Area-derived factors") 
(p1 + p2) &
    geom_boxplot(aes(Label, fill=Label)) &
    labs(x="Cell type", y="Size factor") &
    scale_fill_manual(values=unname(pals::tableau20())) &
    theme_bw() & theme(
        aspect.ratio=1,
        legend.position="none",
        panel.grid.minor=element_blank(),
        plot.title=element_text(hjust=0.5),
        axis.text.x=element_text(angle=45, hjust=1))
```

![](ind-normalization_files/figure-html/plot-size-factors-boxplot-1.png)

We can clearly see that library size factors of tumor cells are systematically higher than those of other cell types, which is a sign of confounding between library size and biology.

We can also compare the area-derived scaling factors with the library size factors:

Code

```
ggplot(df, aes(sizeFactor, library_size)) + 
    geom_point(shape=16, stroke=0, size=0.4) +
    geom_abline(linewidth=0.4, col="blue") + 
    facet_wrap(~Label) +
    labs(x="Area-derived factor", y="Library size factor") + 
    ggtitle("Relation between area-derived and library size factors") +
    coord_equal() + theme_bw() + theme(
        panel.grid.minor=element_blank(),
        plot.title=element_text(hjust=0.5))
```

![](ind-normalization_files/figure-html/plot-size-factors-scatter-1.png)

While there is a general correlation between total counts and cell area, the relationship is cell-type specific, with tumor cells showing systematically higher counts for a given area, and stromal cells showing generally larger areas.

This suggests that the choice of normalization may have a significant impact on downstream analyses, and should be carefully considered in the context of the specific dataset and biological questions at hand.

We next look at the difference in normalized expression values between the two normalization methods, and select the top 10 genes that show the largest differences.

Code

```
mean_libs <- rowMeans(assays(spe)$logcounts)
mean_area <- rowMeans(assays(spe)$normalized_by_area)
diff <- abs(mean_libs - mean_area)
ord <- order(diff, decreasing=TRUE)
top_diff <- diff[ord[1:10]]
names(top_diff)
```

```
##   [1] "EPCAM"   "KRT8"    "TACSTD2" "KRT7"    "GATA3"   "CD9"     "MLPH"   
##   [8] "FASN"    "CDH1"    "FOXA1"
```

By looking at the top genes, we can see several epithelial cell markers, such as EPCAM, KRT8, and KRT7.

Interestingly, looking at the expression of EPCAM across the tissue, we can see that area normalization enhances the contrast between tumor and non-tumor regions, while library size normalization yields a more diffuse expression pattern.

Code

```
cd <- data.frame(colData(spe), spatialCoords(spe))

mx <- logcounts(spe)[ord[1:10], ]
df_ls <- cbind(cd, as.matrix(t(mx)))

mx <- assay(spe, "normalized_by_area")[ord[1:10], ]
df_area <- cbind(cd, as.matrix(t(mx)))

p1 <- ggplot(df_ls) + labs(title="Library size normalization")
p2 <- ggplot(df_area) + labs(title="Area normalization") 

(p1 + p2) & 
    geom_point(
        aes(x_centroid, y_centroid, col=EPCAM), 
        shape=16, stroke=0, size=0.4) &
    scale_color_viridis_c() & coord_equal() & 
    theme_void() & theme(legend.position="bottom")
```

![](ind-normalization_files/figure-html/plot-genes-1.png)

This is even clearer when looking at the expression of EPCAM stratified by cell type.

Code

```
(p1 + p2) &
    geom_boxplot(
        aes(Label, EPCAM, fill=Label), 
        outlier.shape=16, outlier.stroke=0) &
    scale_fill_manual(values=unname(pals::tableau20())) &
    labs(x="Cell type", y="EPCAM expression") &
    theme_bw() & theme(
        aspect.ratio=1,
        legend.position="none",
        panel.grid.minor=element_blank(),
        plot.title=element_text(hjust=0.5),
        axis.text.x=element_text(angle=45, hjust=1))
```

![](ind-normalization_files/figure-html/plot-genes-box-1.png)

## 26.3 Spatially-aware normalization

*[SpaNorm](https://bioconductor.org/packages/3.22/SpaNorm)* ([Salim et al. 2025](https://bioconductor.org/books/release/OSTA/pages/ind-normalization.html#ref-Salim2025-SpaNorm)) is a spatially-aware normalization method that uses spatial information alongside gene expression to decompose spatially-smoothed variation into a technical and biological component. Using generalized linear models and percentile-invariant adjusted counts, `SpaNorm` provides normalized expression values for downstream analyses.

Given the high computational cost of `SpaNorm`, we do not run it here; we refer to the package vignette for details.

## 26.4 Appendix

Back to top