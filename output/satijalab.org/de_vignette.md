---
source_url: https://satijalab.org/seurat/articles/de_vignette
download_date: 2025-11-20
---

# Differential expression testing

#### Compiled: 2024-03-04

Source: [`vignettes/de_vignette.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/de_vignette.Rmd)

`de_vignette.Rmd`

### Load in the data

This vignette highlights some example workflows for performing differential expression in Seurat. For demonstration purposes, we will be using the interferon-beta stimulated human PBMCs dataset ([ifnb](https://www.nature.com/articles/nbt.4042)) that is available via the [SeuratData](https://github.com/satijalab/seurat-data) package.

```
library(Seurat)
library(SeuratData)
library(ggplot2)
ifnb <- LoadData("ifnb")
```

### Perform default differential expression tests

The bulk of Seurat’s differential expression features can be accessed through the `FindMarkers()` function. By default, Seurat performs differential expression (DE) testing based on the non-parametric Wilcoxon rank sum test. To test for DE genes between two specific groups of cells, specify the `ident.1` and `ident.2` parameters.

```
# Normalize the data
ifnb <- NormalizeData(ifnb)

# Find DE features between CD16 Mono and CD1 Mono
Idents(ifnb) <- "seurat_annotations"
monocyte.de.markers <- FindMarkers(ifnb, ident.1 = "CD16 Mono", ident.2 = "CD14 Mono")
# view results
head(monocyte.de.markers)
```

```
##        p_val avg_log2FC pct.1 pct.2 p_val_adj
## VMO1       0   5.700274 0.778 0.084         0
## MS4A4A     0   3.349751 0.748 0.143         0
## FCGR3A     0   3.281942 0.982 0.418         0
## PLAC8      0   3.268470 0.636 0.124         0
## CXCL16     0   2.014758 0.938 0.475         0
## MS4A7      0   2.386436 0.978 0.558         0
```

The results data frame has the following columns :

* p\_val : p-value (unadjusted)
* avg\_log2FC : log fold-change of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group.
* pct.1 : The percentage of cells where the feature is detected in the first group
* pct.2 : The percentage of cells where the feature is detected in the second group
* p\_val\_adj : Adjusted p-value, based on Bonferroni correction using all features in the dataset.

If the `ident.2` parameter is omitted or set to `NULL`, `FindMarkers()` will test for differentially expressed features between the group specified by `ident.1` and all other cells. Additionally, the parameter `only.pos` can be set to `TRUE` to only search for positive markers, i.e. features that are more highly expressed in the `ident.1` group.

```
# Find differentially expressed features between CD14+ Monocytes and all other cells, only
# search for positive markers
monocyte.de.markers <- FindMarkers(ifnb, ident.1 = "CD16 Mono", ident.2 = NULL, only.pos = TRUE)
# view results
head(monocyte.de.markers)
```

```
##        p_val avg_log2FC pct.1 pct.2 p_val_adj
## FCGR3A     0   4.532656 0.982 0.168         0
## MS4A7      0   3.806350 0.978 0.216         0
## CXCL16     0   3.274267 0.938 0.196         0
## VMO1       0   6.254651 0.778 0.044         0
## MS4A4A     0   4.747731 0.748 0.055         0
## LST1       0   2.927351 0.912 0.228         0
```

### Perform DE analysis within the same cell type across conditions

Since this dataset contains treatment information (control versus stimulated with interferon-beta), we can also ask what genes change in different conditions for cells of the same type. First, we create a column in the `meta.data` slot to hold both the cell type and treatment information and switch the current `Idents` to that column. Then we use `FindMarkers()` to find the genes that are different between control and stimulated CD14 monocytes.

```
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
mono.de <- FindMarkers(ifnb, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", verbose = FALSE)
head(mono.de, n = 10)
```

```
##         p_val avg_log2FC pct.1 pct.2 p_val_adj
## IFIT1       0   7.319139 0.985 0.033         0
## CXCL10      0   8.036564 0.984 0.035         0
## RSAD2       0   6.741673 0.988 0.045         0
## TNFSF10     0   6.991279 0.989 0.047         0
## IFIT3       0   6.883785 0.992 0.056         0
## IFIT2       0   7.179929 0.961 0.039         0
## CXCL11      0   8.624208 0.932 0.012         0
## CCL8        0   9.134191 0.918 0.017         0
## IDO1        0   5.455898 0.965 0.089         0
## MX1         0   5.059052 0.960 0.093         0
```

However, the p-values obtained from this analysis should be interpreted with caution, because these tests treat each cell as an independent replicate and ignore inherent correlations between cells originating from the same sample. Such analyses have been shown to find a large number of false positive associations, as has been demonstrated by [Squair et al., 2021](https://www.nature.com/articles/s41467-021-25960-2), [Zimmerman et al., 2021](https://www.nature.com/articles/s41467-021-21038-1), [Junttila et al., 2022](https://academic.oup.com/bib/article/23/5/bbac286/6649780), and others. Below, we show how pseudobulking can be used to account for such within-sample correlation.

#### Perform DE analysis after pseudobulking

To pseudobulk, we will use [AggregateExpression()](https://satijalab.org/seurat/reference/aggregateexpression) to sum together gene counts of all the cells from the same sample for each cell type. This results in one gene expression profile per sample and cell type. We can then perform DE analysis using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) on the sample level. This treats the samples, rather than the individual cells, as independent observations.

First, we need to retrieve the sample information for each cell. This is not loaded in the metadata, so we will load it from the [Github repo](https://github.com/yelabucsf/demuxlet_paper_code/tree/master/) of the source data for the original paper.

**Add sample information to the dataset**

```
# load the inferred sample IDs of each cell
ctrl <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best"), head = T, stringsAsFactors = F)
stim <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye2.stim.8.10.sm.best"), head = T, stringsAsFactors = F)
info <- rbind(ctrl, stim)

# rename the cell IDs by substituting the '-' into '.'
info$BARCODE <- gsub(pattern = "\\-", replacement = "\\.", info$BARCODE)

# only keep the cells with high-confidence sample ID
info <- info[grep(pattern = "SNG", x = info$BEST), ]

# remove cells with duplicated IDs in both ctrl and stim groups
info <- info[!duplicated(info$BARCODE) & !duplicated(info$BARCODE, fromLast = T), ]

# now add the sample IDs to ifnb 
rownames(info) <- info$BARCODE
info <- info[, c("BEST"), drop = F]
names(info) <- c("donor_id")
ifnb <- AddMetaData(ifnb, metadata = info)

# remove cells without donor IDs
ifnb$donor_id[is.na(ifnb$donor_id)] <- "unknown"
ifnb <- subset(ifnb, subset = donor_id != "unknown")
```

We can now perform pseudobulking (`AggregateExpression()`) based on the donor IDs.

```
# pseudobulk the counts based on donor-condition-celltype
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_ifnb))
```

```
## [1] "STIM_SNG-1488_NK"          "STIM_SNG-1488_DC"         
## [3] "STIM_SNG-1488_B Activated" "STIM_SNG-1488_Mk"         
## [5] "STIM_SNG-1488_pDC"         "STIM_SNG-1488_Eryth"
```

```
pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")
```

Next, we perform DE testing on the pseudobulk level for CD14 monocytes, and compare it against the previous single-cell-level DE results.

```
Idents(pseudo_ifnb) <- "celltype.stim"

bulk.mono.de <- FindMarkers(object = pseudo_ifnb, 
                         ident.1 = "CD14 Mono_STIM", 
                         ident.2 = "CD14 Mono_CTRL",
                         test.use = "DESeq2")
head(bulk.mono.de, n = 15)
```

```
##                  p_val avg_log2FC pct.1 pct.2     p_val_adj
## IL1RN    3.701542e-275   6.160156     1     1 5.201777e-271
## IFITM2   1.955626e-250   4.318976     1     1 2.748242e-246
## SSB      2.699554e-203   3.066647     1     1 3.793684e-199
## NT5C3A   2.239898e-198   5.412972     1     1 3.147729e-194
## RTCB     5.700554e-162   3.133362     1     1 8.010989e-158
## RABGAP1L 4.743010e-161   5.562364     1     1 6.665352e-157
## DYNLT1   9.735640e-159   2.402726     1     1 1.368150e-154
## PLSCR1   3.191691e-146   2.676047     1     1 4.485284e-142
## ISG20    9.664488e-145   5.443114     1     1 1.358150e-140
## NAPA     2.858013e-144   1.977719     1     1 4.016365e-140
## DDX58    5.957026e-142   4.640111     1     1 8.371409e-138
## HERC5    6.333722e-133   5.266515     1     1 8.900780e-129
## OASL     3.892853e-130   3.946745     1     1 5.470627e-126
## EIF2AK2  6.636434e-128   3.940167     1     1 9.326180e-124
## TMEM50A  6.731955e-117   1.355947     1     1 9.460417e-113
```

```
# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.mono.de) <- paste0(names(bulk.mono.de), ".bulk")
bulk.mono.de$gene <- rownames(bulk.mono.de)

names(mono.de) <- paste0(names(mono.de), ".sc")
mono.de$gene <- rownames(mono.de)

merge_dat <- merge(mono.de, bulk.mono.de, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                merge_dat$p_val.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val.bulk > 0.05 & 
                                  merge_dat$p_val.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                    merge_dat$p_val.sc > 0.05)]
print(paste0('# Common: ',length(common)))
```

```
## [1] "# Common: 3519"
```

```
print(paste0('# Only in single-cell: ',length(only_sc)))
```

```
## [1] "# Only in single-cell: 1649"
```

```
print(paste0('# Only in bulk: ',length(only_bulk)))
```

```
## [1] "# Only in bulk: 204"
```

We can see that while the p-values are correlated between the single-cell and pseudobulk data, the single-cell p-values are often smaller and suggest higher levels of significance. In particular, there are 3,519 genes with evidence of differential expression (prior to multiple hypothesis testing) in both analyses, 1,649 genes that only appear to be differentially expressed in the single-cell analysis, and just 204 genes that only appear to be differentially expressed in the bulk analysis. We can investigate these discrepancies using `VlnPlot`.

First, we can examine the top genes that are differentially expressed in both analyses.

```
# create a new column to annotate sample-condition-celltype in the single-cell dataset
ifnb$donor_id.stim <- paste0(ifnb$stim, "-", ifnb$donor_id)

# generate violin plot 
Idents(ifnb) <- "celltype.stim"
print(merge_dat[merge_dat$gene%in%common[1:2],c('gene','p_val.sc','p_val.bulk')])
```

```
##        gene p_val.sc    p_val.bulk
## 2785  IL1RN        0 3.701542e-275
## 2739 IFITM2        0 1.955626e-250
```

```
VlnPlot(ifnb, features = common[1:2], idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "stim")
```

![](de_vignette_files/figure-html/unnamed-chunk-8-1.png)

```
VlnPlot(ifnb, features = common[1:2], idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "donor_id.stim", ncol = 1)
```

![](de_vignette_files/figure-html/unnamed-chunk-8-2.png)

In both the pseudobulk and single-cell analyses, the p-values for these two genes are astronomically small. For both of these genes, when just comparing all stimulated CD4 monocytes to all control CD4 monocytes across samples, we see much higher expression in the stimulated cells. When breaking down these cells by sample, we continue to see consistently higher expression levels in the stimulated samples compared to the control samples; in other words, this finding is not driven by just one or two samples. Because of this consistency, we find this signal in both analyses.

By contrast, we can examine examples of genes that are only DE under the single-cell analysis.

```
print(merge_dat[merge_dat$gene%in%c('SRGN','HLA-DRA'),c('gene','p_val.sc','p_val.bulk')])
```

```
##         gene     p_val.sc p_val.bulk
## 5710    SRGN 4.025076e-21  0.1823188
## 2603 HLA-DRA 4.989863e-09  0.1851302
```

```
VlnPlot(ifnb, features <- c('SRGN','HLA-DRA'), idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "stim")
```

![](de_vignette_files/figure-html/unnamed-chunk-9-1.png)

```
VlnPlot(ifnb, features <- c('SRGN','HLA-DRA'), idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "donor_id.stim", ncol = 1)
```

![](de_vignette_files/figure-html/unnamed-chunk-9-2.png)

Here, SRGN and HLA-DRA both have very small p-values in the single-cell analysis (on the orders of \(10^{-21}\) and \(10^{-9}\)), but much larger p-values around 0.18 in the pseudobulk analysis. While there appears to be a difference between control and simulated cells when ignoring sample information, the signal is much weaker on the sample level, and we can see notable variability from sample to sample.

### Perform DE analysis using alternative tests

Finally, we also support many other DE tests using other methods. For completeness, the following tests are currently supported:

* “wilcox” : Wilcoxon rank sum test (default, using ‘[presto](https://github.com/immunogenomics/presto)’ package)
* “wilcox\_limma” : Wilcoxon rank sum test (using ‘[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)’ package)
* “bimod” : Likelihood-ratio test for single cell feature expression, (McDavid et al., Bioinformatics, 2013)
* “roc” : Standard AUC classifier
* “t” : Student’s t-test
* “poisson” : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets
* “negbinom” : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets
* “LR” : Uses a logistic regression framework to determine differentially expressed genes. Constructs a logistic regression model predicting group membership based on each feature individually and compares this to a null model with a likelihood ratio test.
* “MAST” : GLM-framework that treates cellular detection rate as a covariate (Finak et al, Genome Biology, 2015) (Installation instructions)
* “DESeq2” : DE based on a model using the negative binomial distribution (Love et al, Genome Biology, 2014) (Installation instructions) For MAST and DESeq2, please ensure that these packages are installed separately in order to use them as part of Seurat. Once installed, use the test.use parameter can be used to specify which DE test to use.

```
# Test for DE features using the MAST package
Idents(ifnb) <- "seurat_annotations"
head(FindMarkers(ifnb, ident.1 = "CD14 Mono", ident.2 = "CD16 Mono", test.use = "MAST"))
```