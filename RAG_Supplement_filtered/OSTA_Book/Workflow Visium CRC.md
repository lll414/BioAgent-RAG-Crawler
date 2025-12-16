Source URL: https://bioconductor.org/books/release/OSTA/pages/seq-workflow-visium-crc.html
Date Scraped: 2025-12-15

---

1. [Sequencing-based platforms](https://bioconductor.org/books/release/OSTA/pages/seq-introduction.html)
2. [14  Workflow: Visium CRC](https://bioconductor.org/books/release/OSTA/pages/seq-workflow-visium-crc.html)

# 14  Workflow: Visium CRC

## 14.1 Preamble

### 14.1.1 Introduction

In this demo, we will be analyzing Visium data on a human colorectal cancer biopsy from de Oliveira et al. ([2025](https://bioconductor.org/books/release/OSTA/pages/seq-workflow-visium-crc.html#ref-deOliveira2025-high-def)). Rather than recapitulating all possible analyses, our goal is to highlight those that might be of particular interest in the context of these data.

### 14.1.2 Dependencies

Code

```
library(osfr)
library(scran)
library(scater)
library(igraph)
library(AUCell)
library(scuttle)
library(spacexr)
library(msigdbr)
library(VisiumIO)
library(jsonlite)
library(ggspavis)
library(pheatmap)
library(patchwork)
library(OSTA.data)
library(BiocParallel)
library(DropletUtils)
library(SpatialExperiment)
# specify whether/how to 
# perform parallelization
bp <- MulticoreParam(th <- 4)
# set seed for random number generation
# in order to make results reproducible
set.seed(194849)
```

## 14.2 Data import

Code

```
# retrieve dataset from OSF repo
id <- "Visium_HumanColon_Oliveira"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)
```

Code

```
# read into 'SpatialExperiment'
obj <- TENxVisium(
    spacerangerOut=file.path(td, "outs"), 
    format="h5", 
    images="lowres")
(spe <- import(obj))
```

```
##  class: SpatialExperiment 
##  dim: 18085 4269 
##  metadata(2): resources spatialList
##  assays(1): counts
##  rownames(18085): ENSG00000187634 ENSG00000188976 ... ENSG00000198695
##    ENSG00000198727
##  rowData names(3): ID Symbol Type
##  colnames(4269): AACAATGTGCTCCGAG-1 AACACCATTCGCATAC-1 ...
##    TGTTGGTGCGGAATCA-1 TGTTGGTGGACTCAGG-1
##  colData names(4): in_tissue array_row array_col sample_id
##  reducedDimNames(0):
##  mainExpName: Gene Expression
##  altExpNames(0):
##  spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
##  imgData names(4): sample_id image_id data scaleFactor
```

![](seq-workflow-visium-crc_files/figure-html/plt-hne-1.png)

## 14.3 Quality control

Code

```
# use gene symbols as feature names
rownames(spe) <- make.unique(rowData(spe)$Symbol)
# add per-cell quality control metrics
sub <- list(mt=grep("^MT-", rownames(spe)))
spe <- addPerCellQCMetrics(spe, subsets=sub)
```

![](seq-workflow-visium-crc_files/figure-html/plt-xy-qc-1.png)

Code

```
# determine outliers via thresholding on MAD from the median
ol <- perCellQCFilters(spe, sub.fields="subsets_mt_percent")
# add results as cell metadata
colData(spe)[names(ol)] <- ol 
# tabulate # and % of cells that'd 
# be discarded for different reasons
data.frame(
    check.names=FALSE,
    `#`=apply(ol, 2, sum), 
    `%`=round(100*apply(ol, 2, mean), 2))
```

```
##                            #     %
##  low_lib_size              3  0.07
##  low_n_features          636 14.90
##  high_subsets_mt_percent 161  3.77
##  discard                 779 18.25
```

Let’s see which spots would be excluded according to the above criteria:

Code

```
lapply(names(ol), \(.) 
    plotCoords(spe, annotate=.) + ggtitle(.)) |>
    wrap_plots(nrow=1, guides="collect") &
    guides(col=guide_legend(override.aes=list(size=3))) &
    scale_color_manual("discard", values=c("lavender", "purple")) &
    theme(plot.title=element_text(hjust=0.5), legend.key.size=unit(0, "lines"))
```

![](seq-workflow-visium-crc_files/figure-html/plt-ol-xy-1.png)

It seems like low-quality spots are highly spatially organized, so that might might encourage us to not remove them, for now. We will see further below how the quality control metrics used here, and spots deemed to be `discarded`, are distributed across (transcription-based) clusters.

## 14.4 Processing

Code

```
# log-library size normalization
spe <- logNormCounts(spe)
# highly variable feature selection
tbl <- modelGeneVar(spe)
sel <- getTopHVGs(tbl, n=2e3)
# principal component analysis
spe <- runPCA(spe, subset_row=sel)
```

## 14.5 Clustering

As an unsupervised approach, we perform shared nearest-neighbor (SNN) graph-based clustering using the Leiden community detection algorithm.

Code

```
# build shared nearest-neighbor (SNN) graph
g <- buildSNNGraph(spe, use.dimred="PCA", type="jaccard")
# cluster via Leiden community detection algorithm
k <- cluster_leiden(g, objective_function="modularity", resolution=0.5)
table(spe$Leiden <- factor(k$membership))
```

```
##  
##    1   2   3   4   5   6   7   8   9  10 
##  326 504 710 432 732 458 315 171 469 152
```

## 14.6 Deconvolution

In a complementary approach, we deconvolute spot measurements using (annotated) reference single-cell data provided by the authors. Let’s first retrieve these data, alongside corresponding cell metadata, which includes low- (`Level1`) and high-resolution (`Level2`) annotations into 10 and 32 subpopulations, respectively.

Code

```
# retrieve dataset from OSF repo
id <- "Chromium_HumanColon_Oliveira"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)
```

Code

```
# read into 'SingleCellExperiment'
fs <- list.files(td, full.names=TRUE)
h5 <- grep("h5$", fs, value=TRUE)
sce <- read10xCounts(h5, col.names=TRUE)
# add cell metadata
csv <- grep("csv$", fs, value=TRUE)
cd <- read.csv(csv, row.names=1)
colData(sce)[names(cd)] <- cd[colnames(sce), ]
# use gene symbols as feature names
rownames(sce) <- make.unique(rowData(sce)$Symbol)
# exclude cells deemed to be of low-quality
sce <- sce[, sce$QCFilter == "Keep"]
# tabulate subpopulations
table(sce$Level1)
```

```
##  
##                B cells           Endothelial            Fibroblast 
##                  33611                  7969                 28653 
##  Intestinal Epithelial               Myeloid              Neuronal 
##                  22763                 25105                  4199 
##          Smooth Muscle               T cells                 Tumor 
##                  43308                 29272                 65626
```

Next, we perform deconvolution with *[spacexr](https://bioconductor.org/packages/3.22/spacexr)*’s (RCTD) ([Cable et al. 2022](https://bioconductor.org/books/release/OSTA/pages/seq-workflow-visium-crc.html#ref-Cable2022-RCTD)). By default, `runRctd()`’s `rctd_mode="doublet"`, i.e., at most two subpopulations are fit per pixel; here, we set `rctd_mode="full"` in order to allow for an arbitrary number of subpopulations to be fit instead.

Here, we filter the reference data to contain only cells from the same patient, and downsample to retain a limited number of cells per subpopulation. This is not strictly necessary, assuming that clusters are transcriptionally stable across patients, but helps with reducing runtime and memory consumption here.

Code

```
# prep reference data (Chromium);
# subset cells from same patient
.sce <- sce[, grepl("P2", sce$Patient)]
# downsample to at most 2,000 cells per cluster
cs <- split(seq_len(ncol(.sce)), .sce$Level1)
cs <- lapply(cs, \(.) sample(., min(length(.), 2e3)))
.sce <- .sce[, unlist(cs)]
# run 'RCTD' deconvolution
rctd_data <- createRctd(spe, .sce, cell_type_col="Level1")
(res <- runRctd(rctd_data, max_cores=th, rctd_mode="full"))
```

```
##  class: SpatialExperiment 
##  dim: 9 4269 
##  metadata(4): spatial_rna config cell_type_info internal_vars
##  assays(1): weights
##  rownames(9): B cells Endothelial ... T cells Tumor
##  rowData names(0):
##  colnames(4269): AACAATGTGCTCCGAG-1 AACACCATTCGCATAC-1 ...
##    TGTTGGTGCGGAATCA-1 TGTTGGTGGACTCAGG-1
##  colData names(1): sample_id
##  reducedDimNames(0):
##  mainExpName: NULL
##  altExpNames(0):
##  spatialCoords names(2) : x y
##  imgData names(0):
```

Weights inferred by `RCTD` should be normalized such that proportions of cell types sum to 1 in each spot:

Code

```
# scale weights such that they sum to 1
ws <- assay(res)
ws <- sweep(ws, 2, colSums(ws), `/`)
# add proportion estimates as metadata
ws <- data.frame(t(as.matrix(ws)))
colData(spe)[names(ws)] <- ws[colnames(spe), ]
```

For comparison with unsupervised clustering (SNN-based Leiden), we also include assignments we would obtain if we were to assign spots the most frequent label (in terms of deconvolution estimates):

Code

```
ids <- names(ws)[apply(ws, 1, which.max)]
table(spe$RCTD <- factor(ids), spe$Leiden)
```

```
##                         
##                            1   2   3   4   5   6   7   8   9  10
##    B.cells                 0  44   0   8   0   4   0   0   0   2
##    Endothelial             0 114   0   0   0  24   1   0   0   0
##    Fibroblast              0 268   0   0  11  57 311   0   0   0
##    Intestinal.Epithelial   0   5   0   4   0   0   0   0 466  12
##    Myeloid                 0  15   0   0   2   2   0   0   0   0
##    Smooth.Muscle           0   8   0   0   0 362   0   0   0   0
##    T.cells                 0  15   0   0   0   4   0   0   0   0
##    Tumor                 326  35 710 420 719   5   3 171   3 138
```

We can also compartmentalize the tissue into broad biological compartments; here, by grouping `RCTD`-based subpopulation assignments into four classes:

Code

```
lab <- list(
    tum="Tumor",
    epi="Intestinal.Epithelial",
    imm=c("B.cells", "T.cells", "Myeloid"),
    str=c("Endothelial", "Fibroblast", "Smooth.Muscle"))
idx <- match(spe$RCTD, unlist(lab))
lab <- rep.int(names(lab), sapply(lab, length))
table(spe$Domain <- factor(lab[idx]))
```

```
##  
##   epi  imm  str  tum 
##   487   96 1156 2530
```

## 14.7 Exploratory

Let’s visualize deconvolution weights in space, i.e., coloring by the proportion of a given cell type estimated to fall within a given spot:

Code

```
lapply(names(ws), \(.) 
    plotCoords(spe, annotate=.)) |>
    wrap_plots(nrow=3) & theme(
    legend.key.width=unit(0.5, "lines"),
    legend.key.height=unit(1, "lines")) &
    scale_color_gradientn(colors=pals::jet())
```

![](seq-workflow-visium-crc_files/figure-html/plt-dec-xy-1.png)

Code

```
lapply(c("Leiden", "Domain", "RCTD"), 
    \(.) plotCoords(spe, annotate=.)) |>
    wrap_plots(nrow=1) &
    theme(legend.key.size=unit(0, "lines")) &
    scale_color_manual(values=unname(pals::trubetskoy()))
```

![](seq-workflow-visium-crc_files/figure-html/plt-clu-xy-1.png)

To help characterize subpopulations from unsupervised clustering, we can view their distribution across deconvolution-based clusters and broad domains; e.g., tumor spots are quite diverse, while smooth muscle spots and (normal) epithelia map almost completely to a single cluster:

Code

```
cd <- data.frame(colData(spe))
df <- as.data.frame(with(cd, table(RCTD, Leiden)))
fd <- as.data.frame(with(cd, table(Domain, Leiden)))
ggplot(df, aes(Freq, RCTD, fill=Leiden)) + ggtitle("RCTD") +
ggplot(fd, aes(Freq, Domain, fill=Leiden)) + ggtitle("Domain") +
plot_layout(nrow=1, guides="collect") &
    labs(x="Proportion", y=NULL) &
    coord_cartesian(expand=FALSE) &
    geom_col(width=1, col="white", position="fill") &
    scale_fill_manual(values=unname(pals::trubetskoy())) &
    theme_minimal() & theme(aspect.ratio=1,
        legend.key.size=unit(2/3, "lines"),
        plot.title=element_text(hjust=0.5))
```

![](seq-workflow-visium-crc_files/figure-html/plt-clu-fq-1.png)

Let’s inspect the key drivers of (expression) variability in terms of PCs. Considering clustering and deconvolution results from above, we can see that

* PC1 distinguishes stromal from both normal and malignant epithelia;
* PC2 clearly separates (normal) intestinal epithelium from all else;
* PC3 captures a fibroblast-rich region, and normal epithelia;
* PC5 separates fibroblasts and smooth muscle cells; etc.

Code

```
pcs <- reducedDim(spe, "PCA")
colData(spe)[colnames(pcs)] <- pcs
lapply(colnames(pcs)[seq_len(6)], 
    \(.) plotCoords(spe, annotate=.) +
    scale_color_gradientn(., colors=pals::jet())) |>
    wrap_plots(nrow=2) & theme(
        plot.title=element_blank(),
        legend.key.width=unit(0.5, "lines"),
        legend.key.height=unit(1, "lines"))
```

![](seq-workflow-visium-crc_files/figure-html/plt-pcs-xy-1.png)

Quality control metrics tend to be low for specific clusters. Their patch-like pattern, in turn, explains the clustering of low-quality spots seen earlier.

Code

```
lapply(c("detected", "log_sum", "subsets_mt_percent"), \(.)
    plotColData(spe, x=., y="Leiden", color_by="discard", point_size=0.1) +
    scale_x_discrete(limits=names(sort(by(spe[[.]], spe$Leiden, median))))) |>
    wrap_plots(nrow=1, guides="collect") &
    scale_color_manual("discard", values=c("lavender", "purple")) &
    guides(col=guide_legend(override.aes=list(alpha=1, size=3))) &
    theme_minimal() & theme(
        panel.grid.minor=element_blank(), 
        legend.key.size=unit(0, "lines"))
```

![](seq-workflow-visium-crc_files/figure-html/plt-clu-qc-1.png)

## 14.8 Signatures

Rather than investigating single genes, we can also evaluate the expression of sets of genes (e.g., pathway signatures); e.g., malignant tissue may differ in metabolic activity such as glycolysis and fatty acid metabolism, or exhibit increased apoptosis (cell death) etc. Here, we retrieve hallmark gene sets for some biological phenomena from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb), using the *[msigdbr](https://bioconductor.org/packages/3.22/msigdbr)* package:

Note that such gene lists like these may come from many different places - e.g., in-house analyses, other publications matching the research question etc. As such, they may also be read in from a *.csv* file, or stem from upstream analyses. In any case, they should be curated well and interpreted with caution.

Code

```
# retrieve hallmark gene sets from 'MSigDB'
db <- msigdbr(species="Homo sapiens", collection="H")
# get list of gene symbols, one element per set
gs <- split(db$ensembl_gene, db$gs_name)
# simplify set identifiers (drop prefix, use lower case)
names(gs) <- tolower(gsub("HALLMARK_", "", names(gs)))
# how many sets?
length(gs)
```

```
##  [1] 50
```

Code

```
# how many genes in each?
range(sapply(gs, length))
```

```
##  [1]  32 201
```

Next, we will score these using *[AUCell](https://bioconductor.org/packages/3.22/AUCell)* ([Aibar et al. 2017](https://bioconductor.org/books/release/OSTA/pages/seq-workflow-visium-crc.html#ref-Aibar2017-SCENIC)), which works in two steps: (i) rank genes for every observation (here, spots), and (ii) compute AUC values for each gene set. In essence, these represent the fraction of genes (within top-ranked genes; default 5%) that are in a given set; i.e., high values correspond to high activity (in terms of coordinated gene expression).

By definition, `AUCell` yields values in [0,1]. Larger sets (more genes) are more likely to achieve higher scores by chance (e.g., a gene set of *all* genes would score 1 in any dataset). It is thus unfair to compare them directly. However, spatial distribution, correlation, and relative comparisons between subpopulations etc. are still meaningful.

Code

```
# realize (sparse) gene expression matrix
mtx <- as(logcounts(spe), "dgCMatrix") 
# use ensembl identifiers as feature names
rownames(mtx) <- rowData(spe)$ID
# build per-spot gene rankings
rnk <- AUCell_buildRankings(mtx, BPPARAM=bp, plotStats=FALSE, verbose=FALSE)
# calculate AUC for each gene set in each spot
auc <- AUCell_calcAUC(geneSets=gs, rankings=rnk, nCores=th, verbose=FALSE)
# add results as spot metadata
colData(spe)[rownames(auc)] <- res <- t(assay(auc))
```

For simplicity, we’ll continue investigating only those signatures with the highest score variability across spots:

Code

```
var <- colVars(res) # variance across spots
top <- names(tail(sort(var), 8)) # top sets
```

To summarize, MYC signalling is absent in stromal regions; the fibroblast ring surrounding a cancerous patch exhibits EMT, angiogenesis, etc.; INFa response and TNFa signalling is patch-like in both stroma and malignant epithelia.

Code

```
lapply(top, \(.) {
    spe[[.]] <- scale(spe[[.]]) # scaling
    plotCoords(spe, annotate=.) # plotting
}) |> 
    # arrange & prettify
    wrap_plots(nrow=2, guides="collect") & 
    scale_color_gradientn(
        colors=pals::jet(),
        oob=scales::squish, 
        limits=c(-2.5, 2.5)) & 
    theme(
        legend.key.width=unit(0.5, "lines"), 
        legend.key.height=unit(1, "lines"))
```

![](seq-workflow-visium-crc_files/figure-html/plt-auc-xy-1.png)

To ease interpretability, we can stratify `AUCell` scores by spot labels; these may stem from an unsupervised or deconvolution-based approach:

Code

```
for (. in c("Leiden", "RCTD")) {
    # aggregate AUC values by cluster
    mu <- aggregateAcrossCells(auc[top, ], spe[[.]], 
        use.assay.type="AUC", statistics="mean")
    # visualize as (cluster x set) heatmap
    pheatmap(
        mat=t(assay(mu)), scale="column", col=pals::coolwarm(), main=.,
        cellwidth=10, cellheight=10, treeheight_row=5, treeheight_col=5)
}
```

![](seq-workflow-visium-crc_files/figure-html/plt-auc-hm-1.png)

![](seq-workflow-visium-crc_files/figure-html/plt-auc-hm-2.png)

For the latter, we may instead correlate set scores with proportion estimates (rather than discretizing labels according to the dominant subpopulation):

Code

```
# correlate 'AUCell' signature scores with subpopulation
# proportion estimates from deconvolution with 'RCTD'
cm <- cor(as.matrix(ws), t(assay(auc[top, ])))
```

Code

```
pheatmap(cm, 
    col=pals::coolwarm(),
    breaks=seq(-1, 1, length=25),
    cellwidth=10, cellheight=10, 
    treeheight_row=5, treeheight_col=5)
```

![](seq-workflow-visium-crc_files/figure-html/plt-cor-1.png)

## 14.9 Appendix

Back to top