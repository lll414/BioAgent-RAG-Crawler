---
source_url: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
download_date: 2025-11-20
---

# Weighted Nearest Neighbor Analysis

#### Compiled: October 31, 2023

Source: [`vignettes/weighted_nearest_neighbor_analysis.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/weighted_nearest_neighbor_analysis.Rmd)

`weighted_nearest_neighbor_analysis.Rmd`

The simultaneous measurement of multiple modalities, known as multimodal analysis, represents an exciting frontier for single-cell genomics and necessitates new computational methods that can define cellular states based on multiple data types. The varying information content of each modality, even across cells in the same dataset, represents a pressing challenge for the analysis and integration of multimodal datasets. In ([Hao\*, Hao\* et al, Cell 2021](https://doi.org/10.1016/j.cell.2021.04.048)), we introduce ‘weighted-nearest neighbor’ (WNN) analysis, an unsupervised framework to learn the relative utility of each data type in each cell, enabling an integrative analysis of multiple modalities.

This vignette introduces the WNN workflow for the analysis of multimodal single-cell datasets. The workflow consists of three steps

* Independent preprocessing and dimensional reduction of each modality individually
* Learning cell-specific modality ‘weights’, and constructing a WNN graph that integrates the modalities
* Downstream analysis (i.e. visualization, clustering, etc.) of the WNN graph

We demonstrate the use of WNN analysis to two single-cell multimodal technologies: CITE-seq and 10x multiome. We define the cellular states based on both modalities, instead of either individual modality.

## WNN analysis of CITE-seq, RNA + ADT

We use the CITE-seq dataset from ([Stuart\*, Butler\* et al, Cell 2019](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8)), which consists of 30,672 scRNA-seq profiles measured alongside a panel of 25 antibodies from bone marrow. The object contains two assays, RNA and antibody-derived tags (ADT).

To run this vignette please install SeuratData, available on [GitHub](https://github.com/satijalab/seurat-data).

```
library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
```

```
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
```

We first perform pre-processing and dimensional reduction on both assays independently. We use standard normalization, but you can also use SCTransform or any alternative method.

```
DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')
```

For each cell, we calculate its closest neighbors in the dataset based on a weighted combination of RNA and protein similarities. The cell-specific modality weights and multimodal neighbors are calculated in a single function, which takes ~2 minutes to run on this dataset. We specify the dimensionality of each modality (similar to specifying the number of PCs to include in scRNA-seq clustering), but you can vary these settings to see that small changes have minimal effect on the overall results.

```
# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], 
# and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight
bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)
```

We can now use these results for downstream analysis, such as visualization and clustering. For example, we can create a UMAP visualization of the data based on a weighted combination of RNA and protein data We can also perform graph-based clustering and visualize these results on the UMAP, alongside a set of cell annotations.

```
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
```

```
p1 <- DimPlot(bm, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2
```

![](weighted_nearest_neighbor_analysis_files/figure-html/wumap.plot-1.png)

We can also compute UMAP visualization based on only the RNA and protein data and compare. We find that the RNA analysis is more informative than the ADT analysis in identifying progenitor states (the ADT panel contains markers for differentiated cells), while the converse is true of T cell states (where the ADT analysis outperforms RNA).

```
bm <- RunUMAP(bm, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
bm <- RunUMAP(bm, reduction = 'apca', dims = 1:18, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
```

```
p3 <- DimPlot(bm, reduction = 'rna.umap', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(bm, reduction = 'adt.umap', group.by = 'celltype.l2', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4
```

![](weighted_nearest_neighbor_analysis_files/figure-html/umapplot2-1.png)

We can visualize the expression of canonical marker genes and proteins on the multimodal UMAP, which can assist in verifying the provided annotations:

```
p5 <- FeaturePlot(bm, features = c("adt_CD45RA","adt_CD16","adt_CD161"),
                  reduction = 'wnn.umap', max.cutoff = 2, 
                  cols = c("lightgrey","darkgreen"), ncol = 3)
p6 <- FeaturePlot(bm, features = c("rna_TRDC","rna_MPO","rna_AVP"), 
                  reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)
p5 / p6
```

![](weighted_nearest_neighbor_analysis_files/figure-html/ftplot-1.png)

Finally, we can visualize the modality weights that were learned for each cell. Each of the populations with the highest RNA weights represent progenitor cells, while the populations with the highest protein weights represent T cells. This is in line with our biological expectations, as the antibody panel does not contain markers that can distinguish between different progenitor populations.

```
 VlnPlot(bm, features = "RNA.weight", group.by = 'celltype.l2', sort = TRUE, pt.size = 0.1) +
  NoLegend()
```

![](weighted_nearest_neighbor_analysis_files/figure-html/plotwts-1.png)

## WNN analysis of 10x Multiome, RNA + ATAC

Here, we demonstrate the use of WNN analysis to a second multimodal technology, the 10x multiome RNA+ATAC kit. We use a dataset that is publicly available on the 10x website, where paired transcriptomes and ATAC-seq profiles are measured in 10,412 PBMCs.

We use the same WNN methods as we use in the previous tab, where we apply integrated multimodal analysis to a CITE-seq dataset. In this example we will demonstrate how to:

* Create a multimodal Seurat object with paired transcriptome and ATAC-seq profiles
* Perform weighted neighbor clustering on RNA+ATAC data in single cells
* Leverage both modalities to identify putative regulators of different cell types and states

You can download the dataset from the 10x Genomics website [here](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k). Please make sure to download the following files:

* Filtered feature barcode matrix (HDF5)
* ATAC Per fragment information file (TSV.GZ)
* ATAC Per fragment information index (TSV.GZ index)

Finally, in order to run the vignette, make sure the following packages are installed:

* [Seurat](/seurat/articles/install)
* [Signac](https://stuartlab.org/signac/) for the analysis of single-cell chromatin datasets
* [EnsDb.Hsapiens.v86](https://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v86.html) for a set of annotations for hg38
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) to help manipulate data tables

```
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
```

We’ll create a Seurat object based on the gene expression data, and then add in the ATAC-seq data as a second assay. You can explore the [Signac getting started vignette](https://stuartlab.org/signac/articles/pbmc_vignette.html) for more information on the creation and processing of a ChromatinAssay object.

```
# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("../data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "../data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
pbmc[["ATAC"]] <- chrom_assay
```

We perform basic QC based on the number of detected molecules for each modality as well as mitochondrial percentage.

```
VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()
```

![](weighted_nearest_neighbor_analysis_files/figure-html/QCObject-1.png)

```
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)
```

We next perform pre-processing and dimensional reduction on both assays independently, using standard approaches for RNA and ATAC-seq data.

```
# RNA analysis
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
```

We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

```
pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
```

We annotate the clusters below. Note that you could also annotate the dataset using our supervised mapping pipelines, using either our [vignette](/seurat/articles/multimodal_reference_mapping), or [automated web tool, Azimuth](https://satijalab.org/azimuth).

```
# perform sub-clustering on cluster 6 to find additional structure
pbmc <- FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(pbmc) <- "sub.cluster"
```

```
# add annotations
pbmc <- RenameIdents(pbmc, '19' = 'pDC','20' = 'HSPC','15' = 'cDC')
pbmc <- RenameIdents(pbmc, '0' = 'CD14 Mono', '9' ='CD14 Mono', '5' = 'CD16 Mono')
pbmc <- RenameIdents(pbmc, '10' = 'Naive B', '11' = 'Intermediate B', '17' = 'Memory B', '21' = 'Plasma')
pbmc <- RenameIdents(pbmc, '7' = 'NK')
pbmc <- RenameIdents(pbmc, '4' = 'CD4 TCM', '13'= "CD4 TEM", '3' = "CD4 TCM", '16' ="Treg", '1' ="CD4 Naive", '14' = "CD4 Naive")
pbmc <- RenameIdents(pbmc, '2' = 'CD8 Naive', '8'= "CD8 Naive", '12' = 'CD8 TEM_1', '6_0' = 'CD8 TEM_2', '6_1' ='CD8 TEM_2', '6_4' ='CD8 TEM_2')
pbmc <- RenameIdents(pbmc, '18' = 'MAIT')
pbmc <- RenameIdents(pbmc, '6_2' ='gdT', '6_3' = 'gdT')
pbmc$celltype <- Idents(pbmc)
```

We can visualize clustering based on gene expression, ATAC-seq, or WNN analysis. The differences are more subtle than in the previous analysis (you can explore the weights, which are more evenly split than in our CITE-seq example), but we find that WNN provides the clearest separation of cell states.

```
p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

![](weighted_nearest_neighbor_analysis_files/figure-html/UMAPs-1.png)

For example, the ATAC-seq data assists in the separation of CD4 and CD8 T cell states. This is due to the presence of multiple loci that exhibit differential accessibility between different T cell subtypes. For example, we can visualize ‘pseudobulk’ tracks of the CD8A locus alongside violin plots of gene expression levels, using tools in the [Signac visualization vignette](https://stuartlab.org/signac/articles/visualization.html).

```
## to make the visualization easier, subset T cell clusters
celltype.names <- levels(pbmc)
tcell.names <- grep("CD4|CD8|Treg", celltype.names,value = TRUE)
tcells <- subset(pbmc, idents = tcell.names)
CoveragePlot(tcells, region = 'CD8A', features = 'CD8A', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)
```

![](weighted_nearest_neighbor_analysis_files/figure-html/coverageplotcd8-1.png)

Next, we will examine the accessible regions of each cell to determine enriched motifs. As described in the [Signac motifs vignette](https://stuartlab.org/signac/articles/motif_vignette.html), there are a few ways to do this, but we will use the [chromVAR](https://www.nature.com/articles/nmeth.4401) package from the Greenleaf lab. This calculates a per-cell accessibility score for known motifs, and adds these scores as a third assay (`chromvar`) in the Seurat object.

To continue, please make sure you have the following packages installed.

* [chromVAR](https://bioconductor.org/packages/release/bioc/html/chromVAR.html) for the analysis of motif accessibility in scATAC-seq
* [presto](https://github.com/immunogenomics/presto) for fast differential expression analyses.
* [TFBSTools](http://www.bioconductor.org/packages/release/bioc/html/TFBSTools.html) for TFBS analysis
* [JASPAR2020](https://bioconductor.org/packages/release/data/annotation/html/JASPAR2020.html) for JASPAR motif models
* [motifmatchr](https://www.bioconductor.org/packages/release/bioc/html/motifmatchr.html) for motif matching
* [BSgenome.Hsapiens.UCSC.hg38](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html) for chromVAR

**Install command for all dependencies**

```
remotes::install_github("immunogenomics/presto")
BiocManager::install(c("chromVAR", "TFBSTools", "JASPAR2020", "motifmatchr", "BSgenome.Hsapiens.UCSC.hg38"))
```

```
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(pbmc) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(pbmc), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
pbmc <- SetAssayData(pbmc, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
pbmc <- RunChromVAR(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
```

Finally, we explore the multimodal dataset to identify key regulators of each cell state. Paired data provides a unique opportunity to identify transcription factors (TFs) that satisfy multiple criteria, helping to narrow down the list of putative regulators to the most likely candidates. We aim to identify TFs whose expression is enriched in multiple cell types in the RNA measurements, but *also* have enriched accessibility for their motifs in the ATAC measurements.

As an example and positive control, the CCAAT Enhancer Binding Protein (CEBP) family of proteins, including the TF CEBPB, have been repeatedly shown to play important roles in the differentiation and function of myeloid cells including monocytes and dendritic cells. We can see that both the expression of the CEBPB, and the accessibility of the MA0466.2.4 motif (which encodes the binding site for CEBPB), are both enriched in monocytes.

```
#returns MA0466.2
motif.name <- ConvertMotifID(pbmc, name = 'CEBPB')
gene_plot <- FeaturePlot(pbmc, features = "sct_CEBPB", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot
```

![](weighted_nearest_neighbor_analysis_files/figure-html/CEBPB-1.png)

We’d like to quantify this relationship, and search across all cell types to find similar examples. To do so, we will use the `presto` package to perform fast differential expression. We run two tests: one using gene expression data, and the other using chromVAR motif accessibilities. `presto` calculates a p-value based on the Wilcox rank sum test, which is also the default test in Seurat, and we restrict our search to TFs that return significant results in both tests.

`presto` also calculates an “AUC” statistic, which reflects the power of each gene (or motif) to serve as a marker of cell type. A maximum AUC value of 1 indicates a perfect marker. Since the AUC statistic is on the same scale for both genes and motifs, we take the average of the AUC values from the two tests, and use this to rank TFs for each cell type:

```
markers_rna <- presto:::wilcoxauc.Seurat(X = pbmc, group_by = 'celltype', assay = 'data', seurat_assay = 'SCT')
markers_motifs <- presto:::wilcoxauc.Seurat(X = pbmc, group_by = 'celltype', assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(pbmc, id = motif.names)
```

```
# a simple function to implement the procedure above
topTFs <- function(celltype, padj.cutoff = 1e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}
```

We can now compute, and visualize, putative regulators for any cell type. We recover well-established regulators, including [TBX21 for NK cells](https://www.sciencedirect.com/science/article/pii/S1074761304000767), [IRF4 for plasma cells](https://pubmed.ncbi.nlm.nih.gov/16767092/), [SOX4 for hematopoietic progenitors](https://ashpublications.org/blood/article/124/21/1577/88774/Sox4-Is-Required-for-the-Formation-and-Maintenance), [EBF1 and PAX5 for B cells](https://www.nature.com/articles/ni.2641), [IRF8 and TCF4 for pDC](https://www.nature.com/articles/s41590-018-0136-9). We believe that similar strategies can be used to help focus on a set of putative regulators in diverse systems.

```
# identify top markers in NK and visualize
head(topTFs("NK"), 3)
```

```
##   RNA.group  gene   RNA.auc      RNA.pval motif.group motif.feature motif.auc
## 1        NK TBX21 0.7254543  0.000000e+00          NK      MA0690.1 0.9161365
## 2        NK EOMES 0.5889408  5.394245e-99          NK      MA0800.1 0.9239066
## 3        NK RUNX3 0.7705554 1.256004e-119          NK      MA0684.2 0.6813549
##      motif.pval   avg_auc
## 1 2.177841e-204 0.8207954
## 2 5.165634e-212 0.7564237
## 3  2.486287e-40 0.7259552
```

```
motif.name <- ConvertMotifID(pbmc, name = 'TBX21')
gene_plot <- FeaturePlot(pbmc, features = "sct_TBX21", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot
```

![](weighted_nearest_neighbor_analysis_files/figure-html/NK-1.png)

```
# identify top markers in pDC and visualize
head(topTFs("pDC"), 3)
```

```
##   RNA.group gene   RNA.auc      RNA.pval motif.group motif.feature motif.auc
## 1       pDC TCF4 0.9998773 1.881543e-162         pDC      MA0830.2 0.9965481
## 2       pDC IRF8 0.9907541 1.197641e-123         pDC      MA0652.1 0.8833112
## 3       pDC SPIB 0.9113825  0.000000e+00         pDC      MA0081.2 0.9068540
##     motif.pval   avg_auc
## 1 1.785897e-69 0.9982127
## 2 3.982413e-42 0.9370327
## 3 3.083509e-47 0.9091182
```

```
motif.name <- ConvertMotifID(pbmc, name = 'TCF4')
gene_plot <- FeaturePlot(pbmc, features = "sct_TCF4", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot
```

![](weighted_nearest_neighbor_analysis_files/figure-html/pDC-1.png)

```
# identify top markers in HSPC and visualize
head(topTFs("CD16 Mono"),3)
```

```
##   RNA.group  gene   RNA.auc      RNA.pval motif.group motif.feature motif.auc
## 1 CD16 Mono  TCF7 0.6229004  5.572097e-27   CD16 Mono      MA0769.2 0.6956458
## 2 CD16 Mono  LEF1 0.6244160  1.959810e-27   CD16 Mono      MA0768.1 0.6427922
## 3 CD16 Mono GATA3 0.6853318 1.575081e-132   CD16 Mono      MA0037.3 0.5622475
##     motif.pval   avg_auc
## 1 8.605682e-54 0.6592731
## 2 1.845888e-29 0.6336041
## 3 8.972511e-07 0.6237897
```

```
motif.name <- ConvertMotifID(pbmc, name = 'SPI1')
gene_plot <- FeaturePlot(pbmc, features = "sct_SPI1", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot
```

![](weighted_nearest_neighbor_analysis_files/figure-html/CD16Mono-1.png)

```
# identify top markers in other cell types
head(topTFs("Naive B"), 3)
```

```
##   RNA.group   gene   RNA.auc      RNA.pval motif.group motif.feature motif.auc
## 1   Naive B   TCF4 0.8350167 8.961289e-239     Naive B      MA0830.2 0.9066328
## 2   Naive B POU2F2 0.6969289  9.622025e-43     Naive B      MA0507.1 0.9740857
## 3   Naive B   EBF1 0.9114260  0.000000e+00     Naive B      MA0154.4 0.7565424
##      motif.pval   avg_auc
## 1 8.987390e-151 0.8708247
## 2 3.336610e-204 0.8355073
## 3  3.662535e-61 0.8339842
```

```
head(topTFs("HSPC"), 3)
```

```
##   RNA.group  gene   RNA.auc     RNA.pval motif.group motif.feature motif.auc
## 1      HSPC  SOX4 0.9864425 2.831723e-71        HSPC      MA0867.2 0.6830497
## 2      HSPC GATA2 0.7115385 0.000000e+00        HSPC      MA0036.3 0.8275008
## 3      HSPC MEIS1 0.8254177 0.000000e+00        HSPC      MA0498.2 0.6924225
##     motif.pval   avg_auc
## 1 1.241915e-03 0.8347461
## 2 7.591798e-09 0.7695196
## 3 6.877492e-04 0.7589201
```

```
head(topTFs("Plasma"), 3)
```

```
##   RNA.group  gene   RNA.auc     RNA.pval motif.group motif.feature motif.auc
## 1    Plasma  IRF4 0.8189420 5.329976e-35      Plasma      MA1419.1 0.9776046
## 2    Plasma MEF2C 0.9108487 3.135227e-12      Plasma      MA0497.1 0.7596637
## 3    Plasma  TCF4 0.8306796 1.041100e-13      Plasma      MA0830.2 0.7840848
##     motif.pval   avg_auc
## 1 2.334627e-12 0.8982733
## 2 1.374353e-04 0.8352562
## 3 3.028306e-05 0.8073822
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
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.66.3                  
##  [3] rtracklayer_1.58.0                Biostrings_2.66.0                
##  [5] XVector_0.38.0                    motifmatchr_1.20.0               
##  [7] TFBSTools_1.36.0                  JASPAR2020_0.99.10               
##  [9] chromVAR_1.20.2                   EnsDb.Hsapiens.v86_2.99.0        
## [11] ensembldb_2.22.0                  AnnotationFilter_1.22.0          
## [13] GenomicFeatures_1.50.4            AnnotationDbi_1.60.2             
## [15] Biobase_2.58.0                    GenomicRanges_1.50.2             
## [17] GenomeInfoDb_1.34.9               IRanges_2.32.0                   
## [19] S4Vectors_0.36.2                  BiocGenerics_0.44.0              
## [21] Signac_1.12.9000                  ggplot2_3.4.4                    
## [23] dplyr_1.1.3                       cowplot_1.1.1                    
## [25] thp1.eccite.SeuratData_3.1.5      stxBrain.SeuratData_0.1.1        
## [27] ssHippo.SeuratData_3.1.4          pbmcsca.SeuratData_3.0.0         
## [29] pbmcref.SeuratData_1.0.0          pbmcMultiome.SeuratData_0.1.4    
## [31] pbmc3k.SeuratData_3.1.4           panc8.SeuratData_3.0.2           
## [33] ifnb.SeuratData_3.0.0             hcabm40k.SeuratData_3.0.0        
## [35] cbmc.SeuratData_3.1.4             bmcite.SeuratData_0.3.0          
## [37] SeuratData_0.2.2.9001             Seurat_5.0.0                     
## [39] SeuratObject_5.0.0                sp_2.1-1                         
## 
## loaded via a namespace (and not attached):
##   [1] rappdirs_0.3.3              scattermore_1.2            
##   [3] R.methodsS3_1.8.2           ragg_1.2.5                 
##   [5] nabor_0.5.0                 tidyr_1.3.0                
##   [7] bit64_4.0.5                 knitr_1.45                 
##   [9] R.utils_2.12.2              irlba_2.3.5.1              
##  [11] DelayedArray_0.24.0         rpart_4.1.19               
##  [13] data.table_1.14.8           KEGGREST_1.38.0            
##  [15] RCurl_1.98-1.12             generics_0.1.3             
##  [17] RSQLite_2.3.1               RANN_2.6.1                 
##  [19] future_1.33.0               tzdb_0.4.0                 
##  [21] bit_4.0.5                   spatstat.data_3.0-3        
##  [23] xml2_1.3.5                  httpuv_1.6.12              
##  [25] SummarizedExperiment_1.28.0 DirichletMultinomial_1.40.0
##  [27] xfun_0.40                   hms_1.1.3                  
##  [29] jquerylib_0.1.4             evaluate_0.22              
##  [31] promises_1.2.1              fansi_1.0.5                
##  [33] restfulr_0.0.15             progress_1.2.2             
##  [35] caTools_1.18.2              dbplyr_2.3.4               
##  [37] igraph_1.5.1                DBI_1.1.3                  
##  [39] htmlwidgets_1.6.2           spatstat.geom_3.2-7        
##  [41] purrr_1.0.2                 ellipsis_0.3.2             
##  [43] RSpectra_0.16-1             backports_1.4.1            
##  [45] annotate_1.76.0             sparseMatrixStats_1.10.0   
##  [47] biomaRt_2.54.1              deldir_1.0-9               
##  [49] MatrixGenerics_1.10.0       vctrs_0.6.4                
##  [51] ROCR_1.0-11                 abind_1.4-5                
##  [53] cachem_1.0.8                withr_2.5.2                
##  [55] progressr_0.14.0            presto_1.0.0               
##  [57] checkmate_2.2.0             sctransform_0.4.1          
##  [59] GenomicAlignments_1.34.1    prettyunits_1.1.1          
##  [61] goftest_1.2-3               cluster_2.1.4              
##  [63] seqLogo_1.64.0              dotCall64_1.1-0            
##  [65] lazyeval_0.2.2              crayon_1.5.2               
##  [67] hdf5r_1.3.8                 spatstat.explore_3.2-5     
##  [69] pkgconfig_2.0.3             labeling_0.4.3             
##  [71] nlme_3.1-162                vipor_0.4.5                
##  [73] ProtGenerics_1.30.0         nnet_7.3-18                
##  [75] rlang_1.1.1                 globals_0.16.2             
##  [77] lifecycle_1.0.3             miniUI_0.1.1.1             
##  [79] filelock_1.0.2              fastDummies_1.7.3          
##  [81] BiocFileCache_2.6.1         dichromat_2.0-0.1          
##  [83] ggrastr_1.0.1               rprojroot_2.0.3            
##  [85] polyclip_1.10-6             RcppHNSW_0.5.0             
##  [87] matrixStats_1.0.0           lmtest_0.9-40              
##  [89] Matrix_1.6-1.1              zoo_1.8-12                 
##  [91] base64enc_0.1-3             beeswarm_0.4.0             
##  [93] ggridges_0.5.4              png_0.1-8                  
##  [95] viridisLite_0.4.2           rjson_0.2.21               
##  [97] bitops_1.0-7                R.oo_1.25.0                
##  [99] KernSmooth_2.23-22          spam_2.10-0                
## [101] DelayedMatrixStats_1.20.0   blob_1.2.4                 
## [103] stringr_1.5.0               parallelly_1.36.0          
## [105] spatstat.random_3.2-1       readr_2.1.4                
## [107] CNEr_1.34.0                 scales_1.2.1               
## [109] memoise_2.0.1               magrittr_2.0.3             
## [111] plyr_1.8.9                  ica_1.0-3                  
## [113] zlibbioc_1.44.0             compiler_4.2.2             
## [115] BiocIO_1.8.0                RColorBrewer_1.1-3         
## [117] fitdistrplus_1.1-11         Rsamtools_2.14.0           
## [119] cli_3.6.1                   listenv_0.9.0              
## [121] patchwork_1.1.3             pbapply_1.7-2              
## [123] htmlTable_2.4.1             Formula_1.2-5              
## [125] MASS_7.3-58.2               tidyselect_1.2.0           
## [127] stringi_1.7.12              glmGamPoi_1.10.2           
## [129] textshaping_0.3.6           highr_0.10                 
## [131] yaml_2.3.7                  ggrepel_0.9.4              
## [133] grid_4.2.2                  VariantAnnotation_1.44.1   
## [135] sass_0.4.7                  fastmatch_1.1-4            
## [137] tools_4.2.2                 future.apply_1.11.0        
## [139] parallel_4.2.2              rstudioapi_0.14            
## [141] TFMPvalue_0.0.9             foreign_0.8-84             
## [143] gridExtra_2.3               farver_2.1.1               
## [145] Rtsne_0.16                  digest_0.6.33              
## [147] pracma_2.4.2                shiny_1.7.5.1              
## [149] Rcpp_1.0.11                 later_1.3.1                
## [151] RcppAnnoy_0.0.21            httr_1.4.7                 
## [153] biovizBase_1.46.0           colorspace_2.1-0           
## [155] XML_3.99-0.14               fs_1.6.3                   
## [157] tensor_1.5                  reticulate_1.34.0          
## [159] splines_4.2.2               uwot_0.1.16                
## [161] RcppRoll_0.3.0              spatstat.utils_3.0-4       
## [163] pkgdown_2.0.7               plotly_4.10.3              
## [165] systemfonts_1.0.4           xtable_1.8-4               
## [167] poweRlaw_0.70.6             jsonlite_1.8.7             
## [169] R6_2.5.1                    Hmisc_5.1-1                
## [171] pillar_1.9.0                htmltools_0.5.6.1          
## [173] mime_0.12                   DT_0.30                    
## [175] glue_1.6.2                  fastmap_1.1.1              
## [177] BiocParallel_1.32.6         codetools_0.2-19           
## [179] utf8_1.2.4                  lattice_0.21-9             
## [181] bslib_0.5.1                 spatstat.sparse_3.0-3      
## [183] tibble_3.2.1                curl_5.1.0                 
## [185] ggbeeswarm_0.7.1            leiden_0.4.3               
## [187] gtools_3.9.4                GO.db_3.16.0               
## [189] survival_3.5-7              rmarkdown_2.25             
## [191] desc_1.4.2                  munsell_0.5.0              
## [193] GenomeInfoDbData_1.2.9      reshape2_1.4.4             
## [195] gtable_0.3.4
```