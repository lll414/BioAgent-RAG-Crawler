---
source_url: https://satijalab.org/seurat/articles/visiumhd_analysis_vignette
download_date: 2025-11-20
---

# Analysis, visualization, and integration of Visium HD spatial datasets with Seurat

#### Compiled: 2024-05-06

Source: [`vignettes/visiumhd_analysis_vignette.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/visiumhd_analysis_vignette.Rmd)

`visiumhd_analysis_vignette.Rmd`

## Visium HD support in Seurat

We have [previously released support](https://satijalab.org/seurat/articles/spatial_vignette) Seurat for sequencing-based spatial transcriptomic (ST) technologies, including 10x visium and SLIDE-seq. We have now updated Seurat to be compatible with the Visium HD technology, which performs profiling at substantially higher spatial resolution than previous versions.

Users can install the Visium HD-compatible release from Github. Existing Seurat workflows for [clustering, visualization, and downstream analysis](https://satijalab.org/seurat/articles/pbmc3k_tutorial) have been updated to support both Visium and Visium HD data.

We note that Visium HD data is generated from spatially patterned olignocleotides labeled in 2um x 2um bins. However, since the data from this resolution is sparse, adjacent bins are pooled together to create 8um and 16um resolutions. 10x recommends the use of 8um binned data for analysis, but Seurat supports in the simultaneous loading of multiple binnings - and stores them in a single object as multiple assays.

In this vignette, we provide an overview of some of the spatial workflows that Seurat supports for analyzing Visium HD data, in particular:

* Unsupervised clustering
* Identification of spatial tissue domains
* Subsetting spatial regions
* Integration with scRNA-seq data
* Comparing the spatial localization of different cell types

Please note that Visium HD is a new data type, and we expect to update this vignette as we test additional methods for spatial data analysis. We strongly encourage users to explore how different parameter settings affect their results, to analyze data iteratively (and in collaboration with biological experts), and to orthogonally validate unexpected or surprising biological findings.

We focus our analysis on a Visium HD dataset from the mouse brain, available to download [here](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he) but also run clustering workflow on a dataset from the mouse intestine.

### Install Seurat Update

```
# packages required for Visium HD
install.packages("hdf5r")
install.packages("arrow")
```

```
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
```

### Load Visium HD data

* Visium HD mouse brain dataset is available for download [here](https://support.10xgenomics.com/spatial-gene-expression/datasets)
* The Seurat can store multiple binnings/resolutions in different assays
* `bin.size` parameter specifies resolutions to load (8 and 16um are loaded by default)
* Users can switch between resolutions by [changing the assay](https://satijalab.org/seurat/articles/multimodal_vignette)

```
localdir <- "/brahms/lis/visium_hd/mouse/new_mousebrain/"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))

# Setting default assay changes between 8um and 16um binning
Assays(object)
DefaultAssay(object) <- "Spatial.008um"
```

```
vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot
```

![](visiumhd_analysis_vignette_files/figure-html/qc-1.png)

### Normalize datasets

In this vignette we use standard log-normalization for spatial data. We note that the best normalization methods for spatial data are still being developed and evaluated, and encourage users to read manuscripts from the [Phipson/Davis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03241-7) and [Fan](https://www.biorxiv.org/content/10.1101/2023.08.30.555624v2) labs to learn more about potential caveats for spatial normalization.

```
# normalize both 8um and 16um bins
DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object)

DefaultAssay(object) <- "Spatial.016um"
object <- NormalizeData(object)
```

### Visualize gene expression

* Adjusting `pt.size.factor` (set to 1.2 by default) helps to visualize molecular and histological info in this HD dataset
* You can also adjust the `shape`, and `stroke` (outline) parameters for visualization

```
# switch spatial resolution to 16um from 8um
DefaultAssay(object) <- "Spatial.016um"
p1 <- SpatialFeaturePlot(object, features = "Rorb") + ggtitle("Rorb expression (16um)")

# switch back to 8um
DefaultAssay(object) <- "Spatial.008um"
p2 <- SpatialFeaturePlot(object, features = "Hpca") + ggtitle("Hpca expression (8um)")

p1 | p2
```

![](visiumhd_analysis_vignette_files/figure-html/feature.plots-1.png)

### Unsupervised clustering

While the standard scRNA-seq clustering workflow can also be applied to spatial datasets - we have observed that when working with Visium HD datasets, the [Seurat v5 sketch clustering workflow](https://satijalab.org/seurat/articles/seurat5_sketch_analysis) exhibits improved performance, especially for identifying rare and spatially restricted groups.

As described in [Hao et al, Nature Biotechnology 2023](https://www.nature.com/articles/s41587-023-01767-y) and [Hie et al](https://www.sciencedirect.com/science/article/pii/S2405471219301528), sketch-based analyses aim to ‘subsample’ large datasets in a way that preserves rare populations. Here, we sketch the Visium HD dataset, perform clustering on the subsampled cells, and then project the cluster labels back to the full dataset.

Details of the sketching procedure and workflow are described in [Hao et al, Nature Biotechnology 2023](https://www.nature.com/articles/s41587-023-01767-y) and the [Seurat v5 sketch clustering vignette](https://satijalab.org/seurat/articles/seurat5_sketch_analysis). Since the full Visium HD dataset fits in memory, we do not use any of the on-disk capabilities of Seurat v5 in this vignette.

```
# note that data is already normalized
DefaultAssay(object) <- "Spatial.008um"
object <- FindVariableFeatures(object)
object <- ScaleData(object)
# we select 50,0000 cells and create a new 'sketch' assay
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
```

```
# switch analysis to sketched cells
DefaultAssay(object) <- "sketch"

# perform clustering workflow
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 3)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
```

Now we can project the cluster labels, and dimensional reductions (PCA and UMAP) that we learned from the 50,000 sketched cells - to the entire dataset, using the `ProjectData` function.

In the resulting object, for all cells:

* cluster labels will be stored in `object$seurat_cluster.projected`
* Projected PCA embeddings will be stored in `object[["pca.008um"]]`
* Projected UMAP embeddings will be stored in `object[["umap.sketch"]]`

```
object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)
```

We can visualize the clustering results for the sketched cells, as well as the projected clustering results for the full dataset:

```
DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2
```

![](visiumhd_analysis_vignette_files/figure-html/project.plots-1.png)

Of course, we can now also visualize the unsupervised clusters based on their spatial location. Note that running `SpatialDimPlot(object, interactive = TRUE)`, also enables interactive visualization and exploration.

```
SpatialDimPlot(object, label = T, repel = T, label.size = 4)
```

![](visiumhd_analysis_vignette_files/figure-html/dim.plot-1.png)

When there are many different clusters (some of which are spatially restricted and others are mixed), plotting the spatial location of all clusters can be challenging to interpret. We find it helpful to plot the spatial location of different clusters individually. For example, we highlight the spatial localization of a few clusters below, which happen to correspond to different cortical layers:

```
Idents(object) <- "seurat_cluster.projected"
cells <- CellsByIdentities(object, idents = c(0, 4, 32, 34, 35))
p <- SpatialDimPlot(object,
  cells.highlight = cells[setdiff(names(cells), "NA")],
  cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p
```

![](visiumhd_analysis_vignette_files/figure-html/cluster.plot-1.png)

We can also find and visualize the top gene expression markers for each cluster:

```
# Crete downsampled object to make visualization either
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p
```

![](visiumhd_analysis_vignette_files/figure-html/heatmap-1.png)

### Identifying spatially-defined tissue domains

While the previous analyses consider each bin independently, spatial data enables cells to be defined not just by their neighborhood, but also by their broader spatial context.

In [Singhal et al.](https://www.nature.com/articles/s41588-024-01664-3), the authors introduce BANKSY, Building Aggregates with a Neighborhood Kernel and Spatial Yardstick (BANKSY). BANKSY performs multiple tasks, but we find it particularly valuable for identifying and segmenting tissue domains. When performing clustering, BANKSY augments a spot’s expression pattern with both the mean and the gradient of gene expression levels in a spot’s broader neighborhood.

We thank the authors for enabling BANKSY to be compatible with Seurat via the [`SeuratWrappers`](https://github.com/satijalab/seurat-wrappers) framework, which requires separate installation of the BANKSY package:

```
if (!requireNamespace("Banksy", quietly = TRUE)) {
  remotes::install_github("prabhakarlab/Banksy@devel")
}
library(SeuratWrappers)
library(Banksy)
```

Before running BANKSY, there are two important model parameters that users should consider:

* `k_geom` : Local neighborhood size. Larger values will yield larger domains
* `lambda` : Influence of the neighborhood. Larger values yield more spatially coherent domains

The `RunBanksy` function creates a new `BANKSY` assay, which can be used for dimensional reduction and clustering:

```
object <- RunBanksy(object,
  lambda = 0.8, verbose = TRUE,
  assay = "Spatial.008um", slot = "data", features = "variable",
  k_geom = 50
)
```

```
DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)
```

```
Idents(object) <- "banksy_cluster"
p <- SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)
p
```

![](visiumhd_analysis_vignette_files/figure-html/banksy.plot-1.png)

As with unsupervised clustering, we can highlight the spatial location of each tissue domain individually:

```
banksy_cells <- CellsByIdentities(object)
p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
p
```

![](visiumhd_analysis_vignette_files/figure-html/banksy.cluster.plot-1.png)

### Subset out anatomical regions

Users may wish to segment or subset out a restricted region for further downstream analysis. For example, here we create a coordinate-defined segmentation mask marking cortical and hippocampal regions from the entire dataset using the `CreateSegmentation` function, and then identify cells that fall into this region with the `Overlay` function.

The list of coordinates is available for download [here](https://www.dropbox.com/scl/fi/qbs3j1alq33f0qz892ub3/cortex-hippocampus_coordinates.csv?rlkey=lsxglb15jhjdrircy9lb6n0rd&dl=0), and users can identify these boundaries when exploring their own datasets using the `interactive=TRUE` argument to `SpatialDimPlot`.

```
cortex.coordinates <- as.data.frame(read.csv("/brahms/lis/visium_hd/final_mouse/cortex-hippocampus_coordinates.csv"))
cortex <- CreateSegmentation(cortex.coordinates)

object[["cortex"]] <- Overlay(object[["slice1.008um"]], cortex)
cortex <- subset(object, cells = Cells(object[["cortex"]]))
```

### Integration with scRNA-seq data (deconvolution)

Seurat v5 also includes support for [Robust Cell Type Decomposition](https://www.nature.com/articles/s41587-021-00830-w), a computational approach to deconvolve spot-level data from spatial datasets, when provided with an scRNA-seq reference. RCTD has been shown to accurately annotate spatial data from a variety of technologies, including SLIDE-seq, Visium, and the 10x Xenium in-situ spatial platform. We observe good performance with Visium HD as well.

To run RCTD, we first install the `spacexr` package from GitHub which implements RCTD. When running RCTD, we follow the instructions from the [RCTD vignette](https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html).

```
if (!requireNamespace("spacexr", quietly = TRUE)) {
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
}
library(spacexr)
```

RCTD takes an scRNA-seq dataset as a reference, and a spatial dataset as a query. For a reference, we use a mouse scRNA-seq dataset from the Allen Brain Atlas, available for download [here](https://www.dropbox.com/scl/fi/r1mixf4eof2cot891n215/allen_scRNAseq_ref.Rds?rlkey=ynr6s6wu1efqsjsu3h40vitt7&dl=0). The reference scRNAs-eq dataset has been reduced to 200,000 cells (and rare cell types <25 cells have been removed).

We use the cortex Visium HD object as the spatial query. For computational efficiency, we sketch the spatial query dataset, apply RCTD to deconvolute the ‘sketched’ cortical cells and annotate them, and then project these annotations to the full cortical dataset.

```
# sketch the cortical subset of the Visium HD dataset
DefaultAssay(cortex) <- "Spatial.008um"
cortex <- FindVariableFeatures(cortex)
cortex <- SketchData(
  object = cortex,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(cortex) <- "sketch"
cortex <- ScaleData(cortex)
cortex <- RunPCA(cortex, assay = "sketch", reduction.name = "pca.cortex.sketch", verbose = T)
cortex <- FindNeighbors(cortex, reduction = "pca.cortex.sketch", dims = 1:50)
cortex <- RunUMAP(cortex, reduction = "pca.cortex.sketch", reduction.name = "umap.cortex.sketch", return.model = T, dims = 1:50, verbose = T)
```

```
# load in the reference scRNA-seq dataset
ref <- readRDS("/brahms/satijar/allen_scRNAseq_ref.Rds")
```

```
Idents(ref) <- "subclass_label"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$subclass_label)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- cortex[["sketch"]]$counts
cortex_cells_hd <- colnames(cortex[["sketch"]])
coords <- GetTissueCoordinates(cortex)[cortex_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
```

```
# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
cortex <- AddMetaData(cortex, metadata = RCTD@results$results_df)
```

```
# project RCTD labels from sketched cortical cells to all cortical cells
cortex$first_type <- as.character(cortex$first_type)
cortex$first_type[is.na(cortex$first_type)] <- "Unknown"
cortex <- ProjectData(
  object = cortex,
  assay = "Spatial.008um",
  full.reduction = "pca.cortex",
  sketched.assay = "sketch",
  sketched.reduction = "pca.cortex.sketch",
  umap.model = "umap.cortex.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)
```

```
DefaultAssay(object) <- "Spatial.008um"

# we only ran RCTD on the cortical cells
# set labels to all other cells as "Unknown"
object[[]][, "full_first_type"] <- "Unknown"
object$full_first_type[Cells(cortex)] <- cortex$full_first_type[Cells(cortex)]
```

```
Idents(object) <- "full_first_type"

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(object)
excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
p <- SpatialDimPlot(object, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
p
```

![](visiumhd_analysis_vignette_files/figure-html/rctd_results-1.png)

We can now look for associations between the scRNA-seq labels of individual bins, and their tissue domain identity (as assigned by BANKSY). By asking which domains the excitatory neuron cells fall in, we can rename the BANKSY clusters as neuronal layers:

```
plot_cell_types <- function(data, label) {
  p <- ggplot(data, aes(x = get(label), y = n, fill = full_first_type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = ifelse(n >= min_count_to_show_label, full_first_type, "")), position = position_stack(vjust = 0.5), size = 2) +
    xlab(label) +
    ylab("# of Spots") +
    ggtitle(paste0("Distribution of Cell Types across ", label)) +
    theme_minimal()
}

cell_type_banksy_counts <- object[[]] %>%
  dplyr::filter(full_first_type %in% excitatory_names) %>%
  dplyr::count(full_first_type, banksy_cluster)

min_count_to_show_label <- 20

p <- plot_cell_types(cell_type_banksy_counts, "banksy_cluster")
p
```

![](visiumhd_analysis_vignette_files/figure-html/celltype.plot-1.png)

Based on this plot, we can now assign cells (even if they are not excitatory neurons) to individual neuronal layers.

```
Idents(object) <- "banksy_cluster"
object$layer_id <- "Unknown"
object$layer_id[WhichCells(object, idents = c(7))] <- "Layer 2/3"
object$layer_id[WhichCells(object, idents = c(15))] <- "Layer 4"
object$layer_id[WhichCells(object, idents = c(5))] <- "Layer 5"
object$layer_id[WhichCells(object, idents = c(1))] <- "Layer 6"
```

Finally, we can visualize the spatial distribution of other cell types, and ask which cortical layers they fall in. For example, in contrast to excitatory neurons, inhibitory (GABAergic) interneurons in the cortex are not spatially restricted to individual layers - but they do show biases.

In our [previous analysis of STARmap data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687398/), and [consistent with previous work](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6456269/), we found that SST and PV interneuron classes tend to be restricted to layers 4-6, while VIP and Lamp5 interneurons tend to be located in layers 2/3. These results were based on an in-situ imaging technology which captures single-cell profiles - and here we ask whether we can find the same result in Visium HD spot-based data.

```
# set ID to RCTD label
Idents(object) <- "full_first_type"

# Visualize distribution of 4 interneuron subtypes
inhibitory_names <- c("Sst", "Pvalb", "Vip", "Lamp5")
cell_ids <- CellsByIdentities(object, idents = inhibitory_names)
p <- SpatialDimPlot(object, cells.highlight = cell_ids, cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
p
```

![](visiumhd_analysis_vignette_files/figure-html/layer.dim.plot-1.png)

```
# create barplot to show proportions of cell types of interest
layer_table <- table(object$full_first_type, object$layer_id)[inhibitory_names, 1:4]

neuron_props <- reshape2::melt(prop.table(layer_table), margin = 1)
ggplot(neuron_props, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cell type", y = "Proportion", fill = "Layer") +
  theme_classic()
```

![](visiumhd_analysis_vignette_files/figure-html/prop.plot-1.png)

We recapitulate the same findings, previously identified in in-situ imaging data, in the Visium HD dataset. This highlights that the 8um binning of Visium HD, even though it does not represent true single cell resolution, is capable of accurately localizing scRNA-seq-defined cell types, although we strongly encourage users to orthogonally validate unexpected or surprising biological findings.

### Unsupervised clustering: mouse intestine

We briefly demonstrate our sketch-clustering workflow on a second Visium HD dataset, from the Mouse Small Intestine (FFPE), available for download [here](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-intestine). We identify clusters, visualize their spatial locations, and report their top gene expression markers:

```
localdir <- "/brahms/lis/visium_hd/Visium_HD_Public_Data/HD_public_data/Visium_HD_Mouse_Small_Intestine/outs"
object <- Load10X_Spatial(data.dir = localdir, bin.size = 8)
```

```
DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object)
object <- FindVariableFeatures(object)
object <- ScaleData(object)
```

```
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
```

```
DefaultAssay(object) <- "sketch"
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 3)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
```

```
object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)
```

```
Idents(object) <- "seurat_cluster.projected"
DefaultAssay(object) <- "Spatial.008um"

p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + theme(legend.position = "bottom")
p2 <- SpatialDimPlot(object, label = F) + theme(legend.position = "bottom")
p1 | p2
```

![](visiumhd_analysis_vignette_files/figure-html/intestine.project.plots-1.png)

We visualize the location of each cluster individually:

```
Idents(object) <- "seurat_cluster.projected"
cells <- CellsByIdentities(object, idents = c(1, 5, 18, 26))
p <- SpatialDimPlot(object, cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
p
```

![](visiumhd_analysis_vignette_files/figure-html/cluster.plot.intestine-1.png)

```
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)

DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p
```

![](visiumhd_analysis_vignette_files/figure-html/heatmap.compare-1.png)

**Session Info**

```
sessionInfo()
```

```
## R version 4.3.2 (2023-10-31)
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
## time zone: America/New_York
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] dplyr_1.1.4             patchwork_1.2.0         ggplot2_3.5.0          
## [4] Seurat_5.1.0      testthat_3.2.1          SeuratObject_5.0.1.9011
## [7] sp_2.1-4               
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3     rstudioapi_0.15.0      jsonlite_1.8.8        
##   [4] magrittr_2.0.3         ggbeeswarm_0.7.2       spatstat.utils_3.0-4  
##   [7] farver_2.1.1           rmarkdown_2.26         fs_1.6.3              
##  [10] ragg_1.2.7             vctrs_0.6.5            ROCR_1.0-11           
##  [13] memoise_2.0.1          spatstat.explore_3.2-6 htmltools_0.5.7       
##  [16] usethis_2.2.3          sass_0.4.8             sctransform_0.4.1     
##  [19] parallelly_1.37.1      KernSmooth_2.23-22     bslib_0.6.1           
##  [22] htmlwidgets_1.6.4      desc_1.4.3             ica_1.0-3             
##  [25] plyr_1.8.9             plotly_4.10.4          zoo_1.8-12            
##  [28] cachem_1.0.8           igraph_2.0.2           mime_0.12             
##  [31] lifecycle_1.0.4        pkgconfig_2.0.3        Matrix_1.6-5          
##  [34] R6_2.5.1               fastmap_1.1.1          fitdistrplus_1.1-11   
##  [37] future_1.33.1          shiny_1.8.0            digest_0.6.34         
##  [40] colorspace_2.1-0       tensor_1.5             rprojroot_2.0.4       
##  [43] RSpectra_0.16-1        irlba_2.3.5.1          pkgload_1.3.4         
##  [46] textshaping_0.3.7      labeling_0.4.3         progressr_0.14.0      
##  [49] fansi_1.0.6            spatstat.sparse_3.0-3  polyclip_1.10-6       
##  [52] abind_1.4-5            httr_1.4.7             compiler_4.3.2        
##  [55] remotes_2.4.2.1        withr_3.0.0            fastDummies_1.7.3     
##  [58] highr_0.10             pkgbuild_1.4.3         R.utils_2.12.3        
##  [61] MASS_7.3-60.0.1        sessioninfo_1.2.2      tools_4.3.2           
##  [64] vipor_0.4.7            lmtest_0.9-40          ape_5.7-1             
##  [67] beeswarm_0.4.0         httpuv_1.6.14          future.apply_1.11.1   
##  [70] goftest_1.2-3          R.oo_1.26.0            glue_1.7.0            
##  [73] nlme_3.1-164           R.cache_0.16.0         promises_1.2.1        
##  [76] grid_4.3.2             Rtsne_0.17             cluster_2.1.6         
##  [79] reshape2_1.4.4         generics_0.1.3         gtable_0.3.4          
##  [82] spatstat.data_3.0-4    R.methodsS3_1.8.2      tidyr_1.3.1           
##  [85] data.table_1.15.2      utf8_1.2.4             spatstat.geom_3.2-9   
##  [88] RcppAnnoy_0.0.22       ggrepel_0.9.5          RANN_2.6.1            
##  [91] pillar_1.9.0           stringr_1.5.1          limma_3.58.1          
##  [94] spam_2.10-0            RcppHNSW_0.6.0         later_1.3.2           
##  [97] splines_4.3.2          lattice_0.22-5         deldir_2.0-4          
## [100] survival_3.5-8         tidyselect_1.2.0       miniUI_0.1.1.1        
## [103] pbapply_1.7-2          knitr_1.45             gridExtra_2.3         
## [106] scattermore_1.2        xfun_0.42              statmod_1.5.0         
## [109] brio_1.1.4             devtools_2.4.5         matrixStats_1.2.0     
## [112] stringi_1.8.3          lazyeval_0.2.2         yaml_2.3.8            
## [115] evaluate_0.23          codetools_0.2-20       tibble_3.2.1          
## [118] cli_3.6.2              uwot_0.1.16            xtable_1.8-4          
## [121] reticulate_1.35.0      systemfonts_1.0.5      munsell_0.5.0         
## [124] jquerylib_0.1.4        styler_1.10.2          Rcpp_1.0.12           
## [127] spatstat.random_3.2-3  globals_0.16.2         png_0.1-8             
## [130] ggrastr_1.0.2          parallel_4.3.2         ellipsis_0.3.2        
## [133] pkgdown_2.0.7          presto_1.0.0           dotCall64_1.1-1       
## [136] profvis_0.3.8          urlchecker_1.0.1       listenv_0.9.1         
## [139] viridisLite_0.4.2      scales_1.3.0           ggridges_0.5.6        
## [142] leiden_0.4.3.1         purrr_1.0.2            rlang_1.1.3           
## [145] cowplot_1.1.3
```
