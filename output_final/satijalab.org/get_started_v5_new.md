---
source_url: https://satijalab.org/seurat/articles/get_started_v5_new
download_date: 2025-11-20
---

# Getting Started with Seurat

Source: [`vignettes/get_started_v5_new.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/get_started_v5_new.Rmd)

`get_started_v5_new.Rmd`

We provide a series of vignettes, tutorials, and analysis walkthroughs to help users get started with Seurat. You can also check out our [Reference page](/seurat/reference/) which contains a full list of functions available to users.

Our previous Get Started page for Seurat v4 is archived [here](/seurat/articles/get_started).

## Introductory Vignettes

For new users of Seurat, we suggest starting with a guided walk through of a dataset of 2,700 Peripheral Blood Mononuclear Cells (PBMCs) made publicly available by 10X Genomics. This tutorial implements the major components of a standard unsupervised clustering workflow including QC and data filtration, calculation of high-variance genes, dimensional reduction, graph-based clustering, and the identification of cluster markers. We provide additional vignettes introducing visualization techniques in Seurat, the sctransform normalization workflow, and storage/interaction with multimodal datasets. We also provide an ‘essential commands cheatsheet’ as a quick reference.

|  |  |  |
| --- | --- | --- |
| [Guided tutorial — 2,700 PBMCs](/seurat/articles/pbmc3k_tutorial) | [Multimodal analysis](/seurat/articles/multimodal_vignette) | [Visualization](/seurat/articles/visualization_vignette) |
|  |  |  |
| A basic overview of Seurat that includes an introduction to common analytical workflows. | An introduction to working with multi-modal datasets in Seurat. | An overview of the visualization capabilities within Seurat. |
| [GO](/seurat/articles/pbmc3k_tutorial) | [GO](/seurat/articles/multimodal_vignette) | [GO](/seurat/articles/visualization_vignette) |

|  |  |
| --- | --- |
| [SCTransform](/seurat/articles/sctransform_vignette) | [Essential Commands Cheat Sheet](/seurat/articles/essential_commands) |
|  |  |
| Examples of how to perform normalization, feature selection, integration, and differential expression with an updated version of sctransform. | Reference list of commonly used commands to store, access, explore, and analyze datasets. |
| [GO](/seurat/articles/sctransform_vignette) | [GO](/seurat/articles/essential_commands) |

## scRNA Data Integration

We have developed [computational methods](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8) for integrated analysis of single-cell datasets generated across different conditions, technologies, or species. As an example, we provide a guided walk through for integrating and comparing PBMC datasets generated under different stimulation conditions. We provide additional vignettes demonstrating how to leverage an annotated scRNA-seq reference to map and label cells from a query, and to efficiently integrate large datasets.

|  |  |  |
| --- | --- | --- |
| [Introduction to scRNA-seq integration](/seurat/articles/integration_introduction) | [scRNA-seq Integration](/seurat/articles/seurat5_integration) | [Mapping and annotating query datasets](/seurat/articles/integration_mapping) |
|  |  |  |
| An introduction to integrating scRNA-seq datasets in order to identify and compare shared cell types across experiments. | Integrate scRNA-seq datasets using a variety of computational methods. | Learn how to map a query scRNA-seq dataset onto a reference in order to automate the annotation and visualization of query cells. |
| [GO](/seurat/articles/integration_introduction) | [GO](/seurat/articles/seurat5_integration) | [GO](/seurat/articles/integration_mapping) |

## Multi-assay data

Seurat also offers support for a suite of statistical methods for analyzing multimodal single-cell data. These include methods to integrate modalities that are simultaneously measured in the same cells, modalities that are measured in different cells, and techniques to analyze pooled CRISPR screens.

|  |  |  |
| --- | --- | --- |
| [Cross-modality Bridge Integration](/seurat/articles/seurat5_integration_bridge) | [Weighted Nearest Neighbor Analysis](/seurat/articles/weighted_nearest_neighbor_analysis) | [Integrating scRNA-seq and scATAC-seq data](atacseq_integration_vignette.html) |
|  |  |  |
| Map scATAC-seq onto an scRNA-seq reference using a multi-omic bridge dataset in Seurat v5. | Analyze multimodal single-cell data with weighted nearest neighbor analysis in Seurat v4. | Annotate, visualize, and interpret an scATAC-seq experiment using scRNA-seq data from the same biological system in Seurat v3. |
| [GO](/seurat/articles/seurat5_integration_bridge) | [GO](/seurat/articles/weighted_nearest_neighbor_analysis) | [GO](atacseq_integration_vignette.html) |

|  |  |
| --- | --- |
| [Reference Mapping for Multimodal Data](/seurat/articles/multimodal_reference_mapping) | [Mixscape](/seurat/articles/mixscape_vignette) |
|  |  |
| Analyze query data in the context of multimodal reference atlases. | Explore new methods to analyze pooled single-celled perturbation screens. |
| [GO](/seurat/articles/multimodal_reference_mapping) | [GO](/seurat/articles/mixscape_vignette) |

## Flexible analysis of massively scalable datasets

In Seurat v5, we introduce new infrastructure and methods to analyze, interpret, and explore datasets that extend to millions of cells. We introduce support for ‘sketch-based’ techniques, where a subset of representative cells are stored in memory to enable rapid and iterative exploration, while the remaining cells are stored on-disk. Users can flexibly switch between both data representations, and we leverage the [BPCells package](https://bnprks.github.io/BPCells/) from Ben Parks in the Greenleaf lab to enable high-performance analysis of disk-backed data.

The vignettes below demonstrate three scalable analyses in Seurat v5: Unsupervised clustering analysis of a large dataset (1.3M neurons), Unsupervised integration and comparison of 1M PBMC from healthy and diabetic patients, and Supervised mapping of 1.5M immune cells from healthy and COVID donors. In all cases, the vignettes perform these analyses without ever loading the full datasets into memory.

|  |  |  |
| --- | --- | --- |
| [Unsupervised clustering of 1.3M neurons](/seurat/articles/seurat5_sketch_analysis) | [Integrating/comparing healthy and diabetic samples](/seurat/articles/parsebio_sketch_integration) | [Supervised mapping of 1.5M immune cells](/seurat/articles/covid_sctmapping) |
|  |  |  |
| Analyze a 1.3 million cell mouse brain dataset using on-disk capabilities powered by BPCells. | Perform sketch integration on a large dataset from Parse Biosciences. | Map PBMC datasets from COVID-19 patients to a healthy PBMC reference. |
| [GO](/seurat/articles/seurat5_sketch_analysis) | [GO](/seurat/articles/parsebio_sketch_integration) | [GO](/seurat/articles/covid_sctmapping) |

|  |
| --- |
| [BPCells Interaction](/seurat/articles/seurat5_bpcells_interaction_vignette) |
|  |
| Load and save large on-disk matrices using BPCells. |
| [GO](/seurat/articles/seurat5_bpcells_interaction_vignette) |

## Spatial analysis

These vignettes will help introduce users to the analysis of spatial datasets in Seurat v5, including technologies that leverage sequencing-based readouts, as well as technologies that leverage in-situ imaging-based readouts. The vignettes introduce data from multiple platforms including 10x Visium, SLIDE-seq, Vizgen MERSCOPE, 10x Xenium, Nanostring CosMx, and Akoya CODEX.

|  |  |  |
| --- | --- | --- |
| [Analysis of spatial datasets (Imaging-based)](/seurat/articles/seurat5_spatial_vignette_2) | [Analysis of spatial datasets (Sequencing-based)](/seurat/articles/spatial_vignette) | [Analysis of Visium HD spatial datasets](/seurat/articles/visiumhd_analysis_vignette) |
|  |  |  |
| Learn to explore spatially-resolved data from multiplexed imaging technologies, including MERSCOPE, Xenium, CosMx SMI, and CODEX. | Learn to explore spatially-resolved transcriptomic data with examples from 10x Visium and Slide-seq v2. | Learn to explore spatially-resolved transcriptomic data in high-definition from 10x Visium HD. |
| [GO](/seurat/articles/seurat5_spatial_vignette_2) | [GO](/seurat/articles/spatial_vignette) | [GO](/seurat/articles/visiumhd_analysis_vignette) |

## Other

Here we provide a series of short vignettes to demonstrate a number of features that are commonly used in Seurat. We’ve focused the vignettes around questions that we frequently receive from users.

|  |  |  |
| --- | --- | --- |
| [Cell Cycle Regression](/seurat/articles/cell_cycle_vignette) | [Differential Expression Testing](/seurat/articles/de_vignette) | [Demultiplex Cell Hashing data](/seurat/articles/hashing_vignette) |
|  |  |  |
| Mitigate the effects of cell cycle heterogeneity by computing cell cycle phase scores based on marker genes. | Perform differential expression (DE) testing in Seurat using a number of frameworks. | Learn how to work with data produced with Cell Hashing. |
| [GO](/seurat/articles/cell_cycle_vignette) | [GO](/seurat/articles/de_vignette) | [GO](/seurat/articles/hashing_vignette) |

## SeuratWrappers

In order to facilitate the use of community tools with Seurat, we provide the Seurat Wrappers package, which contains code to run other analysis tools on Seurat objects. For the initial release, we provide wrappers for a few packages in the table below but would encourage other package developers interested in interfacing with Seurat to check out our contributor guide [here](https://github.com/satijalab/seurat.wrappers/wiki/Submission-Process).

| Package | Vignette | Reference | Source |
| --- | --- | --- | --- |
