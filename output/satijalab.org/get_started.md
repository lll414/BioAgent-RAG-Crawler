---
source_url: https://satijalab.org/seurat/articles/get_started
download_date: 2025-11-20
---

# Getting Started with Seurat v4

Source: [`vignettes/get_started.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/get_started.Rmd)

`get_started.Rmd`

We provide a series of vignettes, tutorials, and analysis walkthroughs to help users get started with Seurat. You can also check out our [Reference page](/seurat/reference/) which contains a full list of functions available to users.

## Introductory Vignettes

For new users of Seurat, we suggest starting with a guided walk through of a dataset of 2,700 Peripheral Blood Mononuclear Cells (PBMCs) made publicly available by 10X Genomics. This tutorial implements the major components of a standard unsupervised clustering workflow including QC and data filtration, calculation of high-variance genes, dimensional reduction, graph-based clustering, and the identification of cluster markers.

We provide additional introductory vignettes for users who are interested in analyzing multimodal single-cell datasets (e.g. from CITE-seq, or the 10x multiome kit), or spatial datasets (e.g. 10x Visium or Vizgen MERFISH).

|  |  |  |
| --- | --- | --- |
| [Guided tutorial — 2,700 PBMCs](/seurat/articles/pbmc3k_tutorial) | [Multimodal analysis](/seurat/articles/multimodal_vignette) | [Analysis of spatial datasets (Sequencing-based)](/seurat/articles/spatial_vignette) |
|  |  |  |
| A basic overview of Seurat that includes an introduction to common analytical workflows. | An introduction to working with multi-modal datasets in Seurat. | Learn to explore spatially-resolved transcriptomic data with examples from 10x Visium and Slide-seq v2. |
| [GO](/seurat/articles/pbmc3k_tutorial) | [GO](/seurat/articles/multimodal_vignette) | [GO](/seurat/articles/spatial_vignette) |

|  |
| --- |
| [Analysis of spatial datasets (Imaging-based)](spatial_vignette_2.html) |
|  |
| Learn to explore spatially-resolved data from multiplexed imaging technologies, including MERFISH, Xenium, CosMx SMI, and CODEX. |
| [GO](spatial_vignette_2.html) |

## Data Integration

Recently, we have developed [computational methods](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8) for integrated analysis of single-cell datasets generated across different conditions, technologies, or species. As an example, we provide a guided walk through for integrating and comparing PBMC datasets generated under different stimulation conditions. We provide additional vignettes demonstrating how to leverage an annotated scRNA-seq reference to map and label cells from a query, and to efficiently integrate large datasets.

|  |  |  |
| --- | --- | --- |
| [Introduction to scRNA-seq integration](/seurat/articles/integration_introduction) | [Mapping and annotating query datasets](/seurat/articles/integration_mapping) | [Fast integration using reciprocal PCA (RPCA)](/seurat/articles/integration_rpca) |
|  |  |  |
| An introduction to integrating scRNA-seq datasets in order to identify and compare shared cell types across experiments. | Learn how to map a query scRNA-seq dataset onto a reference in order to automate the annotation and visualization of query cells. | Identify anchors using the reciprocal PCA (rPCA) workflow, which performs a faster and more conservative integration. |
| [GO](/seurat/articles/integration_introduction) | [GO](/seurat/articles/integration_mapping) | [GO](/seurat/articles/integration_rpca) |

|  |  |  |
| --- | --- | --- |
| [Tips for integrating large datasets](integration_large_datasets.html) | [Integrating scRNA-seq and scATAC-seq data](atacseq_integration_vignette.html) | [Multimodal Reference Mapping](/seurat/articles/multimodal_reference_mapping) |
|  |  |  |
| Tips and examples for integrating very large scRNA-seq datasets (including >200,000 cells). | Annotate, visualize, and interpret an scATAC-seq experiment using scRNA-seq data from the same biological system. | Analyze query data in the context of multimodal reference atlases. |
| [GO](integration_large_datasets.html) | [GO](atacseq_integration_vignette.html) | [GO](/seurat/articles/multimodal_reference_mapping) |

## Additional New Methods

Seurat also offers additional novel statistical methods for analyzing single-cell data. These include:

* Weighted-nearest neighbor (WNN) analysis: to define cell state based on multiple modalities [[paper](https://doi.org/10.1016/j.cell.2021.04.048)]
* Mixscape: to analyze data from pooled single-cell CRISPR screens [[paper](https://doi.org/10.1038/s41588-021-00778-2)]
* SCTransform: Improved normalization for single-cell RNA-seq data [[paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)]]
* SCTransform, v2 regularization [[paper](https://www.biorxiv.org/content/10.1101/2021.07.07.451498v1.full)]]

|  |  |  |
| --- | --- | --- |
| [Weighted Nearest Neighbor Analysis](/seurat/articles/weighted_nearest_neighbor_analysis) | [Mixscape](/seurat/articles/mixscape_vignette) | [SCTransform](/seurat/articles/sctransform_vignette) |
|  |  |  |
| Analyze multimodal single-cell data with weighted nearest neighbor analysis in Seurat v4. | Explore new methods to analyze pooled single-celled perturbation screens. | Examples of how to use the SCTransform wrapper in Seurat. |
| [GO](/seurat/articles/weighted_nearest_neighbor_analysis) | [GO](/seurat/articles/mixscape_vignette) | [GO](/seurat/articles/sctransform_vignette) |

|  |
| --- |
| [SCTransform, v2 regularization](sctransform_v2_vignette.html) |
|  |
| Examples of how to perform normalization, feature selection, integration, and differential expression with an updated version of sctransform. |
| [GO](sctransform_v2_vignette.html) |

## Other

Here we provide a series of short vignettes to demonstrate a number of features that are commonly used in Seurat. We’ve focused the vignettes around questions that we frequently receive from users. Click on a vignette to get started.

|  |  |  |
| --- | --- | --- |
| [Visualization](/seurat/articles/visualization_vignette) | [Cell Cycle Regression](/seurat/articles/cell_cycle_vignette) | [Differential Expression Testing](/seurat/articles/de_vignette) |
|  |  |  |
| An overview of the major visualization functionality within Seurat. | Mitigate the effects of cell cycle heterogeneity by computing cell cycle phase scores based on marker genes. | Perform differential expression (DE) testing in Seurat using a number of frameworks. |
| [GO](/seurat/articles/visualization_vignette) | [GO](/seurat/articles/cell_cycle_vignette) | [GO](/seurat/articles/de_vignette) |

|  |  |  |
| --- | --- | --- |
| [Demultiplex Cell Hashing data](/seurat/articles/hashing_vignette) | [Interoperability with Other Analysis Tools](conversion_vignette.html) | [Parallelization](future_vignette.html) |
|  |  |  |
| Learn how to work with data produced with Cell Hashing. | Convert data between formats for different analysis tools. | Speed up compute-intensive functions with parallelization. |
| [GO](/seurat/articles/hashing_vignette) | [GO](conversion_vignette.html) | [GO](future_vignette.html) |

## SeuratWrappers

In order to facilitate the use of community tools with Seurat, we provide the Seurat Wrappers package, which contains code to run other analysis tools on Seurat objects. For the initial release, we provide wrappers for a few packages in the table below but would encourage other package developers interested in interfacing with Seurat to check out our contributor guide [here](https://github.com/satijalab/seurat.wrappers/wiki/Submission-Process).

| Package | Vignette | Reference | Source |
| --- | --- | --- | --- |
| alevin | [Import alevin counts into Seurat](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/alevin.html) | [Srivastava et. al., Genome Biology 2019](https://doi.org/10.1186/s13059-019-1670-y) | <https://github.com/k3yavi/alevin-Rtools> |
| ALRA | [Zero-preserving imputation with ALRA](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/alra.html) | [Linderman et al, bioRxiv 2018](https://doi.org/10.1101/397588) | <https://github.com/KlugerLab/ALRA> |
| CoGAPS | [Running CoGAPS on Seurat Objects](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/cogaps.html) | [Stein-O’Brien et al, Cell Systems 2019](https://doi.org/10.1016/j.cels.2019.04.004) | <https://www.bioconductor.org/packages/release/bioc/html/CoGAPS.html> |
| Conos | [Integration of datasets using Conos](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/conos.html) | [Barkas et al, Nature Methods 2019](https://doi.org/10.1038/s41592-019-0466-z) | <https://github.com/hms-dbmi/conos> |
| fastMNN | [Running fastMNN on Seurat Objects](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/fast_mnn.html) | [Haghverdi et al, Nature Biotechnology 2018](https://doi.org/10.1038/nbt.4091) | <https://bioconductor.org/packages/release/bioc/html/scran.html> |
| glmpca | [Running GLM-PCA on a Seurat Object](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/glmpca.html) | [Townes et al, Genome Biology 2019](https://doi.org/10.1186/s13059-019-1861-6) | <https://github.com/willtownes/glmpca> |
| Harmony | [Integration of datasets using Harmony](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/harmony.html) | [Korsunsky et al, Nature Methods 2019](https://doi.org/10.1038/s41592-019-0619-0) | <https://github.com/immunogenomics/harmony> |
| LIGER | [Integrating Seurat objects using LIGER](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/liger.html) | [Welch et al, Cell 2019](https://doi.org/10.1016/j.cell.2019.05.006) | <https://github.com/MacoskoLab/liger> |
| Monocle3 | [Calculating Trajectories with Monocle 3 and Seurat](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html) | [Cao et al, Nature 2019](https://doi.org/10.1038/s41586-019-0969-x) | <https://cole-trapnell-lab.github.io/monocle3> |
| Nebulosa | [Visualization of gene expression with Nebulosa](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/nebulosa.html) | [Jose Alquicira-Hernandez and Joseph E. Powell, Under Review](https://github.com/powellgenomicslab/Nebulosa) | <https://github.com/powellgenomicslab/Nebulosa> |
| schex | [Using schex with Seurat](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/schex.html) | [Freytag, R package 2019](https://doi.org/0.1242/dev.173807) | <https://github.com/SaskiaFreytag/schex> |
| scVelo | [Estimating RNA Velocity using Seurat and scVelo](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html) | [Bergen et al, bioRxiv 2019](https://doi.org/10.1101/820936) | <https://scvelo.readthedocs.io/> |
| Velocity | [Estimating RNA Velocity using Seurat](https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/velocity.html) | [La Manno et al, Nature 2018](10.1038/s41586-018-0414-6) | <https://velocyto.org> |
| CIPR | [Using CIPR with human PBMC data](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/cipr.html) | [Ekiz et. al., BMC Bioinformatics 2020](https://doi.org/10.1186/s12859-020-3538-2) | <https://github.com/atakanekiz/CIPR-Package> |
| miQC | [Running miQC on Seurat objects](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/miQC.html) | [Hippen et. al., bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.03.03.433798v1) | <https://github.com/greenelab/miQC> |
| tricycle | [Running estimate\_cycle\_position from tricycle on Seurat Objects](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/tricycle.html) | [Zheng et. al., bioRxiv 2021](https://doi.org/10.1101/2021.04.06.438463) | <https://www.bioconductor.org/packages/release/bioc/html/tricycle.html> |