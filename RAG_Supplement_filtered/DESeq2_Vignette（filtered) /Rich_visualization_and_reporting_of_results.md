# Rich visualization and reporting of results

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

**regionReport** An HTML and PDF summary of the results with plots can also be generated using the [regionReport](http://bioconductor.org/packages/regionReport) package. The *DESeq2Report* function should be run on a *DESeqDataSet* that has been processed by the *DESeq* function. For more details see the manual page for *DESeq2Report* and an example vignette in the [regionReport](http://bioconductor.org/packages/regionReport) package.

**Glimma** Interactive visualization of DESeq2 output, including MA-plots (also called MD-plots) can be generated using the [Glimma](http://bioconductor.org/packages/Glimma) package. See the manual page for *glMDPlot.DESeqResults*.

**pcaExplorer** Interactive visualization of DESeq2 output, including PCA plots, boxplots of counts and other useful summaries can be generated using the [pcaExplorer](http://bioconductor.org/packages/pcaExplorer) package. See the *Launching the application* section of the package vignette.

**iSEE** Provides functions for creating an interactive Shiny-based graphical user interface for exploring data stored in SummarizedExperiment objects, including row- and column-level metadata. Particular attention is given to single-cell data in a SingleCellExperiment object with visualization of dimensionality reduction results. [iSEE](https://bioconductor.org/packages/iSEE) is on Bioconductor. An example wrapper function for converting a *DESeqDataSet* to a SingleCellExperiment object for use with *iSEE* can be found at the following gist, written by Federico Marini:

* <https://gist.github.com/federicomarini/4a543eebc7e7091d9169111f76d59de1>

The [iSEEde](https://bioconductor.org/packages/iSEEde) package provides additional panels that facilitate the interactive visualisation of differential expression results in iSEE applications.

**DEvis** DEvis is a powerful, integrated solution for the analysis of differential expression data. This package includes an array of tools for manipulating and aggregating data, as well as a wide range of customizable visualizations, and project management functionality that simplify RNA-Seq analysis and provide a variety of ways of exploring and analyzing data. *DEvis* can be found on [CRAN](https://cran.r-project.org/package=DEVis) and [GitHub](https://github.com/price0416/DEvis).