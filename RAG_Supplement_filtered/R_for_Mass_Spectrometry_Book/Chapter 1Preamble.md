Source URL: https://rformassspectrometry.github.io/book/index.html
Date Scraped: 2025-12-15

---

# Chapter 1 Preamble

The aim of the [R for Mass
Spectrometry](https://www.rformassspectrometry.org/) initiative is to
provide efficient, thoroughly documented, tested and flexible R
software for the analysis and interpretation of high throughput mass
spectrometry assays, including proteomics and metabolomics
experiments. The project formalises the longtime collaborative
development efforts of its core members under the RforMassSpectrometry
organisation to facilitate dissemination and accessibility of their
work.

This material introduces participants to the analysis and exploration
of mass spectrometry (MS) based proteomics data using R and
Bioconductor. The course will cover all levels of MS data, from raw
data to identification and quantitation data, up to the statistical
interpretation of a typical shotgun MS experiment and will focus on
hands-on tutorials. At the end of this course, the participants will
be able to manipulate MS data in R and use existing packages for their
exploratory and statistical proteomics data analysis.

## 1.1 Targeted audience and assumed background

The course material is targeted to either proteomics practitioners or
data analysts/bioinformaticians that would like to learn how to use R
and Bioconductor to analyse proteomics data. Familiarity with MS or
proteomics in general is desirable, but not essential as we will walk
through and describe a typical MS data as part of learning about the
tools. For approachable introductions to sample preparation, mass
spectrometry, data interpretation and analysis, readers are redirected
to:

* *A beginner’s guide to mass spectrometry–based proteomics* (Sinha and Mann 2020)
* *The ABC’s (and XYZ’s) of peptide sequencing* (Steen and Mann 2004)
* *How do shotgun proteomics algorithms identify proteins?* (Marcotte 2007)
* *An Introduction to Mass Spectrometry-Based Proteomics* (Shuken 2023)

A working knowledge of R (R syntax, commonly used functions, basic
data structures such as data frames, vectors, matrices, … and their
manipulation) is required. Familiarity with other Bioconductor omics
data classes and the tidyverse syntax is useful, but not necessary.

## 1.2 Setup

This material uses the latest version of the R for Mass Spectrometry
package and their dependencies. It might thus be possible that even
the latest Bioconductor stable version isn’t recent enough.

To install all the necessary package, please use the latest release of
R and execute:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("remotes")
BiocManager::install("tidyverse")
BiocManager::install("factoextra")
BiocManager::install("MsDataHub")
BiocManager::install("mzR")
BiocManager::install("rhdf5")
BiocManager::install("rpx")
BiocManager::install("MsCoreUtils")
BiocManager::install("QFeatures")
BiocManager::install("Spectra")
BiocManager::install("ProtGenerics")
BiocManager::install("PSMatch")
BiocManager::install("pheatmap")
BiocManager::install("limma")
BiocManager::install("MSnID")
BiocManager::install("Biostrings")
BiocManager::install("cleaver")
BiocManager::install("RforMassSpectrometry/SpectraVis")
```

After installation, you can download some data that will be used in
the latter chapter running the following:

```
library(rpx)
px <- PXDataset("PXD000001") ## answer yes if asked to create a cache directory
fn <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML"
mzf <- pxget(px, fn)
px <- PXDataset("PXD022816")
pxget(px, grep("mzID", pxfiles(px))[1:3])
pxget(px, grep("mzML", pxfiles(px))[1:3])
```

All software versions used to generate this document are recoded at
the end of the book in chapter [8](sec-si.html#sec-si).

## 1.3 Questions and help

For questions about specific software or their usage, please refer to
the software’s github issue page, or use the [Bioconductor support
site](http://support.bioconductor.org/).

## 1.4 Citation

If you need to cite this book, please use the following reference:

[![DOI](https://zenodo.org/badge/349528091.svg)](https://doi.org/10.5281/zenodo.15180829)

Laurent Gatto, Sebastian Gibb and Johannes Rainer, *R for Mass
Spectrometry* (2025)
[DOI:10.5281/zenodo.15180830](https://doi.org/10.5281/zenodo.15180829).

## 1.5 Acknowledgments

Thank you to [Charlotte Soneson](https://github.com/csoneson) for
fixing many typos in a previous version of this book.

## 1.6 License

[![Creative Commons Licence](https://i.creativecommons.org/l/by-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-sa/4.0/)  
This material is licensed under a [Creative Commons
Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/). You are free to
**share** (copy and redistribute the material in any medium or format)
and **adapt** (remix, transform, and build upon the material) for any
purpose, even commercially, as long as you give appropriate credit and
distribute your contributions under the same license as the original.