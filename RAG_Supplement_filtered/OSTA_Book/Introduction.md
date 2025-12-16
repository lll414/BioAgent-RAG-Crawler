Source URL: https://bioconductor.org/books/release/OSTA/pages/bkg-introduction.html
Date Scraped: 2025-12-15

---

1. [Background](https://bioconductor.org/books/release/OSTA/pages/bkg-introduction.html)
2. [1  Introduction](https://bioconductor.org/books/release/OSTA/pages/bkg-introduction.html)

# 1  Introduction

## 1.1 Introduction

This book provides reproducible examples and discussion on computational analysis workflows for spatial omics data using [Bioconductor](https://bioconductor.org/) in R. The book contains chapters describing individual analysis steps as well as extended workflows, each with examples including R code and datasets. In some examples, R code is also integrated with Python tools.

## 1.2 Contents

Book chapters are organized into several parts:

* **Background:** introduction and background on spatial omics, data representations, related R/Bioconductor infrastructure, and interoperability with Python
* **Sequencing-based platforms:** analysis steps and workflows for data from sequencing-based platforms
* **Imaging-based platforms:** analysis steps and workflows for data from imaging-based platforms
* **Platform-independent analyses:** downstream analyses and workflows that are applicable to data from both types of platforms
* **Multiple-sample analyses:** analyses and workflows applicable to datasets consisting of multiple samples (e.g. multiple tissue sections)
* **Cross-platform analyses:** downstream analyses and workflows to integrate or combine information across platforms
* **Appendices:** acknowledgments, related resources, and session information

## 1.3 Scope and who this book is for

The aim of this book is to demonstrate key principles of computational analysis workflows for spatial omics data through examples and discussion, including reproducible R code and a variety of datasets. We assume some familiarity with R programming and an understanding of the types of biological questions that single-cell and spatial omics can be used to answer. Previous experience with Bioconductor is not required.

The book covers both preprocessing and downstream analyses, so the starting point can be either raw spatial omics data, or processed spot/cell-level expression matrices and sets of spatial coordinates as the main inputs. Since preprocessing procedures vary from platform to platform, we cover some key examples only.

For most analysis steps, multiple methods are available to choose from. In general, we showcase methods that we have found to work well and are computationally scalable, with a preference for methods available through Bioconductor. The book is not intended to provide a comprehensive listing of all available methods.

The code examples will only include methods available as software packages from either Bioconductor or CRAN (in R) or PyPI (in Python). This restriction helps ensure long-term stability and maintainability, enables regular testing via the Bioconductor build system, and makes it easier for readers to adapt the examples to integrate new methods or build extended Bioconductor-based workflows. Methods available from GitHub are also discussed in the text, but are not included in the code examples.

## 1.4 Bioconductor

[Bioconductor](https://bioconductor.org/) is an “open source and open development” project providing a cohesive and flexible framework for rigorous and reproducible analyses of high-throughput genomic data in R ([Carey 2025](https://bioconductor.org/books/release/OSTA/pages/bkg-introduction.html#ref-Carey2025-Bioconductor); [Huber et al. 2015](https://bioconductor.org/books/release/OSTA/pages/bkg-introduction.html#ref-Huber2015-Bioconductor); [Gentleman et al. 2004](https://bioconductor.org/books/release/OSTA/pages/bkg-introduction.html#ref-Gentleman2004-Bioconductor)). Bioconductor provides access to more than 2,000 contributed R packages, as well as infrastructure maintained by the Bioconductor Core Team, providing a rich analysis environment for users.

A key strength of the Bioconductor framework is the modularity and open development philosophy. Packages are contributed by numerous research groups, with the Bioconductor Core Team coordinating the overall project and maintaining infrastructure, build testing, and development guidelines. Contributed packages use consistent data structures, enabling users to connect packages developed by different research groups to build analysis workflows that include the latest state-of-the-art methods. Bioconductor packages also include comprehensive documentation, including extended tutorials and package vignettes.

## 1.5 Additional introductory resources

For readers who are new to R and/or Bioconductor, additional useful resources include:

* The [Orchestrating Single-Cell Analysis with Bioconductor (OSCA)](https://bioconductor.org/books/OSCA/) online book ([Amezquita et al. 2020](https://bioconductor.org/books/release/OSTA/pages/bkg-introduction.html#ref-Amezquita2020-OSCA)), which contains comprehensive materials on analysis workflows for single-cell (non-spatial) data, as well as further introductory materials on R and Bioconductor.
* [R for Data Science (2e)](https://r4ds.hadley.nz/) (online book) provides an excellent introduction to R. Additional related books include [Advanced R](https://adv-r.hadley.nz/) and [ggplot2: Elegant Graphics for Data Analysis (3e)](https://ggplot2-book.org/).
* [Data Carpentry](https://datacarpentry.org/) and [Software Carpentry](https://software-carpentry.org/) provide online lesson materials on scientific computing skills including R programming, Python programming, the Unix shell, and version control.
* The R/Bioconductor Data Science Team at the Lieber Institute for Brain Development has a [detailed guide](https://lcolladotor.github.io/bioc_team_ds/rbioconductor-data-science-bootcamps.html) of free resources and videos to learn more about R and Bioconductor, as well as [YouTube videos](https://www.youtube.com/@lcolladotor/videos), including some of the basics of Bioconductor and infrastructure for storing gene expression data, and a wide range of topics in genomics.

## 1.6 Feedback and suggestions

We welcome feedback, suggestions, and contributions from readers in the research community. These may be provided as [GitHub issues](https://github.com/lmweber/OSTA) for further discussion with the developers. Note that all methods used within code examples must be available as packages from either Bioconductor or CRAN (in R) or PyPI (in Python), as discussed above.

## 1.7 Appendix

Back to top