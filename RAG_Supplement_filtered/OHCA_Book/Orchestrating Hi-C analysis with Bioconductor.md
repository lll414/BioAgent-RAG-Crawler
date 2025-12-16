Source URL: https://bioconductor.org/books/release/OHCA/
Date Scraped: 2025-12-15

---

# Orchestrating Hi-C analysis with Bioconductor

![](assets/cover.jpg "Orchestrating Hi-C analysis with Bioconductor")

**Package:** OHCA  
 **Authors:** Jacques Serizay [aut, cre]  
 **Compiled:** 2025-08-22  
 **Package version:** 1.5.0  
 **R version:** **R version 4.5.1 (2025-06-13)**  
 **BioC version:** **3.22**  
 **License:** MIT + file LICENSE

This is the landing page of the **“Orchestrating Hi-C analysis with Bioconductor”** book. **The primary aim of this book is to introduce the `R` user to Hi-C analysis**. This book starts with key concepts important for the analysis of chromatin conformation capture and then presents `Bioconductor` tools that can be leveraged to process, analyze, explore and visualize Hi-C data.

## Table of contents

This book is divided in three parts:

**Part I: Introduction to Hi-C analysis**

* [Chapter 1: General principles and Hi-C data pre-processing](https://bioconductor.org/books/release/OHCA/pages/principles.html)
* [Chapter 2: The different R classes implemented to analyze Hi-C](https://bioconductor.org/books/release/OHCA/pages/data-representation.html)
* [Chapter 3: Manipulating Hi-C data in R](https://bioconductor.org/books/release/OHCA/pages/parsing.html)
* [Chapter 4: Hi-C data visualization](https://bioconductor.org/books/release/OHCA/pages/visualization.html)

**Part II: In-depth Hi-C analysis**

* [Chapter 5: Matrix-centric analysis](https://bioconductor.org/books/release/OHCA/pages/matrix-centric.html)
* [Chapter 6: Interactions-centric analysis](https://bioconductor.org/books/release/OHCA/pages/interactions-centric.html)
* [Chapter 7: Finding topological features from a Hi-C contact matrix](https://bioconductor.org/books/release/OHCA/pages/topological-features.html)

**Part III: Hi-C analysis workflows**

* [Data gateways: accessing public Hi-C data portals](https://bioconductor.org/books/release/OHCA/pages/disseminating.html)
* [Interoperability: Using HiCExperiment with other R packages](https://bioconductor.org/books/release/OHCA/pages/interoperability.html)
* [Workflow 1: Distance-dependent interactions across yeast mutants](https://bioconductor.org/books/release/OHCA/pages/workflow-yeast.html)
* [Workflow 2: Chromosome compartment cohesion upon mitosis entry](https://bioconductor.org/books/release/OHCA/pages/workflow-chicken.html)
* [Workflow 3: Inter-centromere interactions in yeast](https://bioconductor.org/books/release/OHCA/pages/workflow-centros.html)

## General audience

This books aims to demonstrate how to pre-process, parse and investigate Hi-C data in `R`. For this reason, a significant portion of this book consists of executable R code chunks. To be able to reproduce the examples demonstrated in this book and go further in the analysis of *your **real** datasets*, you will need to rely on several dependencies.

* `R >= 4.3` is required. You can check R version by typing `version` in an R console or in RStudio. If you do not have `R >= 4.3` installed, you will need to update your `R` version, as most extra dependencies will require `R >= 4.3`.

Installing R 4.3 👇

Detailed instructions are available [here](https://github.com/js2264/setup_ubuntu) to install `R 4.3` on a Linux machine (Ubuntu 22.04).

Briefly, to install pre-compiled version of `R 4.3.0`:

```
# This is adapted from Posit (https://docs.posit.co/resources/install-r/)
export R_VERSION=4.3.0

# Install curl and gdebi-core
sudo apt update -qq
sudo apt install curl gdebi-core -y

# Fetching the `.deb` install file from Posit repository
curl -O https://cdn.rstudio.com/r/ubuntu-2204/pkgs/r-${R_VERSION}_1_amd64.deb

# Install R
sudo gdebi r-${R_VERSION}_1_amd64.deb --non-interactive -q

# Optional: create a symlink to add R to your PATH
sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R
```

If you have some issues when installing the Hi-C packages listed below, you may need to install the following system libraries:

```
sudo apt update -qq
sudo apt install -y \
    automake make cmake fort77 gfortran \
    bzip2 unzip ftp build-essential \
    libc6 libreadline-dev \
    libpng-dev libjpeg-dev libtiff-dev \
    libx11-dev libxt-dev x11-common \
    libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libfontconfig1-dev \
    libbz2-dev liblzma-dev libtool \
    libxml2 libxml2-dev \
    libzstd-dev zlib1g-dev \
    libdb-dev libglu1-mesa-dev \
    libncurses5-dev libghc-zlib-dev libncurses-dev \
    libpcre3-dev libxml2-dev libblas-dev libzmq3-dev \
    libssl-dev libcurl4-openssl-dev \
    libgsl-dev libeigen3-dev libboost-all-dev \
    libgtk2.0-dev xvfb xauth xfonts-base apt-transport-https \
    libhdf5-dev libudunits2-dev libgdal-dev libgeos-dev \
    libproj-dev libnode-dev libmagick++-dev
```

* `Bioconductor >= 3.18` is also required. You can check whether `Bioconductor` is available and its version in `R` by typing `BiocManager::version()`. If you do not have `BiocManager` >= 3.18 installed, you will need to update it as follows:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")
```

* You will also need important packages, which will be described in length in this book. The following `R` code will set up most of the extra dependencies:

```
BiocManager::install("HiCExperiment", ask = FALSE)
BiocManager::install("HiCool", ask = FALSE)
BiocManager::install("HiContacts", ask = FALSE)
BiocManager::install("HiContactsData", ask = FALSE)
BiocManager::install("fourDNData", ask = FALSE)
BiocManager::install("DNAZooData", ask = FALSE)
```

## Developers

For developers or advanced R users, the `devel` versions of these packages can be installed by installing Bioc `devel` version prior to package installation:

```
BiocManager::install(version = "devel")
BiocManager::install("HiCExperiment", ask = FALSE)
BiocManager::install("HiCool", ask = FALSE)
BiocManager::install("HiContacts", ask = FALSE)
BiocManager::install("HiContactsData", ask = FALSE)
BiocManager::install("fourDNData", ask = FALSE)
BiocManager::install("DNAZooData", ask = FALSE)
```

## Docker image

If you have `docker` installed, the easiest approach would be to run the following command in a `shell` terminal:

```
docker run -it ghcr.io/js2264/ohca:latest R
```

This will fetch a `docker` image with the latest development versions of the aforementioned packages pre-installed, and initiate an interactive R session.

## Building book

The OHCA book has been rendered in R thanks to a number of packages, including but not only:

* **`BiocBook`**
* `devtools`
* `quarto`
* `rebook`

To build this book locally, you can run:

```
git clone git@github.com:js2264/OHCA.git && cd OHCA
quarto render
```

Warning

All dependencies listed above will be required!

The actual rendering of this book is done by GitHub Actions, and the rendered static website is hosted by GitHub Pages.

Back to top