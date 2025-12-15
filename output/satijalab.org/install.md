---
source_url: https://satijalab.org/seurat/articles/install
download_date: 2025-11-20
---

# Installation Instructions for Seurat

Source: [`vignettes/install.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/install.Rmd)

`install.Rmd`

To install Seurat, [R](https://www.r-project.org/) version 4.0 or greater is required. We also recommend installing [R Studio](https://www.rstudio.com/).

## Seurat v5: Seurat 5: Install from GitHub

Copy the code below to install Seurat v5:

```
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
```

The following packages are not required but are used in many Seurat v5 vignettes:

* SeuratData: automatically load datasets pre-packaged as Seurat objects
* Azimuth: local annotation of scRNA-seq and scATAC-seq queries across multiple organs and tissues
* SeuratWrappers: enables use of additional integration and differential expression methods
* Signac: analysis of single-cell chromatin data

```
remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
```

Seurat v5 utilizes BPCells to support analysis of extremely large datasets:

```
remotes::install_github("bnprks/BPCells/r")
```

For more information on BPCells installation, please see the [installation instructions](https://bnprks.github.io/BPCells/#installation). For macOS users, the following GitHub issues concerning [M1 chip installation](https://github.com/bnprks/BPCells/issues/6) and [compiler compatibility](https://github.com/bnprks/BPCells/issues/3) may be of use.

## Install from CRAN

Seurat is available on [CRAN](https://cran.r-project.org/package=Seurat) for all platforms. To install, run:

```
# Enter commands in R (or R studio, if installed)
install.packages('Seurat')
library(Seurat)
```

If you see the warning message below, enter `y`:

```
package which is only available in source form, and may need compilation of C/C++/Fortran: 'Seurat'
Do you want to attempt to install these from sources?
y/n:
```

## Install previous versions of Seurat

### Install any version 3 release

Any of the Seurat version 3 releases can be installed with the following command:

```
remotes::install_version("Seurat", version = "3.X.X")
```

### Install the last version 2 release (2.3.4)

To facilitate easy re-installation of the last version 2 release, we are hosting the binaries on our website. These can be installed with the following command:

```
source("https://z.umn.edu/archived-seurat")
```

View the script

### Older versions of Seurat

Old versions of Seurat, from Seurat v2.0.1 and up, are hosted in CRAN’s archive. To install an old version of Seurat, run:

```
# Enter commands in R (or R studio, if installed)
# Install the remotes package 
install.packages('remotes')
# Replace '2.3.0' with your desired version
remotes::install_version(package = 'Seurat', version = package_version('2.3.0'))
library(Seurat)
```

For versions of Seurat older than those not hosted on CRAN (versions 1.3.0 and 1.4.0), please download the packaged source code from our [releases page](https://github.com/satijalab/seurat/releases) and [install from the tarball](https://stackoverflow.com/questions/4739837/how-do-i-install-an-r-package-from-the-source-tarball-on-windows).

## Install the development version of Seurat

Install the development version of Seurat - directly from [GitHub](https://github.com/satijalab/seurat/tree/develop).

```
# Enter commands in R (or R studio, if installed)
# Install the remotes package
install.packages('remotes')
remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
```

## Docker

We provide docker images for Seurat via [dockerhub](https://hub.docker.com/r/satijalab/seurat).

To pull the latest image from the command line:

```
docker pull satijalab/seurat:latest
```

To use as a base image in a new Dockerfile:

```
FROM satijalab/seurat:latest
```