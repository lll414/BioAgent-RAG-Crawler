# I have trouble installing DESeq2 on Ubuntu/Linuxâ¦

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

â*I try to install DESeq2, but I get an error trying to install the R packages XML and/or RCurl:*â

`ERROR: configuration failed for package XML`

`ERROR: configuration failed for package RCurl`

You need to install the following devel versions of packages using your standard package manager, e.g.Â `sudo apt-get install` or `sudo apt install`

* libxml2-dev
* libcurl4-openssl-dev