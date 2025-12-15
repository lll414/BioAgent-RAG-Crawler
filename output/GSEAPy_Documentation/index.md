# Welcome to GSEAPY’s documentation!

Source URL: https://gseapy.readthedocs.io/en/latest/index.html
Date Scraped: 2025-12-11

---

## GSEAPY: Gene Set Enrichment Analysis in Python.[](#gseapy-gene-set-enrichment-analysis-in-python "Link to this heading")

[![https://badge.fury.io/py/gseapy.svg](https://badge.fury.io/py/gseapy.svg)](https://badge.fury.io/py/gseapy)
[![https://img.shields.io/conda/vn/bioconda/GSEApy.svg?style=plastic](https://img.shields.io/conda/vn/bioconda/GSEApy.svg?style=plastic)](http://bioconda.github.io)
[![Action Status](https://github.com/zqfang/GSEApy/workflows/GSEApy/badge.svg?branch=master)](https://github.com/zqfang/GSEApy/actions)
[![Documentation Status](http://readthedocs.org/projects/gseapy/badge/?version=master)](http://gseapy.readthedocs.io/en/master/?badge=master)
[![https://img.shields.io/badge/license-MIT-blue.svg](https://img.shields.io/badge/license-MIT-blue.svg)](https://img.shields.io/badge/license-MIT-blue.svg)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/gseapy.svg)

**Release notes** : <https://github.com/zqfang/GSEApy/releases>

## Citation[](#citation "Link to this heading")

```
Zhuoqing Fang, Xinyuan Liu, Gary Peltz, GSEApy: a comprehensive package for performing gene set enrichment analysis in Python,
Bioinformatics, 2022;, btac757, https://doi.org/10.1093/bioinformatics/btac757
```

## Installation[](#installation "Link to this heading")

Install gseapy package from bioconda or pypi.

```
# if you have conda (MacOS_x86-64 and Linux only)
$ conda install -c bioconda gseapy

# or use pip to install the latest release
$ pip install gseapy
```

## GSEApy is a Python/Rust implementation of **GSEA** and wrapper for **Enrichr**.[](#gseapy-is-a-python-rust-implementation-of-gsea-and-wrapper-for-enrichr "Link to this heading")

GSEApy has multiple subcommands: `gsea`, `prerank`, `ssgsea`, `gsva`, `replot` `enrichr`, `biomart`.

1. The `gsea` module produces **GSEA** results.
The input requries a txt file(FPKM, Expected Counts, TPM, et.al), a cls file, and gene\_sets file in gmt format.

2. The `prerank` module produces **Prerank tool** results.
The input expects a pre-ranked gene list dataset with correlation values, which in .rnk format, and gene\_sets file in gmt format. `prerank` module is an API to GSEA pre-rank tools.

3. The `ssgsea` module performs **single sample GSEA(ssGSEA)** analysis.
The input expects a gene list with expression values(same with `.rnk` file, and gene\_sets file in gmt format. ssGSEA enrichment score for the gene set as described by [D. Barbie et al 2009](http://www.nature.com/nature/journal/v462/n7269/abs/nature08460.html).

4. The `gsva` module performs **GSVA** analysis, which described by [Hänzelmann et al](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7).

5. The `replot` module reproduces GSEA desktop version results.
The only input for GSEAPY is the location to GSEA Desktop output results.

6. The `enrichr` module enables you to perform gene set enrichment analysis using `Enrichr` API.
Enrichr is open source and freely available online at: <http://amp.pharm.mssm.edu/Enrichr> . It runs very fast and generates results in txt format.

7. The `biomart` module helps you convert gene ids using BioMart API.

GSEApy could be used for **RNA-seq, ChIP-seq, Microarry** data. It’s used for convenient GO enrichments and produce **publishable quality figures** in python.

The full `GSEA` is far too extensive to describe here; see
[GSEA](http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page) documentation for more information. All files’ formats for GSEApy are identical to `GSEA` desktop version.

## Why GSEAPY[](#why-gseapy "Link to this heading")

I would like to use Pandas to explore my data, but I did not find a convenient tool to
do gene set enrichment analysis in python. So, here are my reasons:

* **Ability to run inside python interactive console without having to switch to R!!!**
* User friendly for both wet and dry lab users.
* Produce or reproduce publishable figures.
* Perform batch jobs easy.
* Easy to use in bash shell or your data analysis workflow, e.g. snakemake.

Table of Contents

* [1. Welcome to GSEAPY’s documentation!](introduction.html)
* [2. GSEAPY Example](gseapy_example.html)
* [3. GSVA example](gseapy_example.html#GSVA-example)
* [4. scRNA-seq Example](singlecell_example.html)
* [5. A Protocol to Prepare files for GSEApy](gseapy_tutorial.html)
* [6. Developmental Guide](run.html)
* [7. Frequently Asked Questions](faq.html)

# Indices and tables[](#indices-and-tables "Link to this heading")

* [Index](genindex.html)
* [Module Index](py-modindex.html)
* [Search Page](search.html)