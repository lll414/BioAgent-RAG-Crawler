Source URL: https://bioconductor.org/books/release/csawBook/
Date Scraped: 2025-12-15

---

# Chapter 1 Welcome



## 1.1 Introduction

Chromatin immunoprecipitation with sequencing (ChIP-seq) is a widely used technique for identifying the genomic binding sites of a target protein.
Conventional analyses of ChIP-seq data aim to detect absolute binding (i.e., the presence or absence of a binding site) based on peaks in the read coverage.
An alternative analysis strategy is to detect of changes in the binding profile between conditions (Ross-Innes et al. [2012](#ref-rossinnes2012differential); Pal et al. [2013](#ref-pal2013)).
These differential binding (DB) analyses involve counting reads into genomic intervals and testing those counts for significant differences between conditions.
This defines a set of putative DB regions for further examination.
DB analyses are statistically easier to perform than their conventional counterparts,
as the effect of genomic biases is largely mitigated when counts for different libraries are compared at the same genomic region.
DB regions may also be more relevant as the change in binding can be associated with the biological difference between conditions.

This book describes the use of the *[csaw](https://bioconductor.org/packages/3.22/csaw)* Bioconductor package to detect differential binding (DB) in ChIP-seq experiments with sliding windows (Lun and Smyth [2016](#ref-lun2016csaw)).
In these analyses, we detect and summarize DB regions between conditions in a *de novo* manner, i.e., without making any prior assumptions about the location or width of bound regions.
We demonstrate on data from a variety of real studies focusing on changes in transcription factor binding and histone mark enrichment.
Our aim is to facilitate the practical implementation of window-based DB analyses by providing detailed code and expected output.
The code here can be adapted to any dataset with multiple experimental conditions and with multiple biological samples within one or more of the conditions;
it is similarly straightforward to accommodate batch effects, covariates and additional experimental factors.
Indeed, though the book focuses on ChIP-seq, the same software can be adapted to data from any sequencing technique where reads represent coverage of enriched genomic regions.

## 1.2 How to read this book

The descriptions in this book explore the theoretical and practical motivations behind each step of a *[csaw](https://bioconductor.org/packages/3.22/csaw)* analysis.
While all users are welcome to read it from start to finish, new users may prefer to examine the case studies presented in the later sections (Lun and Smyth [2015](#ref-lun2015from)),
which provides the important information in a more concise format.
Experienced users (or those looking for some nighttime reading!) are more likely to benefit from the in-depth discussions in this document.

All of the workflows described here start from sorted and indexed BAM files in the *[chipseqDBData](https://bioconductor.org/packages/3.22/chipseqDBData)* package.
For application to user-specified data, the raw read sequences have to be aligned to the appropriate reference genome beforehand.
Most aligners can be used for this purpose, but we have used *[Rsubread](https://bioconductor.org/packages/3.22/Rsubread)* (Liao, Smyth, and Shi [2013](#ref-liao2013)) due to the convenience of its R interface.
It is also recommended to mark duplicate reads using tools like `Picard` prior to starting the workflow.

The statistical methods described here are based upon those in the *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* package (Robinson, McCarthy, and Smyth [2010](#ref-robinson2010)).
Knowledge of *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* is useful but not a prerequesite for reading this guide.

## 1.3 How to get help

Most questions about *[csaw](https://bioconductor.org/packages/3.22/csaw)* should be answered by the documentation.
Every function mentioned in this guide has its own help page.
For example, a detailed description of the arguments and output of the `windowCounts()` function can be obtained by typing `?windowCounts` or `help(windowCounts)` at the R prompt.
Further detail on the methods or the underlying theory can be found in the references at the bottom of each help page.

The authors of the package always appreciate receiving reports of bugs in the package functions or in the documentation.
The same goes for well-considered suggestions for improvements.
Other questions about how to use *[csaw](https://bioconductor.org/packages/3.22/csaw)* are best sent to the [Bioconductor support site](https://support.bioconductor.org).
Please send requests for general assistance and advice to the support site, rather than to the individual authors.
Users posting to the support site for the first time may find it helpful to read the [posting guide](http://www.bioconductor.org/help/support/posting-guide).

## 1.4 How to cite this book

Most users of *[csaw](https://bioconductor.org/packages/3.22/csaw)* should cite the following in any publications:

> A. T. Lun and G. K. Smyth.
> csaw: a Bioconductor package for differential binding analysis of ChIP-seq data using sliding windows.
> *Nucleic Acids Res.*, 44(5):e45, Mar 2016

To cite the workflows specifically, we can use:

> A. T. L. Lun and G. K. Smyth.
> From reads to regions: a Bioconductor workflow to detect differential binding in ChIP-seq data.
> *F1000Research*, 4, 2015

For people interested in combined p-values, their use in DB analyses was proposed in:

> A. T. Lun and G. K. Smyth.
> De novo detection of differentially bound regions for ChIP-seq data using peaks and windows: controlling error rates correctly.
> *Nucleic Acids Res.*, 42(11):e95, Jul 2014

The DB analyses shown here use methods from the *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* package, which has its own citation recommendations.
See the appropriate section of the *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* user’s guide for more details.

## 1.5 Quick start

A typical ChIP-seq analysis in *[csaw](https://bioconductor.org/packages/3.22/csaw)* would look something like that described below.
This assumes that a vector of file paths to sorted and indexed BAM files is provided in and a design matrix in supplied in .
The code is split across several steps:

```
library(chipseqDBData)
tf.data <- NFYAData()
tf.data <- head(tf.data, -1) # skip the input.
bam.files <- tf.data$Path

cell.type <- sub("NF-YA ([^ ]+) .*", "\\1", tf.data$Description)
design <- model.matrix(~factor(cell.type))
colnames(design) <- c("intercept", "cell.type")
```

1. Loading in data from BAM files.

   ```
   library(csaw)
   param <- readParam(minq=20)
   data <- windowCounts(bam.files, ext=110, width=10, param=param)
   ```
2. Filtering out uninteresting regions.

   ```
   binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
   keep <- filterWindowsGlobal(data, binned)$filter > log2(5)
   data <- data[keep,]
   ```
3. Calculating normalization factors.

   ```
   data <- normFactors(binned, se.out=data)
   ```
4. Identifying DB windows.

   ```
   library(edgeR)
   y <- asDGEList(data)
   y <- estimateDisp(y, design)
   fit <- glmQLFit(y, design, robust=TRUE)
   results <- glmQLFTest(fit)
   ```
5. Correcting for multiple testing.

   ```
   merged <- mergeResults(data, results$table, tol=1000L)
   ```

## Session information

View session info

### Bibliography

Liao, Y., G. K. Smyth, and W. Shi. 2013. “The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote.” *Nucleic Acids Res.* 41 (10): e108.

Lun, A. T. L., and G. K. Smyth. 2015. “From Reads to Regions: A Bioconductor Workflow to Detect Differential Binding in ChIP-Seq Data.” *F1000Research* 4.

Lun, A. 2016. “csaw: a Bioconductor package for differential binding analysis of ChIP-seq data using sliding windows.” *Nucleic Acids Res.* 44 (5): e45.

Pal, B., T. Bouras, W. Shi, F. Vaillant, J. M. Sheridan, N. Fu, K. Breslin, et al. 2013. “Global changes in the mammary epigenome are induced by hormonal cues and coordinated by Ezh2.” *Cell Rep.* 3 (2): 411–26.

Robinson, M. D., D. J. McCarthy, and G. K. Smyth. 2010. “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” *Bioinformatics* 26 (1): 139–40.

Ross-Innes, C. S., R. Stark, A. E. Teschendorff, K. A. Holmes, H. R. Ali, M. J. Dunning, G. D. Brown, et al. 2012. “Differential oestrogen receptor binding is associated with clinical outcome in breast cancer.” *Nature* 481 (7381): 389–93.