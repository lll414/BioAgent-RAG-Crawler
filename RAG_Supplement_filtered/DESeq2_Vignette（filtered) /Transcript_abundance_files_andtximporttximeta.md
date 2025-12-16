# Transcript abundance files andtximport/tximeta

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Our recommended pipeline for *DESeq2* is to use fast transcript abundance quantifiers upstream of DESeq2, and then to create gene-level count matrices for use with DESeq2 by importing the quantification data using [tximport](http://bioconductor.org/packages/tximport) (Soneson, Love, and Robinson 2015). This workflow allows users to import transcript abundance estimates from a variety of external software, including the following methods:

* [Salmon](http://combine-lab.github.io/salmon/) (Patro et al. 2017)
* [Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish/) (Patro, Mount, and Kingsford 2014)
* [kallisto](https://pachterlab.github.io/kallisto/about.html) (Bray et al. 2016)
* [RSEM](http://deweylab.github.io/RSEM/) (Li and Dewey 2011)

Some advantages of using the above methods for transcript abundance estimation are: (i) this approach corrects for potential changes in gene length across samples (e.g.Â from differential isoform usage) (Trapnell et al. 2013), (ii) some of these methods (*Salmon*, *Sailfish*, *kallisto*) are substantially faster and require less memory and disk usage compared to alignment-based methods that require creation and storage of BAM files, and (iii) it is possible to avoid discarding those fragments that can align to multiple genes with homologous sequence, thus increasing sensitivity (Robert and Watson 2015).

Full details on the motivation and methods for importing transcript level abundance and count estimates, summarizing to gene-level count matrices and producing an offset which corrects for potential changes in average transcript length across samples are described in (Soneson, Love, and Robinson 2015). Note that the tximport-to-DESeq2 approach uses *estimated* gene counts from the transcript abundance quantifiers, but not *normalized* counts.

A tutorial on how to use the *Salmon* software for quantifying transcript abundance can be found [here](https://combine-lab.github.io/salmon/getting_started/). We recommend using the `--gcBias` [flag](http://salmon.readthedocs.io/en/latest/salmon.html#gcbias) which estimates a correction factor for systematic biases commonly present in RNA-seq data (Love, Hogenesch, and Irizarry 2016; Patro et al. 2017), unless you are certain that your data do not contain such bias.

Here, we demonstrate how to import transcript abundances and construct a gene-level *DESeqDataSet* object from *Salmon* `quant.sf` files, which are stored in the [tximportData](http://bioconductor.org/packages/tximportData) package. You do not need the `tximportData` package for your analysis, it is only used here for demonstration.

Note that, instead of locating `dir` using *system.file*, a user would typically just provide a path, e.g.Â `/path/to/quant/files`. For a typical use, the `condition` information should already be present as a column of the sample table `samples`, while here we construct artificial condition labels for demonstration.

```
library("tximport")
library("readr")
library("tximportData")
dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
samples$condition <- factor(rep(c("A","B"),each=3))
rownames(samples) <- samples$run
samples[,c("pop","center","run","condition")]
```

```
##           pop center       run condition
## ERR188297 TSI  UNIGE ERR188297         A
## ERR188088 TSI  UNIGE ERR188088         A
## ERR188329 TSI  UNIGE ERR188329         A
## ERR188288 TSI  UNIGE ERR188288         B
## ERR188021 TSI  UNIGE ERR188021         B
## ERR188356 TSI  UNIGE ERR188356         B
```

Next we specify the path to the files using the appropriate columns of `samples`, and we read in a table that links transcripts to genes for this dataset.

```
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- samples$run
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
```

We import the necessary quantification data for DESeq2 using the *tximport* function. For further details on use of *tximport*, including the construction of the `tx2gene` table for linking transcripts to genes in your dataset, please refer to the [tximport](http://bioconductor.org/packages/tximport) package vignette.

```
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
```

Finally, we can construct a *DESeqDataSet* from the `txi` object and sample information in `samples`.

```
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
```

The `ddsTxi` object here can then be used as `dds` in the following analysis steps.