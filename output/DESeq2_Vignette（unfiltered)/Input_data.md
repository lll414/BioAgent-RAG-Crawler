# Input data

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

### Why un-normalized counts?

As input, the DESeq2 package expects count data as obtained, e.g., from RNA-seq or another high-throughput sequencing experiment, in the form of a matrix of integer values. The value in the *i*-th row and the *j*-th column of the matrix tells how many reads can be assigned to gene *i* in sample *j*. Analogously, for other types of assays, the rows of the matrix might correspond e.g.Â to binding regions (with ChIP-Seq) or peptide sequences (with quantitative mass spectrometry). We will list method for obtaining count matrices in sections below.

The values in the matrix should be un-normalized counts or estimated counts of sequencing reads (for single-end RNA-seq) or fragments (for paired-end RNA-seq). The [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/) describes multiple techniques for preparing such count matrices. It is important to provide count matrices as input for DESeq2âs statistical model (Love, Huber, and Anders 2014) to hold, as only the count values allow assessing the measurement precision correctly. The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input.

### The DESeqDataSet

The object class used by the DESeq2 package to store the read counts and the intermediate estimated quantities during statistical analysis is the *DESeqDataSet*, which will usually be represented in the code here as an object `dds`.

A technical detail is that the *DESeqDataSet* class extends the *RangedSummarizedExperiment* class of the [SummarizedExperiment](http://bioconductor.org/packages/SummarizedExperiment) package. The âRangedâ part refers to the fact that the rows of the assay data (here, the counts) can be associated with genomic ranges (the exons of genes). This association facilitates downstream exploration of results, making use of other Bioconductor packagesâ range-based functionality (e.g.Â find the closest ChIP-seq peaks to the differentially expressed genes).

A *DESeqDataSet* object must have an associated *design formula*. The design formula expresses the variables which will be used in modeling. The formula should be a tilde (~) followed by the variables with plus signs between them (it will be coerced into an *formula* if it is not already). The design can be changed later, however then all differential analysis steps should be repeated, as the design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model.

*Note*: In order to benefit from the default settings of the package, you should put the variable of interest at the end of the formula and make sure the control level is the first level.

We will now show 4 ways of constructing a *DESeqDataSet*, depending on what pipeline was used upstream of DESeq2 to generated counts or estimated counts:

1. From [transcript abundance files and tximport](#tximport)
2. From a [count matrix](#countmat)
3. From [htseq-count files](#htseq)
4. From a [SummarizedExperiment](#se) object

### Transcript abundance files and *tximport* / *tximeta*

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

### Tximeta for import with automatic metadata

Another Bioconductor package, [tximeta](https://bioconductor.org/packages/tximeta) (Love et al. 2020), extends *tximport*, offering the same functionality, plus the additional benefit of automatic addition of annotation metadata for commonly used transcriptomes (GENCODE, Ensembl, RefSeq for human and mouse). See the [tximeta](https://bioconductor.org/packages/tximeta) package vignette for more details. *tximeta* produces a *SummarizedExperiment* that can be loaded easily into *DESeq2* using the `DESeqDataSet` function, with an example in the *tximeta* package vignette, and below:

```
coldata <- samples
coldata$files <- files
coldata$names <- coldata$run
```

```
library("tximeta")
se <- tximeta(coldata)
ddsTxi <- DESeqDataSet(se, design = ~ condition)
```

The `ddsTxi` object here can then be used as `dds` in the following analysis steps. If *tximeta* recognized the reference transcriptome as one of those with a pre-computed hashed checksum, the `rowRanges` of the `dds` object will be pre-populated. Again, see the *tximeta* vignette for full details.

### Count matrix input

Alternatively, the function *DESeqDataSetFromMatrix* can be used if you already have a matrix of read counts prepared from another source. Another method for quickly producing count matrices from alignment files is the *featureCounts* function (Liao, Smyth, and Shi 2013) in the [Rsubread](http://bioconductor.org/packages/Rsubread) package. To use *DESeqDataSetFromMatrix*, the user should provide the counts matrix, the information about the samples (the columns of the count matrix) as a *DataFrame* or *data.frame*, and the design formula.

To demonstrate the use of *DESeqDataSetFromMatrix*, we will read in count data from the [pasilla](http://bioconductor.org/packages/pasilla) package, which have been copied into the *DESeq2* package in the `inst/extdata` directory, or the `extdata` directory for the installed package. We read in a count matrix, which we will name `cts`, and the sample information table, which we will name `coldata`. Further below we describe how to extract these objects from, e.g.Â *featureCounts* output.

```
library("DESeq2")
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv.gz",
                      package="DESeq2", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="DESeq2", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
```

We examine the count matrix and column data to see if they are consistent in terms of sample order.

```
head(cts,2)
```

```
##             untreated1 untreated2 untreated3 untreated4 treated1 treated2
## FBgn0000003          0          0          0          0        0        0
## FBgn0000008         92        161         76         70      140       88
##             treated3
## FBgn0000003        1
## FBgn0000008       70
```

```
coldata
```

```
##              condition        type
## treated1fb     treated single-read
## treated2fb     treated  paired-end
## treated3fb     treated  paired-end
## untreated1fb untreated single-read
## untreated2fb untreated single-read
## untreated3fb untreated  paired-end
## untreated4fb untreated  paired-end
```

Note that these are not in the same order with respect to samples!

It is absolutely critical that the columns of the count matrix and the rows of the column data (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order.

As they are not in the correct order as given, we need to re-arrange one or the other so that they are consistent in terms of sample order (if we do not, later functions would produce an error). We additionally need to chop off the `"fb"` of the row names of `coldata`, so the naming is consistent.

```
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
```

```
## [1] TRUE
```

```
all(rownames(coldata) == colnames(cts))
```

```
## [1] FALSE
```

```
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
```

```
## [1] TRUE
```

If you have used the *featureCounts* function (Liao, Smyth, and Shi 2013) in the [Rsubread](http://bioconductor.org/packages/Rsubread) package, the matrix of read counts can be directly provided from the `"counts"` element in the list output. The count matrix and column data can typically be read into R from flat files using base R functions such as *read.csv* or *read.delim*. For *htseq-count* files, see the dedicated input function below.

With the count matrix, `cts`, and the sample information, `coldata`, we can construct a *DESeqDataSet*:

```
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
```

```
## class: DESeqDataSet 
## dim: 14599 7 
## metadata(1): version
## assays(1): counts
## rownames(14599): FBgn0000003 FBgn0000008 ... FBgn0261574 FBgn0261575
## rowData names(0):
## colnames(7): treated1 treated2 ... untreated3 untreated4
## colData names(2): condition type
```

If you have additional feature data, it can be added to the *DESeqDataSet* by adding to the metadata columns of a newly constructed object. (Here we add redundant data just for demonstration, as the gene names are already the rownames of the `dds`.)

```
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
```

```
## DataFrame with 14599 rows and 1 column
##                    gene
##             <character>
## FBgn0000003 FBgn0000003
## FBgn0000008 FBgn0000008
## FBgn0000014 FBgn0000014
## FBgn0000015 FBgn0000015
## FBgn0000017 FBgn0000017
## ...                 ...
## FBgn0261571 FBgn0261571
## FBgn0261572 FBgn0261572
## FBgn0261573 FBgn0261573
## FBgn0261574 FBgn0261574
## FBgn0261575 FBgn0261575
```

### *htseq-count* input

You can use the function *DESeqDataSetFromHTSeqCount* if you have used *htseq-count* from the [HTSeq](http://www-huber.embl.de/users/anders/HTSeq) python package (Anders, Pyl, and Huber 2014). For an example of using the python scripts, see the [pasilla](http://bioconductor.org/packages/pasilla) data package. First you will want to specify a variable which points to the directory in which the *htseq-count* output files are located.

```
directory <- "/path/to/your/files/"
```

However, for demonstration purposes only, the following code points to the directory for the demo *htseq-count* output files packages for the [pasilla](http://bioconductor.org/packages/pasilla) package.

We specify which files to read in using *list.files*, and select those files which contain the string `"treated"` using *grep*. The *sub* function is used to chop up the sample filename to obtain the condition status, or you might alternatively read in a phenotypic table using *read.table*.

```
# un-evaluated
directory <- system.file("extdata", package="pasilla",
                         mustWork=TRUE)
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)
library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq
```

### *SummarizedExperiment* input

If one has already created or obtained a *SummarizedExperiment*, it can be easily input into DESeq2 as follows. First we load the package containing the `airway` dataset.

```
library("airway")
data("airway")
se <- airway
```

The constructor function below shows the generation of a *DESeqDataSet* from a *RangedSummarizedExperiment* `se`.

```
library("DESeq2")
ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
ddsSE
```

```
## class: DESeqDataSet 
## dim: 63677 8 
## metadata(2): '' version
## assays(1): counts
## rownames(63677): ENSG00000000003 ENSG00000000005 ... ENSG00000273492
##   ENSG00000273493
## rowData names(10): gene_id gene_name ... seq_coord_system symbol
## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
## colData names(9): SampleName cell ... Sample BioSample
```

### Pre-filtering

While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the `dds` data object, and we increase the speed of count modeling within DESeq2. It can also improve visualizations, as features with no information for differential expression are not plotted in dispersion plots or MA-plots.

Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. The count of 10 is a reasonable choice for bulk RNA-seq. A recommendation for the minimal number of samples is to specify the smallest group size, e.g.Â here there are 3 treated samples. If there are not discrete groups, one can use the minimal number of samples where non-zero counts would be considered interesting. One can also omit this step entirely and just rely on the independent filtering procedures available in `results()`, either *IHW* or *genefilter*. See [independent filtering](#indfilt) section.

```
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
```

### Note on factor levels

By default, R will choose a *reference level* for factors based on alphabetical order. Then, if you never tell the DESeq2 functions which level you want to compare against (e.g.Â which level represents the control group), the comparisons will be based on the alphabetical order of the levels. There are two solutions: you can either explicitly tell *results* which comparison to make using the `contrast` argument (this will be shown later), or you can explicitly set the factors levels. In order to see the change of reference levels reflected in the results names, you need to either run `DESeq` or `nbinomWaldTest`/`nbinomLRT` after the re-leveling operation. Setting the factor levels can be done in two ways, either using factor:

```
dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
```

â¦or using *relevel*, just specifying the reference level:

```
dds$condition <- relevel(dds$condition, ref = "untreated")
```

If you need to subset the columns of a *DESeqDataSet*, i.e., when removing certain samples from the analysis, it is possible that all the samples for one or more levels of a variable in the design formula would be removed. In this case, the *droplevels* function can be used to remove those levels which do not have samples in the current *DESeqDataSet*:

```
dds$condition <- droplevels(dds$condition)
```

### Collapsing technical replicates

DESeq2 provides a function *collapseReplicates* which can assist in combining the counts from technical replicates into single columns of the count matrix. The term *technical replicate* implies multiple sequencing runs of the same library. You should not collapse biological replicates using this function. See the manual page for an example of the use of *collapseReplicates*.

### About the pasilla dataset

We continue with the [pasilla](http://bioconductor.org/packages/pasilla) data constructed from the count matrix method above. This data set is from an experiment on *Drosophila melanogaster* cell cultures and investigated the effect of RNAi knock-down of the splicing factor *pasilla* (Brooks et al. 2011). The detailed transcript of the production of the [pasilla](http://bioconductor.org/packages/pasilla) data is provided in the vignette of the data package [pasilla](http://bioconductor.org/packages/pasilla).