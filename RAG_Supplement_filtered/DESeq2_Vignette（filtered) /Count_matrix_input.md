# Count matrix input

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

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