# htseq-countinput

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

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