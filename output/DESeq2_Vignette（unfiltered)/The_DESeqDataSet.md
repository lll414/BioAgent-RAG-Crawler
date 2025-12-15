# The DESeqDataSet

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

The object class used by the DESeq2 package to store the read counts and the intermediate estimated quantities during statistical analysis is the *DESeqDataSet*, which will usually be represented in the code here as an object `dds`.

A technical detail is that the *DESeqDataSet* class extends the *RangedSummarizedExperiment* class of the [SummarizedExperiment](http://bioconductor.org/packages/SummarizedExperiment) package. The âRangedâ part refers to the fact that the rows of the assay data (here, the counts) can be associated with genomic ranges (the exons of genes). This association facilitates downstream exploration of results, making use of other Bioconductor packagesâ range-based functionality (e.g.Â find the closest ChIP-seq peaks to the differentially expressed genes).

A *DESeqDataSet* object must have an associated *design formula*. The design formula expresses the variables which will be used in modeling. The formula should be a tilde (~) followed by the variables with plus signs between them (it will be coerced into an *formula* if it is not already). The design can be changed later, however then all differential analysis steps should be repeated, as the design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model.

*Note*: In order to benefit from the default settings of the package, you should put the variable of interest at the end of the formula and make sure the control level is the first level.

We will now show 4 ways of constructing a *DESeqDataSet*, depending on what pipeline was used upstream of DESeq2 to generated counts or estimated counts:

1. From [transcript abundance files and tximport](#tximport)
2. From a [count matrix](#countmat)
3. From [htseq-count files](#htseq)
4. From a [SummarizedExperiment](#se) object