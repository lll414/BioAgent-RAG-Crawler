# Extracting transformed values

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

These transformation functions return an object of class *DESeqTransform* which is a subclass of *RangedSummarizedExperiment*. For ~20 samples, running on a newly created `DESeqDataSet`, *rlog* may take 30 seconds, while *vst* takes less than 1 second. The running times are shorter when using `blind=FALSE` and if the function *DESeq* has already been run, because then it is not necessary to re-estimate the dispersion values. The *assay* function is used to extract the matrix of normalized values.

```
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
```

```
##              treated1  treated2  treated3 untreated1 untreated2 untreated3
## FBgn0000008  7.777746  7.984491  7.765448   7.735026   7.807720   7.997378
## FBgn0000017 11.954934 12.035700 12.028610  12.049838  12.295662  12.471723
## FBgn0000018  9.219162  9.089390  9.028791   9.374577   9.173538   9.052143
##             untreated4
## FBgn0000008   7.832458
## FBgn0000017  12.089675
## FBgn0000018   9.142848
```