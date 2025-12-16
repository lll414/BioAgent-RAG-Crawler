# Note on factor levels

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

By default, R will choose a *reference level* for factors based on alphabetical order. Then, if you never tell the DESeq2 functions which level you want to compare against (e.g.Â which level represents the control group), the comparisons will be based on the alphabetical order of the levels. There are two solutions: you can either explicitly tell *results* which comparison to make using the `contrast` argument (this will be shown later), or you can explicitly set the factors levels. In order to see the change of reference levels reflected in the results names, you need to either run `DESeq` or `nbinomWaldTest`/`nbinomLRT` after the re-leveling operation. Setting the factor levels can be done in two ways, either using factor:

```
dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
```

â¦or using *relevel*, just specifying the reference level:

```
dds$condition <- relevel(dds$condition, ref = "untreated")
```

If you need to subset the columns of a *DESeqDataSet*, i.e., when removing certain samples from the analysis, it is possible that all the samples for one or more levels of a variable in the design formula would be removed. In this case, the *droplevels* function can be used to remove those levels which do not have samples in the current *DESeqDataSet*:

```
dds$condition <- droplevels(dds$condition)
```