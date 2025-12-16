# Multi-factor designs

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Experiments with more than one factor influencing the counts can be analyzed using design formulas that include the additional variables. In fact, DESeq2 can analyze any possible experimental design that can be expressed with fixed effects terms (multiple factors, designs with interactions, designs with continuous variables, splines, and so on are all possible).

By adding variables to the design, one can control for additional variation in the counts. For example, if the condition samples are balanced across experimental batches, by including the `batch` factor to the design, one can increase the sensitivity for finding differences due to `condition`. There are multiple ways to analyze experiments when the additional variables are of interest and not just controlling factors (see [section on interactions](#interactions)).

**Experiments with many samples**: in experiments with many samples (e.g.Â 50, 100, etc.) it is highly likely that there will be technical variation affecting the observed counts. Failing to model this additional technical variation will lead to spurious results. Many methods exist that can be used to model technical variation, which can be easily included in the DESeq2 design to control for technical variation while estimating effects of interest. See the [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene) for examples of using RUV or SVA in combination with DESeq2. For more details on why it is important to control for technical variation in large sample experiments, see the following [thread](https://twitter.com/mikelove/status/1513468597288452097), also archived [here](https://htmlpreview.github.io/?https://github.com/frederikziebell/science_tweetorials/blob/master/DESeq2_many_samples.html) by Frederik Ziebell.

The data in the [pasilla](http://bioconductor.org/packages/pasilla) package have a condition of interest (the column `condition`), as well as information on the type of sequencing which was performed (the column `type`), as we can see below:

```
colData(dds)
```

```
## DataFrame with 7 rows and 3 columns
##            condition        type sizeFactor
##             <factor>    <factor>  <numeric>
## treated1   treated   single-read   1.629707
## treated2   treated   paired-end    0.761162
## treated3   treated   paired-end    0.830312
## untreated1 untreated single-read   1.143904
## untreated2 untreated single-read   1.791281
## untreated3 untreated paired-end    0.645994
## untreated4 untreated paired-end    0.750728
```

We create a copy of the *DESeqDataSet*, so that we can rerun the analysis using a multi-factor design.

```
ddsMF <- dds
```

We change the levels of `type` so it only contains letters (numbers, underscore and period are also allowed in design factor levels). Be careful when changing level names to use the same order as the current levels.

```
levels(ddsMF$type)
```

```
## [1] "paired-end"  "single-read"
```

```
levels(ddsMF$type) <- sub("-.*", "", levels(ddsMF$type))
levels(ddsMF$type)
```

```
## [1] "paired" "single"
```

We can account for the different types of sequencing, and get a clearer picture of the differences attributable to the treatment. As `condition` is the variable of interest, we put it at the end of the formula. Thus the *results* function will by default pull the `condition` results unless `contrast` or `name` arguments are specified.

Then we can re-run *DESeq*:

```
design(ddsMF) <- formula(~ type + condition)
ddsMF <- DESeq(ddsMF)
```

Again, we access the results using the *results* function.

```
resMF <- results(ddsMF)
head(resMF)
```

```
## log2 fold change (MLE): condition treated vs untreated 
## Wald test p-value: condition treated vs untreated 
## DataFrame with 6 rows and 6 columns
##               baseMean log2FoldChange     lfcSE      stat    pvalue      padj
##              <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
## FBgn0000008   95.28865     -0.0390130  0.218997 -0.178144 0.8586100  0.947833
## FBgn0000017 4359.09632     -0.2548984  0.113535 -2.245099 0.0247617  0.131475
## FBgn0000018  419.06811     -0.0625571  0.129956 -0.481372 0.6302523  0.852180
## FBgn0000024    6.41105      0.3097331  0.750231  0.412850 0.6797164  0.877741
## FBgn0000032  990.79225     -0.0465134  0.120215 -0.386918 0.6988171  0.886082
## FBgn0000037   14.11443      0.4541562  0.523436  0.867644 0.3855893  0.691941
```

It is also possible to retrieve the log2 fold changes, *p* values and adjusted *p* values of variables other than the last one in the design. While in this case, `type` is not biologically interesting as it indicates differences across sequencing protocol, for other hypothetical designs, such as `~genotype + condition + genotype:condition`, we may actually be interested in the difference in baseline expression across genotype, which is not the last variable in the design.

In any case, the `contrast` argument of the function *results* takes a character vector of length three: the name of the variable, the name of the factor level for the numerator of the log2 ratio, and the name of the factor level for the denominator. The `contrast` argument can also take other forms, as described in the help page for *results* and [below](#contrasts)

```
resMFType <- results(ddsMF,
                     contrast=c("type", "single", "paired"))
head(resMFType)
```

```
## log2 fold change (MLE): type single vs paired 
## Wald test p-value: type single vs paired 
## DataFrame with 6 rows and 6 columns
##               baseMean log2FoldChange     lfcSE      stat    pvalue      padj
##              <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
## FBgn0000008   95.28865      -0.265123  0.217459 -1.219189 0.2227725  0.492729
## FBgn0000017 4359.09632      -0.103203  0.113399 -0.910090 0.3627748  0.642817
## FBgn0000018  419.06811       0.225857  0.128864  1.752669 0.0796589  0.271610
## FBgn0000024    6.41105       0.302083  0.745703  0.405099 0.6854049        NA
## FBgn0000032  990.79225       0.233891  0.119646  1.954855 0.0506002  0.206081
## FBgn0000037   14.11443      -0.053260  0.521939 -0.102043 0.9187228  0.969866
```

If the variable is continuous or an interaction term (see [section on interactions](#interactions)) then the results can be extracted using the `name` argument to *results*, where the name is one of elements returned by `resultsNames(dds)`.