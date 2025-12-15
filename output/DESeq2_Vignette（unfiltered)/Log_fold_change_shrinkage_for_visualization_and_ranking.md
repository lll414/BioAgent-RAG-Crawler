# Log fold change shrinkage for visualization and ranking

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. To shrink the LFC, we pass the `dds` object to the function `lfcShrink`. Below we specify to use the *apeglm* method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.

We provide the `dds` object and the name or number of the coefficient we want to shrink, where the number refers to the order of the coefficient as it appears in `resultsNames(dds)`.

```
resultsNames(dds)
```

```
## [1] "Intercept"                      "condition_treated_vs_untreated"
```

```
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC
```

```
## log2 fold change (MAP): condition treated vs untreated 
## Wald test p-value: condition treated vs untreated 
## DataFrame with 8148 rows and 5 columns
##               baseMean log2FoldChange     lfcSE    pvalue      padj
##              <numeric>      <numeric> <numeric> <numeric> <numeric>
## FBgn0000008   95.28865     0.00195376  0.152654 0.9858470  0.996699
## FBgn0000017 4359.09632    -0.18810628  0.120870 0.0606585  0.289604
## FBgn0000018  419.06811    -0.06893831  0.122805 0.4870968  0.822681
## FBgn0000024    6.41105     0.01786546  0.199499 0.7566555  0.939146
## FBgn0000032  990.79225    -0.06001511  0.121962 0.5430003  0.848881
## ...                ...            ...       ...       ...       ...
## FBgn0261564   1160.028     -0.0669829 0.0976567 0.4288481  0.789246
## FBgn0261565    620.388     -0.2284564 0.1362122 0.0361772  0.206423
## FBgn0261570   3212.969      0.2395981 0.1237304 0.0190379  0.133380
## FBgn0261573   2243.936      0.0115395 0.0981689 0.8952617  0.977565
## FBgn0261574   4863.807      0.0101618 0.1417667 0.9262385  0.986726
```

Shrinkage estimation is discussed more in a [later section](#moreshrink).