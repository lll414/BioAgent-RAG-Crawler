# Extended section on shrinkage estimators

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Here we extend the [discussion of shrinkage estimators](#shrink). Below is a summary table of differences between methods available in `lfcShrink` via the `type` argument (and for further technical reference on use of arguments please see `?lfcShrink`):

| method: | `apeglm`1 | `ashr`2 | `normal`3 |
| --- | --- | --- | --- |
| Good for ranking by LFC | â | â | â |
| Preserves size of large LFC | â | â |  |
| Can compute *s-values* (Stephens 2016) | â | â |  |
| Allows use of `coef` | â | â | â |
| Allows use of `lfcThreshold` | â | â | â |
| Allows use of `contrast` |  | â | â |
| Can shrink interaction terms | â | â |  |

**References:** 1. Zhu, Ibrahim, and Love (2018); 2. Stephens (2016); 3. Love, Huber, and Anders (2014)

Beginning with the first row, all shrinkage methods provided by DESeq2 are good for ranking genes by âeffect sizeâ, that is the log2 fold change (LFC) across groups, or associated with an interaction term. It is useful to contrast ranking by effect size with ranking by a p-value or adjusted p-value associated with a null hypothesis: while increasing the number of samples will tend to decrease the associated p-value for a gene that is differentially expressed, the estimated effect size or LFC becomes more precise. Also, a gene can have a small p-value although the change in expression is not great, as long as the standard error associated with the estimated LFC is small.

The next two rows point out that `apeglm` and `ashr` shrinkage methods help to preserve the size of large LFC, and can be used to compute *s-values*. These properties are related. As noted in the [previous section](#moreshrink), the original DESeq2 shrinkage estimator used a Normal distribution, with a scale that adapts to the spread of the observed LFCs. Because the tails of the Normal distribution become thin relatively quickly, it was important when we designed the method that the prior scaling is sensitive to the very largest observed LFCs. As you can read in the DESeq2 paper, under the section, â*Empirical prior estimate*â, we used the top 5% of the LFCs by absolute value to set the scale of the Normal prior (we later added weighting the quantile by precision). `ashr`, published in 2016, and `apeglm` use wide-tailed priors to avoid shrinking large LFCs. While a typical RNA-seq experiment may have many LFCs between -1 and 1, we might consider a LFC of >4 to be very large, as they represent 16-fold increases or decreases in expression. `ashr` and `apeglm` can adapt to the scale of the entirety of LFCs, while not over-shrinking the few largest LFCs. The potential for over-shrinking LFC is also why DESeq2âs shrinkage estimator is not recommended for designs with interaction terms.

What are *s-values*? This quantity proposed by Stephens (2016) gives the estimated rate of *false sign* among genes with equal or smaller s-value. Stephens (2016) points out they are analogous to the *q*-value of Storey (2003). The s-value has a desirable property relative to the adjusted p-value or *q*-value, in that it does not require supposing there to be a set of null genes with LFC = 0 (the most commonly used null hypothesis). Therefore, it can be benchmarked by comparing estimated LFC and s-value to the âtrue LFCâ in a setting where this can be reasonably defined. For these estimated probabilities to be accurate, the scale of the prior needs to match the scale of the distribution of effect sizes, and so the original DESeq2 shrinkage method is not really compatible with computing s-values.

The last four rows explain differences in whether coefficients or contrasts can have shrinkage applied by the various methods. All three methods can use `coef` with either the name or numeric index from `resultsNames(dds)` to specify which coefficient to shrink. All three methods allow for a positive `lfcThreshold` to be specified, in which case, they will return p-values and adjusted p-values or s-values for the LFC being greater in absolute value than the threshold (see [this section](#thresh) for `normal`). For `apeglm` and `ashr`, setting a threshold means that the s-values will give the âfalse sign or smallâ rate (FSOS) among genes with equal or small s-value. We found FSOS to be a useful description for when the LFC is either the wrong sign or less than the threshold distance from 0.

```
resApeT <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)
```


```
resAshT <- lfcShrink(dds, coef=2, type="ashr", lfcThreshold=1)
plotMA(resAshT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)
```


Finally, `normal` and `ashr` can be used with arbitrary specified `contrast` because `normal` shrinks multiple coefficients simultaneously (`apeglm` does not), and because `ashr` does not estimate a vector of coefficients but models estimated coefficients and their standard errors from upstream methods (here, DESeq2âs MLE). Although `apeglm` cannot be used with `contrast`, we note that many designs can be easily rearranged such that what was a contrast becomes its own coefficient. In this case, the dispersion does not have to be estimated again, as the designs are equivalent, up to the meaning of the coefficients. Instead, one need only run `nbinomWaldTest` to re-estimate MLE coefficients â these are necessary for `apeglm` â and then run `lfcShrink` specifying the coefficient of interest in `resultsNames(dds)`.

We give some examples below of producing equivalent designs for use with `coef`. We show how the coefficients change with `model.matrix`, but the user would, for example, either change the levels of `dds$condition` or replace the design using `design(dds)<-`, then run `nbinomWaldTest` followed by `lfcShrink`.

Three groups:

```
condition <- factor(rep(c("A","B","C"),each=2))
model.matrix(~ condition)
```

```
##   (Intercept) conditionB conditionC
## 1           1          0          0
## 2           1          0          0
## 3           1          1          0
## 4           1          1          0
## 5           1          0          1
## 6           1          0          1
## attr(,"assign")
## [1] 0 1 1
## attr(,"contrasts")
## attr(,"contrasts")$condition
## [1] "contr.treatment"
```

```
# to compare C vs B, make B the reference level,
# and select the last coefficient
condition <- relevel(condition, "B")
model.matrix(~ condition)
```

```
##   (Intercept) conditionA conditionC
## 1           1          1          0
## 2           1          1          0
## 3           1          0          0
## 4           1          0          0
## 5           1          0          1
## 6           1          0          1
## attr(,"assign")
## [1] 0 1 1
## attr(,"contrasts")
## attr(,"contrasts")$condition
## [1] "contr.treatment"
```

Three groups, compare condition effects:

```
grp <- factor(rep(1:3,each=4))
cnd <- factor(rep(rep(c("A","B"),each=2),3))
model.matrix(~ grp + cnd + grp:cnd)
```

```
##    (Intercept) grp2 grp3 cndB grp2:cndB grp3:cndB
## 1            1    0    0    0         0         0
## 2            1    0    0    0         0         0
## 3            1    0    0    1         0         0
## 4            1    0    0    1         0         0
## 5            1    1    0    0         0         0
## 6            1    1    0    0         0         0
## 7            1    1    0    1         1         0
## 8            1    1    0    1         1         0
## 9            1    0    1    0         0         0
## 10           1    0    1    0         0         0
## 11           1    0    1    1         0         1
## 12           1    0    1    1         0         1
## attr(,"assign")
## [1] 0 1 1 2 3 3
## attr(,"contrasts")
## attr(,"contrasts")$grp
## [1] "contr.treatment"
## 
## attr(,"contrasts")$cnd
## [1] "contr.treatment"
```

```
# to compare condition effect in group 3 vs 2,
# make group 2 the reference level,
# and select the last coefficient
grp <- relevel(grp, "2")
model.matrix(~ grp + cnd + grp:cnd)
```

```
##    (Intercept) grp1 grp3 cndB grp1:cndB grp3:cndB
## 1            1    1    0    0         0         0
## 2            1    1    0    0         0         0
## 3            1    1    0    1         1         0
## 4            1    1    0    1         1         0
## 5            1    0    0    0         0         0
## 6            1    0    0    0         0         0
## 7            1    0    0    1         0         0
## 8            1    0    0    1         0         0
## 9            1    0    1    0         0         0
## 10           1    0    1    0         0         0
## 11           1    0    1    1         0         1
## 12           1    0    1    1         0         1
## attr(,"assign")
## [1] 0 1 1 2 3 3
## attr(,"contrasts")
## attr(,"contrasts")$grp
## [1] "contr.treatment"
## 
## attr(,"contrasts")$cnd
## [1] "contr.treatment"
```

Two groups, two individuals per group, compare within-individual condition effects:

```
grp <- factor(rep(1:2,each=4))
ind <- factor(rep(rep(1:2,each=2),2))
cnd <- factor(rep(c("A","B"),4))
model.matrix(~grp + grp:ind + grp:cnd)
```

```
##   (Intercept) grp2 grp1:ind2 grp2:ind2 grp1:cndB grp2:cndB
## 1           1    0         0         0         0         0
## 2           1    0         0         0         1         0
## 3           1    0         1         0         0         0
## 4           1    0         1         0         1         0
## 5           1    1         0         0         0         0
## 6           1    1         0         0         0         1
## 7           1    1         0         1         0         0
## 8           1    1         0         1         0         1
## attr(,"assign")
## [1] 0 1 2 2 3 3
## attr(,"contrasts")
## attr(,"contrasts")$grp
## [1] "contr.treatment"
## 
## attr(,"contrasts")$ind
## [1] "contr.treatment"
## 
## attr(,"contrasts")$cnd
## [1] "contr.treatment"
```

```
# to compare condition effect across group,
# add a main effect for 'cnd',
# and select the last coefficient
model.matrix(~grp + cnd + grp:ind + grp:cnd)
```

```
##   (Intercept) grp2 cndB grp1:ind2 grp2:ind2 grp2:cndB
## 1           1    0    0         0         0         0
## 2           1    0    1         0         0         0
## 3           1    0    0         1         0         0
## 4           1    0    1         1         0         0
## 5           1    1    0         0         0         0
## 6           1    1    1         0         0         1
## 7           1    1    0         0         1         0
## 8           1    1    1         0         1         1
## attr(,"assign")
## [1] 0 1 2 3 3 4
## attr(,"contrasts")
## attr(,"contrasts")$grp
## [1] "contr.treatment"
## 
## attr(,"contrasts")$cnd
## [1] "contr.treatment"
## 
## attr(,"contrasts")$ind
## [1] "contr.treatment"
```