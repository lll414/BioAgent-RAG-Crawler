# Levels without samples

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

The base R function for creating model matrices will produce a column of zeros if a level is missing from a factor or a combination of levels is missing from an interaction of factors. The solution to the first case is to call *droplevels* on the column, which will remove levels without samples. This was shown in the beginning of this vignette.

The second case is also solvable, by manually editing the model matrix, and then providing this to *DESeq*. Here we construct an example dataset to illustrate:

```
group <- factor(rep(1:3,each=6))
condition <- factor(rep(rep(c("A","B","C"),each=2),3))
d <- DataFrame(group, condition)[-c(17,18),]
as.data.frame(d)
```

```
##    group condition
## 1      1         A
## 2      1         A
## 3      1         B
## 4      1         B
## 5      1         C
## 6      1         C
## 7      2         A
## 8      2         A
## 9      2         B
## 10     2         B
## 11     2         C
## 12     2         C
## 13     3         A
## 14     3         A
## 15     3         B
## 16     3         B
```

Note that if we try to estimate all interaction terms, we introduce a column with all zeros, as there are no condition C samples for group 3. (Here, *unname* is used to display the matrix concisely.)

```
m1 <- model.matrix(~ condition*group, d)
colnames(m1)
```

```
## [1] "(Intercept)"       "conditionB"        "conditionC"       
## [4] "group2"            "group3"            "conditionB:group2"
## [7] "conditionC:group2" "conditionB:group3" "conditionC:group3"
```

```
unname(m1)
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
##  [1,]    1    0    0    0    0    0    0    0    0
##  [2,]    1    0    0    0    0    0    0    0    0
##  [3,]    1    1    0    0    0    0    0    0    0
##  [4,]    1    1    0    0    0    0    0    0    0
##  [5,]    1    0    1    0    0    0    0    0    0
##  [6,]    1    0    1    0    0    0    0    0    0
##  [7,]    1    0    0    1    0    0    0    0    0
##  [8,]    1    0    0    1    0    0    0    0    0
##  [9,]    1    1    0    1    0    1    0    0    0
## [10,]    1    1    0    1    0    1    0    0    0
## [11,]    1    0    1    1    0    0    1    0    0
## [12,]    1    0    1    1    0    0    1    0    0
## [13,]    1    0    0    0    1    0    0    0    0
## [14,]    1    0    0    0    1    0    0    0    0
## [15,]    1    1    0    0    1    0    0    1    0
## [16,]    1    1    0    0    1    0    0    1    0
## attr(,"assign")
## [1] 0 1 1 2 2 3 3 3 3
## attr(,"contrasts")
## attr(,"contrasts")$condition
## [1] "contr.treatment"
## 
## attr(,"contrasts")$group
## [1] "contr.treatment"
```

```
all.zero <- apply(m1, 2, function(x) all(x==0))
all.zero
```

```
##       (Intercept)        conditionB        conditionC            group2 
##             FALSE             FALSE             FALSE             FALSE 
##            group3 conditionB:group2 conditionC:group2 conditionB:group3 
##             FALSE             FALSE             FALSE             FALSE 
## conditionC:group3 
##              TRUE
```

We can remove this column like so:

```
idx <- which(all.zero)
m1 <- m1[,-idx]
unname(m1)
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
##  [1,]    1    0    0    0    0    0    0    0
##  [2,]    1    0    0    0    0    0    0    0
##  [3,]    1    1    0    0    0    0    0    0
##  [4,]    1    1    0    0    0    0    0    0
##  [5,]    1    0    1    0    0    0    0    0
##  [6,]    1    0    1    0    0    0    0    0
##  [7,]    1    0    0    1    0    0    0    0
##  [8,]    1    0    0    1    0    0    0    0
##  [9,]    1    1    0    1    0    1    0    0
## [10,]    1    1    0    1    0    1    0    0
## [11,]    1    0    1    1    0    0    1    0
## [12,]    1    0    1    1    0    0    1    0
## [13,]    1    0    0    0    1    0    0    0
## [14,]    1    0    0    0    1    0    0    0
## [15,]    1    1    0    0    1    0    0    1
## [16,]    1    1    0    0    1    0    0    1
```

Now this matrix `m1` can be provided to the `full` argument of *DESeq*. For a likelihood ratio test of interactions, a model matrix using a reduced design such as `~ condition + group` can be given to the `reduced` argument. Wald tests can also be generated instead of the likelihood ratio test, but for user-supplied model matrices, the argument `betaPrior` must be set to `FALSE`.