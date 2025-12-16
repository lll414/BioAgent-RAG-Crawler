# âModel matrix not full rankâ

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

While most experimental designs run easily using design formula, some design formulas can cause problems and result in the *DESeq* function returning an error with the text: âthe model matrix is not full rank, so the model cannot be fit as specified.â There are two main reasons for this problem: either one or more columns in the model matrix are linear combinations of other columns, or there are levels of factors or combinations of levels of multiple factors which are missing samples. We address these two problems below and discuss possible solutions:

### Linear combinations

The simplest case is the linear combination, or linear dependency problem, when two variables contain exactly the same information, such as in the following sample table. The software cannot fit an effect for `batch` and `condition`, because they produce identical columns in the model matrix. This is also referred to as *perfect confounding*. A unique solution of coefficients (the \(\beta\_i\) in the formula [below](#theory)) is not possible.

```
## DataFrame with 4 rows and 2 columns
##      batch condition
##   <factor>  <factor>
## 1        1         A
## 2        1         A
## 3        2         B
## 4        2         B
```

Another situation which will cause problems is when the variables are not identical, but one variable can be formed by the combination of other factor levels. In the following example, the effect of batch 2 vs 1 cannot be fit because it is identical to a column in the model matrix which represents the condition C vs A effect.

```
## DataFrame with 6 rows and 2 columns
##      batch condition
##   <factor>  <factor>
## 1        1         A
## 2        1         A
## 3        1         B
## 4        1         B
## 5        2         C
## 6        2         C
```

In both of these cases above, the batch effect cannot be fit and must be removed from the model formula. There is just no way to tell apart the condition effects and the batch effects. The options are either to assume there is no batch effect (which we know is highly unlikely given the literature on batch effects in sequencing datasets) or to repeat the experiment and properly balance the conditions across batches. A balanced design would look like:

```
## DataFrame with 6 rows and 2 columns
##      batch condition
##   <factor>  <factor>
## 1        1         A
## 2        1         B
## 3        1         C
## 4        2         A
## 5        2         B
## 6        2         C
```

### Group-specific condition effects, individuals nested within groups

Finally, there is a case where we *can* in fact perform inference, but we may need to re-arrange terms to do so. Consider an experiment with grouped individuals, where we seek to test the group-specific effect of a condition or treatment, while controlling for individual effects. The individuals are nested within the groups: an individual can only be in one of the groups, although each individual has one or more observations across condition.

An example of such an experiment is below:

```
coldata <- DataFrame(grp=factor(rep(c("X","Y"),each=6)),
                     ind=factor(rep(1:6,each=2)),
                     cnd=factor(rep(c("A","B"),6)))
coldata
```

```
## DataFrame with 12 rows and 3 columns
##          grp      ind      cnd
##     <factor> <factor> <factor>
## 1          X        1        A
## 2          X        1        B
## 3          X        2        A
## 4          X        2        B
## 5          X        3        A
## ...      ...      ...      ...
## 8          Y        4        B
## 9          Y        5        A
## 10         Y        5        B
## 11         Y        6        A
## 12         Y        6        B
```

Note that individual (`ind`) is a *factor* not a numeric. This is very important.

To make R display all the rows, we can do:

```
as.data.frame(coldata)
```

```
##    grp ind cnd
## 1    X   1   A
## 2    X   1   B
## 3    X   2   A
## 4    X   2   B
## 5    X   3   A
## 6    X   3   B
## 7    Y   4   A
## 8    Y   4   B
## 9    Y   5   A
## 10   Y   5   B
## 11   Y   6   A
## 12   Y   6   B
```

We have two groups of samples X and Y, each with three distinct individuals (labeled here 1-6). For each individual, we have conditions A and B (for example, this could be control and treated).

This design can be analyzed by DESeq2 but requires a bit of refactoring in order to fit the model terms. Here we will use a trick described in the [edgeR](http://bioconductor.org/packages/edgeR) user guide, from the section *Comparisons Both Between and Within Subjects*. If we try to analyze with a formula such as, `~ ind + grp*cnd`, we will obtain an error, because the effect for group is a linear combination of the individuals.

However, the following steps allow for an analysis of group-specific condition effects, while controlling for differences in individual. For object construction, you can use a simple design, such as `~ ind + cnd`, as long as you remember to replace it before running *DESeq*. Then add a column `ind.n` which distinguishes the individuals nested within a group. Here, we add this column to coldata, but in practice you would add this column to `dds`.

```
coldata$ind.n <- factor(rep(rep(1:3,each=2),2))
as.data.frame(coldata)
```

```
##    grp ind cnd ind.n
## 1    X   1   A     1
## 2    X   1   B     1
## 3    X   2   A     2
## 4    X   2   B     2
## 5    X   3   A     3
## 6    X   3   B     3
## 7    Y   4   A     1
## 8    Y   4   B     1
## 9    Y   5   A     2
## 10   Y   5   B     2
## 11   Y   6   A     3
## 12   Y   6   B     3
```

Now we can reassign our *DESeqDataSet* a design of `~ grp + grp:ind.n + grp:cnd`, before we call *DESeq*. This new design will result in the following model matrix:

```
model.matrix(~ grp + grp:ind.n + grp:cnd, coldata)
```

```
##    (Intercept) grpY grpX:ind.n2 grpY:ind.n2 grpX:ind.n3 grpY:ind.n3 grpX:cndB
## 1            1    0           0           0           0           0         0
## 2            1    0           0           0           0           0         1
## 3            1    0           1           0           0           0         0
## 4            1    0           1           0           0           0         1
## 5            1    0           0           0           1           0         0
## 6            1    0           0           0           1           0         1
## 7            1    1           0           0           0           0         0
## 8            1    1           0           0           0           0         0
## 9            1    1           0           1           0           0         0
## 10           1    1           0           1           0           0         0
## 11           1    1           0           0           0           1         0
## 12           1    1           0           0           0           1         0
##    grpY:cndB
## 1          0
## 2          0
## 3          0
## 4          0
## 5          0
## 6          0
## 7          0
## 8          1
## 9          0
## 10         1
## 11         0
## 12         1
## attr(,"assign")
## [1] 0 1 2 2 2 2 3 3
## attr(,"contrasts")
## attr(,"contrasts")$grp
## [1] "contr.treatment"
## 
## attr(,"contrasts")$ind.n
## [1] "contr.treatment"
## 
## attr(,"contrasts")$cnd
## [1] "contr.treatment"
```

Note that, if you have unbalanced numbers of individuals in the two groups, you will have zeros for some of the interactions between `grp` and `ind.n`. You can remove these columns manually from the model matrix and pass the corrected model matrix to the `full` argument of the *DESeq* function. See example code in the next section. Note that, in this case, you will not be able to create the *DESeqDataSet* with the design that leads to less than full rank model matrix. You can either use `design=~1` when creating the dataset object, or you can provide the corrected model matrix to the `design` slot of the dataset from the start.

Above, the terms `grpX.cndB` and `grpY.cndB` give the group-specific condition effects, in other words, the condition B vs A effect for group X samples, and likewise for group Y samples. These terms control for all of the six individual effects. These group-specific condition effects can be extracted using *results* with the `name` argument.

Furthermore, `grpX.cndB` and `grpY.cndB` can be contrasted using the `contrast` argument, in order to test if the condition effect is different across group:

```
results(dds, contrast=list("grpY.cndB","grpX.cndB"))
```

### Levels without samples

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