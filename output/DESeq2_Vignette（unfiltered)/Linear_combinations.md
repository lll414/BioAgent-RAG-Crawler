# Linear combinations

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

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