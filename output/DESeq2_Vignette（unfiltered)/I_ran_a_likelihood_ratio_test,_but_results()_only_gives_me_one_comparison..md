# I ran a likelihood ratio test, but results() only gives me one comparison.

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

ââ¦ How do I get the *p* values for all of the variables/levels that were removed in the reduced design?â

This is explained in the help page for `?results` in the section about likelihood ratio test p-values, but we will restate the answer here. When one performs a likelihood ratio test, the *p* values and the test statistic (the `stat` column) are values for the test that removes all of the variables which are present in the full design and not in the reduced design. This tests the null hypothesis that all the coefficients from these variables and levels of these factors are equal to zero.

The likelihood ratio test *p* values therefore represent a test of *all the variables and all the levels of factors* which are among these variables. However, the results table only has space for one column of log fold change, so a single variable and a single comparison is shown (among the potentially multiple log fold changes which were tested in the likelihood ratio test). This is indicated at the top of the results table with the text, e.g., log2 fold change (MLE): condition C vs A, followed by, LRT p-value: â~ batch + conditionâ vs â~ batchâ. This indicates that the *p* value is for the likelihood ratio test of *all the variables and all the levels*, while the log fold change is a single comparison from among those variables and levels. See the help page for *results* for more details.