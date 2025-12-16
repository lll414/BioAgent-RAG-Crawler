# Contrasts

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Contrasts can be calculated for a *DESeqDataSet* object for which the GLM coefficients have already been fit using the Wald test steps (*DESeq* with `test="Wald"` or using *nbinomWaldTest*). The vector of coefficients \(\beta\) is left multiplied by the contrast vector \(c\) to form the numerator of the test statistic. The denominator is formed by multiplying the covariance matrix \(\Sigma\) for the coefficients on either side by the contrast vector \(c\). The square root of this product is an estimate of the standard error for the contrast. The contrast statistic is then compared to a Normal distribution as are the Wald statistics for the DESeq2 package.

\[ W = \frac{c^t \beta}{\sqrt{c^t \Sigma c}} \]