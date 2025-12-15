# Regularized log transformation

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

The function *rlog*, stands for *regularized log*, transforming the original count data to the log2 scale by fitting a model with a term for each sample and a prior distribution on the coefficients which is estimated from the data. This is the same kind of shrinkage (sometimes referred to as regularization, or moderation) of log fold changes used by *DESeq* and *nbinomWaldTest*. The resulting data contains elements defined as:

\[ \log\_2(q\_{ij}) = \beta\_{i0} + \beta\_{ij} \]

where \(q\_{ij}\) is a parameter proportional to the expected true concentration of fragments for gene *i* and sample *j* (see formula [below](#theory)), \(\beta\_{i0}\) is an intercept which does not undergo shrinkage, and \(\beta\_{ij}\) is the sample-specific effect which is shrunk toward zero based on the dispersion-mean trend over the entire dataset. The trend typically captures high dispersions for low counts, and therefore these genes exhibit higher shrinkage from the *rlog*.

Note that, as \(q\_{ij}\) represents the part of the mean value \(\mu\_{ij}\) after the size factor \(s\_j\) has been divided out, it is clear that the rlog transformation inherently accounts for differences in sequencing depth. Without priors, this design matrix would lead to a non-unique solution, however the addition of a prior on non-intercept betas allows for a unique solution to be found.