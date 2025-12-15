# Can I run DESeq2 to contrast the levels of many groups?

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

DESeq2 will work with any kind of design specified using the R formula. We enourage users to consider exploratory data analysis such as principal components analysis rather than performing statistical testing of all pairs of many groups of samples. Statistical testing is one of many ways of describing differences between samples.

As a speed concern with fitting very large models, note that each additional level of a factor in the design formula adds another parameter to the GLM which is fit by DESeq2. Users might consider first removing genes with very few reads, as this will speed up the fitting procedure.