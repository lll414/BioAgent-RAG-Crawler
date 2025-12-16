# Do normalized counts correct for variables in the design?

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

No.Â The design variables are not used when estimating the size factors, and `counts(dds, normalized=TRUE)` is providing counts scaled by size or normalization factors. The design is only used when estimating dispersion and log2 fold changes.

The only case in which there is more than size factor scaling on the counts is when either normalization factors have been provided (e.g.Â from `cqn` or `EDASeq`), or if `tximport` is used and the upstream software corrected for various technical biases (e.g.Â *Salmon* quantification with GC bias correction). In this case, the average transcript length is taken into account when scaling the counts with `counts(dds, normalized=TRUE)`. For details, see the *tximport* package vignette and citation (Soneson, Love, and Robinson 2015).