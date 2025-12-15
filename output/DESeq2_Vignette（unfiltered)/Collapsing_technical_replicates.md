# Collapsing technical replicates

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

DESeq2 provides a function *collapseReplicates* which can assist in combining the counts from technical replicates into single columns of the count matrix. The term *technical replicate* implies multiple sequencing runs of the same library. You should not collapse biological replicates using this function. See the manual page for an example of the use of *collapseReplicates*.