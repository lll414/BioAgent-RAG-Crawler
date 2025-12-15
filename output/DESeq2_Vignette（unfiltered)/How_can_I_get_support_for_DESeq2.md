# How can I get support for DESeq2?

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

We welcome questions about our software, and want to ensure that we eliminate issues if and when they appear. We have a few requests to optimize the process:

* all questions should take place on the Bioconductor support site: <https://support.bioconductor.org>, which serves as a repository of questions and answers. This helps to save the developersâ time in responding to similar questions. Make sure to tag your post with `deseq2`. It is often very helpful in addition to describe the aim of your experiment.
* before posting, first search the Bioconductor support site mentioned above for past threads which might have answered your question.
* if you have a question about the behavior of a function, read the sections of the manual page for this function by typing a question mark and the function name, e.g.Â `?results`. We spend a lot of time documenting individual functions and the exact steps that the software is performing.
* include all of your R code, especially the creation of the *DESeqDataSet* and the design formula. Include complete warning or error messages, and conclude your message with the full output of `sessionInfo()`.
* if possible, include the output of `as.data.frame(colData(dds))`, so that we can have a sense of the experimental setup. If this contains confidential information, you can replace the levels of those factors using *levels()*.