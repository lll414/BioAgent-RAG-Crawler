# If I have multiple groups, should I run all together or split into pairs of groups?

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Typically, we recommend users to run samples from all groups together, and then use the `contrast` argument of the *results* function to extract comparisons of interest after fitting the model using *DESeq*.

The model fit by *DESeq* estimates a single dispersion parameter for each gene, which defines how far we expect the observed count for a sample will be from the mean value from the model given its size factor and its condition group. See the section [above](#theory) and the DESeq2 paper for full details. Having a single dispersion parameter for each gene is usually sufficient for analyzing multi-group data, as the final dispersion value will incorporate the within-group variability across all groups.

However, for some datasets, exploratory data analysis (EDA) plots could reveal that one or more groups has much higher within-group variability than the others. A simulated example of such a set of samples is shown below. This is case where, by comparing groups A and B separately â subsetting a *DESeqDataSet* to only samples from those two groups and then running *DESeq* on this subset â will be more sensitive than a model including all samples together. It should be noted that such an extreme range of within-group variability is not common, although it could arise if certain treatments produce an extreme reaction (e.g.Â cell death). Again, this can be easily detected from the EDA plots such as PCA described in this vignette.

Here we diagram an extreme range of within-group variability with a simulated dataset. Typically, it is recommended to run *DESeq* across samples from all groups, for datasets with multiple groups. However, this simulated dataset shows a case where it would be preferable to compare groups A and B by creating a smaller dataset without the C samples. Group C has much higher within-group variability, which would inflate the per-gene dispersion estimate for groups A and B as well:
