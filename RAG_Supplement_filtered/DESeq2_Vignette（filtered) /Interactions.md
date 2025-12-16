# Interactions

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

Interaction terms can be added to the design formula, in order to test, for example, if the log2 fold change attributable to a given condition is *different* based on another factor, for example if the condition effect differs across genotype.

**Initial note:** Many users begin to add interaction terms to the design formula, when in fact a much simpler approach would give all the results tables that are desired. We will explain this approach first, because it is much simpler to perform. If the comparisons of interest are, for example, the effect of a condition for different sets of samples, a simpler approach than adding interaction terms explicitly to the design formula is to perform the following steps:

* combine the factors of interest into a single factor with all combinations of the original factors
* change the design to include just this factor, e.g.Â ~ group

Using this design is similar to adding an interaction term, in that it models multiple condition effects which can be easily extracted with *results*. Suppose we have two factors `genotype` (with values I, II, and III) and `condition` (with values A and B), and we want to extract the condition effect specifically for each genotype. We could use the following approach to obtain, e.g.Â the condition effect for genotype I:

```
dds$group <- factor(paste0(dds$genotype, dds$condition))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)
results(dds, contrast=c("group", "IB", "IA"))
```

**Adding interactions to the design:** The following two plots diagram genotype-specific condition effects, which could be modeled with interaction terms by using a design of `~genotype + condition + genotype:condition`.

In the first plot (Gene 1), note that the condition effect is consistent across genotypes. Although condition A has a different baseline for I,II, and III, the condition effect is a log2 fold change of about 2 for each genotype. Using a model with an interaction term `genotype:condition`, the interaction terms for genotype II and genotype III will be nearly 0.

Here, the y-axis represents log2(n+1), and each group has 20 samples (black dots). A red line connects the mean of the groups within each genotype.


In the second plot (Gene 2), we can see that the condition effect is not consistent across genotype. Here the main condition effect (the effect for the reference genotype I) is again 2. However, this time the interaction terms will be around 1 for genotype II and -4 for genotype III. This is because the condition effect is higher by 1 for genotype II compared to genotype I, and lower by 4 for genotype III compared to genotype I. The condition effect for genotype II (or III) is obtained by adding the main condition effect and the interaction term for that genotype. Such a plot can be made using the *plotCounts* function as shown above.


Now we will continue to explain the use of interactions in order to test for *differences* in condition effects. We continue with the example of condition effects across three genotypes (I, II, and III).

The key point to remember about designs with interaction terms is that, unlike for a design `~genotype + condition`, where the condition effect represents the *overall* effect controlling for differences due to genotype, by adding `genotype:condition`, the main condition effect only represents the effect of condition for the *reference level* of genotype (I, or whichever level was defined by the user as the reference level). The interaction terms `genotypeII.conditionB` and `genotypeIII.conditionB` give the *difference* between the condition effect for a given genotype and the condition effect for the reference genotype.

This genotype-condition interaction example is examined in further detail in Example 3 in the help page for *results*, which can be found by typing `?results`. In particular, we show how to test for differences in the condition effect across genotype, and we show how to obtain the condition effect for non-reference genotypes.