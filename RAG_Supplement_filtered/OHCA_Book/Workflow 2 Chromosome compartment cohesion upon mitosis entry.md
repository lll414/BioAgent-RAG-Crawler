Source URL: https://bioconductor.org/books/release/OHCA/pages/workflow-chicken.html
Date Scraped: 2025-12-15

---

1. [Advanced Hi-C topics](https://bioconductor.org/books/release/OHCA/pages/disseminating.html)
2. [Workflow 2: Chromosome compartment cohesion upon mitosis entry](https://bioconductor.org/books/release/OHCA/pages/workflow-chicken.html)

# Workflow 2: Chromosome compartment cohesion upon mitosis entry

Pre-loading packages and objects 📦

```
library(ggplot2)
library(cowplot)
library(purrr)
library(HiCExperiment)
##  Consider using the `HiContacts` package to perform advanced genomic operations 
##  on `HiCExperiment` objects.
##  
##  Read "Orchestrating Hi-C analysis with Bioconductor" online book to learn more:
##  https://js2264.github.io/OHCA/
##  
##  Attaching package: 'HiCExperiment'
##  The following object is masked from 'package:ggplot2':
##  
##      resolution
library(HiContactsData)
##  Loading required package: ExperimentHub
##  Loading required package: BiocGenerics
##  Loading required package: generics
##  
##  Attaching package: 'generics'
##  The following objects are masked from 'package:base':
##  
##      as.difftime, as.factor, as.ordered, intersect, is.element,
##      setdiff, setequal, union
##  
##  Attaching package: 'BiocGenerics'
##  The following object is masked from 'package:HiCExperiment':
##  
##      as.data.frame
##  The following objects are masked from 'package:stats':
##  
##      IQR, mad, sd, var, xtabs
##  The following objects are masked from 'package:base':
##  
##      Filter, Find, Map, Position, Reduce, anyDuplicated, aperm,
##      append, as.data.frame, basename, cbind, colnames, dirname,
##      do.call, duplicated, eval, evalq, get, grep, grepl, is.unsorted,
##      lapply, mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##      pmin.int, rank, rbind, rownames, sapply, saveRDS, table, tapply,
##      unique, unsplit, which.max, which.min
##  Loading required package: AnnotationHub
##  Loading required package: BiocFileCache
##  Loading required package: dbplyr
library(fourDNData)
```

Aims

This chapter illustrates how to:

* Annotate compartments for a list of HiC experiments
* Generate saddle plots for a list of HiC experiments
* Quantify changes in interactions between compartments between different timepoints

Datasets

We leverage five chicken datasets in this notebook, published in Gibcus et al. ([2018](https://bioconductor.org/books/release/OHCA/pages/workflow-chicken.html#ref-Gibcus_2018)). They are all available from the 4DN data portal using the `fourDNData` package.

* `4DNES9LEZXN7`: chicken cell culture blocked in G2
* `4DNESNWWIFZU`: chicken cell culture released from G2 block (5min)
* `4DNESGDXKM2I`: chicken cell culture released from G2 block (10min)
* `4DNESIR416OW`: chicken cell culture released from G2 block (15min)
* `4DNESS8PTK6F`: chicken cell culture released from G2 block (30min)

## Importing data

The 4DN consortium provides access to the datasets published in Gibcus et al. ([2018](https://bioconductor.org/books/release/OHCA/pages/workflow-chicken.html#ref-Gibcus_2018)). in `R`, they can be obtained thanks to the `fourDNData` gateway package.

Beware

The first time the following chunk of code is executed, it will cache a large amount of data (mostly consisting of contact matrices stored in `.mcool` files).

```
library(HiCExperiment)
library(fourDNData)
library(BiocParallel)
samples <- list(
    '4DNES9LEZXN7' = 'G2 block', 
    '4DNESNWWIFZU' = 'prophase (5m)', 
    '4DNESGDXKM2I' = 'prophase (10m)', 
    '4DNESIR416OW' = 'prometaphase (15m)', 
    '4DNESS8PTK6F' = 'prometaphase (30m)' 
)
bpparam <- MulticoreParam(workers = 5, progressbar = TRUE)
##  Warning:   'IS_BIOC_BUILD_MACHINE' environment variable detected, setting
##    BiocParallel workers to 4 (was 5)
hics <- bplapply(names(samples), fourDNHiCExperiment, BPPARAM = bpparam)
##  
  |                                                                         
  |                                                                   |   0%
  |                                                                         
  |=============                                                      |  20%
  |                                                                         
  |===========================                                        |  40%
  |                                                                         
  |========================================                           |  60%
  |                                                                         
  |======================================================             |  80%
  |                                                                         
  |===================================================================| 100%
```

```
names(hics) <- samples

hics[["G2 block"]]
##  `HiCExperiment` object with 150,494,008 contacts over 4,109 regions 
##  -------
##  fileName: "/home/biocbuild/.cache/R/fourDNData/bf1d554dfc19b_4DNFIT479GDR.mcool" 
##  focus: "whole genome" 
##  resolutions(13): 1000 2000 ... 5000000 10000000
##  active resolution: 250000 
##  interactions: 7262748 
##  scores(2): count balanced 
##  topologicalFeatures: compartments(891) borders(3465) 
##  pairsFile: N/A 
##  metadata(3): 4DN_info eigens diamond_insulation
```

## Plotting whole chromosome matrices

We can visualize the five different Hi-C maps on the entire chromosome `3` with `HiContacts` by iterating over each of the `HiCExperiment` objects.

```
library(purrr)
library(HiContacts)
##  Registered S3 methods overwritten by 'readr':
##    method                    from 
##    as.data.frame.spec_tbl_df vroom
##    as_tibble.spec_tbl_df     vroom
##    format.col_spec           vroom
##    print.col_spec            vroom
##    print.collector           vroom
##    print.date_names          vroom
##    print.locale              vroom
##    str.col_spec              vroom
library(ggplot2)
pl <- imap(hics, ~ .x['chr3'] |> 
    zoom(100000) |> 
    plotMatrix(use.scores = 'balanced', limits = c(-4, -1), caption = FALSE) + 
    ggtitle(.y)
)
library(cowplot)
plot_grid(plotlist = pl, nrow = 1)
```

![](workflow-chicken_files/figure-html/unnamed-chunk-4-1.png)

This highlights the progressive remodeling of chromatin into condensed chromosomes, starting as soon as 5’ after release from G2 phase.

## Zooming on a chromosome section

Zooming on a chromosome section, we can plot the Hi-C autocorrelation matrix for each timepoint. These matrices are generally used to highlight the overall correlation of interaction profiles between different segments of a chromosome section (see [Chapter 5](https://bioconductor.org/books/release/OHCA/pages/matrix-centric.html#computing-autocorrelated-map) for more details).

```
## --- Format compartment positions of chr. 4 segment
.chr <- 'chr4'
.start <- 59000000L
.stop <- 75000000L
library(GenomicRanges)
##  Loading required package: stats4
##  Loading required package: S4Vectors
##  
##  Attaching package: 'S4Vectors'
##  The following object is masked from 'package:HiCExperiment':
##  
##      metadata<-
##  The following object is masked from 'package:utils':
##  
##      findMatches
##  The following objects are masked from 'package:base':
##  
##      I, expand.grid, unname
##  Loading required package: IRanges
##  
##  Attaching package: 'IRanges'
##  The following object is masked from 'package:purrr':
##  
##      reduce
##  Loading required package: Seqinfo
coords <- GRanges(paste0(.chr, ':', .start, '-', .stop))
compts_df <- topologicalFeatures(hics[["G2 block"]], "compartments") |> 
    subsetByOverlaps(coords, type = 'within') |> 
    as.data.frame()
compts_gg <- geom_rect(
    data = compts_df, 
    mapping = aes(xmin = start, xmax = end, ymin = -500000, ymax = 0, alpha = compartment), 
    col = 'black', inherit.aes = FALSE
)

## --- Subset contact matrices to chr. 4 segment and computing autocorrelation scores
g2 <- hics[["G2 block"]] |> 
    zoom(100000) |> 
    subsetByOverlaps(coords) |>
    autocorrelate()
##  
pro5 <- hics[["prophase (5m)"]] |> 
    zoom(100000) |> 
    subsetByOverlaps(coords) |>
    autocorrelate()
pro30 <- hics[["prometaphase (30m)"]] |> 
    zoom(100000) |> 
    subsetByOverlaps(coords) |>
    autocorrelate()

## --- Plot autocorrelation matrices
plot_grid(
    plotMatrix(
        subsetByOverlaps(g2, coords),
        use.scores = 'autocorrelated', 
        scale = 'linear', 
        limits = c(-1, 1), 
        cmap = bwrColors(), 
        maxDistance = 10000000, 
        caption = FALSE
    ) + ggtitle('G2') + compts_gg,
    plotMatrix(
        subsetByOverlaps(pro5, coords),
        use.scores = 'autocorrelated', 
        scale = 'linear', 
        limits = c(-1, 1), 
        cmap = bwrColors(), 
        maxDistance = 10000000, 
        caption = FALSE
    ) + ggtitle('Prophase 5min') + compts_gg,
    plotMatrix(
        subsetByOverlaps(pro30, coords),
        use.scores = 'autocorrelated', 
        scale = 'linear', 
        limits = c(-1, 1), 
        cmap = bwrColors(), 
        maxDistance = 10000000, 
        caption = FALSE
    ) + ggtitle('Prometaphase 30min') + compts_gg,
    nrow = 1
)
##  Warning: Using alpha for a discrete variable is not advised.
##  Using alpha for a discrete variable is not advised.
##  Using alpha for a discrete variable is not advised.
```

![](workflow-chicken_files/figure-html/unnamed-chunk-5-1.png)

These correlation matrices suggest that there are two different regimes of chromatin compartment remodeling in this chromosome section:

1. Correlation scores between genomic bins within the compartment A remain positive 5’ after G2 release (albeit reduced compared to G2 block) and eventually become null 30’ after G2 release.
2. Correlation scores between genomic bins within the compartment B are overall null as soon as 5’ after G2 release.

## Generating saddle plots

Saddle plots are typically used to measure the `observed` vs. `expected` interaction scores within or between genomic loci belonging to A and B compartments. Here, they can be used to check whether the two regimes of chromatin compartment remodeling are observed genome-wide.

Non-overlapping genomic windows are grouped by `nbins` quantiles (typically between 10 and 50 bins) according to their A/B compartment eigenvector value, from lowest eigenvector values (i.e. strongest B compartments) to highest eigenvector values (i.e. strongest A compartments). The average `observed` vs. `expected` interaction scores are computed for pairwise eigenvector quantiles and plotted in a 2D heatmap.

```
pl <- imap(hics, ~ plotSaddle(.x, nbins = 38, BPPARAM = bpparam) + ggtitle(.y)) 
##  
  |                                                                         
  |                                                                   |   0%
  |                                                                         
  |==                                                                 |   3%
  |                                                                         
  |=====                                                              |   7%
  |                                                                         
  |=======                                                            |  10%
  |                                                                         
  |=========                                                          |  14%
  |                                                                         
  |============                                                       |  17%
  |                                                                         
  |==============                                                     |  21%
  |                                                                         
  |================                                                   |  24%
  |                                                                         
  |==================                                                 |  28%
  |                                                                         
  |=====================                                              |  31%
  |                                                                         
  |=======================                                            |  34%
  |                                                                         
  |=========================                                          |  38%
  |                                                                         
  |============================                                       |  41%
  |                                                                         
  |==============================                                     |  45%
  |                                                                         
  |================================                                   |  48%
  |                                                                         
  |===================================                                |  52%
  |                                                                         
  |=====================================                              |  55%
  |                                                                         
  |=======================================                            |  59%
  |                                                                         
  |==========================================                         |  62%
  |                                                                         
  |============================================                       |  66%
  |                                                                         
  |==============================================                     |  69%
  |                                                                         
  |=================================================                  |  72%
  |                                                                         
  |===================================================                |  76%
  |                                                                         
  |=====================================================              |  79%
  |                                                                         
  |=======================================================            |  83%
  |                                                                         
  |==========================================================         |  86%
  |                                                                         
  |============================================================       |  90%
  |                                                                         
  |==============================================================     |  93%
  |                                                                         
  |=================================================================  |  97%
  |                                                                         
  |===================================================================| 100%
##  
##  
  |                                                                         
  |                                                                   |   0%
  |                                                                         
  |==                                                                 |   3%
  |                                                                         
  |=====                                                              |   7%
  |                                                                         
  |=======                                                            |  10%
  |                                                                         
  |=========                                                          |  14%
  |                                                                         
  |============                                                       |  17%
  |                                                                         
  |==============                                                     |  21%
  |                                                                         
  |================                                                   |  24%
  |                                                                         
  |==================                                                 |  28%
  |                                                                         
  |=====================                                              |  31%
  |                                                                         
  |=======================                                            |  34%
  |                                                                         
  |=========================                                          |  38%
  |                                                                         
  |============================                                       |  41%
  |                                                                         
  |==============================                                     |  45%
  |                                                                         
  |================================                                   |  48%
  |                                                                         
  |===================================                                |  52%
  |                                                                         
  |=====================================                              |  55%
  |                                                                         
  |=======================================                            |  59%
  |                                                                         
  |==========================================                         |  62%
  |                                                                         
  |============================================                       |  66%
  |                                                                         
  |==============================================                     |  69%
  |                                                                         
  |=================================================                  |  72%
  |                                                                         
  |===================================================                |  76%
  |                                                                         
  |=====================================================              |  79%
  |                                                                         
  |=======================================================            |  83%
  |                                                                         
  |==========================================================         |  86%
  |                                                                         
  |============================================================       |  90%
  |                                                                         
  |==============================================================     |  93%
  |                                                                         
  |=================================================================  |  97%
  |                                                                         
  |===================================================================| 100%
##  
##  
  |                                                                         
  |                                                                   |   0%
  |                                                                         
  |==                                                                 |   3%
  |                                                                         
  |=====                                                              |   7%
  |                                                                         
  |=======                                                            |  10%
  |                                                                         
  |=========                                                          |  14%
  |                                                                         
  |============                                                       |  17%
  |                                                                         
  |==============                                                     |  21%
  |                                                                         
  |================                                                   |  24%
  |                                                                         
  |==================                                                 |  28%
  |                                                                         
  |=====================                                              |  31%
  |                                                                         
  |=======================                                            |  34%
  |                                                                         
  |=========================                                          |  38%
  |                                                                         
  |============================                                       |  41%
  |                                                                         
  |==============================                                     |  45%
  |                                                                         
  |================================                                   |  48%
  |                                                                         
  |===================================                                |  52%
  |                                                                         
  |=====================================                              |  55%
  |                                                                         
  |=======================================                            |  59%
  |                                                                         
  |==========================================                         |  62%
  |                                                                         
  |============================================                       |  66%
  |                                                                         
  |==============================================                     |  69%
  |                                                                         
  |=================================================                  |  72%
  |                                                                         
  |===================================================                |  76%
  |                                                                         
  |=====================================================              |  79%
  |                                                                         
  |=======================================================            |  83%
  |                                                                         
  |==========================================================         |  86%
  |                                                                         
  |============================================================       |  90%
  |                                                                         
  |==============================================================     |  93%
  |                                                                         
  |=================================================================  |  97%
  |                                                                         
  |===================================================================| 100%
##  
##  
  |                                                                         
  |                                                                   |   0%
  |                                                                         
  |==                                                                 |   3%
  |                                                                         
  |=====                                                              |   7%
  |                                                                         
  |=======                                                            |  10%
  |                                                                         
  |=========                                                          |  14%
  |                                                                         
  |============                                                       |  17%
  |                                                                         
  |==============                                                     |  21%
  |                                                                         
  |================                                                   |  24%
  |                                                                         
  |==================                                                 |  28%
  |                                                                         
  |=====================                                              |  31%
  |                                                                         
  |=======================                                            |  34%
  |                                                                         
  |=========================                                          |  38%
  |                                                                         
  |============================                                       |  41%
  |                                                                         
  |==============================                                     |  45%
  |                                                                         
  |================================                                   |  48%
  |                                                                         
  |===================================                                |  52%
  |                                                                         
  |=====================================                              |  55%
  |                                                                         
  |=======================================                            |  59%
  |                                                                         
  |==========================================                         |  62%
  |                                                                         
  |============================================                       |  66%
  |                                                                         
  |==============================================                     |  69%
  |                                                                         
  |=================================================                  |  72%
  |                                                                         
  |===================================================                |  76%
  |                                                                         
  |=====================================================              |  79%
  |                                                                         
  |=======================================================            |  83%
  |                                                                         
  |==========================================================         |  86%
  |                                                                         
  |============================================================       |  90%
  |                                                                         
  |==============================================================     |  93%
  |                                                                         
  |=================================================================  |  97%
  |                                                                         
  |===================================================================| 100%
##  
##  
  |                                                                         
  |                                                                   |   0%
  |                                                                         
  |==                                                                 |   3%
  |                                                                         
  |====                                                               |   7%
  |                                                                         
  |=======                                                            |  10%
  |                                                                         
  |=========                                                          |  13%
  |                                                                         
  |===========                                                        |  17%
  |                                                                         
  |=============                                                      |  20%
  |                                                                         
  |================                                                   |  23%
  |                                                                         
  |==================                                                 |  27%
  |                                                                         
  |====================                                               |  30%
  |                                                                         
  |======================                                             |  33%
  |                                                                         
  |=========================                                          |  37%
  |                                                                         
  |===========================                                        |  40%
  |                                                                         
  |=============================                                      |  43%
  |                                                                         
  |===============================                                    |  47%
  |                                                                         
  |==================================                                 |  50%
  |                                                                         
  |====================================                               |  53%
  |                                                                         
  |======================================                             |  57%
  |                                                                         
  |========================================                           |  60%
  |                                                                         
  |==========================================                         |  63%
  |                                                                         
  |=============================================                      |  67%
  |                                                                         
  |===============================================                    |  70%
  |                                                                         
  |=================================================                  |  73%
  |                                                                         
  |===================================================                |  77%
  |                                                                         
  |======================================================             |  80%
  |                                                                         
  |========================================================           |  83%
  |                                                                         
  |==========================================================         |  87%
  |                                                                         
  |============================================================       |  90%
  |                                                                         
  |===============================================================    |  93%
  |                                                                         
  |=================================================================  |  97%
  |                                                                         
  |===================================================================| 100%
plot_grid(plotlist = pl, nrow = 1)
```

![](workflow-chicken_files/figure-html/unnamed-chunk-6-1.png)

These plots confirm the previous observation made on chr. `4` and reveal that intra-B compartment interactions are generally lost 5’ after G2 release, while intra-A interactions take up to 15’ after G2 release to disappear.

Beware

The `plotSaddle()` function requires an eigenvector corresponding to A/B compartments. In this example, this eigenvector is recovered from the 4DN data portal. If not already available, this eigenvector can be computed from the contact matrix using the `getCompartments()` function.

## Quantifying interactions within and between compartments

We can leverage the replicate-merged contact matrices to quantify the interaction frequencies within A or B compartments or between A and B compartments, at different timepoints.

We can use the A/B compartment annotations obtained at the `G2 block` timepoint and extract `O/E` (observed vs expected) scores for interactions within A or B compartments or between A and B compartments, at different timepoints.

```
## --- Extract the A/B compartments identified in G2 block
compts <- topologicalFeatures(hics[["G2 block"]], "compartments")
compts$ID <- paste0(compts$compartment, seq_along(compts))

## --- Iterate over timepoints to extract `detrended` (O/E) scores and 
##     compartment annotations
library(tibble)
library(plyranges)
##  
##  Attaching package: 'plyranges'
##  The following object is masked from 'package:IRanges':
##  
##      slice
##  The following object is masked from 'package:stats':
##  
##      filter
df <- imap(hics[c(1, 2, 5)], ~ {
    ints <- cis(.x) |> ## Filter out trans interactions
        detrend() |> ## Compute O/E scores
        interactions() ## Recover interactions 
    ints$comp_first <- join_overlap_left(anchors(ints, "first"), compts)$ID
    ints$comp_second <- join_overlap_left(anchors(ints, "second"), compts)$ID
    tibble(
        sample = .y, 
        bin1 = ints$comp_first, 
        bin2 = ints$comp_second, 
        dist = InteractionSet::pairdist(ints), 
        OE = ints$detrended 
    ) |> 
        filter(dist > 5e6) |>
        mutate(type = dplyr::case_when(
            grepl('A', bin1) & grepl('A', bin2) ~ 'AA',
            grepl('B', bin1) & grepl('B', bin2) ~ 'BB',
            grepl('A', bin1) & grepl('B', bin2) ~ 'AB',
            grepl('B', bin1) & grepl('A', bin2) ~ 'BA'
        )) |> 
        filter(bin1 != bin2)
}) |> list_rbind() |> mutate(
    sample = factor(sample, names(hics)[c(1, 2, 5)])
)
```

We can now plot the changes in O/E scores for intra-A, intra-B, A-B or B-A interactions, splitting boxplots by timepoint.

```
ggplot(df, aes(x = type, y = OE, group = type, fill = type)) + 
    geom_boxplot(outlier.shape = NA) + 
    facet_grid(~sample) + 
    theme_bw() + 
    ylim(c(-2, 2))
##  Warning: Removed 66307 rows containing non-finite outside the scale range
##  (`stat_boxplot()`).
```

![](workflow-chicken_files/figure-html/unnamed-chunk-8-1.png)

This visualization suggests that interactions between genomic loci belonging to the B compartment are lost more rapidly than those between genomic loci belonging to the A compartment, when cells are released from G2 to enter mitosis.

Back to top