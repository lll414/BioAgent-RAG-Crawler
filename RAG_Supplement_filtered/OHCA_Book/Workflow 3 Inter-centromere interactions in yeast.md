Source URL: https://bioconductor.org/books/release/OHCA/pages/workflow-centros.html
Date Scraped: 2025-12-15

---

1. [Advanced Hi-C topics](https://bioconductor.org/books/release/OHCA/pages/disseminating.html)
2. [Workflow 3: Inter-centromere interactions in yeast](https://bioconductor.org/books/release/OHCA/pages/workflow-centros.html)

# Workflow 3: Inter-centromere interactions in yeast

Pre-loading packages and objects 📦

```
library(ggplot2)
library(purrr)
library(GenomicRanges)
##  Loading required package: stats4
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
##  Loading required package: S4Vectors
##  
##  Attaching package: 'S4Vectors'
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
library(InteractionSet)
##  Loading required package: SummarizedExperiment
##  Loading required package: MatrixGenerics
##  Loading required package: matrixStats
##  
##  Attaching package: 'MatrixGenerics'
##  The following objects are masked from 'package:matrixStats':
##  
##      colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##      colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##      colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##      colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##      colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##      colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##      colWeightedMeans, colWeightedMedians, colWeightedSds,
##      colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##      rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##      rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##      rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##      rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##      rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs,
##      rowVars, rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##      rowWeightedSds, rowWeightedVars
##  Loading required package: Biobase
##  Welcome to Bioconductor
##  
##      Vignettes contain introductory material; view with
##      'browseVignettes()'. To cite Bioconductor, see
##      'citation("Biobase")', and for packages 'citation("pkgname")'.
##  
##  Attaching package: 'Biobase'
##  The following object is masked from 'package:MatrixGenerics':
##  
##      rowMedians
##  The following objects are masked from 'package:matrixStats':
##  
##      anyMissing, rowMedians
library(HiCExperiment)
##  Consider using the `HiContacts` package to perform advanced genomic operations 
##  on `HiCExperiment` objects.
##  
##  Read "Orchestrating Hi-C analysis with Bioconductor" online book to learn more:
##  https://js2264.github.io/OHCA/
##  
##  Attaching package: 'HiCExperiment'
##  The following object is masked from 'package:SummarizedExperiment':
##  
##      metadata<-
##  The following object is masked from 'package:S4Vectors':
##  
##      metadata<-
##  The following object is masked from 'package:ggplot2':
##  
##      resolution
library(HiContactsData)
##  Loading required package: ExperimentHub
##  Loading required package: AnnotationHub
##  Loading required package: BiocFileCache
##  Loading required package: dbplyr
##  
##  Attaching package: 'AnnotationHub'
##  The following object is masked from 'package:Biobase':
##  
##      cache
library(multiHiCcompare)
##  
##  Attaching package: 'multiHiCcompare'
##  The following object is masked from 'package:HiCExperiment':
##  
##      resolution
##  The following object is masked from 'package:ggplot2':
##  
##      resolution
```

Aims

This chapter illustrates how to plot the aggregate signal over pairs of genomic ranges, in this case pairs of yeast centromeres.

Datasets

We leverage two yeast datasets in this notebook.

* One from a WT yeast strain in G1 phase
* One from a WT yeast strain in G2/M phase

## Importing Hi-C data and plotting contact matrices

```
library(HiContactsData)
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
library(purrr)
library(ggplot2)
hics <- list(
    'G1' = import(HiContactsData('yeast_g1', 'mcool'), format = 'cool', resolution = 4000),
    'G2M' = import(HiContactsData('yeast_g2m', 'mcool'), format = 'cool', resolution = 4000)
)
##  see ?HiContactsData and browseVignettes('HiContactsData') for documentation
##  loading from cache
##  see ?HiContactsData and browseVignettes('HiContactsData') for documentation
##  loading from cache
imap(hics, ~ plotMatrix(
    .x, use.scores = 'balanced', limits = c(-4, -1), caption = FALSE
) + ggtitle(.y))
##  $G1
```

![](workflow-centros_files/figure-html/unnamed-chunk-2-1.png)

```
##  
##  $G2M
```

![](workflow-centros_files/figure-html/unnamed-chunk-2-2.png)

We can visually appreciate that inter-chromosomal interactions, notably between centromeres, are less prominent in G2/M.

## Checking P(s) and cis/trans interactions ratio

```
library(dplyr)
##  
##  Attaching package: 'dplyr'
##  The following objects are masked from 'package:dbplyr':
##  
##      ident, sql
##  The following object is masked from 'package:Biobase':
##  
##      combine
##  The following object is masked from 'package:matrixStats':
##  
##      count
##  The following objects are masked from 'package:GenomicRanges':
##  
##      intersect, setdiff, union
##  The following object is masked from 'package:Seqinfo':
##  
##      intersect
##  The following objects are masked from 'package:IRanges':
##  
##      collapse, desc, intersect, setdiff, slice, union
##  The following objects are masked from 'package:S4Vectors':
##  
##      first, intersect, rename, setdiff, setequal, union
##  The following objects are masked from 'package:BiocGenerics':
##  
##      combine, intersect, setdiff, setequal, union
##  The following object is masked from 'package:generics':
##  
##      explain
##  The following objects are masked from 'package:stats':
##  
##      filter, lag
##  The following objects are masked from 'package:base':
##  
##      intersect, setdiff, setequal, union
pairs <- list(
    'G1' = PairsFile(HiContactsData('yeast_g1', 'pairs')),
    'G2M' = PairsFile(HiContactsData('yeast_g2m', 'pairs')) 
)
##  see ?HiContactsData and browseVignettes('HiContactsData') for documentation
##  loading from cache
##  see ?HiContactsData and browseVignettes('HiContactsData') for documentation
##  loading from cache
ps <- imap_dfr(pairs, ~ distanceLaw(.x, by_chr = TRUE) |> 
    mutate(sample = .y) 
)
##  Importing pairs file /home/biocbuild/.cache/R/ExperimentHub/c040417bddb78_8630 in memory. This may take a while...
##  Importing pairs file /home/biocbuild/.cache/R/ExperimentHub/c040430b74c78_8631 in memory. This may take a while...
plotPs(ps, aes(x = binned_distance, y = norm_p, group = interaction(sample, chr), color = sample)) + 
    scale_color_manual(values = c('black', 'red'))
##  Warning: Removed 2133 rows containing missing values or values outside the scale
##  range (`geom_line()`).
```

![](workflow-centros_files/figure-html/unnamed-chunk-3-1.png)

```
plotPsSlope(ps, ggplot2::aes(x = binned_distance, y = slope, group = interaction(sample, chr), color = sample)) + 
    scale_color_manual(values = c('black', 'red'))
##  Warning: Removed 2183 rows containing missing values or values outside the scale
##  range (`geom_line()`).
```

![](workflow-centros_files/figure-html/unnamed-chunk-3-2.png)

This confirms that interactions in cells synchronized in G2/M are enriched for 10-30kb-long interactions.

```
ratios <- imap_dfr(hics, ~ cisTransRatio(.x) |> mutate(sample = .y))
ggplot(ratios, aes(x = chr, y = trans_pct, fill = sample)) + 
    geom_col() + 
    labs(x = 'Chromosomes', y = "% of trans interactions") + 
    scale_y_continuous(labels = scales::percent) + 
    facet_grid(~sample)
```

![](workflow-centros_files/figure-html/unnamed-chunk-4-1.png)

We can also highlight that trans (inter-chromosomal) interactions are proportionally decreasing in G2/M-synchronized cells.

## Centromere virtual 4C profiles

```
data(centros_yeast)
v4c_centro <- imap_dfr(hics, ~ virtual4C(.x, GenomicRanges::resize(centros_yeast[2], 8000)) |> 
    as_tibble() |> 
    mutate(sample = .y) |> 
    filter(seqnames == 'IV')
) 
ggplot(v4c_centro, aes(x = start, y = score, fill = sample)) +
    geom_area() +
    theme_bw() +
    labs(
        x = "chrIV position", 
        y = "Contacts with chrII centromere", 
        title = "Interaction profile of chrII centromere"
    ) + 
    coord_cartesian(ylim = c(0, 0.015))
```

![](workflow-centros_files/figure-html/unnamed-chunk-5-1.png)

## Aggregated 2D signal over all pairs of centromeres

We can start by computing all possible pairs of centromeres.

```
centros_pairs <- lapply(1:length(centros_yeast), function(i) {
    lapply(1:length(centros_yeast), function(j) {
        S4Vectors::Pairs(centros_yeast[i], centros_yeast[j])
    })
}) |> 
    do.call(c, args = _) |>
    do.call(c, args = _) |> 
    InteractionSet::makeGInteractionsFromGRangesPairs()
centros_pairs <- centros_pairs[anchors(centros_pairs, 'first') != anchors(centros_pairs, 'second')]

centros_pairs
##  GInteractions object with 240 interactions and 0 metadata columns:
##          seqnames1       ranges1     seqnames2       ranges2
##              <Rle>     <IRanges>         <Rle>     <IRanges>
##      [1]         I 151583-151641 ---        II 238361-238419
##      [2]         I 151583-151641 ---       III 114322-114380
##      [3]         I 151583-151641 ---        IV 449879-449937
##      [4]         I 151583-151641 ---         V 152522-152580
##      [5]         I 151583-151641 ---        VI 147981-148039
##      ...       ...           ... ...       ...           ...
##    [236]       XVI 556255-556313 ---        XI 440229-440287
##    [237]       XVI 556255-556313 ---       XII 151366-151424
##    [238]       XVI 556255-556313 ---      XIII 268222-268280
##    [239]       XVI 556255-556313 ---       XIV 628588-628646
##    [240]       XVI 556255-556313 ---        XV 326897-326955
##    -------
##    regions: 16 ranges and 0 metadata columns
##    seqinfo: 17 sequences (1 circular) from R64-1-1 genome
```

Then we can aggregate the Hi-C signal over each pair of centromeres.

```
aggr_maps <- purrr::imap(hics, ~ {
    aggr <- aggregate(.x, centros_pairs, maxDistance = 1e999)
    plotMatrix(
        aggr, use.scores = 'balanced', limits = c(-5, -1), 
        cmap = HiContacts::rainbowColors(), 
        caption = FALSE
    ) + ggtitle(.y)
})
##  Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 120 out-of-bound ranges located on sequences I,
##    III, V, VI, VIII, IX, XII, and XIV. Note that ranges located on a sequence
##    whose length is unknown (NA) or on a circular sequence are not considered
##    out-of-bound (use seqlengths() and isCircular() to get the lengths and
##    circularity flags of the underlying sequences). You can use trim() to trim
##    these ranges. See ?`trim,GenomicRanges-method` for more information.
##  Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 120 out-of-bound ranges located on sequences III,
##    V, VI, VIII, IX, XII, XIV, and I. Note that ranges located on a sequence
##    whose length is unknown (NA) or on a circular sequence are not considered
##    out-of-bound (use seqlengths() and isCircular() to get the lengths and
##    circularity flags of the underlying sequences). You can use trim() to trim
##    these ranges. See ?`trim,GenomicRanges-method` for more information.
##  Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 240 out-of-bound ranges located on sequences I,
##    III, V, VI, VIII, IX, XII, and XIV. Note that ranges located on a sequence
##    whose length is unknown (NA) or on a circular sequence are not considered
##    out-of-bound (use seqlengths() and isCircular() to get the lengths and
##    circularity flags of the underlying sequences). You can use trim() to trim
##    these ranges. See ?`trim,GenomicRanges-method` for more information.
##  Going through preflight checklist...
##  Parsing the entire contact matrice as a sparse matrix...
##  Modeling distance decay...
##  Filtering for contacts within provided targets...
##  Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 120 out-of-bound ranges located on sequences I,
##    III, V, VI, VIII, IX, XII, and XIV. Note that ranges located on a sequence
##    whose length is unknown (NA) or on a circular sequence are not considered
##    out-of-bound (use seqlengths() and isCircular() to get the lengths and
##    circularity flags of the underlying sequences). You can use trim() to trim
##    these ranges. See ?`trim,GenomicRanges-method` for more information.
##  Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 120 out-of-bound ranges located on sequences III,
##    V, VI, VIII, IX, XII, XIV, and I. Note that ranges located on a sequence
##    whose length is unknown (NA) or on a circular sequence are not considered
##    out-of-bound (use seqlengths() and isCircular() to get the lengths and
##    circularity flags of the underlying sequences). You can use trim() to trim
##    these ranges. See ?`trim,GenomicRanges-method` for more information.
##  Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 240 out-of-bound ranges located on sequences I,
##    III, V, VI, VIII, IX, XII, and XIV. Note that ranges located on a sequence
##    whose length is unknown (NA) or on a circular sequence are not considered
##    out-of-bound (use seqlengths() and isCircular() to get the lengths and
##    circularity flags of the underlying sequences). You can use trim() to trim
##    these ranges. See ?`trim,GenomicRanges-method` for more information.
##  Going through preflight checklist...
##  Parsing the entire contact matrice as a sparse matrix...
##  Modeling distance decay...
##  Filtering for contacts within provided targets...

cowplot::plot_grid(plotlist = aggr_maps, nrow = 1)
```

![](workflow-centros_files/figure-html/unnamed-chunk-7-1.png)

## Aggregated 1D interaction profile of centromeres

One can generalize the previous virtual 4C plot, by extracting the interaction profile between all possible pairs of centromeres in each dataset.

```
df <- map_dfr(1:{length(centros_yeast)-1}, function(i) {
    centro1 <- GenomicRanges::resize(centros_yeast[i], fix = 'center', 8000)
    map_dfr({i+1}:length(centros_yeast), function(j) {
        centro2 <- GenomicRanges::resize(centros_yeast[j], fix = 'center', 80000)
        gi <- InteractionSet::GInteractions(centro1, centro2)
        imap_dfr(hics, ~ .x[gi] |> 
            interactions() |> 
            as_tibble() |>
            mutate(
                sample = .y, 
                center = center2 - start(GenomicRanges::resize(centro2, fix = 'center', 1))
            ) |> 
            select(sample, seqnames1, seqnames2, center, balanced)
        )
    })
}) 
ggplot(df, aes(x = center/1e3, y = balanced)) + 
    geom_line(aes(group = interaction(seqnames1, seqnames2)), alpha = 0.03, col = "black") + 
    geom_smooth(col = "red", fill = "red") + 
    theme_bw() + 
    theme(legend.position = 'none') + 
    labs(
        x = "Distance from centromere (kb)", y = "Normalized interaction frequency", 
        title = "Centromere pairwise interaction profiles"
    ) +
    facet_grid(~sample)
##  `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'
##  Warning: Removed 25 rows containing non-finite outside the scale range
##  (`stat_smooth()`).
##  Warning: Removed 3 rows containing missing values or values outside the scale range
##  (`geom_line()`).
```

![](workflow-centros_files/figure-html/unnamed-chunk-8-1.png)

Back to top