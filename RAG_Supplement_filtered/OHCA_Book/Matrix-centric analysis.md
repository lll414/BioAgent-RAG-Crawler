Source URL: https://bioconductor.org/books/release/OHCA/pages/matrix-centric.html
Date Scraped: 2025-12-15

---

1. [In-depth Hi-C analysis](https://bioconductor.org/books/release/OHCA/pages/matrix-centric.html)
2. [5  Matrix-centric analysis](https://bioconductor.org/books/release/OHCA/pages/matrix-centric.html)

# 5  Matrix-centric analysis

Pre-loading packages and objects 📦

```
library(ggplot2)
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
library(rtracklayer)
##  
##  Attaching package: 'rtracklayer'
##  The following object is masked from 'package:AnnotationHub':
##  
##      hubUrl
```

Aims

This chapter focuses on the various analytical tools offered by `HiContacts` to compute matrix-related metrics from a `HiCExperiment` object.

In the first part of this book, we have seen how to query parts or all of the data contained in Hi-C contact matrices using the `HiCExperiment` object ([Chapter 2](https://bioconductor.org/books/release/OHCA/pages/data-representation.html#hicexperiment-class)), how to manipulate `HiCExperiment` objects ([Chapter 3](https://bioconductor.org/books/release/OHCA/pages/parsing.html)) and how to visualize Hi-C contact matrices as heatmaps ([Chapter 4](https://bioconductor.org/books/release/OHCA/pages/visualization.html)).

The `HiContacts` package directly operates on `HiCExperiment` objects and extends its usability by providing a comprehensive toolkit to analyze Hi-C data, focusing on four main topics:

* Contact matrix-centric analyses (this chapter)
* Interactions-centric analyses ([Chapter 6](https://bioconductor.org/books/release/OHCA/pages/interactions-centric.html))
* Structural feature annotations ([Chapter 7](https://bioconductor.org/books/release/OHCA/pages/topological-features.html))
* Hi-C visualization (see [previous chapter](https://bioconductor.org/books/release/OHCA/pages/visualization.html))

**Matrix-centric** analyses consider a `HiCExperiment` object from the “matrix” perspective to perform a range of matrix-based operations. This encompasses:

* Computing observed/expected (O/E) map
* Computing auto-correlation map
* Smoothing out a contact map
* Merging multiple Hi-C maps together
* Comparing two Hi-C maps to each other

![](images/20230421134800.jpg)

Note

* All the functions described in this chapter are **endomorphisms**: they take `HiCExperiment` objects as input and return modified `HiCExperiment` objects.
* Internally, most of the functions presented in this chapter make a call to `as.matrix(<HiCExperiment>)` to coerce it into a `matrix`.

Generating the example `hic` object 👇

To demonstrate `HiContacts` functionalities, we will create an `HiCExperiment` object from an example `.cool` file provided in the `HiContactsData` package.

```
library(HiCExperiment)
library(HiContactsData)

# ---- This downloads an example `.mcool` file and caches it locally 
coolf <- HiContactsData('yeast_wt', 'mcool')
##  see ?HiContactsData and browseVignettes('HiContactsData') for documentation
##  loading from cache

# ---- This creates a connection to the disk-stored `.mcool` file
cf <- CoolFile(coolf)
cf
##  CoolFile object
##  .mcool file: /home/biocbuild/.cache/R/ExperimentHub/15c3e2641b8180_7752 
##  resolution: 1000 
##  pairs file: 
##  metadata(0):

# ---- This imports contacts from the chromosome `II` at resolution `2000`
hic <- import(cf, focus = 'II', resolution = 2000)
```

```
hic
##  `HiCExperiment` object with 471,364 contacts over 407 regions 
##  -------
##  fileName: "/home/biocbuild/.cache/R/ExperimentHub/15c3e2641b8180_7752" 
##  focus: "II" 
##  resolutions(5): 1000 2000 4000 8000 16000
##  active resolution: 2000 
##  interactions: 34063 
##  scores(2): count balanced 
##  topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
##  pairsFile: N/A 
##  metadata(0):
```

## 5.1 Operations in an individual matrix

### 5.1.1 Balancing a raw interaction count map

Hi-C sequencing coverage is systematically affected by multiple confounding factors, e.g.  density of restriction sites, GC%, genome mappability, etc.. Overall, it generally ends up not homogenous throughout the entire genome and this leads to artifacts in un-normalized `count` matrices.

To correct for sequencing coverage heterogeneity of raw `count` maps, Hi-C data can be normalized using matrix balancing approaches (Cournac et al. ([2012](https://bioconductor.org/books/release/OHCA/pages/matrix-centric.html#ref-Cournac_2012)), Imakaev et al. ([2012](https://bioconductor.org/books/release/OHCA/pages/matrix-centric.html#ref-Imakaev_2012))). This is generally done directly on the disk-stored matrices using out-of-memory strategies (e.g. with `cooler balance <.cool>`). However, if contact matrix files are imported into a `HiCExperiment` object but no `balanced` scores are available, in-memory balancing can be performed using the `normalize` function. This adds an extra `ICE` element in `scores` list (while the `interactions` themselves are unmodified).

```
normalized_hic <- normalize(hic)
normalized_hic
##  `HiCExperiment` object with 471,364 contacts over 407 regions 
##  -------
##  fileName: "/home/biocbuild/.cache/R/ExperimentHub/15c3e2641b8180_7752" 
##  focus: "II" 
##  resolutions(5): 1000 2000 4000 8000 16000
##  active resolution: 2000 
##  interactions: 34063 
##  scores(3): count balanced ICE 
##  topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
##  pairsFile: N/A 
##  metadata(0):
```

It is possible to plot the different `scores` of the resulting object to visualize the newly computed `scores`. In this example, `ICE` scores should be nearly identical to `balanced` scores, which were originally imported from the disk-stored contact matrix.

```
cowplot::plot_grid(
    plotMatrix(normalized_hic, use.scores = 'count', caption = FALSE),
    plotMatrix(normalized_hic, use.scores = 'balanced', caption = FALSE),
    plotMatrix(normalized_hic, use.scores = 'ICE', caption = FALSE), 
    nrow = 1
)
```

![](matrix-centric_files/figure-html/unnamed-chunk-5-1.png)

### 5.1.2 Computing observed/expected (O/E) map

The most prominent feature of a balanced Hi-C matrix is the strong main diagonal. This main diagonal is observed because interactions between immediate adjacent genomic loci are more prone to happen than interactions spanning longer genomic distances. This “expected” behavior is due to the polymer nature of the chromosomes being studied, and can be locally estimated using the distance-dependent interaction frequency (a.k.a. the “distance law”, or P(s)). It can be used to compute an `expected` matrix on interactions.

When it is desirable to “mask” this polymer behavior to emphasize topological structures formed by chromosomes, one can divide a given balanced matrix by its `expected` matrix, i.e. calculate the observed/expected (O/E) map. This is sometimes called “detrending”, as it effectively removes the average polymer behavior from the balanced matrix.

The `detrend` function performs this operation on a given `HiCExperiment` object. It adds two extra elements in `scores` list: `expected` and `detrended` metrics (while the `interactions` themselves are unmodified).

```
detrended_hic <- detrend(hic)
detrended_hic
##  `HiCExperiment` object with 471,364 contacts over 407 regions 
##  -------
##  fileName: "/home/biocbuild/.cache/R/ExperimentHub/15c3e2641b8180_7752" 
##  focus: "II" 
##  resolutions(5): 1000 2000 4000 8000 16000
##  active resolution: 2000 
##  interactions: 34063 
##  scores(4): count balanced expected detrended 
##  topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
##  pairsFile: N/A 
##  metadata(0):
```

Topological features will be visually more prominent in the O/E `detrended` Hi-C map.

```
cowplot::plot_grid(
    plotMatrix(detrended_hic, use.scores = 'balanced', scale = 'log10', limits = c(-3.5, -1.2), caption = FALSE),
    plotMatrix(detrended_hic, use.scores = 'expected', scale = 'log10', limits = c(-3.5, -1.2), caption = FALSE),
    plotMatrix(detrended_hic, use.scores = 'detrended', scale = 'linear', limits = c(-1, 1), cmap = bwrColors(), caption = FALSE), 
    nrow = 1
)
```

![](matrix-centric_files/figure-html/unnamed-chunk-7-1.png)

Scale for `detrended` scores

* `expected` scores are in `linear` scale and ± in the same amplitude than `balanced` scores;
* `detrended` scores are in `log2` scale, in general approximately centered around 0. When plotting `detrended` scores, `scale = linear` should be set to prevent the default `log10` scaling.

### 5.1.3 Computing autocorrelated map

Correlation matrices are often calculated from balanced Hi-C matrices. For instance, in genomes composed of eu- and heterochromatin, a correlation matrix can be used to reveal a checkerboard pattern emphasizing the segregation of chromatin into two A/B compartments (Lieberman-Aiden et al. ([2009](https://bioconductor.org/books/release/OHCA/pages/matrix-centric.html#ref-Lieberman_Aiden_2009))).

The `autocorrelate` function is used to compute a correlation matrix of a `HiCExperiment` object. For each pair of interacting loci, the `autocorrelated` score represents the correlation between their respective interaction profiles with the rest of the genome.

```
autocorr_hic <- autocorrelate(hic)
##  
autocorr_hic
##  `HiCExperiment` object with 471,364 contacts over 407 regions 
##  -------
##  fileName: "/home/biocbuild/.cache/R/ExperimentHub/15c3e2641b8180_7752" 
##  focus: "II" 
##  resolutions(5): 1000 2000 4000 8000 16000
##  active resolution: 2000 
##  interactions: 34063 
##  scores(5): count balanced expected detrended autocorrelated 
##  topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
##  pairsFile: N/A 
##  metadata(0):
```

Since these metrics represent correlation scores, they range between `-1` and `1`. Two loci with an `autocorrelated` score close to `-1` have anti-correlated interaction profiles, while two loci with a `autocorrelated` score close to `1` are likely to interact with shared targets.

```
summary(scores(autocorr_hic, 'autocorrelated'))
##      Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
##  -0.41561  0.00249  0.05040  0.06447  0.10360  1.00000      564
```

Correlated and anti-correlated loci will be visually represented in the `autocorrelated` Hi-C map in red and blue pixels, respectively.

Note

Here we have illustrated how to compute an autocorrelation matrix from a `HiCExperiment` object using the example **yeast** Hi-C experiment. Bear in mind that this is unusual and not very useful, as yeast chromatin is not segregated in two compartments but rather follows a Rabl conformation (Duan et al. ([2010](https://bioconductor.org/books/release/OHCA/pages/matrix-centric.html#ref-Duan_2010))). An example of autocorrelation map from a vertebrate Hi-C experiment (for which chromatin is segregated in A/B compartments) is shown in [Chapter 10](https://bioconductor.org/books/release/OHCA/pages/workflow-chicken.html).

```
plotMatrix(
    autocorr_hic, 
    use.scores = 'autocorrelated', 
    scale = 'linear', 
    limits = c(-0.4, 0.4), 
    cmap = bgrColors()
)
```

![](matrix-centric_files/figure-html/unnamed-chunk-10-1.png)

Scale for `autocorrelated` scores

* `autocorrelated` scores are in `linear` scale, in general approximately centered around 0. When plotting `autocorrelated` scores, `scale = linear` should be set to prevent the default `log10` scaling.
* `limits` should be manually set to `c(-x, x)` (`0 < x <= 1`) to ensure that the color range is effectively centered on `0`.

### 5.1.4 Despeckling (smoothing out) a contact map

Shallow-sequenced Hi-C libraries or matrices binned with an overly small bin size sometimes produce “grainy” Hi-C maps with noisy backgrounds. A grainy map may also be obtained when dividing two matrices, e.g. when computing the O/E ratio with `detrend`. This is particularly true for sparser long-range interactions. To overcome such limitations, `HiCExperiment` objects can be “`despeckle`d” to smooth out focal speckles.

```
hic2 <- detrend(hic['II:400000-700000'])
hic2 <- despeckle(hic2, use.scores = 'detrended', focal.size = 2)
hic2
##  `HiCExperiment` object with 168,785 contacts over 150 regions 
##  -------
##  fileName: "/home/biocbuild/.cache/R/ExperimentHub/15c3e2641b8180_7752" 
##  focus: "II:400,000-700,000" 
##  resolutions(5): 1000 2000 4000 8000 16000
##  active resolution: 2000 
##  interactions: 11325 
##  scores(5): count balanced expected detrended detrended.despeckled 
##  topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
##  pairsFile: N/A 
##  metadata(0):
```

The added `<use.scores>.despeckled` scores correspond to scores averaged using a window, whose width is provided with the `focal.size` argument. This results in a smoother Hi-C heatmap, effectively removing the “speckles” observed at longer range.

```
library(InteractionSet)
loops <- system.file('extdata', 'S288C-loops.bedpe', package = 'HiCExperiment') |> 
    import() |> 
    makeGInteractionsFromGRangesPairs()
borders <- system.file('extdata', 'S288C-borders.bed', package = 'HiCExperiment') |> 
    import()
cowplot::plot_grid(
    plotMatrix(hic2, caption = FALSE),
    plotMatrix(hic2, use.scores = 'detrended', scale = 'linear', limits = c(-1, 1), caption = FALSE),
    plotMatrix(
        hic2, 
        use.scores = 'detrended.despeckled', 
        scale = 'linear', 
        limits = c(-1, 1), 
        caption = FALSE, 
        loops = loops, 
        borders = borders
    ),
    nrow = 1
)
```

![](matrix-centric_files/figure-html/unnamed-chunk-12-1.png)

Scale for `despeckled` scores

`despeckled` scores are in the same scale than the `scores` they were computed from.

## 5.2 Operations between multiple matrices

### 5.2.1 Merging maps

Hi-C libraries are often sequenced in multiple rounds, for example when high genome coverage is required. This results in multiple contact matrix files being generated. The `merge` function can be used to bind several `HiCExperiment` objects into a single one.

The different `HiCExperiment` objects do not need to all have identical `regions`, as shown in the following example.

```
hic_sub1 <- subsetByOverlaps(hic, GRanges("II:100001-200000"))
hic_sub2 <- subsetByOverlaps(hic, GRanges("II:300001-400000"))
bound_hic <- merge(hic_sub1, hic_sub2)
plotMatrix(bound_hic)
```

![](matrix-centric_files/figure-html/unnamed-chunk-13-1.png)

### 5.2.2 Computing ratio between two maps

Comparing two Hi-C maps can be useful to infer which genomic loci are differentially interacting between experimental conditions. Comparing two `HiCExperiment` objects can be done in `R` using the `divide` function.

For example, we can divide the *eco1* mutant Hi-C data by wild-type Hi-C dataset using the `divide` function.

```
hic_eco1 <- import(
    CoolFile(HiContactsData('yeast_eco1', 'mcool')), 
    focus = 'II', 
    resolution = 2000
)
##  see ?HiContactsData and browseVignettes('HiContactsData') for documentation
##  loading from cache
```

```
div_contacts <- divide(hic_eco1, by = hic) 
div_contacts
##  `HiCExperiment` object with 996,154 contacts over 407 regions 
##  -------
##  fileName: N/A 
##  focus: "II" 
##  resolutions(1): 2000
##  active resolution: 2000 
##  interactions: 60894 
##  scores(6): count.x balanced.x count.by balanced.by balanced.fc balanced.l2fc 
##  topologicalFeatures: () 
##  pairsFile: N/A 
##  metadata(2): hce_list operation
```

We can visually compare wild-type and *eco1* maps side by side (left) and their ratio map (right). This highlights the depletion of short-range and increase of long-range interactions in the *eco1* dataset.

```
cowplot::plot_grid(
    plotMatrix(hic_eco1, compare.to = hic, limits = c(-4, -1)), 
    plotMatrix(
        div_contacts, 
        use.scores = 'balanced.fc', 
        scale = 'log2', 
        limits = c(-1, 1),
        cmap = bwrColors()
    )
)
##  [1] "/home/biocbuild/.cache/R/ExperimentHub/15c47362e415a0_7754 | /home/biocbuild/.cache/R/ExperimentHub/15c3e2641b8180_7752"
```

![](matrix-centric_files/figure-html/unnamed-chunk-16-1.png)

Back to top