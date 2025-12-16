Source URL: https://bioconductor.org/books/release/csawBook/correction-for-multiple-testing.html
Date Scraped: 2025-12-15

---

# Chapter 7 Correction for multiple testing



## 7.1 Overview

The false discovery rate (FDR) is usually the most appropriate measure of error for high-throughput experiments.
Control of the FDR can be provided by applying the Benjamini-Hochberg (BH) method (Benjamini and Hochberg [1995](#ref-benjamini1995)) to a set of p-values.
This is less conservative than the alternatives (e.g., Bonferroni) yet still provides some measure of error control.
The most obvious approach is to apply the BH method to the set of p-values across all windows.
This will control the FDR across the set of putative DB windows.

However, the FDR across all detected windows is not necessarily the most relevant error rate.
Interpretation of ChIP-seq experiments is more concerned with regions of the genome in which (differential) protein binding is found, rather than the individual windows.
In other words, the FDR across all detected DB regions is usually desired.
This is not equivalent to that across all DB windows as each region will often consist of multiple overlapping windows.
Control of one will not guarantee control of the other (Lun and Smyth [2014](#ref-lun2014)).

To illustrate this difference, consider an analysis where the FDR across all window positions is controlled at 10%.
In the results, there are 18 adjacent window positions in one region and 2 windows in a separate region.
The first set of windows is a truly DB region whereas the second set is a false positive.
A window-based interpretation of the FDR is correct as only 2 of the 20 window positions are false positives.
However, a region-based interpretation results in an actual FDR of 50%.

To avoid misinterpretation of the FDR, *[csaw](https://bioconductor.org/packages/3.22/csaw)* provides a number of strategies to obtain region-level results.
This involves defining the regions of interest - possibly from the windows themselves -
and converting per-window statistics into a p-value for each region.
Application of the BH method to the per-region p-values will then control the relevant FDR across regions.
These strategies are demonstrated below using the NF-YA data.

## 7.2 Grouping windows into regions

### 7.2.1 Quick and dirty clustering

The `mergeWindows()` function provides a simple single-linkage algorithm to cluster windows into regions.
Windows that are less than `tol` apart are considered to be adjacent and are grouped into the same cluster.
The chosen `tol` represents the minimum distance at which two binding events are treated as separate sites.
Large values (500 - 1000 bp) reduce redundancy and favor a region-based interpretation of the results,
while smaller values (< 200 bp) allow resolution of individual binding sites.

View set-up code

```
library(csaw)
merged <- mergeWindows(filtered.data, tol=1000L)
merged$regions
```

```
## GRanges object with 3577 ranges and 0 metadata columns:
##                      seqnames            ranges strand
##                         <Rle>         <IRanges>  <Rle>
##      [1]                 chr1   7397901-7398110      *
##      [2]                 chr1   9541401-9541510      *
##      [3]                 chr1   9545301-9545360      *
##      [4]                 chr1 10007401-10007460      *
##      [5]                 chr1 13134451-13134510      *
##      ...                  ...               ...    ...
##   [3573] chrX_GL456233_random     336801-336910      *
##   [3574]                 chrY     143051-143060      *
##   [3575]                 chrY     259151-259210      *
##   [3576]                 chrY 90808851-90808860      *
##   [3577]                 chrY 90812851-90812910      *
##   -------
##   seqinfo: 66 sequences from an unspecified genome
```

If many adjacent windows are present, very large clusters may be formed that are difficult to interpret.
We perform a simple check below to determine whether most clusters are of an acceptable size.
Huge clusters indicate that more aggressive filtering from Chapter [4](chap-filter.html#chap-filter) is required.  
This mitigates chaining effects by reducing the density of windows in the genome.

```
summary(width(merged$regions))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##      10      60     110     165     160   15660
```

Alternatively, chaining can be limited by setting `max.width` to restrict the size of the merged intervals.
Clusters substantially larger than `max.width` are split into several smaller subclusters of roughly equal size.
The chosen value should be small enough so as to separate DB regions from unchanged neighbors,
yet large enough to avoid misinterpretation of the FDR.
Any value from 2000 to 10000 bp is recommended.
This paramater can also interpreted as the maximum distance at which two binding sites are considered part of the same event.

```
merged.max <- mergeWindows(filtered.data, tol=1000L, max.width=5000L)
summary(width(merged.max$regions))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##      10      60     110     164     160    4860
```

### 7.2.2 Using external information

Another approach is to group together windows that overlap with a pre-specified region of interest.
The most obvious source of pre-specified regions is that of annotated features such as promoters or gene bodies.
Alternatively, called peaks can be used provided that sufficient care has been taken to avoid loss of error control from data snooping (Lun and Smyth [2014](#ref-lun2014)).
Regardless of how they are specified, each region of interest corresponds to a group that contains all overlapping windows,  
as identified by the `findOverlaps` function from the *[GenomicRanges](https://bioconductor.org/packages/3.22/GenomicRanges)* package.

```
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
broads <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
broads <- resize(broads, width(broads)+3000, fix="end")

olap <- findOverlaps(broads, rowRanges(filtered.data))
olap
```

```
## Hits object with 12867 hits and 0 metadata columns:
##           queryHits subjectHits
##           <integer>   <integer>
##       [1]         7        6995
##       [2]        18        8323
##       [3]        18        8324
##       [4]        18        8325
##       [5]        18        8326
##       ...       ...         ...
##   [12863]     24521        6840
##   [12864]     24521        6841
##   [12865]     24524        6601
##   [12866]     24524        6602
##   [12867]     24524        6603
##   -------
##   queryLength: 24528 / subjectLength: 12352
```

At this point, one might imagine that it would be simpler to just collect and analyze counts over the pre-specified regions.
This is a valid strategy but will yield different results.
Consider a promoter containing two separate sites that are identically DB in opposite directions.
Counting reads across the promoter will give equal counts for each condition so changes within the promoter will not be detected.
Similarly, imprecise peak boundaries can lead to loss of detection power due to “contamination” by reads in background regions.
Window-based methods may be more robust as each interval of the promoter/peak region is examined separately (Lun and Smyth [2014](#ref-lun2014)),
avoiding potential problems with peak-calling errors and incorrect/incomplete annotation.

## 7.3 Obtaining per-region p-value

### 7.3.1 Combining window-level p-values

We compute a combined p-value for each region based on the p-values of the constituent windows (Simes [1986](#ref-simes1986)).
This tests the joint null hypothesis for each region, i.e., that no enrichment is observed across any of its windows.
Any DB within the region will reject the joint null and yield a low p-value for the entire region.
The combined p-values are then adjusted using the BH method to control the region-level FDR.

```
tabcom <- combineTests(merged$ids, rowData(filtered.data))
is.sig.region <- tabcom$FDR <= 0.05
summary(is.sig.region)
```

```
##    Mode   FALSE    TRUE 
## logical    1593    1984
```

Summarizing the direction of DB for each cluster requires some care as the direction of DB can differ between constituent windows.
The `num.up.tests` and `num.down.tests` fields contain the number of windows that change in each direction,
and can be used to gauge whether binding increases or decreases across the cluster.
A complex DB event may be present if both `num.up.tests` and `num.down.tests` are non-zero
(i.e., opposing changes within the region) or if the total number of windows is much larger than either number
(e.g., interval of constant binding adjacent to the DB interval).

Alternatively, the `direction` field specifies which DB direction contributes to the combined p-value.
If `"up"`, the combined p-value for this cluster is driven by p-values of windows with positive log-fold changes.
If `"down"`, the combined p-value is driven by windows with negative log-fold changes.
If `"mixed"`, windows with both positive and negative log-fold changes are involved.
This allows the dominant DB in significant clusters to be quickly summarized, as shown below.

```
table(tabcom$direction[is.sig.region])
```

```
## 
## down   up 
##  178 1806
```

For pre-specified regions, the `combineOverlaps()` function will combine the p-values for all windows in each region.
This is a wrapper around `combineTests()` for `Hits` objects.
It returns a single combined p-value (and its BH-adjusted value) for each region.
Regions that do not overlap any windows have values of `NA` in all fields for the corresponding rows.

```
tabbroad <- combineOverlaps(olap, rowData(filtered.data))
head(tabbroad[!is.na(tabbroad$PValue),])
```

```
## DataFrame with 6 rows and 8 columns
##    num.tests num.up.logFC num.down.logFC      PValue         FDR   direction
##    <integer>    <integer>      <integer>   <numeric>   <numeric> <character>
## 7          1            1              0 2.55503e-05 0.000479036          up
## 18         4            4              0 1.58040e-04 0.001420838          up
## 23         3            3              0 2.56903e-02 0.045240657          up
## 25         2            0              0 7.33967e-01 0.756264158          up
## 28         3            3              0 1.06179e-04 0.001128648          up
## 36         4            4              0 5.21069e-03 0.014422161          up
##     rep.test rep.logFC
##    <integer> <numeric>
## 7       6995  3.396598
## 18      8326  3.445923
## 23       315  1.331098
## 25      9977  0.201879
## 28      8774  3.503338
## 36      2716  2.105379
```

```
is.sig.gene <- tabcom$FDR <= 0.05
table(tabbroad$direction[is.sig.gene])
```

```
## 
##  down mixed    up 
##    94    30  1991
```

### 7.3.2 Based on the most significant window

Another approach is to use the single window with the strongest DB as a representative of the entire region.
This is useful when a log-fold change is required for each cluster, e.g., for plotting.
(In contrast, taking the average log-fold change across all windows in a region will understate the magnitude of DB,
especially if the region includes some non-DB background intervals of the genome.)
Identification of the most significant (i.e., “best”) window is performed using the `getBestTest()` function.
This reports the index of the window with the lowest p-value in each cluster as well as the associated statistics.

```
tab.best <- getBestTest(merged$ids, rowData(filtered.data))
head(tab.best)
```

```
## DataFrame with 6 rows and 8 columns
##   num.tests num.up.logFC num.down.logFC    PValue       FDR   direction
##   <integer>    <integer>      <integer> <numeric> <numeric> <character>
## 1         5            2              0 0.0172290 0.0361031          up
## 2         3            2              0 0.0121320 0.0277825          up
## 3         2            2              0 0.0275741 0.0518027          up
## 4         2            0              0 0.1882341 0.2401260          up
## 5         2            0              0 0.1627876 0.2124375          up
## 6         5            1              0 0.0362304 0.0637777          up
##    rep.test rep.logFC
##   <integer> <numeric>
## 1         3   1.57358
## 2         8   1.98351
## 3        10   1.69016
## 4        11   1.12022
## 5        14   1.08339
## 6        17   1.32711
```

A Bonferroni correction is applied to the p-value of the best window in each region,
based on the number of constituent windows in that region.
This is necessary to account for the implicit multiple testing across all windows in each region.
The corrected p-value is reported as `PValue` in `tab.best`,
and can be used for correction across regions using the BH method to control the region-level FDR.

In addition, it is often useful to report the start location of the best window within each cluster.
This allows users to easily identify a relevant DB subinterval in large regions.
For example, the sequence of the DB subinterval can be extracted for motif discovery.

```
tabcom$rep.start <- start(rowRanges(filtered.data))[tab.best$rep.test]
head(tabcom[,c("rep.logFC", "rep.start")])
```

```
## DataFrame with 6 rows and 2 columns
##   rep.logFC rep.start
##   <numeric> <integer>
## 1   1.57358   7398001
## 2   1.81157   9541501
## 3   1.57463   9545351
## 4   1.03337  10007401
## 5   1.08339  13134501
## 6   1.32711  13372551
```

The same approach can be applied to the overlaps between windows and pre-specified regions,
using the `getBestOverlaps()` wrapper function.
This is demonstrated below for the broad gene body example.
As with `combineOverlaps()`, regions with no windows are assigned `NA` in the output table,
but these are removed here to show some actual results.

```
tab.best.broad <- getBestOverlaps(olap, rowData(filtered.data))
tabbroad$rep.start <- start(rowRanges(filtered.data))[tab.best.broad$rep.test]
head(tabbroad[!is.na(tabbroad$PValue),c("rep.logFC", "rep.start")])
```

```
## DataFrame with 6 rows and 2 columns
##    rep.logFC rep.start
##    <numeric> <integer>
## 7   3.396598  32657101
## 18  3.445923   8259301
## 23  1.331098  92934601
## 25  0.201879  71596101
## 28  3.503338   4137001
## 36  2.105379 100187601
```

### 7.3.3 Wrapper functions

For convenience, the steps of merging windows and computing statistics are implemented in a single wrapper function.
This simply calls `mergeWindows()` followed by `combineTests()` and `getBestTest()`.

```
merge.res <- mergeResults(filtered.data, rowData(filtered.data), tol=100,
    merge.args=list(max.width=5000))
names(merge.res)
```

```
## [1] "regions"  "combined" "best"
```

An equivalent wrapper function is also available for handling overlaps to pre-specified regions.
This simply calls `findOverlaps()` followed by `combineOverlaps()` and `getBestOverlaps()`.

```
broad.res <- overlapResults(filtered.data, regions=broads,
    tab=rowData(filtered.data))
names(broad.res)
```

```
## [1] "regions"  "combined" "best"
```

## 7.4 Squeezing out more detection power

### 7.4.1 Integrating across multiple window sizes

Repeating the analysis with different window sizes may uncover new DB events at different resolutions.
Multiple sets of DB results are integrated by clustering adjacent windows together (even if they differ in size) and combining p-values within each of the resulting clusters.
The example below uses the H3 acetylation data from Chapter [5](chap-norm.html#chap-norm).
Some filtering is performed to avoid excessive chaining in this demonstration.
Corresponding tables of DB results should also be obtained – for brevity, mock results are used here.

```
library(chipseqDBData)
ac.files <- H3K9acData()$Path
ac.small <- windowCounts(ac.files, width=150L, spacing=100L, 
    filter=25, param=param)
ac.large <- windowCounts(ac.files, width=1000L, spacing=500L, 
    filter=35, param=param)

# TODO: actually do the analysis here.
# In the meantime, mocking up results for demonstration purposes.
ns <- nrow(ac.small)
mock.small <- data.frame(logFC=rnorm(ns), logCPM=0, PValue=runif(ns)) 
nl <- nrow(ac.large)
mock.large <- data.frame(logFC=rnorm(nl), logCPM=0, PValue=runif(nl))
```

The `mergeResultsList()` function merges windows of all sizes into a single set of regions,
and computes a combined p-value from the associated p-values for each region.
Equal contributions from each window size are enforced by setting `equiweight=TRUE`,
which uses a weighted version of Simes’ method (Benjamini and Hochberg [1997](#ref-benjamini1997)).
The weight assigned to each window is inversely proportional to the number of windows of that size in the same cluster.
This avoids the situation where, if a cluster contains many small windows,
the DB results for the analysis with the small window size contribute most to the combined p-value.
This is not ideal when results from all window sizes are of equal interest.

```
cons.res <- mergeResultsList(list(ac.small, ac.large), 
    tab.list=list(mock.small, mock.large), 
    equiweight=TRUE, tol=1000)
cons.res$regions
```

```
## GRanges object with 30486 ranges and 0 metadata columns:
##           seqnames            ranges strand
##              <Rle>         <IRanges>  <Rle>
##       [1]     chr1   4774001-4776500      *
##       [2]     chr1   4784501-4787000      *
##       [3]     chr1   4806501-4809000      *
##       [4]     chr1   4856501-4860000      *
##       [5]     chr1   5082501-5084500      *
##       ...      ...               ...    ...
##   [30482]     chrY 38230601-38230750      *
##   [30483]     chrY 73037501-73039000      *
##   [30484]     chrY 75445901-75446150      *
##   [30485]     chrY 88935501-88937000      *
##   [30486]     chrY 90812501-90814000      *
##   -------
##   seqinfo: 66 sequences from an unspecified genome
```

```
cons.res$combined
```

```
## DataFrame with 30486 rows and 8 columns
##       num.tests num.up.logFC num.down.logFC    PValue       FDR   direction
##       <integer>    <integer>      <integer> <numeric> <numeric> <character>
## 1             4            0              0  0.943272  0.999837       mixed
## 2            13            0              0  0.141228  0.991002       mixed
## 3            11            0              0  0.102143  0.991002       mixed
## 4            22            0              0  0.621215  0.996897       mixed
## 5             9            0              0  0.458395  0.996622       mixed
## ...         ...          ...            ...       ...       ...         ...
## 30482         1            0              0  0.948568  0.999837          up
## 30483         5            0              0  0.580198  0.996897       mixed
## 30484         2            0              0  0.539253  0.996897       mixed
## 30485         5            0              0  0.290521  0.993682          up
## 30486         2            0              0  0.528792  0.996897        down
##        rep.test rep.logFC
##       <integer> <numeric>
## 1        238289 -0.312926
## 2        238294 -0.953913
## 3        238299  0.806627
## 4        238304 -0.750927
## 5            35  0.234733
## ...         ...       ...
## 30482    238280  0.493052
## 30483    238282  0.690223
## 30484    238284  0.800469
## 30485    238286  1.045164
## 30486    391253 -0.517564
```

Similarly, the `overlapResultsList()` function is used to merge windows of varying size that overlap pre-specified regions.

```
cons.broad <- overlapResultsList(list(ac.small, ac.large),
    tab.list=list(mock.small, mock.large), 
    equiweight=TRUE, region=broads)
cons.broad$regions
```

```
## GRanges object with 24528 ranges and 1 metadata column:
##             seqnames              ranges strand |     gene_id
##                <Rle>           <IRanges>  <Rle> | <character>
##   100009600     chr9   21062393-21076096      - |   100009600
##   100009609     chr7   84935565-84967115      - |   100009609
##   100009614    chr10   77708457-77712009      + |   100009614
##   100009664    chr11   45805087-45841171      + |   100009664
##      100012     chr4 144157557-144165663      - |      100012
##         ...      ...                 ...    ... .         ...
##       99889     chr3   84496093-85890516      - |       99889
##       99890     chr3 110246109-110253998      - |       99890
##       99899     chr3 151730922-151752960      - |       99899
##       99929     chr3   65525410-65555518      + |       99929
##       99982     chr4 136550540-136605723      - |       99982
##   -------
##   seqinfo: 66 sequences (1 circular) from mm10 genome
```

```
cons.res$combined
```

```
## DataFrame with 30486 rows and 8 columns
##       num.tests num.up.logFC num.down.logFC    PValue       FDR   direction
##       <integer>    <integer>      <integer> <numeric> <numeric> <character>
## 1             4            0              0  0.943272  0.999837       mixed
## 2            13            0              0  0.141228  0.991002       mixed
## 3            11            0              0  0.102143  0.991002       mixed
## 4            22            0              0  0.621215  0.996897       mixed
## 5             9            0              0  0.458395  0.996622       mixed
## ...         ...          ...            ...       ...       ...         ...
## 30482         1            0              0  0.948568  0.999837          up
## 30483         5            0              0  0.580198  0.996897       mixed
## 30484         2            0              0  0.539253  0.996897       mixed
## 30485         5            0              0  0.290521  0.993682          up
## 30486         2            0              0  0.528792  0.996897        down
##        rep.test rep.logFC
##       <integer> <numeric>
## 1        238289 -0.312926
## 2        238294 -0.953913
## 3        238299  0.806627
## 4        238304 -0.750927
## 5            35  0.234733
## ...         ...       ...
## 30482    238280  0.493052
## 30483    238282  0.690223
## 30484    238284  0.800469
## 30485    238286  1.045164
## 30486    391253 -0.517564
```

In this manner, DB results from multiple window widths can be gathered together and reported as a single set of regions.
Consolidation is most useful for histone marks and other analyses involving diffuse regions of enrichment.
For such studies, the ideal window size is not known or may not even exist,
e.g., if the widths of the enriched regions or DB subintervals are variable.

### 7.4.2 Weighting windows on abundance

Windows that are more likely to be DB can be upweighted to improve detection power.
For example, in TF ChIP-seq data, the window of highest abundance within each enriched region probably contains the binding site.
It is reasonable to assume that this window will also have the strongest DB.
To improve power, the weight assigned to the most abundant window is increased relative to that of other windows in the same cluster.
This means that the p-value of this window will have a greater influence on the final combined p-value.

Weights are computed in a manner to minimize conservativeness relative to the optimal unweighted approaches in each possible scenario.
If the strongest DB event is at the most abundant window, the weighted approach will yield a combined p-value that is no larger than twice the p-value of the most abundant window.
(Here, the optimal approach would be to use the p-value of the most abundance window directly as a proxy for the p-value of the cluster.)
If the strongest DB event is *not* at the most abundant window, the weighted approach will yield a combined p-value that is no larger than twice the combined p-value without wweighting (which is optimal as all windows have equal probabilities of containing the strongest DB).
All windows have non-zero weights, which ensures that any DB events in the other windows will still be considered when the p-values are combined.

The application of this weighting scheme is demonstrated in the example below.
First, the `getBestTest} function with \Rcode{by.pval=FALSE()` is used to identify the most abundant window in each cluster.
Window-specific weights are then computed using the `upweightSummits} function, and supplied to \Rcode{combineTests()` to use in computing combined p-values.

```
tab.ave <- getBestTest(merged$id, rowData(filtered.data), by.pval=FALSE)
weights <- upweightSummit(merged$id, tab.ave$rep.test)
head(weights)
```

```
## [1] 1 5 1 1 1 1
```

```
tabcom.w <- combineTests(merged$id, rowData(filtered.data), weight=weights)
head(tabcom.w)
```

```
## DataFrame with 6 rows and 8 columns
##   num.tests num.up.logFC num.down.logFC     PValue       FDR   direction
##   <integer>    <integer>      <integer>  <numeric> <numeric> <character>
## 1         5            3              0 0.01035174 0.0242018          up
## 2         3            2              0 0.00753833 0.0196020          up
## 3         2            2              0 0.02509503 0.0461279          up
## 4         2            0              0 0.11664146 0.1541857          up
## 5         2            0              0 0.12209068 0.1602636          up
## 6         5            2              0 0.01304296 0.0282927          up
##    rep.test rep.logFC
##   <integer> <numeric>
## 1         2   1.40804
## 2         7   1.81157
## 3         9   1.57463
## 4        12   1.03337
## 5        14   1.08339
## 6        17   1.32711
```

The weighting approach can also be applied to the clusters from the broad gene body example.
This is done by replacing the call to `getBestTest} with one to \Rfunction{getBestOverlaps()`, as before.
Similarly, `upweightSummit} can be replaced with \Rfunction{summitOverlaps()`.
These wrappers are designed to minimize book-keeping problems when one window overlaps multiple regions.

```
broad.best <- getBestOverlaps(olap, rowData(filtered.data), by.pval=FALSE)
head(broad.best[!is.na(broad.best$PValue),])
```

```
## DataFrame with 6 rows and 8 columns
##    num.tests num.up.logFC num.down.logFC      PValue         FDR   direction
##    <integer>    <integer>      <integer>   <numeric>   <numeric> <character>
## 7          1            1              0 2.55503e-05 0.000483900          up
## 18         4            1              0 6.00178e-03 0.016513735          up
## 23         3            1              0 2.56903e-02 0.047031665          up
## 25         2            0              0 7.33967e-01 0.754181359          up
## 28         3            1              0 3.53931e-05 0.000581726          up
## 36         4            1              0 2.40124e-03 0.008896761          up
##     rep.test rep.logFC
##    <integer> <numeric>
## 7       6995  3.396598
## 18      8324  1.690223
## 23       315  1.331098
## 25      9977  0.201879
## 28      8774  3.503338
## 36      2717  1.957661
```

```
broad.weights <- summitOverlaps(olap, region.best=broad.best$rep.test)
tabbroad.w <- combineOverlaps(olap, rowData(filtered.data), o.weight=broad.weights)
```

### 7.4.3 Filtering after testing but before correction

Most of the filters in Chapter~[4](chap-filter.html#chap-filter) are applied before the statistical analysis.
However, some of the approaches may be too aggressive, e.g., filtering to retain only local maxima or based on pre-defined regions.
In such cases, it may be preferable to initially apply one of the other, milder filters.
This ensures that sufficient windows are retained for stable normalization and/or EB shrinkage.
The aggressive filters can then be applied after the window-level statistics have been calculated, but before clustering into regions and calculation of cluster-level statistics.
This is still beneficial as it removes irrelevant windows that would increase the severity of the BH correction.
It may also reduce chaining effects during clustering.

## 7.5 FDR control in difficult situations

### 7.5.1 Clustering only on DB windows for diffuse marks

The clustering procedures described above rely on independent filtering to remove irrelevant windows.
This ensures that the regions of interest are reasonably narrow and can be easily interpreted,
which is typically the case for most protein targets, e.g., TFs, narrow histone marks.
However, enriched regions may be very large for more diffuse marks.
Such regions may be difficult to interpret when only the DB subinterval is of interest.
To overcome this, a post-hoc analysis can be performed whereby only significant windows are used for clustering.

```
postclust <- clusterWindows(rowRanges(filtered.data), rowData(filtered.data),
                            target=0.05, tol=100, max.width=1000)
postclust$FDR
```

```
## [1] 0.04978
```

```
postclust$region
```

```
## GRanges object with 2069 ranges and 0 metadata columns:
##                      seqnames              ranges strand
##                         <Rle>           <IRanges>  <Rle>
##      [1]                 chr1     7397951-7398010      *
##      [2]                 chr1     9541451-9541510      *
##      [3]                 chr1   13372551-13372560      *
##      [4]                 chr1   13590001-13590010      *
##      [5]                 chr1   15805551-15805660      *
##      ...                  ...                 ...    ...
##   [2065]                 chrX 104482701-104482760      *
##   [2066]                 chrX 106187051-106187060      *
##   [2067]                 chrX 136741351-136741360      *
##   [2068]                 chrX 140456551-140456560      *
##   [2069] chrX_GL456233_random       336801-336910      *
##   -------
##   seqinfo: 66 sequences from an unspecified genome
```

This will define and cluster significant windows in a manner that controls the cluster-level FDR at 5%.
The clustering step itself is performed using `mergeWindows()` with the specified parameters.
Each cluster consists entirely of DB windows and can be directly interpreted as a DB region or a DB subinterval of a larger enriched region.
This reduces the pressure on abundance filtering to obtain well-separated regions prior to clustering, e.g., for diffuse marks or in data sets with weak IP signal.
That said, users should be aware that calculation of the cluster-level FDR is not entirely rigorous.
As such, independent clustering and FDR control via Simes’ method should be considered as the default for routine analyses.

### 7.5.2 Using the empirical FDR for noisy data

Some analyses involve comparisons of ChIP samples to negative controls.
In such cases, any region exhibiting enrichment in the negative control over the ChIP samples must be a false positive.
The number of significant regions that change in the “wrong” direction can be used as an estimate of the number of false positives at any given p-value threshold.
Division by the number of discoveries changing in the “right” direction yields an estimate of the FDR, i.e., the empirical FDR (Zhang et al. [2008](#ref-zhang2008)).
This strategy is implemented in the `empiricalFDR()` function, which controls the empirical FDR across clusters based on their combined p-values.
Its use is demonstrated below, though the output is not meaningful in this situation as genuine changes in binding can be present in both directions.

```
empres <- empiricalFDR(merged$id, rowData(filtered.data))
```

The empirical FDR is useful for analyses of noisy data with high levels of non-specific binding.
This is because the estimate of the number of false positives adapts to the observed number of regions exhibiting enrichment in the negative controls.
In contrast, the standard BH method in `combineTests()` relies on proper type I error control during hypothesis testing.
As non-specific binding events tend to be condition-specific, they are indistinguishable from DB events and assigned low p-values, resulting in loss of FDR control.
Thus, for noisy data, use of the empirical FDR may be more appropriate to control the proportion of “experimental” false positives.
However, calculation of the empirical FDR is not as statistically rigorous as that of the BH method, so users are advised to only apply it when necessary.

### 7.5.3 Detecting complex DB

Complex DB events involve changes to the shape of the binding profile, not just a scaling increase/decrease to binding intensity.
Such regions may contain multiple sites that change in binding strength in opposite directions,
or peaks that change in width or position between conditions.
This often manifests as DB in opposite directions in different subintervals of a region.
Some of these events can be identified using the `mixedTests()` function.

```
tab.mixed <- mixedTests(merged$ids, rowData(filtered.data))
tab.mixed
```

```
## DataFrame with 3577 rows and 10 columns
##      num.tests num.up.logFC num.down.logFC    PValue       FDR   direction
##      <integer>    <integer>      <integer> <numeric> <numeric> <character>
## 1            5            5              0  0.998277         1       mixed
## 2            3            3              0  0.997978         1       mixed
## 3            2            2              0  0.993106         1       mixed
## 4            2            0              0  0.952941         1       mixed
## 5            2            0              0  0.959303         1       mixed
## ...        ...          ...            ...       ...       ...         ...
## 3573         3            0              3  1.000000         1       mixed
## 3574         1            0              0  0.936583         1       mixed
## 3575         2            0              0  0.856997         1       mixed
## 3576         1            1              0  0.971223         1       mixed
## 3577         2            0              0  0.754199         1       mixed
##      rep.up.test rep.up.logFC rep.down.test rep.down.logFC
##        <integer>    <numeric>     <integer>      <numeric>
## 1              3      1.57358             3        1.57358
## 2              7      1.81157             8        1.98351
## 3              9      1.57463            10        1.69016
## 4             12      1.03337            11        1.12022
## 5             14      1.08339            14        1.08339
## ...          ...          ...           ...            ...
## 3573       12345    -2.607150         12345      -2.607150
## 3574       12347    -0.987682         12347      -0.987682
## 3575       12349    -0.636873         12348      -0.440188
## 3576       12350     1.246123         12350       1.246123
## 3577       12351    -0.430191         12352      -0.163265
```

`mixedTests()` converts the p-value for each window into two one-sided p-values.
The one-sided p-values in each direction are combined using Simes’ method,
and the two one-sided combined p-values are themselves combined using an intersection-union test (Berger and Hsu [1996](#ref-berger1996bioequivalence)).
The resulting p-value is only low if a region contains strong DB in both directions.

`combineTests()` also computes some statistics for informal detection of complex DB.
For example, the `num.up.tests` and `num.down.tests` fields can be used to identify regions with changes in both directions.
The `direction` field will also label some regions as `"mixed"`, though this is not comprehensive.
Indeed, regions labelled as `"up"` or `"down"` in the `direction` field may also correspond to complex DB events,
but will not be labelled as `"mixed"` if the significance calculations are dominated by windows changing in only one direction.

### 7.5.4 Enforcing a minimal number of DB windows

On occasion, we may be interested in genomic regions that contain at least a minimal number or proportion of DB windows.
This is motivated by the desire to avoid detecting DB regions where only a small subinterval exhibits a change,
instead favoring more systematic changes throughout the region that are easier to interpret.
We can identify these regions using the `minimalTests()` function.

```
tab.min <- minimalTests(merged$ids, rowData(filtered.data),
    min.sig.n=3, min.sig.prop=0.5)
tab.min
```

```
## DataFrame with 3577 rows and 8 columns
##      num.tests num.up.logFC num.down.logFC      PValue         FDR   direction
##      <integer>    <integer>      <integer>   <numeric>   <numeric> <character>
## 1            5            2              0   0.0835097   0.1513244          up
## 2            3            2              0   0.0988417   0.1707179          up
## 3            2            2              0   0.0275741   0.0719419          up
## 4            2            0              0   0.1882341   0.2739601          up
## 5            2            0              0   0.2330766   0.3209065          up
## ...        ...          ...            ...         ...         ...         ...
## 3573         3            0              3 7.97835e-06 0.000636499        down
## 3574         1            0              0 1.26834e-01 0.204823467        down
## 3575         2            0              0 5.72013e-01 0.670039730        down
## 3576         1            0              0 5.75535e-02 0.116873491          up
## 3577         2            0              0 9.83204e-01 1.000000000        down
##       rep.test rep.logFC
##      <integer> <numeric>
## 1            1  1.238149
## 2            6  1.094827
## 3            9  1.574630
## 4           12  1.033368
## 5           13  0.739037
## ...        ...       ...
## 3573     12346 -2.534975
## 3574     12347 -0.987682
## 3575     12348 -0.440188
## 3576     12350  1.246123
## 3577     12352 -0.163265
```

`minimalTests()` applies a Holm-Bonferroni correction to all windows in the same cluster and picks the xth-smallest adjusted p-value (where x is defined from `min.sig.n` and `min.sig.prop`).
This tests the joint null hypothesis that the per-window null hypothesis is false for fewer than x windows in the cluster.
If the xth-smallest p-value is low, this provides strong evidence against the joint null for that cluster.

As an aside, this function also has some utility outside of ChIP-seq contexts.
For example, we might want to obtain a single p-value for a gene set based on the presence of a minimal percentage of differentially expressed genes.
Alternatively, we may be interested in ranking genes in loss-of-function screens based on a minimal number of shRNA/CRISPR guides that exhibit a significant effect.
These problems are equivalent to that of identifying a genomic region with a minimal number of DB windows.

## Session information

View session info

### Bibliography

Benjamini, Y., and Y. Hochberg. 1995. “Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.” *J. R. Stat. Soc. Series B* 57: 289–300.

Benjamini, Y., and Y. Hochberg. 1995. “Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.” *J. R. Stat. Soc. Series B* 57: 289–300.

1997. “Multiple Hypotheses Testing with Weights.” *Scand. J. Stat.* 24: 407–18.

Berger, R. L., and J. C. Hsu. 1996. “Bioequivalence Trials, Intersection-Union Tests and Equivalence Confidence Sets.” *Statist. Sci.* 11 (4): 283–319.

Lun, A. T., and G. K. Smyth. 2014. “De novo detection of differentially bound regions for ChIP-seq data using peaks and windows: controlling error rates correctly.” *Nucleic Acids Res.* 42 (11): e95.

Simes, R. J. 1986. “An Improved Bonferroni Procedure for Multiple Tests of Significance.” *Biometrika* 73 (3): 751–54.

Zhang, Y., T. Liu, C. A. Meyer, J. Eeckhoute, D. S. Johnson, B. E. Bernstein, C. Nusbaum, et al. 2008. “Model-based analysis of ChIP-Seq (MACS).” *Genome Biol.* 9 (9): R137.