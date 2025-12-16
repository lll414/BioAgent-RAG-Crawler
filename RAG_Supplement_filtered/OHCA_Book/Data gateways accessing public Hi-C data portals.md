Source URL: https://bioconductor.org/books/release/OHCA/pages/disseminating.html
Date Scraped: 2025-12-15

---

1. [Advanced Hi-C topics](https://bioconductor.org/books/release/OHCA/pages/disseminating.html)
2. [8  Data gateways: accessing public Hi-C data portals](https://bioconductor.org/books/release/OHCA/pages/disseminating.html)

# 8  Data gateways: accessing public Hi-C data portals

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
library(fourDNData)
library(DNAZooData)
library(rtracklayer)
```

Aims

This chapter focuses on introducing two important portals hosting public Hi-C datasets: the [4DN Consortium](https://www.4dnucleome.org/) and the [DNA Zoo project](https://www.dnazoo.org/). Two `R` packages provide a programmatic access to these portals:

* `fourDNData`
* `DNAZooData`

The Hi-C experimental approach has gained significant traction across multiple fields related to genome biology, and several consortia developed large-scale programs based on this technique. The `fourDNData` and `DNAZooData` R packages were designed to accelerate the investigation of chromatin structure using these public resources.

## 8.1 4DN data portal

The 4D Nucleome Data Coordination and Integration Center ([DCIC](http://dcic.4dnucleome.org/)) has developed and actively maintains a [data portal](https://data.4dnucleome.org/browse) providing public access to a wealth of resources to investigate 3D chromatin architecture. Notably, 3D chromatin conformation libraries relying on different technologies (“in situ” or “dilution” Hi-C, Capture Hi-C, Micro-C, DNase Hi-C, …), generated by 50+ collaborating labs, were homogeneously processed, yielding more than 350 sets of processed files.

`fourDNData` (read *4DN-Data*) is a package giving programmatic access to these uniformly processed Hi-C contact files.

The `fourDNData()` function provides a gateway to 4DN-hosted Hi-C files, including contact matrices (in `.hic` or `.mcool`) and other Hi-C derived files such as annotated compartments, domains, insulation scores, or `.pairs` files.

```
library(fourDNData)
head(fourDNData())
##    experimentSetAccession     fileType     size organism experimentType
##  1           4DNES18BMU79        pairs 10151.53    mouse   in situ Hi-C
##  3           4DNES18BMU79          hic  5285.82    mouse   in situ Hi-C
##  4           4DNES18BMU79        mcool  6110.75    mouse   in situ Hi-C
##  5           4DNES18BMU79   boundaries     0.12    mouse   in situ Hi-C
##  6           4DNES18BMU79   insulation     7.18    mouse   in situ Hi-C
##  7           4DNES18BMU79 compartments     0.18    mouse   in situ Hi-C
##    details                              dataset
##  1   DpnII Hi-C on Mouse Olfactory System cells
##  3   DpnII Hi-C on Mouse Olfactory System cells
##  4   DpnII Hi-C on Mouse Olfactory System cells
##  5   DpnII Hi-C on Mouse Olfactory System cells
##  6   DpnII Hi-C on Mouse Olfactory System cells
##  7   DpnII Hi-C on Mouse Olfactory System cells
##                                                          condition
##  1 Mature olfactory sensory neurons with conditional Ldb1 knockout
##  3 Mature olfactory sensory neurons with conditional Ldb1 knockout
##  4 Mature olfactory sensory neurons with conditional Ldb1 knockout
##  5 Mature olfactory sensory neurons with conditional Ldb1 knockout
##  6 Mature olfactory sensory neurons with conditional Ldb1 knockout
##  7 Mature olfactory sensory neurons with conditional Ldb1 knockout
##                  biosource biosourceType             publication
##  1 olfactory receptor cell  primary cell Monahan K et al. (2019)
##  3 olfactory receptor cell  primary cell Monahan K et al. (2019)
##  4 olfactory receptor cell  primary cell Monahan K et al. (2019)
##  5 olfactory receptor cell  primary cell Monahan K et al. (2019)
##  6 olfactory receptor cell  primary cell Monahan K et al. (2019)
##  7 olfactory receptor cell  primary cell Monahan K et al. (2019)
##                                                                                                                                    URL
##  1 https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/49504f97-904e-48c1-8c20-1033680b66da/4DNFIC5AHBPV.pairs.gz
##  3      https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/6cd4378a-8f51-4e65-99eb-15f5c80abf8d/4DNFIT4I5C6Z.hic
##  4    https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/01fb704f-2fd7-48c6-91af-c5f4584529ed/4DNFIVPAXJO8.mcool
##  5   https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/5c07cdee-53e2-43e0-8853-cfe5f057b3f1/4DNFIR3XCIMA.bed.gz
##  6       https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/d1f4beb9-701f-4188-abe2-6271fe658770/4DNFIXKKNMS7.bw
##  7       https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/3d429647-51c8-4e3a-a18b-eec0b1480905/4DNFIN13N8C1.bw
```

### 8.1.1 Querying individual files

The `fourDNData()` function can be used to directly fetch specific files from the 4DN data portal:

```
cf <- fourDNData(experimentSetAccession = '4DNESJNPEKZD', type = 'mcool')
```

This effectively downloads and caches the queried file locally.

```
cf
##  [1] "/home/biocbuild/.cache/R/fourDNData/bf00639ad4c1d_4DNFIZL8OZE1.mcool"

availableChromosomes(cf)
##  Seqinfo object with 24 sequences from an unspecified genome:
##    seqnames seqlengths isCircular genome
##    chr1      248956422       <NA>   <NA>
##    chr2      242193529       <NA>   <NA>
##    chr3      198295559       <NA>   <NA>
##    chr4      190214555       <NA>   <NA>
##    chr5      181538259       <NA>   <NA>
##    ...             ...        ...    ...
##    chr20      64444167       <NA>   <NA>
##    chr21      46709983       <NA>   <NA>
##    chr22      50818468       <NA>   <NA>
##    chrX      156040895       <NA>   <NA>
##    chrY       57227415       <NA>   <NA>

availableResolutions(cf)
##  resolutions(13): 1000 2000 ... 5e+06 1e+07
##  

import(cf, focus = "chr4:10000001-20000000", resolution = 5000)
##  `HiCExperiment` object with 656 contacts over 2,000 regions 
##  -------
##  fileName: "/home/biocbuild/.cache/R/fourDNData/bf00639ad4c1d_4DNFIZL8OZE1.mcool" 
##  focus: "chr4:10,000,001-20,000,000" 
##  resolutions(13): 1000 2000 ... 5000000 10000000
##  active resolution: 5000 
##  interactions: 614 
##  scores(2): count balanced 
##  topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
##  pairsFile: N/A 
##  metadata(0):
```

Different Hi-C related genomic files are provided by the 4DN consortium. The type of file to fetch can be specified with the `type` argument:

* `type = 'pairs'` will fetch the pairs file which was generated by the 4DN consortium and binned into a contact matrix. Once fetched from the 4DN data portal, the local file can be imported in `R` using the `import` function, which will generate a `GInteractions` object.

```
## Not evaluated for now
pairs_f <- fourDNData(experimentSetAccession = '4DNESJNPEKZD', type = 'pairs') 
print(pairs_f)
import(pairs_f)
```

Watch out

`.pairs` files can be particularly large and therefore will take both and long time to download and a large storage footprint.

* `type = 'insulation'` will fetch a `.bigwig` track file precomputed by the 4DN consortium. This track corresponds to the genome-wide insulation score computed by `cooltools` as described in Crane et al. ([2015](https://bioconductor.org/books/release/OHCA/pages/disseminating.html#ref-Crane_2015)). To know more about this, read [the excerpt from 4DN data portal](https://data.4dnucleome.org/resources/data-analysis/insulation_compartment_scores#insulation_scores_and_boundaries_page_all). Once fetched from the 4DN data portal, the local file can be imported in `R` using the `import` function, which will generate a `RleList` object.

```
library(rtracklayer)
fourDNData(experimentSetAccession = '4DNES25ABNZ1', type = 'insulation') |> 
    import(as = 'Rle')
##  RleList of length 21
##  $chr1
##  numeric-Rle of length 195471971 with 38145 runs
##    Lengths:      3065000         5000 ...         5000       171971
##    Values :  0.00000e+00  1.01441e-01 ...     0.807009     0.000000
##  
##  $chr10
##  numeric-Rle of length 130694993 with 25100 runs
##    Lengths:     3175000        5000        5000 ...        5000      169993
##    Values :  0.00000000  0.37584546  0.33597839 ...    0.628601    0.000000
##  
##  $chr11
##  numeric-Rle of length 122082543 with 23536 runs
##    Lengths:    3165000       5000       5000 ...       5000     162543
##    Values :  0.0000000 -0.7906257 -0.7930040 ...   0.515919   0.000000
##  
##  $chr12
##  numeric-Rle of length 120129022 with 22578 runs
##    Lengths:   3075000      5000      5000 ...      5000      5000    164022
##    Values :  0.000000  0.411216  0.400357 ... 0.1650951 0.2175749 0.0000000
##  
##  $chr13
##  numeric-Rle of length 120421639 with 22807 runs
##    Lengths:     3080000        5000        5000 ...        5000      171639
##    Values :  0.00000000  0.17005745  0.10652249 ...  1.14856148  0.00000000
##  
##  ...
##  <16 more elements>
```

* `type = 'boundaries'` will fetch a `.bed` file precomputed by the 4DN consortium, listing the annotated borders between topological domains. These borders correspond to local minima identified from the genome-wide insulation track. It can also be imported in `R` using the `import` function, which will generate a `GRanges` object.

```
fourDNData(experimentSetAccession = '4DNES25ABNZ1', type = 'boundaries') |> 
    import()
##  GRanges object with 6103 ranges and 2 metadata columns:
##           seqnames            ranges strand |        name     score
##              <Rle>         <IRanges>  <Rle> | <character> <numeric>
##       [1]     chr1   4380001-4385000      * |      Strong  0.695274
##       [2]     chr1   4760001-4765000      * |        Weak  0.444476
##       [3]     chr1   4910001-4915000      * |        Weak  0.353184
##       [4]     chr1   5180001-5185000      * |      Strong  0.565763
##       [5]     chr1   6170001-6175000      * |      Strong  1.644911
##       ...      ...               ...    ... .         ...       ...
##    [6099]     chrY 89725001-89730000      * |        Weak  0.258094
##    [6100]     chrY 89790001-89795000      * |        Weak  0.442186
##    [6101]     chrY 89895001-89900000      * |        Weak  0.279879
##    [6102]     chrY 90025001-90030000      * |      Strong  0.660699
##    [6103]     chrY 90410001-90415000      * |      Strong  1.160018
##    -------
##    seqinfo: 21 sequences from an unspecified genome; no seqlengths
```

* `type = 'compartments'` will fetch a `.bigwig` track file precomputed by the 4DN consortium. This track corresponds to the selected genome-wide eigenvector computed by `cooltools` and representing A/B compartments. To know more about this, read [the excerpt from 4DN data portal](https://data.4dnucleome.org/resources/data-analysis/insulation_compartment_scores#compartments_page_all). Once fetched from the 4DN data portal, the local file can be imported in `R` using the `import` function, which will generate a `RleList` object. The score represents the eigenvector values, and by convention a genomic bin with a positive score is associated with the A compartment whereas a genomic bin with a negative score is associated with the B compartment.

```
fourDNData(experimentSetAccession = '4DNES25ABNZ1', type = 'compartments') |> 
    import()
##  GRanges object with 10911 ranges and 1 metadata column:
##            seqnames            ranges strand |     score
##               <Rle>         <IRanges>  <Rle> | <numeric>
##        [1]     chr1          1-250000      * |       NaN
##        [2]     chr1     250001-500000      * |       NaN
##        [3]     chr1     500001-750000      * |       NaN
##        [4]     chr1    750001-1000000      * |       NaN
##        [5]     chr1   1000001-1250000      * |       NaN
##        ...      ...               ...    ... .       ...
##    [10907]     chrY 90500001-90750000      * | 0.0237907
##    [10908]     chrY 90750001-91000000      * |       NaN
##    [10909]     chrY 91000001-91250000      * |       NaN
##    [10910]     chrY 91250001-91500000      * |       NaN
##    [10911]     chrY 91500001-91744698      * |       NaN
##    -------
##    seqinfo: 21 sequences from an unspecified genome
```

### 8.1.2 Querying a complete experiment dataset

Rather than importing multiple files corresponding to a single experimentSet accession ID one by one, one can import all the available files associated with a experimentSet accession ID into a `HiCExperiment` object by using the `fourDNHiCExperiment()` function.

```
hic <- fourDNHiCExperiment('4DNESJNPEKZD')
##  Fetching local Hi-C contact map from Bioc cache
##  Fetching local compartments bigwig file from Bioc cache
##  Insulation not found for the provided experimentSet accession.
##  Borders not found for the provided experimentSet accession.
##  Importing contacts in memory
```

This is a more efficient way to import datasets, as it aggregates the different bits together into a single `HiCExperiment` object with populated `topologicalFeatures` and `metadata` slots.

```
hic
##  `HiCExperiment` object with 453,301 contacts over 12,366 regions 
##  -------
##  fileName: "/home/biocbuild/.cache/R/fourDNData/bf00639ad4c1d_4DNFIZL8OZE1.mcool" 
##  focus: "whole genome" 
##  resolutions(13): 1000 2000 ... 5000000 10000000
##  active resolution: 250000 
##  interactions: 289086 
##  scores(2): count balanced 
##  topologicalFeatures: compartments(5437) borders(0) 
##  pairsFile: N/A 
##  metadata(2): 4DN_info eigens
```

```
metadata(hic)
##  $`4DN_info`
##       experimentSetAccession     fileType   size organism experimentType
##  1376           4DNESJNPEKZD        pairs   6.67    human   in situ Hi-C
##  1378           4DNESJNPEKZD          hic 179.51    human   in situ Hi-C
##  1379           4DNESJNPEKZD        mcool  30.17    human   in situ Hi-C
##  1380           4DNESJNPEKZD compartments   0.21    human   in situ Hi-C
##       details                                     dataset
##  1376    MboI Hi-C on GM12878 cells - protocol variations
##  1378    MboI Hi-C on GM12878 cells - protocol variations
##  1379    MboI Hi-C on GM12878 cells - protocol variations
##  1380    MboI Hi-C on GM12878 cells - protocol variations
##                                                       condition biosource
##  1376 in situ Hi-C on GM12878 with MboI and bio-dUTP (Tri-Link)   GM12878
##  1378 in situ Hi-C on GM12878 with MboI and bio-dUTP (Tri-Link)   GM12878
##  1379 in situ Hi-C on GM12878 with MboI and bio-dUTP (Tri-Link)   GM12878
##  1380 in situ Hi-C on GM12878 with MboI and bio-dUTP (Tri-Link)   GM12878
##                biosourceType          publication
##  1376 immortalized cell line Rao SS et al. (2014)
##  1378 immortalized cell line Rao SS et al. (2014)
##  1379 immortalized cell line Rao SS et al. (2014)
##  1380 immortalized cell line Rao SS et al. (2014)
##                                                                                                                                       URL
##  1376 https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/0bdd4745-7203-49d0-adf6-291cef1a96b7/4DNFIOZ7D1OQ.pairs.gz
##  1378      https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/1201682a-a223-482d-913d-3c3972b8eb65/4DNFIIRIHBR2.hic
##  1379    https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/356fab42-5562-4cfd-a3f8-592aa060b992/4DNFIZL8OZE1.mcool
##  1380       https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/333aabfd-b747-447c-b93a-8138f9488fad/4DNFIO9V5G93.bw
##  
##  $eigens
##  GRanges object with 11280 ranges and 2 metadata columns:
##            seqnames              ranges strand |       score       eigen
##               <Rle>           <IRanges>  <Rle> |   <numeric>   <numeric>
##        [1]     chr1      750001-1000000      * |   1.6911879   1.6911879
##        [2]     chr1     1000001-1250000      * |   0.0809129   0.0809129
##        [3]     chr1     1250001-1500000      * |   0.0690173   0.0690173
##        [4]     chr1     1500001-1750000      * |  -0.1903324  -0.1903324
##        [5]     chr1     1750001-2000000      * |   0.3283633   0.3283633
##        ...      ...                 ...    ... .         ...         ...
##    [11276]     chrX 154750001-155000000      * | -0.10909061 -0.10909061
##    [11277]     chrX 155000001-155250000      * | -1.39655280 -1.39655280
##    [11278]     chrX 155250001-155500000      * |  0.00264734  0.00264734
##    [11279]     chrX 155500001-155750000      * | -0.15279847 -0.15279847
##    [11280]     chrX 155750001-156000000      * | -1.41699576 -1.41699576
##    -------
##    seqinfo: 24 sequences from an unspecified genome
```

## 8.2 DNA Zoo

The [DNA Zoo Consortium](https://www.dnazoo.org/) is a collaborative group whose aim is to correct and refine genome assemblies across the tree of life using Hi-C approaches. As of 2023, they have performed Hi-C across more than 300 animal, plant and fungi species.

`DNAZooData` is a package giving programmatic access to these uniformly processed Hi-C contact files, as well as the refined genome assemblies.

The `DNAZooData()` function provides a gateway to DNA Zoo-hosted Hi-C files, fetching and caching relevant contact matrices in `.hic` format It returns a `HicFile` object, which can then be imported in memory using `import()`.

```
library(DNAZooData)
head(DNAZooData())
##                    species                              readme
##  1        Acinonyx_jubatus        Acinonyx_jubatus/README.json
##  2      Acropora_millepora      Acropora_millepora/README.json
##  3     Addax_nasomaculatus     Addax_nasomaculatus/README.json
##  4           Aedes_aegypti           Aedes_aegypti/README.json
##  5   Aedes_aegypti__AaegL4   Aedes_aegypti__AaegL4/README.json
##  6 Aedes_aegypti__AaegL5.0 Aedes_aegypti__AaegL5.0/README.json
##                                                            readme_link
##  1        https://dnazoo.s3.wasabisys.com/Acinonyx_jubatus/README.json
##  2      https://dnazoo.s3.wasabisys.com/Acropora_millepora/README.json
##  3     https://dnazoo.s3.wasabisys.com/Addax_nasomaculatus/README.json
##  4           https://dnazoo.s3.wasabisys.com/Aedes_aegypti/README.json
##  5   https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL4/README.json
##  6 https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL5.0/README.json
##    original_assembly     new_assembly
##  1           aciJub1      aciJub1_HiC
##  2       amil_sf_1.1  amil_sf_1.1_HiC
##  3      ASM1959352v1 ASM1959352v1_HiC
##  4        AGWG.draft         AaegL5.0
##  5            AaegL3           AaegL4
##  6        AGWG.draft         AaegL5.0
##                                                                new_assembly_link
##  1         https://dnazoo.s3.wasabisys.com/Acinonyx_jubatus/aciJub1_HiC.fasta.gz
##  2   https://dnazoo.s3.wasabisys.com/Acropora_millepora/amil_sf_1.1_HiC.fasta.gz
##  3 https://dnazoo.s3.wasabisys.com/Addax_nasomaculatus/ASM1959352v1_HiC.fasta.gz
##  4               https://dnazoo.s3.wasabisys.com/Aedes_aegypti/AaegL5.0.fasta.gz
##  5         https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL4/AaegL4.fasta.gz
##  6     https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL5.0/AaegL5.0.fasta.gz
##    new_assembly_link_status
##  1                      200
##  2                      200
##  3                      200
##  4                      404
##  5                      200
##  6                      200
##                                                                    hic_link
##  1    https://dnazoo.s3.wasabisys.com/Acinonyx_jubatus/aciJub1.rawchrom.hic
##  2   https://dnazoo.s3.wasabisys.com/Acropora_millepora/amil_sf_1.1_HiC.hic
##  3 https://dnazoo.s3.wasabisys.com/Addax_nasomaculatus/ASM1959352v1_HiC.hic
##  4                                                                     <NA>
##  5         https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL4/AaegL4.hic
##  6     https://dnazoo.s3.wasabisys.com/Aedes_aegypti__AaegL5.0/AaegL5.0.hic
```

For example, we can directly fetch a Hi-C dataset generated from a tardigrade sample by specifying the right `species` argument.

```
hicfile <- DNAZooData(species = 'Hypsibius_dujardini')

hicfile
##  HicFile object
##  .hic file: /home/biocbuild/.cache/R/DNAZooData/404991a1daed0_nHd_3.1_HiC.hic 
##  resolution: 5000 
##  pairs file: 
##  metadata(6): organism draftAssembly ... credits assemblyURL
```

Here again, the resulting `HicFile` is populated with metadata parsed from the DNA Zoo data portal.

```
metadata(hicfile)$organism
##  $vernacular
##  [1] "Tardigrade"
##  
##  $binomial
##  [1] "Hypsibius dujardini"
##  
##  $funFact
##  [1] "<i>Hypsibius dujardini</i> is a species of tardigrade, a tiny microscopic organism. They are also commonly called water bears. This species is found world-wide!"
##  
##  $extraInfo
##  [1] "on BioKIDS website"
##  
##  $extraInfoLink
##  [1] "http://www.biokids.umich.edu/critters/Hypsibius_dujardini/"
##  
##  $image
##  [1] "https://static.wixstatic.com/media/2b9330_82db39c219f24b20a75cb38943aad1fb~mv2.jpg"
##  
##  $imageCredit
##  [1] "By Willow Gabriel, Goldstein Lab - https://www.flickr.com/photos/waterbears/1614095719/ Template:Uploader Transferred from en.wikipedia to Commons., CC BY-SA 2.5, https://commons.wikimedia.org/w/index.php?curid=2261992"
##  
##  $isChromognomes
##  [1] "FALSE"
##  
##  $taxonomy
##  [1] "Species:202423-914154-914155-914158-155166-155362-710171-710179-710192-155390-155420"
```

`HiCFile` metadata also points to a URL to directly fetch the genome assembly corrected by the DNA Zoo consortium.

```
metadata(hicfile)$assemblyURL
##  [1] "https://dnazoo.s3.wasabisys.com/Hypsibius_dujardini/nHd_3.1_HiC.fasta.gz"
```

Back to top