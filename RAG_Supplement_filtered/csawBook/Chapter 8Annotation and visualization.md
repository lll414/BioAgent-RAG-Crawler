Source URL: https://bioconductor.org/books/release/csawBook/annotation-and-visualization.html
Date Scraped: 2025-12-15

---

# Chapter 8 Annotation and visualization



## 8.1 Adding gene-based annotation

Annotation can be added to a given set of regions using the `detailRanges()` function.
This will identify overlaps between the regions and annotated genomic features such as exons, introns and promoters.
Here, the promoter region of each gene is defined as some interval 3 kbp up- and 1 kbp downstream of the TSS for that gene.
Any exonic features within `dist` on the left or right side of each supplied region will also be reported.

View set-up code

```
library(csaw)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

anno <- detailRanges(merged$regions, txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
    orgdb=org.Mm.eg.db, promoter=c(3000, 1000), dist=5000)
head(anno$overlap)
```

```
## [1] ""                    ""                    "Rrs1:+:P,Adhfe1:+:P"
## [4] "Ppp1r42:-:I"         ""                    "Ncoa2:-:PI"
```

```
head(anno$left)
```

```
## [1] ""               ""               ""               "Ppp1r42:-:3948"
## [5] ""               "Ncoa2:-:12"
```

```
head(anno$right)
```

```
## [1] ""                        "Rrs1:+:3898"            
## [3] "Rrs1:+:48,Adhfe1:+:2588" "Ppp1r42:-:1612"         
## [5] "Ncoa2:-:4595"            "Ncoa2:-:1278"
```

Character vectors of compact string representations are provided to summarize the features overlapped by each supplied region.
Each pattern contains `GENE|STRAND|TYPE` to describe the strand and overlapped features of that gene.
Exons are labelled as `E`, promoters are `P` and introns are `I`.
For `left` and `right`, `TYPE` is replaced by `DISTANCE`.
This indicates the gap (in base pairs) between the supplied region and the closest non-overlapping exon of the annotated feature.
All of this annotation can be stored in the metadata of the `GRanges` object for later use.

```
merged$regions$overlap <- anno$overlap
merged$regions$left <- anno$left
merged$regions$right <- anno$right
```

While the string representation saves space in the output, it is not easy to work with.
If the annotation needs to manipulated directly, users can obtain it from the `detailRanges()` command by not specifying the regions of interest.
This can then be used for interactive manipulation, e.g., to identify all genes where the promoter contains DB sites.

```
anno.ranges <- detailRanges(txdb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
    orgdb=org.Mm.eg.db)
anno.ranges
```

```
## GRanges object with 505738 ranges and 2 metadata columns:
##             seqnames              ranges strand |      symbol        type
##                <Rle>           <IRanges>  <Rle> | <character> <character>
##   100009600     chr9   21062393-21062717      - |       Zglp1           E
##   100009600     chr9   21062400-21062717      - |       Zglp1           E
##   100009600     chr9   21062894-21062987      - |       Zglp1           E
##   100009600     chr9   21063314-21063396      - |       Zglp1           E
##   100009600     chr9   21066024-21066377      - |       Zglp1           E
##         ...      ...                 ...    ... .         ...         ...
##       99982     chr4 136554325-136558324      - |       Kdm1a           P
##       99982     chr4 136560502-136564501      - |       Kdm1a           P
##       99982     chr4 136567608-136571607      - |       Kdm1a           P
##       99982     chr4 136576159-136580158      - |       Kdm1a           P
##       99982     chr4 136550540-136602723      - |       Kdm1a           G
##   -------
##   seqinfo: 66 sequences (1 circular) from mm10 genome
```

## 8.2 Checking bimodality for TF studies

For TF experiments, a simple measure of strand bimodality can be reported as a diagnostic.
Given a set of regions, the `checkBimodality()` function will return the maximum bimodality score across all base positions in each region.
The bimodality score at each base position is defined as the minimum of the ratio of the number of forward- to reverse-stranded reads to the left of that position, and the ratio of the reverse- to forward-stranded reads to the right.
A high score is only possible if both ratios are large, i.e., strand bimodality is present.

```
# TODO: make this less weird.
spacing <- metadata(data)$spacing
expanded <- resize(merged$regions, fix="center", 
    width=width(merged$regions)+spacing)
sbm.score <- checkBimodality(bam.files, expanded, width=frag.len)
head(sbm.score)
```

```
## [1] 1.456 2.263 1.364 5.333 1.933 3.405
```

In the above code, all regions are expanded by `spacing`, i.e., 50 bp.
This ensures that the optimal bimodality score can be computed for the centre of the binding site, even if that position is not captured by a window.
The `width` argument specifies the span with which to count reads for the score calculation.
This should be set to the average fragment length.
If multiple `bam.files` are provided, they will be pooled during counting.

For typical TF binding sites, bimodality scores can be considered to be “high” if they are larger than 4.
This allows users to distinguish between genuine binding sites and high-abundance artifacts such as repeats or read stacks.
However, caution is still required as some high scores may be driven by the stochastic distribution of reads.
Obviously, the concept of strand bimodality is less relevant for diffuse targets like histone marks.

## 8.3 Saving the results to file

It is a simple matter to save the results for later perusal, e.g., to a tab-separated file.

```
ofile <- gzfile("clusters.tsv.gz", open="w")
write.table(as.data.frame(merged), file=ofile, 
    row.names=FALSE, quote=FALSE, sep="\t")
close(ofile)
```

Of course, other formats can be used depending on the purpose of the file.
For example, significantly DB regions can be exported to BED files through the *[rtracklayer](https://bioconductor.org/packages/3.22/rtracklayer)* package for visual inspection with genomic browsers.
A transformed FDR is used here for the score field.

```
is.sig <- merged$combined$FDR <= 0.05
test <- merged$regions[is.sig]
test$score <- -10*log10(merged$combined$FDR[is.sig])
names(test) <- paste0("region", 1:sum(is.sig))

library(rtracklayer)
export(test, "clusters.bed")
head(read.table("clusters.bed"))
```

```
##     V1       V2       V3      V4    V5 V6
## 1 chr1  7397900  7398110 region1 14.63  .
## 2 chr1  9541400  9541510 region2 16.60  .
## 3 chr1  9545300  9545360 region3 13.35  .
## 4 chr1 13589950 13590010 region4 15.36  .
## 5 chr1 15805500 15805660 region5 24.04  .
## 6 chr1 16657100 16657160 region6 13.60  .
```

Alternatively, the `GRanges` object can be directly saved to file and reloaded later for direct manipulation in the R environment, e.g., to find overlaps with other regions of interest.

```
saveRDS(merged$regions, "ranges.rds")
```

## 8.4 Simple visualization of genomic coverage

Visualization of the read depth around interesting features is often desired.
This is facilitated by the `extractReads()` function, which pulls out the reads from the BAM file.
The returned `GRanges` object can then be used to plot the sequencing coverage or any other statistic of interest.
Note that the `extractReads()` function also accepts a `readParam` object.
This ensures that the same reads used in the analysis will be pulled out during visualization.

```
cur.region <- GRanges("chr18", IRanges(77806807, 77807165))
extractReads(bam.files[[1]], cur.region, param=param)
```

```
## GRanges object with 55 ranges and 0 metadata columns:
##        seqnames            ranges strand
##           <Rle>         <IRanges>  <Rle>
##    [1]    chr18 77806886-77806922      +
##    [2]    chr18 77806887-77806923      +
##    [3]    chr18 77806887-77806923      +
##    [4]    chr18 77806887-77806923      +
##    [5]    chr18 77806890-77806926      +
##    ...      ...               ...    ...
##   [51]    chr18 77807063-77807095      -
##   [52]    chr18 77807068-77807104      -
##   [53]    chr18 77807082-77807119      -
##   [54]    chr18 77807084-77807120      -
##   [55]    chr18 77807087-77807123      -
##   -------
##   seqinfo: 1 sequence from an unspecified genome
```

Here, coverage is visualized as the number of reads covering each base pair in the interval of interest.
Specifically, the reads-per-million is shown to allow comparisons between libraries of different size.
The plots themselves are constructed using methods from the *[Gviz](https://bioconductor.org/packages/3.22/Gviz)* package.
The blue and red tracks represent the coverage on the forward and reverse strands, respectively.
Strong strand bimodality is consistent with a genuine TF binding site.
For paired-end data, coverage can be similarly plotted for fragments, i.e., proper read pairs.

```
library(Gviz)
collected <- vector("list", length(bam.files))
for (i in seq_along(bam.files)) { 
    reads <- extractReads(bam.files[[i]], cur.region, param=param)
    adj.total <- data$totals[i]/1e6
    pcov <- as(coverage(reads[strand(reads)=="+"])/adj.total, "GRanges")
    ncov <- as(coverage(reads[strand(reads)=="-"])/adj.total, "GRanges")
    ptrack <- DataTrack(pcov, type="histogram", lwd=0, fill=rgb(0,0,1,.4), 
        ylim=c(0,1.1), name=tf.data$Name[i], col.axis="black", 
        col.title="black")
    ntrack <- DataTrack(ncov, type="histogram", lwd=0, fill=rgb(1,0,0,.4), 
        ylim=c(0,1.1))
    collected[[i]] <- OverlayTrack(trackList=list(ptrack,ntrack))
}
gax <- GenomeAxisTrack(col="black")
plotTracks(c(gax, collected), from=start(cur.region), to=end(cur.region))
```

![](postprocess_files/figure-html/nfya-regionplot-1.png)

## Session information

View session info