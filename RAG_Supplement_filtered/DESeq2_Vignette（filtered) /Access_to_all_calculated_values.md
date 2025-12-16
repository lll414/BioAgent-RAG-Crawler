# Access to all calculated values

Source URL: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
Date Scraped: 2025-12-11

---

All row-wise calculated values (intermediate dispersion calculations, coefficients, standard errors, etc.) are stored in the *DESeqDataSet* object, e.g.Â `dds` in this vignette. These values are accessible by calling *mcols* on `dds`. Descriptions of the columns are accessible by two calls to *mcols*. Note that the call to `substr` below is only for display purposes.

```
mcols(dds,use.names=TRUE)[1:4,1:4]
```

```
## DataFrame with 4 rows and 4 columns
##                    gene   baseMean     baseVar   allZero
##             <character>  <numeric>   <numeric> <logical>
## FBgn0000008 FBgn0000008   95.28865 2.29337e+02     FALSE
## FBgn0000017 FBgn0000017 4359.09632 3.71585e+05     FALSE
## FBgn0000018 FBgn0000018  419.06811 2.24013e+03     FALSE
## FBgn0000024 FBgn0000024    6.41105 3.74104e+00     FALSE
```

```
substr(names(mcols(dds)),1,10)
```

```
##  [1] "gene"       "baseMean"   "baseVar"    "allZero"    "dispGeneEs"
##  [6] "dispGeneIt" "dispFit"    "dispersion" "dispIter"   "dispOutlie"
## [11] "dispMAP"    "Intercept"  "condition_" "SE_Interce" "SE_conditi"
## [16] "WaldStatis" "WaldStatis" "WaldPvalue" "WaldPvalue" "betaConv"  
## [21] "betaIter"   "deviance"   "maxCooks"
```

```
mcols(mcols(dds), use.names=TRUE)[1:4,]
```

```
## DataFrame with 4 rows and 2 columns
##                  type            description
##           <character>            <character>
## gene                                        
## baseMean intermediate mean of normalized c..
## baseVar  intermediate variance of normaliz..
## allZero  intermediate all counts for a gen..
```

The mean values \(\mu\_{ij} = s\_j q\_{ij}\) and the Cookâs distances for each gene and sample are stored as matrices in the assays slot:

```
head(assays(dds)[["mu"]])
```

```
##               treated1    treated2    treated3 untreated1 untreated2
## FBgn0000008  154.18297   72.011890   78.553988  107.92325  169.00095
## FBgn0000017 6441.09688 3008.344869 3281.645427 5333.50912 8351.93629
## FBgn0000018  657.53926  307.106832  335.006715  495.29423  775.59928
## FBgn0000024   11.43421    5.340404    5.825567    6.91794   10.83305
## FBgn0000032 1559.80292  728.513359  794.696964 1164.47564 1823.49484
## FBgn0000037   27.23930   12.722244   13.878028   13.84125   21.67451
##              untreated3  untreated4
## FBgn0000008   60.947217   70.828452
## FBgn0000017 3011.978798 3500.304091
## FBgn0000018  279.706228  325.054365
## FBgn0000024    3.906750    4.540143
## FBgn0000032  657.611313  764.228344
## FBgn0000037    7.816532    9.083808
```

```
head(assays(dds)[["cooks"]])
```

```
##               treated1    treated2    treated3 untreated1   untreated2
## FBgn0000008 0.08564565 0.297964299 0.077021521 0.10556148 0.0135846198
## FBgn0000017 0.01275474 0.004120273 0.002353591 0.08760807 0.0105562769
## FBgn0000018 0.09813824 0.005595703 0.054059461 0.17277269 0.0021395029
## FBgn0000024 0.06482596 0.129406896 0.030982961 0.26445639 0.0005741831
## FBgn0000032 0.07625620 0.017481712 0.020056352 0.32377428 0.0211262017
## FBgn0000037 0.45819676 0.026766117 0.151000426 0.01530776 0.0859061169
##             untreated3   untreated4
## FBgn0000008 0.19761575 0.0004998745
## FBgn0000017 0.18275047 0.0548693390
## FBgn0000018 0.07090752 0.0104402584
## FBgn0000024 0.02968291 0.0809265573
## FBgn0000032 0.02151283 0.0764805449
## FBgn0000037 0.02257542 0.2489127110
```

The dispersions \(\alpha\_i\) can be accessed with the *dispersions* function.

```
head(dispersions(dds))
```

```
## [1] 0.03082228 0.01305447 0.01518649 0.23983770 0.01654718 0.12963629
```

```
head(mcols(dds)$dispersion)
```

```
## [1] 0.03082228 0.01305447 0.01518649 0.23983770 0.01654718 0.12963629
```

The size factors \(s\_j\) are accessible via *sizeFactors*:

```
sizeFactors(dds)
```

```
##   treated1   treated2   treated3 untreated1 untreated2 untreated3 untreated4 
##  1.6297067  0.7611622  0.8303119  1.1439041  1.7912811  0.6459940  0.7507275
```

For advanced users, we also include a convenience function *coef* for extracting the matrix \([\beta\_{ir}]\) for all genes *i* and model coefficients \(r\). This function can also return a matrix of standard errors, see `?coef`. The columns of this matrix correspond to the effects returned by *resultsNames*. Note that the *results* function is best for building results tables with *p* values and adjusted *p* values.

```
head(coef(dds))
```

```
##             Intercept condition_treated_vs_untreated
## FBgn0000008  6.559896                     0.00399148
## FBgn0000017 12.186903                    -0.23842494
## FBgn0000018  8.758176                    -0.10185506
## FBgn0000024  2.596376                     0.21429657
## FBgn0000032  9.991499                    -0.08896298
## FBgn0000037  3.596936                     0.46606939
```

The beta prior variance \(\sigma\_r^2\) is stored as an attribute of the *DESeqDataSet*:

```
attr(dds, "betaPriorVar")
```

```
## [1] 1e+06 1e+06
```

General information about the prior used for log fold change shrinkage is also stored in a slot of the *DESeqResults* object. This would also contain information about what other packages were used for log2 fold change shrinkage.

```
priorInfo(resLFC)
```

```
## $type
## [1] "apeglm"
## 
## $package
## [1] "apeglm"
## 
## $version
## [1] '1.32.0'
## 
## $prior.control
## $prior.control$no.shrink
## [1] 1
## 
## $prior.control$prior.mean
## [1] 0
## 
## $prior.control$prior.scale
## [1] 0.2031588
## 
## $prior.control$prior.df
## [1] 1
## 
## $prior.control$prior.no.shrink.mean
## [1] 0
## 
## $prior.control$prior.no.shrink.scale
## [1] 15
## 
## $prior.control$prior.var
## [1] 0.0412735
```

```
priorInfo(resNorm)
```

```
## $type
## [1] "normal"
## 
## $package
## [1] "DESeq2"
## 
## $version
## [1] '1.50.2'
## 
## $betaPriorVar
##        Intercept conditiontreated 
##     1.000000e+06     1.049206e-01
```

```
priorInfo(resAsh)
```

```
## $type
## [1] "ashr"
## 
## $package
## [1] "ashr"
## 
## $version
## [1] '2.2.63'
## 
## $fitted_g
## $pi
##  [1] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.128392406
##  [7] 0.000000000 0.452702361 0.057818455 0.000000000 0.176425328 0.000000000
## [13] 0.120967208 0.027145716 0.005428328 0.018807790 0.012312407 0.000000000
## [19] 0.000000000 0.000000000 0.000000000 0.000000000
## 
## $mean
##  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
## 
## $sd
##  [1] 0.006752731 0.009549803 0.013505461 0.019099607 0.027010923 0.038199213
##  [7] 0.054021845 0.076398426 0.108043690 0.152796852 0.216087381 0.305593704
## [13] 0.432174761 0.611187408 0.864349522 1.222374817 1.728699044 2.444749633
## [19] 3.457398088 4.889499267 6.914796176 9.778998533
## 
## attr(,"class")
## [1] "normalmix"
## attr(,"row.names")
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
```

The dispersion prior variance \(\sigma\_d^2\) is stored as an attribute of the dispersion function:

```
dispersionFunction(dds)
```

```
## function (q) 
## coefs[1] + coefs[2]/q
## <bytecode: 0x60ae19294228>
## <environment: 0x60ae19293bd0>
## attr(,"coefficients")
## asymptDisp  extraPois 
## 0.01377513 2.88413302 
## attr(,"fitType")
## [1] "parametric"
## attr(,"varLogDispEsts")
## [1] 1.025797
## attr(,"dispPriorVar")
## [1] 0.5354389
```

```
attr(dispersionFunction(dds), "dispPriorVar")
```

```
## [1] 0.5354389
```

The version of DESeq2 which was used to construct the *DESeqDataSet* object, or the version used when *DESeq* was run, is stored here:

```
metadata(dds)[["version"]]
```

```
## [1] '1.50.2'
```