---
source_url: https://satijalab.org/seurat/articles/cell_cycle_vignette
download_date: 2025-11-20
---

# Cell-Cycle Scoring and Regression

#### Compiled: 2023-10-31

Source: [`vignettes/cell_cycle_vignette.Rmd`](https://github.com/satijalab/seurat/blob/HEAD/vignettes/cell_cycle_vignette.Rmd)

`cell_cycle_vignette.Rmd`

We demonstrate how to mitigate the effects of cell cycle heterogeneity in scRNA-seq data by calculating cell cycle phase scores based on canonical markers, and regressing these out of the data during pre-processing. We demonstrate this on a dataset of murine hematopoietic progenitors ([Nestorowa *et al*., Blood 2016](http://www.bloodjournal.org/content/early/2016/06/30/blood-2016-05-716480?sso-checked=true)).You can download the files needed to run this vignette [here](https://www.dropbox.com/s/3dby3bjsaf5arrw/cell_cycle_vignette_files.zip?dl=1).

```
library(Seurat)

# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "/brahms/shared/vignette-data/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt",
    header = TRUE, as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create our Seurat object and complete the initalization steps
marrow <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(exp.mat), sparse = T))
marrow <- NormalizeData(marrow)
marrow <- FindVariableFeatures(marrow, selection.method = "vst")
marrow <- ScaleData(marrow, features = rownames(marrow))
```

If we run a PCA on our object, using the variable genes we found in `FindVariableFeatures()` above, we see that while most of the variance can be explained by lineage, PC8 and PC10 are split on cell-cycle genes including *TOP2A* and *MKI67*. We will attempt to regress this signal from the data, so that cell-cycle heterogeneity does not contribute to PCA or downstream analysis.

```
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 6:10, nfeatures.print = 10)
```

```
## PC_ 6 
## Positive:  SELL, ARL6IP1, CCL9, CD34, ADGRL4, BPIFC, NUSAP1, FAM64A, CD244, C030034L19RIK 
## Negative:  LY6C2, AA467197, CYBB, MGST2, ITGB2, PF4, CD74, ATP1B1, GP1BB, TREM3 
## PC_ 7 
## Positive:  F13A1, LY86, CFP, IRF8, CSF1R, TIFAB, IFI209, CCR2, TNS4, MS4A6C 
## Negative:  HDC, CPA3, PGLYRP1, MS4A3, NKG7, UBE2C, CCNB1, NUSAP1, PLK1, FUT8 
## PC_ 8 
## Positive:  NUSAP1, UBE2C, KIF23, PLK1, CENPF, FAM64A, CCNB1, H2AFX, ID2, CDC20 
## Negative:  WFDC17, SLC35D3, ADGRL4, VLDLR, CD33, H2AFY, P2RY14, IFI206, CCL9, CD34 
## PC_ 9 
## Positive:  IGKC, JCHAIN, LY6D, MZB1, CD74, IGLC2, FCRLA, IGKV4-50, IGHM, IGHV9-1 
## Negative:  SLC2A6, HBA-A1, HBA-A2, IGHV8-7, FCER1G, F13A1, HBB-BS, PLD4, HBB-BT, IGFBP4 
## PC_ 10 
## Positive:  H2AFX, FAM64A, ZFP383, NUSAP1, CDC25B, CENPF, GBP10, TOP2A, GBP6, GFRA1 
## Negative:  CTSW, XKRX, PRR5L, RORA, MBOAT4, A630014C17RIK, ZFP105, COL9A3, CLEC2I, TRAT1
```

```
DimHeatmap(marrow, dims = c(8, 10))
```

![](cell_cycle_vignette_files/figure-html/justification-1.png)

## Assign Cell-Cycle Scores

First, we assign each cell a score, based on its expression of G2/M and S phase markers. These marker sets should be anticorrelated in their expression levels, and cells expressing neither are likely not cycling and in G1 phase.

We assign scores in the `CellCycleScoring()` function, which stores S and G2/M scores in object meta data, along with the predicted classification of each cell in either G2M, S or G1 phase. `CellCycleScoring()` can also set the identity of the Seurat object to the cell-cycle phase by passing `set.ident = TRUE` (the original identities are stored as `old.ident`). Please note that Seurat does not use the discrete classifications (G2M/G1/S) in downstream cell cycle regression. Instead, it uses the quantitative scores for G2M and S phase. However, we provide our predicted classifications in case they are of interest.

```
marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(marrow[[]])
```

```
##          orig.ident nCount_RNA nFeature_RNA     S.Score  G2M.Score Phase
## Prog_013       Prog    2563089        10211 -0.14248691 -0.4680395    G1
## Prog_019       Prog    3030620         9991 -0.16915786  0.5851766   G2M
## Prog_031       Prog    1293487        10192 -0.34627038 -0.3971879    G1
## Prog_037       Prog    1357987         9599 -0.44270212  0.6820229   G2M
## Prog_008       Prog    4079891        10540  0.55854051  0.1284359     S
## Prog_014       Prog    2569783        10788  0.07116218  0.3166073   G2M
##          old.ident
## Prog_013      Prog
## Prog_019      Prog
## Prog_031      Prog
## Prog_037      Prog
## Prog_008      Prog
## Prog_014      Prog
```

```
# Visualize the distribution of cell cycle markers across
RidgePlot(marrow, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
```

![](cell_cycle_vignette_files/figure-html/cc_score-1.png)

```
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)
```

![](cell_cycle_vignette_files/figure-html/cc_score-2.png)

```
library(ggplot2)
plot <- DimPlot(marrow) + theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) +
    guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../output/images/cell_cycle_vignette.jpg", height = 7, width = 12, plot = plot,
    quality = 50)
```

We score single cells based on the scoring strategy described in [Tirosh *et al*. 2016](http://science.sciencemag.org/content/352/6282/189). See `?AddModuleScore()` in Seurat for more information, this function can be used to calculate supervised module scores for any gene list.

## Regress out cell cycle scores during data scaling

We now attempt to subtract (‘regress out’) this source of heterogeneity from the data. For users of Seurat v1.4, this was implemented in `RegressOut`. However, as the results of this procedure are stored in the scaled data slot (therefore overwriting the output of `ScaleData()`), we now merge this functionality into the `ScaleData()` function itself.

For each gene, Seurat models the relationship between gene expression and the S and G2M cell cycle scores. The scaled residuals of this model represent a ‘corrected’ expression matrix, that can be used downstream for dimensional reduction.

```
marrow <- ScaleData(marrow, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(marrow))
```

```
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)
```

```
## PC_ 1 
## Positive:  BLVRB, CAR2, KLF1, AQP1, CES2G, ERMAP, CAR1, FAM132A, RHD, SPHK1 
## Negative:  TMSB4X, H2AFY, CORO1A, PLAC8, EMB, MPO, PRTN3, CD34, LCP1, BC035044 
## PC_ 2 
## Positive:  ANGPT1, ADGRG1, MEIS1, ITGA2B, MPL, DAPP1, APOE, RAB37, GATA2, F2R 
## Negative:  LY6C2, ELANE, HP, IGSF6, ANXA3, CTSG, CLEC12A, TIFAB, SLPI, ALAS1 
## PC_ 3 
## Positive:  APOE, GATA2, NKG7, MUC13, MS4A3, RAB44, HDC, CPA3, FCGR3, TUBA8 
## Negative:  FLT3, DNTT, LSP1, WFDC17, MYL10, GIMAP6, LAX1, GPR171, TBXA2R, SATB1 
## PC_ 4 
## Positive:  CSRP3, ST8SIA6, DNTT, MPEG1, SCIN, LGALS1, CMAH, RGL1, APOE, MFSD2B 
## Negative:  PROCR, MPL, HLF, MMRN1, SERPINA3G, ESAM, GSTM1, D630039A03RIK, MYL10, LY6A 
## PC_ 5 
## Positive:  CPA3, LMO4, IKZF2, IFITM1, FUT8, MS4A2, SIGLECF, CSRP3, HDC, RAB44 
## Negative:  PF4, GP1BB, SDPR, F2RL2, RAB27B, SLC14A1, TREML1, PBX1, F2R, TUBA8
```

```
# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)
```

![](cell_cycle_vignette_files/figure-html/pca3-1.png)

As the best cell cycle markers are extremely well conserved across tissues and species, we have found this procedure to work robustly and reliably on diverse datasets.

## Alternate Workflow

The procedure above removes all signal associated with cell cycle. In some cases, we’ve found that this can negatively impact downstream analysis, particularly in differentiating processes (like murine hematopoiesis), where stem cells are quiescent and differentiated cells are proliferating (or vice versa). In this case, regressing out all cell cycle effects can blur the distinction between stem and progenitor cells as well.

As an alternative, we suggest regressing out the **difference** between the G2M and S phase scores. This means that signals separating non-cycling cells and cycling cells will be maintained, but differences in cell cycle phase among proliferating cells (which are often uninteresting), will be regressed out of the data

```
marrow$CC.Difference <- marrow$S.Score - marrow$G2M.Score
marrow <- ScaleData(marrow, vars.to.regress = "CC.Difference", features = rownames(marrow))
```

```
# cell cycle effects strongly mitigated in PCA
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)
```

```
## PC_ 1 
## Positive:  BLVRB, KLF1, ERMAP, FAM132A, CAR2, RHD, CES2G, SPHK1, AQP1, SLC38A5 
## Negative:  TMSB4X, CORO1A, PLAC8, H2AFY, LAPTM5, CD34, LCP1, TMEM176B, IGFBP4, EMB 
## PC_ 2 
## Positive:  APOE, GATA2, RAB37, ANGPT1, ADGRG1, MEIS1, MPL, F2R, PDZK1IP1, DAPP1 
## Negative:  CTSG, ELANE, LY6C2, HP, CLEC12A, ANXA3, IGSF6, TIFAB, SLPI, MPO 
## PC_ 3 
## Positive:  APOE, GATA2, NKG7, MUC13, ITGA2B, TUBA8, CPA3, RAB44, SLC18A2, CD9 
## Negative:  DNTT, FLT3, WFDC17, LSP1, MYL10, LAX1, GIMAP6, IGHM, CD24A, MN1 
## PC_ 4 
## Positive:  CSRP3, ST8SIA6, SCIN, LGALS1, APOE, ITGB7, MFSD2B, RGL1, DNTT, IGHV1-23 
## Negative:  MPL, MMRN1, PROCR, HLF, SERPINA3G, ESAM, PTGS1, D630039A03RIK, NDN, PPIC 
## PC_ 5 
## Positive:  GP1BB, PF4, SDPR, F2RL2, TREML1, RAB27B, SLC14A1, PBX1, PLEK, TUBA8 
## Negative:  HDC, LMO4, CSRP3, IFITM1, FCGR3, HLF, CPA3, PROCR, PGLYRP1, IKZF2
```

```
# when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1
# cells however, within actively proliferating cells, G2M and S phase cells group together
marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
DimPlot(marrow)
```

![](cell_cycle_vignette_files/figure-html/pca5-1.png)

**Session Info**

```
sessionInfo()
```

```
## R version 4.2.2 Patched (2022-11-10 r83330)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.6 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_3.4.4      Seurat_5.0.0       SeuratObject_5.0.0 sp_2.1-1          
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.7.1       Rtsne_0.16             colorspace_2.1-0      
##   [4] deldir_1.0-9           ellipsis_0.3.2         ggridges_0.5.4        
##   [7] rprojroot_2.0.3        RcppHNSW_0.5.0         fs_1.6.3              
##  [10] spatstat.data_3.0-3    farver_2.1.1           leiden_0.4.3          
##  [13] listenv_0.9.0          ggrepel_0.9.4          RSpectra_0.16-1       
##  [16] fansi_1.0.5            codetools_0.2-19       splines_4.2.2         
##  [19] cachem_1.0.8           knitr_1.45             polyclip_1.10-6       
##  [22] spam_2.10-0            jsonlite_1.8.7         ica_1.0-3             
##  [25] cluster_2.1.4          png_0.1-8              uwot_0.1.16           
##  [28] spatstat.sparse_3.0-3  sctransform_0.4.1      shiny_1.7.5.1         
##  [31] compiler_4.2.2         httr_1.4.7             Matrix_1.6-1.1        
##  [34] fastmap_1.1.1          lazyeval_0.2.2         cli_3.6.1             
##  [37] later_1.3.1            formatR_1.14           htmltools_0.5.6.1     
##  [40] tools_4.2.2            igraph_1.5.1           dotCall64_1.1-0       
##  [43] gtable_0.3.4           glue_1.6.2             reshape2_1.4.4        
##  [46] RANN_2.6.1             dplyr_1.1.3            Rcpp_1.0.11           
##  [49] scattermore_1.2        jquerylib_0.1.4        pkgdown_2.0.7         
##  [52] vctrs_0.6.4            nlme_3.1-162           spatstat.explore_3.2-5
##  [55] progressr_0.14.0       lmtest_0.9-40          spatstat.random_3.2-1 
##  [58] xfun_0.40              stringr_1.5.0          globals_0.16.2        
##  [61] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3       
##  [64] irlba_2.3.5.1          goftest_1.2-3          future_1.33.0         
##  [67] MASS_7.3-58.2          zoo_1.8-12             scales_1.2.1          
##  [70] spatstat.utils_3.0-4   ragg_1.2.5             promises_1.2.1        
##  [73] parallel_4.2.2         RColorBrewer_1.1-3     yaml_2.3.7            
##  [76] gridExtra_2.3          memoise_2.0.1          reticulate_1.34.0     
##  [79] pbapply_1.7-2          ggrastr_1.0.1          sass_0.4.7            
##  [82] stringi_1.7.12         highr_0.10             desc_1.4.2            
##  [85] fastDummies_1.7.3      rlang_1.1.1            pkgconfig_2.0.3       
##  [88] systemfonts_1.0.4      matrixStats_1.0.0      evaluate_0.22         
##  [91] lattice_0.21-9         tensor_1.5             ROCR_1.0-11           
##  [94] purrr_1.0.2            labeling_0.4.3         patchwork_1.1.3       
##  [97] htmlwidgets_1.6.2      cowplot_1.1.1          tidyselect_1.2.0      
## [100] parallelly_1.36.0      RcppAnnoy_0.0.21       plyr_1.8.9            
## [103] magrittr_2.0.3         R6_2.5.1               generics_0.1.3        
## [106] withr_2.5.2            pillar_1.9.0           fitdistrplus_1.1-11   
## [109] abind_1.4-5            survival_3.5-7         tibble_3.2.1          
## [112] future.apply_1.11.0    KernSmooth_2.23-22     utf8_1.2.4            
## [115] spatstat.geom_3.2-7    plotly_4.10.3          rmarkdown_2.25        
## [118] grid_4.2.2             data.table_1.14.8      digest_0.6.33         
## [121] xtable_1.8-4           tidyr_1.3.0            httpuv_1.6.12         
## [124] textshaping_0.3.6      munsell_0.5.0          beeswarm_0.4.0        
## [127] viridisLite_0.4.2      vipor_0.4.5            bslib_0.5.1
```