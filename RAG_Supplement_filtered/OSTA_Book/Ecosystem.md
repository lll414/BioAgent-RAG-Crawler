Source URL: https://bioconductor.org/books/release/OSTA/pages/bkg-ecosystem.html
Date Scraped: 2025-12-15

---

1. [Background](https://bioconductor.org/books/release/OSTA/pages/bkg-introduction.html)
2. [4  Ecosystem](https://bioconductor.org/books/release/OSTA/pages/bkg-ecosystem.html)

# 4  Ecosystem

## 4.1 Introduction

This chapter provides additional background about the Bioconductor ecosystem by exploring the space of packages that revolve around spatial and single cell omics data. We deliberately include the latter, since many tools developed for single cell data are either directly applicable to, or lay the foundation for tools developed for spatial data.

Here, we will rely on *[biocViews](https://bioconductor.org/packages/3.22/biocViews)* ([Carey et al. 2025](https://bioconductor.org/books/release/OSTA/pages/bkg-ecosystem.html#ref-Carey2006-BiocViews)) that are provided in the `DESCRIPTION` of each package. These may be browsed online – in a non-programmatic manner – on the [Bioconductor BiocViews website](https://www.bioconductor.org/packages/release/BiocViews.html).

Note that the data reported herein may be inaccurate, since packages might provide inaccurate/incomplete `biocViews`. Especially the “Spatial” term relied upon here was added fairly recently and might be missing from older packages.

### 4.1.1 Exploratory

The *[BiocPkgTools](https://bioconductor.org/packages/3.22/BiocPkgTools)* package ([Su et al. 2025](https://bioconductor.org/books/release/OSTA/pages/bkg-ecosystem.html#ref-Su2018-BiocPkgTools)) provides tools to access and explore Bioconductor package metadata in R, including their details (e.g., authors) and (monthly) download statistics.

`biocExplore()` provides an interactive way to explore the package space:

Code

```
library(BiocPkgTools)
biocExplore()
```

GenomeInfoDb: 10357461GenoBiocGenerics: 11569707BiocIRanges: 11288318IRanS4Vectors: 10774534S4VeBiocVersion: 9593759BiocXVector: 9053660XVecBiostrings: 9035741BiosBiobase: 10861435BiobGenomicRanges: 8932856GenoSparseArray: 2216077SpaDelayedArray: 7477514DelaSummarizedExperiment: 7386091SummMatrixGenerics: 5162605MatBiocParallel: 8375546BiocS4Arrays: 2419023S4AUCSC.utils: 1311343UCAnnotationDbi: 8100683AnnoKEGGREST: 4421864KEGlimma: 6739471limmfgsea: 2101697fgsrhdf5: 3109237rhdRsamtools: 5240406RsaRhtslib: 3700551RhtGenomicAlignments: 4755870GenRhdf5lib: 2961881RhdBiocFileCache: 3026853Bioenrichplot: 2064513enredgeR: 4102339edgrtracklayer: 4967375rtrDESeq2: 3965233DESclusterProfiler: 2233888clubiomaRt: 4968221biorhdf5filters: 2018274rhdgraph: 3350321graggtree: 2065590ggtDOSE: 2069657DOStreeio: 1764286trSingleCellExperiment: 1892409SinGenomicFeatures: 4164501GenDelayedMatrixStats: 1844768DelBiocIO: 1820840Bioannotate: 4605794annGOSemSim: 1905707GOSqvalue: 2265568qvabeachmat: 1655543beComplexHeatmap: 1672635CoHDF5Array: 1616815HDsparseMatrixStats: 1499572spBSgenome: 2689013BSgpreprocessCore: 2759235preBiocSingular: 1273221BiScaledMatrix: 1044827Scgenefilter: 4066778genBiocNeighbors: 1120182Bimulttest: 2553572mulAnnotationHub: 1414622AnAnnotationFilter: 1273597AnProtGenerics: 1454596Prscuttle: 1002488scensembldb: 1494924enimpute: 1902997impassorthead: 256441aGEOquery: 2262560GEORBGL: 1850214RBGVariantAnnotation: 1998626VarGSVA: 688399GSGSEABase: 1325119GSGenomicDataCommons: 221403GTCGAutils: 177762Taffy: 2189145affscater: 993270scaffyio: 1860286affRgraphviz: 1858267RgrExperimentHub: 715110Exbiomformat: 806424biphyloseq: 1029885phpwalign: 208787pglmSparseNet: 105965gBiocBaseUtils: 343742BShortRead: 1921814Shosva: 1113017svbluster: 530917blMultiAssayExperiment: 476963MuglmGamPoi: 348851gDirichletMultinomial: 469682DibiovizBase: 1292797biBiocStyle: 997026BiSpatialExperiment: 277450Sscran: 641254scEnhancedVolcano: 461988EnpcaMethods: 803945pcmetapod: 442230meTFBSTools: 451191TFvsn: 1122315vsgeneplotter: 2727194genbatchelor: 413280batxdbmaker: 149153tbiocViews: 594348biDNAcopy: 756361DNTCGAbiolinks: 541007TCSingleR: 363809SiCNEr: 330236CKEGGgraph: 727314KEseqLogo: 483334seResidualMatrix: 309562Rapeglm: 428270appathview: 664158pamonocle: 451852momzR: 1123894mzMSnbase: 558746MStximport: 517158txcytolib: 349580cGviz: 991657GvMsCoreUtils: 314184MinteractiveDisplayBase: 991096inalabaster.base: 118097aalabaster.matrix: 110379adada2: 868164daflowCore: 4921345floEBImage: 628906EBmixOmics: 332604mdir.expiry: 248827dChIPseeker: 337185Cbasilisk: 270058balabaster.ranges: 103794agraphite: 360797gralabaster.se: 102492aReactomePA: 567396RemzID: 403608mzbasilisk.utils: 241025bregioneR: 307787rChIPpeakAnno: 319107Calabaster.schemas: 91999aQFeatures: 134865QDECIPHER: 332807DAUCell: 225931Abumphunter: 533516bugypsum: 86933gRProtoBufLib: 309178RConsensusClusterPlus: 420936CoDropletUtils: 301073DtopGO: 576994tominfi: 917068mimaftools: 314725msiggenes: 570613siMassSpecWavelet: 295972MRsubread: 401045Rsilluminaio: 561534ilOrganismDbi: 554956OrTreeSummarizedExperiment: 141845TMAST: 251167MTrajectoryUtils: 157337Tggbio: 542583gggdsfmt: 309277gsnpStats: 474323snmetagenomeSeq: 225227mPSMatch: 78128Pxcms: 951930xcBiocCheck: 284412Bslingshot: 175840suniversalmotif: 71207uscDblFinder: 144253schromVAR: 148594cmsa: 264140mmarray: 632166mamotifmatchr: 155257mzellkonverter: 116802zAnnotationForge: 570241AnDiffBind: 292889Daroma.light: 403214armicrobiome: 211945mgoseq: 290251gSTRINGdb: 166965SInteractionSet: 187360IANCOMBC: 113661Ainfercnv: 149204iSNPRelate: 254294SsystemPipeR: 265814sDEXSeq: 315801DEDASeq: 388919EDoligo: 489932oldecontam: 127292dUCell: 78350Umia: 166888mlpsymphony: 184505lflowWorkspace: 339626fMetaboCoreUtils: 100103Maffxparser: 391929afWrench: 142838Wmethylumi: 375047mesingscore: 114842stkWidgets: 257627tALDEx2: 132117ASpectra: 129706SMfuzz: 169703MwidgetTools: 260179wCategory: 584399CadecoupleR: 71484dsimplifyEnrichment: 60516svariancePartition: 114193vncdfFlow: 287537nplyranges: 137306pIHW: 143149IoligoClasses: 424133olDynDoc: 448705DydittoSeq: 105691dGreyListChIP: 118855Gbsseq: 272493bGOstats: 569881GOGENIE3: 121581Gtximeta: 161397tglobaltest: 285839gMsFeatures: 126555MNebulosa: 74960NFlowSOM: 183499FMsExperiment: 58042Mlumi: 394190lusurvcomp: 175559sSeqArray: 137591SDMRcate: 188337DRcisTarget: 121026Rropls: 128283rflowViz: 347517fPCAtools: 118111Pgcrma: 673069gcROC: 299780RTCseq: 37520TMultiDataSet: 106609Mcelda: 69121cggtreeExtra: 79130gmotifStack: 141572mChemmineR: 521116Chgenomation: 115686gGlobalAncova: 139670GmissMethyl: 174482mOmnipathR: 55926Oggcyto: 173090gRUVSeq: 191193REnrichedHeatmap: 85371Ebamsignals: 142461balabaster.sce: 36325agage: 219827gGenomicScores: 83998GseqPattern: 93898saffyPLM: 618794afMaaslin2: 75421Mdensvis: 116222dzinbwave: 125229zGenomicFiles: 199268GGlimma: 178713Gfastseg: 141339fwateRmelon: 173722wflowClust: 217009fReportingTools: 229025RmethylKit: 141603mrGREAT: 74864rkaryoploteR: 126427kMotifDb: 134414MopenCyto: 195918ochipseq: 141308ctrackViewer: 79564tGeneOverlap: 64245GtradeSeq: 63176tbaySeq: 150681bRCy3: 112854RarrayQualityMetrics: 208329aMSstats: 89627MATACseqQC: 62680Aquantsmooth: 136410qrGADEM: 105397rRaggedExperiment: 106752RsangerseqR: 90762sscRepertoire: 44687sdestiny: 177985dminet: 154292mLEA: 78101LEpiDISH: 51695EAnVIL: 55168AGWASTools: 132881GCytoML: 133333Csesame: 81361slefser: 57630lDSS: 189835Dbeadarray: 226469bviper: 92537vflowAI: 72339fMungeSumstats: 47500MRdisop: 132586RCATALYST: 78524CChAMP: 151895CEBSeq: 110218EarrayQuality: 73468aDEP: 75194DflowStats: 230669fR4RNA: 52403RDEGreport: 84179Dballgown: 198143bgwascat: 88892gmiloR: 38087mcummeRbund: 591293cuderfinder: 95343dTOAST: 32910Tquantiseqr: 31713qBeadDataPackR: 173157Bmbkmeans: 64005mChemmineOB: 236068Cprogeny: 48425pNOISeq: 127932NUniProt.ws: 77028Umuscat: 46522mclusterExperiment: 74476cRTCGA: 87505RderfinderHelper: 84499dBioNERO: 22551BMSstatsConvert: 39479MsoGGi: 41362sChIPQC: 74374CLoomExperiment: 46718LCGHbase: 98409Cggmsa: 41741gSPIA: 129935SInteractiveComplexHeatmap: 35155IHeatplus: 149553HRnBeads: 75739RArrayExpress: 149728AiClusterPlus: 66118iheatmaps: 37905hMOFA2: 58087MTSCAN: 64734TRbcBook1: 46206RcmapR: 40392cQuasR: 109018QBumpyMatrix: 47063Bsafe: 104821sdiffcyt: 47869dM3C: 82283MRbowtie2: 43300RRTCGAToolbox: 88059Raffycoretools: 136311aannotatr: 84319aMLInterfaces: 119285MbioDist: 101507bbambu: 29459bCAMERA: 135255Csimona: 16293Organism.dplyr: 54911OAnnotationHubData: 50365Acsaw: 100909cSeqVarTools: 89094SHMMcopy: 79907HRbowtie: 87736RGeomxTools: 31588Gpiano: 98895pscds: 49539srWikiPathways: 74794rDeconRNASeq: 51150DM3Drop: 66694Mescape: 33178eAIMS: 66670AGenomicInteractions: 64771GQDNAseq: 73547Qaffycomp: 68051aIcens: 138113IBioNet: 79839BDRIMSeq: 59628DCGHcall: 91718Cscmap: 45921shpar: 64230htracktables: 27593tpowerTCR: 42055pGLAD: 96164Gannaffy: 391238anbiobroom: 45651bBayesSpace: 28193BHybridMTest: 33967HBiocSet: 30697BCelliD: 20438MutationalPatterns: 63688Mqusage: 64900qEnrichmentBrowser: 70980ESC3: 86772Sggkegg: 18316iSEE: 57725imsmsTests: 48069mmuscle: 58214mrrvgo: 38823rmade4: 101321mscde: 72092sRankProd: 122914Rsplatter: 70307sAgiMicroRna: 57626AENmix: 45211EMicrobiotaProcess: 36289Mmosaics: 273178morthogene: 25321oROTS: 37571Rfishpond: 46993fcytomapper: 31820cSGSeq: 57163Sctc: 101168cdmrseq: 31815dadSplit: 46559aGENESIS: 79409Glfa: 52836lfmcsR: 64024fscry: 24829spcaExplorer: 64720pa4Core: 51608aannotationTools: 56713aGOfuncR: 35126GHIBAG: 28848Hquantro: 45021qSRAdb: 127452Srecount: 70416rmakecdfenv: 87841mASICS: 20599satuRn: 19667Rqc: 44930RmaSigPro: 88126mbiocGraph: 61328bACME: 51884AstageR: 31146sExperimentHubData: 39379EBaalChIP: 19912BicARE: 50311BHiCBricks: 14930EBarrays: 65926EAGDEX: 35455AapComplex: 50762adeepSNV: 45007dplgem: 43840pzFPKM: 23037ztilingArray: 58294thopach: 78338hYAPSA: 26700YGenVisR: 62255Gcicero: 36355chypergraph: 59764hmygene: 74714mCoreGx: 38113Crols: 71816rmsmsEDA: 48265mABSSeq: 32844AaffylmGUI: 75686agtrellis: 32397ga4: 48193ainteractiveDisplay: 50863iHiTC: 51918HMVCClass: 39201MBiRewire: 35835Bbiocthis: 29423bLinnorm: 33177Lbacon: 28445bBioMVCClass: 41816BBiocPkgTools: 25242Bcancerclass: 34961cGuitar: 28621Ganimalcules: 19435BgeeDB: 25952BRBioFormats: 13076Repitools: 84215RDEsingle: 34884DflowDensity: 49918fscDD: 30334srsbml: 62135raCGH: 75662aBaseSpaceR: 31433Btricycle: 18460NanoStringNCTools: 28898NRedeR: 61288RsingleCellTK: 34134sa4Classif: 45963aHiCcompare: 33134HiCOBRA: 34336ia4Base: 49235aABarray: 51872AADaCGH2: 38667ASpatialFeatureExperiment: 14693StructuralVariantAnnotation: 28748SBiocWorkflowTools: 27346BmiaViz: 21047MMUPHin: 17170BAGS: 27285Bcn.mops: 65707cOmicCircos: 53836Oparody: 47265pASSET: 30452APharmacoGx: 44245Panota2seq: 21558scMerge: 44268scleaver: 44761ccqn: 67557cmicroRNA: 90510mMsBackendMgf: 39220MsigFeature: 22575smdqc: 41264maffyContam: 45923aMethylAid: 26931MbiomvRCNS: 28998byarn: 22557y

Filter:

NoneATACSeqAlignmentAlternativeSplicingAnnotationAssayDomainBatchEffectBayesianBioinformaticsBiologicalQuestionBiomedicalInformaticsCRISPRCellBasedAssaysCellBiologyChIPSeqCheminformaticsChipOnChipClassificationClusteringComparativeGenomicsCopyNumberVariationCoverageCpGIslandDNA3DStructureDNAMethylationDNASeqDataImportDataRepresentationDifferentialExpressionDifferentialMethylationDifferentialPeakCallingDifferentialSplicingDimensionReductionEpigeneticsExonArrayExperimentDataExperimentHubSoftwareExperimentalDesignFeatureExtractionFlowCytometryFunctionalGenomicsFunctionalPredictionGOGUIGeneExpressionGenePredictionGeneRegulationGeneSetEnrichmentGeneSignalingGeneTargetGeneticVariabilityGeneticsGenomeAnnotationGenomeAssemblyGenomeWideAssociationGenomicVariationGraphAndNetworkHiCHiddenMarkovModelImmunoOncologyInfrastructureKEGGLipidomicsMassSpectrometryMetabolomicsMetagenomicsMethylSeqMethylationArrayMicroarrayMicrobiomeMicrotitrePlateAssayMotifAnnotationMotifDiscoveryMultipleComparisonMultipleSequenceAlignmentNetworkNetworkEnrichmentNetworkInferenceNormalizationOneChannelPathwaysPeakDetectionPharmacogeneticsPharmacogenomicsPhylogeneticsPreprocessingPrincipalComponentProprietaryPlatformsProteomicsQualityControlRNASeqReactomeRegressionReportWritingResearchFieldSNPSequenceMatchingSequencingShinyAppsSingleCellSoftwareSomaticMutationSpatialStatisticalMethodStructuralPredictionSupportVectorMachineSurvivalSystemsBiologyTechnologyThirdPartyClientTimeCourseTranscriptionTranscriptomicsTwoChannelVariantAnnotationVariantDetectionVisualizationWholeGenomeWorkflowStepmRNAMicroarraymiRNAqPCR

`biocPkgList()` retrieves the full Bioconductor software package listing (at the time of calling the function), including associated metadata. These include `biocViews`, which we can use to identify packages of interest:

Code

```
# retrieve current package record
df <- biocPkgList()
# helper function to get indices of packages
# that contain 'biocViews' specified by 'x'
.f <- \(x, y=df) vapply(y$biocViews, \(.) all(x %in% .), logical(1))
# view "Spatial" packages
df$Package[.f("Spatial")]
```

```
##   [1] "alabaster.sfe"            "Banksy"                   "BatchSVG"                
##   [4] "betaHMM"                  "BulkSignalR"              "CARDspa"                 
##   [7] "CatsCradle"               "clustSIGNAL"              "concordexR"              
##  [10] "CTSV"                     "cytoviewer"               "DESpace"                 
##  [13] "escheR"                   "FuseSOM"                  "GeomxTools"              
##  [16] "ggsc"                     "ggspavis"                 "HiCPotts"                
##  [19] "hoodscanR"                "HuBMAPR"                  "imcRtools"               
##  [22] "jazzPanda"                "knowYourCG"               "lisaClust"               
##  [25] "miRspongeR"               "mistyR"                   "mitology"                
##  [28] "MoleculeExperiment"       "nnSVG"                    "OSTA.data"               
##  [31] "pengls"                   "poem"                     "RegionalST"              
##  [34] "retrofit"                 "scatterHatch"             "sccomp"                  
##  [37] "scDesign3"                "scider"                   "SEraster"                
##  [40] "shinyDSP"                 "signifinder"              "simpleSeg"               
##  [43] "smoothclust"              "smoppix"                  "sosta"                   
##  [46] "SpaceMarkers"             "SpaceTrooper"             "spacexr"                 
##  [49] "SpaNorm"                  "spARI"                    "spaSim"                  
##  [52] "SpatialDecon"             "SpatialExperiment"        "SpatialExperimentIO"     
##  [55] "spatialFDA"               "SpatialFeatureExperiment" "SpatialOmicsOverlay"     
##  [58] "spatialSimGP"             "SPIAT"                    "spicyR"                  
##  [61] "spoon"                    "SpotClean"                "SPOTlight"               
##  [64] "SpotSweeper"              "standR"                   "Statial"                 
##  [67] "stJoincount"              "stPipe"                   "SVP"                     
##  [70] "tidySpatialExperiment"    "tomoda"                   "tomoseqr"                
##  [73] "tpSVG"                    "VisiumIO"                 "visiumStitched"          
##  [76] "Voyager"                  "XeniumIO"                 "scFeatures"              
##  [79] "spatialHeatmap"
```

We can also browse for packages that also include more specific terms, for example:

Code

```
# view "Spatial Clustering" packages
df$Package[.f(c("Spatial", "Clustering"))]
```

```
##   [1] "Banksy"         "clustSIGNAL"    "concordexR"     "FuseSOM"       
##   [5] "hoodscanR"      "imcRtools"      "poem"           "smoothclust"   
##   [9] "spARI"          "SPIAT"          "stJoincount"    "stPipe"        
##  [13] "tomoda"         "spatialHeatmap"
```

### 4.1.2 Metrics

We will now have a look at the package ecosystem in a more quantitative manner. Specifically, we will investigate how the number of available packages evolves over time, and how long individual packages remain available (their “lifetime”). Lastly, we will quantify the number of times more specific subterms (co-)occur.

#### 4.1.2.1 Number of packages

Let’s first summarize the number of packages available over time, focusing specifically on “SingleCell” and “Spatial” packages.

Note that we here rely on the first and last date at which a package was available through Bioconductor; some packages might have been deprecated at one point or another.

Code

```
# dependencies
library(dplyr)
library(ggplot2)
# specify 'biocViews' of interest
names(ids) <- ids <- c("SingleCell", "Spatial")
now <- as.Date(format(Sys.Date(), "%Y-%m-%d"))
gg <- lapply(ids, \(id) {
    # get metadata & simplify naming
    nm <- df$Package[.f(id)]
    ys <- getPkgYearsInBioc(nm) |>
        mutate(first=first_version_release_date) |>
        mutate(last=last_version_release_date) |>
        mutate(last=case_when(is.na(last)~now, TRUE~last)) |>
        filter(!is.na(first))
    # complete months between first/last dates
    ys <- lapply(split(ys, ys$package), \(.) {
        data.frame(
            package=.$package, 
            date=seq(.$first, .$last))
    }) |> do.call(what=rbind) 
    # get cumulative number of packages available each month
    ys |> group_by(date) |> count()
}) |> bind_rows(.id="biocViews")
```

Code

```
ggplot(gg, aes(date, n, col=biocViews)) + 
    geom_line(linewidth=0.8) +
    geom_smooth(data=filter(gg, n >= 5), 
        method="lm", se=FALSE, linewidth=1) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    labs(x=NULL, y="# packages") +
    theme_bw() + theme(
        aspect.ratio=1, 
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))
```

![](bkg-ecosystem_files/figure-html/unnamed-chunk-6-1.png)

To get an idea of the project’s “growth” in this space, we can regress the number of packages against time using a linear model (LM). We here skip dates at which very few packages were availabile/being developed:

Code

```
# fit linear model where x = 'date', y = # packages
# (starting at 'date' where at least 5 packages exist)
dy <- round(as.integer(diff(slice_min(group_by(.gg, biocViews), date)$date))/365, 2)
(bs <- .gg |>
    group_by(biocViews) |>
    group_split() |>
    # return coefficients = # packages 
    # added each month (on average)
    sapply(\(.) coef(lm(n~date, .))[[2]]) |>
    setNames(sort(unique(gg$biocViews))))
```

```
##  SingleCell    Spatial 
##  0.09524586 0.04039224
```

The resulting coefficients tell us that about 1.14 single cell and 0.48 spatial packages are being added each year (on average), and there is a delay of about 4.49 years between them.

#### 4.1.2.2 Package lifetimes

We can further inspect the “lifetime” of packages, i.e., how long they remain (installable) on Bioconductor:

A package will be deprecated if it fails to build and/or pass checks on the Bioconductor Build System (BBS), provided its maintainer is unresponsive and does not take action to fix the package before the next (six-monthly) release.

Code

```
gg <- lapply(ids, \(id) {
    nm <- df$Package[.f(id)]
    ys <- BiocPkgTools:::getPkgYearsInBioc(nm)
}) |> bind_rows(.id="biocViews") |>
    mutate(years=approx_years_in) |>
    filter(!is.na(years))
mu <- gg |>
    group_by(biocViews) |>
    summarise_at("years", mean)
# print
cat("years in Bioconductor:\n")
summary(gg$approx_years_in)
# plot
ggplot(gg, aes(years, fill=biocViews)) + 
    geom_histogram(alpha=1/3, binwidth=0.5) +
    geom_vline(
        linewidth=1, data=mu,
        aes(xintercept=years, col=biocViews)) +
    scale_x_continuous(breaks=seq(0, 100, 2)) +
    labs(x="lifetime (years)", y="# packages") +
    theme_bw() + theme(
        aspect.ratio=1, 
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))
```

![](bkg-ecosystem_files/figure-html/unnamed-chunk-8-1.png)

```
##  years in Bioconductor:
##     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    0.000   1.500   3.500   3.755   5.500  17.000
```

#### 4.1.2.3 Subterms

Finally, let’s investigate more specific `biocViews`, e.g., the type of work tasks different packages cover. We will first count the number of packages that list specific terms that might be of interest in the context of single cell and spatial omics data analysis. Secondly, we will visualize their co-occurence.

Note that packages can list an arbitary number of terms so that packages can be counted in more than one category; e.g., “Visualization” appears in *most*.

Code

```
i <- c("SingleCell", "Spatial")
j <- c(
    "BatchEffect", "Normalization", "QualityControl", "Visualization", # WorkflowStep
    "Clustering", "DimensionReduction", "FeatureExtraction", # StatisticalMethod
    "DifferentialExpression", "GeneSetEnrichment") # BiologicalQuestion
names(i) <- i; names(j) <- j
# count packages for each pair of terms
gg <- mapply(
    i=i, j=rep(j, each=2), 
    SIMPLIFY=FALSE, \(i, j) {
        n <- sum(.f(c(i, j)))
        data.frame(i, j, n)
    }) |> do.call(what=rbind)
# order x-axis by total
xo <- gg |>
    group_by(j) |>
    summarise_at("n", sum) |>
    arrange(n) |>
    pull(j)
ggplot(gg, aes(j, n, fill=i)) + 
    scale_x_discrete(limits=xo) +
    geom_col(position="dodge", alpha=2/3) +
    labs(x=NULL, y="# packages", fill="biocViews") +
    theme_bw() + theme(
        aspect.ratio=1, 
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))
```

![](bkg-ecosystem_files/figure-html/unnamed-chunk-9-1.png)

Because there are yet very few “Spatial” packages, we will pool “Spatial” and “SingleCell” tools when counting the number of times different terms appear together; i.e., `biocViews` contain “Spatial” OR “SingleCell”, together with “Clustering” AND “Visualization”, etc. (For concise labeling, we abbreviate `biocViews` to capital letters only; e.g., “BatchEffects” becomes “BE”.)

Code

```
gg <- lapply(i, \(i) {
    lapply(seq_along(j), \(n) {
        js <- combn(j, n, simplify=FALSE)
        lapply(js, \(j) {
            n <- sum(.f(c(i, j)))
            j <- gsub("[a-z]", "", j)
            j <- paste(j, collapse="+")
            data.frame(i, j, n)
        }) |> do.call(what=rbind)
    }) |> do.call(what=rbind)
}) |> do.call(what=rbind) |>
    group_by(j) |>
    summarize_at("n", sum) |>
    slice_max(n, n=30)
yo <- gg$j[order(gg$n)]
ggplot(gg, aes(n, j)) + 
    geom_col(alpha=1/3, fill="blue") +
    labs(y=NULL, x="# packages") +
    scale_y_discrete(limits=yo) +
    theme_bw() + theme(
        aspect.ratio=1, 
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))
```

![](bkg-ecosystem_files/figure-html/unnamed-chunk-10-1.png)

## 4.2 Appendix

Back to top