Source URL: https://bioconductor.org/books/release/OSTA/pages/img-introduction.html
Date Scraped: 2025-12-15

---

1. [Imaging-based platforms](https://bioconductor.org/books/release/OSTA/pages/img-introduction.html)
2. [17  Introduction](https://bioconductor.org/books/release/OSTA/pages/img-introduction.html)

# 17  Introduction

## 17.1 Introduction to imaging-based ST

Imaging-based ST assays have gone from resolving 100s to 1000s of features. Their commercialization by Vizgen, NanoString, and 10x Genomics has made these data increasingly popular. In general, imaging-based ST platforms rely on multiplexed error-robust fluorescence *in situ* hybridization (FISH), as originally proposed by Chen et al. ([2015](https://bioconductor.org/books/release/OSTA/pages/img-introduction.html#ref-Chen2015-MERFISH)). Briefly, targets are preassigned a barcode of (multi-color) ‘on’ and (no-color) ‘off’ bits, whereas barcodes are optimized to differ from each other by some bits as to minimize readout errors. After iterative imaging-bleaching, individual reporter binding events are identified computationally (spot calling). Aligning spots along the image z-stack then allows to identify individual targets based on a predefined barcoding scheme.

high-plex barcoding

1,000-plex CosMx, for example, employs 64-bit barcodes (4 immunofluorescent reports over 16 cycles of imaging-bleaching) with a Hamming weight and distance of 4 (i.e., every target is ‘on’ in 4 rounds and ‘off’ in 12 rounds; each barcode differs from all others by 4 bits as to minimize RNA decoding errors). Combinatorially, a 64-bit barcoding scheme is sufficient to encode a much larger number of targets, namely, ~20,000 protein-coding genes.

Commercially available panels can typically include three types of barcodes:

* **RNA targets** determine the ‘plexity’ of the panel.
* **Negative probes** serve to quantify non-specific binding.
* **System controls** serve to quantify spot calling errors.

non-RNA targets

Modeled after synthetic sequences from the External RNA Controls Consortium (ERCC), **negative probes** contain hybridization regions that are not complementary to the genome or transcriptome of the organism under study; their detection thus corresponds to non-specific ISH probe hybridization events (e.g., in ‘sticky’ regions of the tissue).

Depending on the barcoding scheme (number of bits, Hamming weight/distance), a number of barcodes will be left unoccupied; these **system controls** (also referred to as blank or false codes) serve to quantify misidentification of reporter readout.

## 17.2 Platforms

### 17.2.1 Xenium (10x Genomics)

The [10x Genomics Xenium](https://www.10xgenomics.com/platforms/xenium) platform employs rolling circle amplification (RCA) to boost fluorescent signals from reporter-binding events, thereby improving overall spot calling sensitivity.

Xenium gives an imageable area of 12mm x 24 mm, although slightly less to avoid fiducial markers. Several pre-designed panels are available, including a “Pan Tissue and Pathways” for both mouse and human, which both include over 5000 genes and there is a possibility to custom design an additional 100 genes. More focused panels are available for human brain, lung, colon and breast tissue as well as cancer immunological subpopulations; for mouse, focused brain and multi-tissue panels are available.

### 17.2.2 CosMx (NanoString)

Unlike Xenium, [CosMx](https://nanostring.com/products/cosmx-spatial-molecular-imager/single-cell-imaging-overview/) ([He et al. 2022](https://bioconductor.org/books/release/OSTA/pages/img-introduction.html#ref-He2022-CosMx)) does not rely on RCA. As a result, CosMx tends to be less sensitive than Xenium (albeit higher-plex at present). Furthermore, segmentation is performed separately for each FOV (without stitching). This results in a variety of technical artefacts, such as fractured and possibly duplicated cells near FOV borders.

The NanoString CosMx platform gives an imageable area of 15x20 mm and allows up to 6000 RNAs (and now up to 68 proteins) to be assayed in a spatial context. Their human “discovery” 6000 gene panel is currently being expanded to full-transcriptome (>18,000 targets), which should become available in 2025. Other more targeted panels are available, including a 1000-gene “Universal Cell Characterization”, a focused 100-gene immuno-oncology panel or a 1000-gene mouse neuroscience panel. In addition, custom panels are available in either an add-on (up to 50 genes) or standalone (up to 300 genes) format. (Antibody-based) Protein panels are currently available for human immuno-oncology and mouse neuroscience contexts.

### 17.2.3 MERSCOPE (Vizgen)

The [Vizgen MERSCOPE](https://vizgen.com/products/) platform, among the fluorescence-based targeted spatial transcriptomics platforms, allows the highest resolution (≤20 nm), while also providing a high capacity for multiplexing, with (customizable) panels containing up to 1000 genes, and an imageable area of 1 cm2. MERSCOPE was recently shown to be amongst the most sensitive in a recent head-to-head comparison ([Hartman and Satija 2024](https://bioconductor.org/books/release/OSTA/pages/img-introduction.html#ref-Hartman2024-comparative)) using mouse brain tissue.

## 17.3 Appendix

### Resources

* Available datasets
  + [Xenium (10x Genomics)](https://www.10xgenomics.com/datasets)
  + [CosMx (NanoString)](https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/)
  + [MERSCOPE (Vizgen)](https://vizgen.com/data-release-program)
* Platform-specific methods
  + [CosMx analysis](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space)

Back to top