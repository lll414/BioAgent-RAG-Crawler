Source URL: https://rformassspectrometry.github.io/book/sec-ex.html
Date Scraped: 2025-12-15

---

# Chapter 6 Supplementary exercises

## 6.1 Raw data and identification results

* Download the 3 first mzML and mzID files from the
  [PXD022816](https://www.ebi.ac.uk/pride/archive/projects/PXD022816)
  project (Morgenstern, Barzilay, and Levin 2021).
* Generate a `Spectra` object and a table of filtered PSMs. Visualise
  the total ion chromatograms and check the quality of the
  identification data by comparing the density of the decoy and target
  PSMs id scores for each file.
* Join the raw and identification data. Beware though that the joining
  must now be performed by spectrum ids and by files.
* Extract the PSMs that have been matched to peptides from protein
  `O43175` and compare and cluster the scans. Hint: once you have
  created the smaller `Spectra` object with the scans of interest,
  switch to an in-memory backend to seed up the calculations.
* Generate total ion chromatograms for each acquisition and annotate
  the MS1 scans with the number of PSMs using the
  `countIdentifications()` function, as shown above. The function will
  automatically perform the counts in parallel for each acquisition.

## 6.2 Search engine

Download the [spectra and protein
database](https://dataverse.uclouvain.be/dataset.xhtml?persistentId=doi:10.14428/DVN/QR0OYG&faces-redirect=true#)
needed for the exercise (here is a [direct
link](https://dataverse.uclouvain.be/api/access/dataset/:persistentId/?persistentId=doi:10.14428/DVN/QR0OYG%5D). The
protein database is in fasta format and can be processed as described
in the section [4.8](sec-id.html#sec-id-seq) *Reading and processing protein
sequences*. The MS2 spectra are provided in the Mascot Generic Format
(MGF) format, that can be loaded using the dedicated
[MsBackendMgf](https://rformassspectrometry.github.io/MsBackendMgf/)
backend as `Spectra` objects.

You are asked to write code to identify the spectra, following the
principles defined in the *Identification data* chapter, include ways
to provide confidence in your identification results, beyond a single
identification score.

**Hints:**

* Focus on expected peptides sequences that longer than 6 and shorter
  than 28 amino acids to reduce the search space.
* Do not search each MS2 scan against the whole database, but focus on
  peptides that have a mass that is close to the scan’s precuror mass.
* To calculate the mass of a peptides, use `m/z * c - proton_mass * c`, where `m/z` and `c` is the mass-over-charge and the charge of
  the precursor and `proton_mass` is the mass of a proton (available
  with `PSMatch::getAtomicMass()[["p"]]`).
* The `PSMatch::getAminoAcids()` function returns a `data.frame` of
  amino acid properties.
* Consider using `spectrapply` to iterate of the individual scans of a
  `Spectra` object.

## 6.3 Quantitative data processing

Following up from the quantitative data analysis seen on chapter
[5](sec-quant.html#sec-quant), the following file includes a third condition C and
a two additional lab, tallying now 27 samples.

```
f <- MsDataHub::cptac_a_b_c_peptides.txt()
```

|  | LTQ-Orbitrap\_86 | LTQ-OrbitrapO\_65 | LTQ-OrbitrapW\_56 |
| --- | --- | --- | --- |
| 6A | 3 | 3 | 3 |
| 6B | 3 | 3 | 3 |
| 6C | 3 | 3 | 3 |

The full design is shown below.

| TRUE | id | condition | lab | previous |
| --- | --- | --- | --- | --- |
| 6A\_1 | 1 | 6A | LTQ-Orbitrap\_86 | new |
| 6A\_2 | 2 | 6A | LTQ-Orbitrap\_86 | new |
| 6A\_3 | 3 | 6A | LTQ-Orbitrap\_86 | new |
| 6A\_4 | 4 | 6A | LTQ-OrbitrapO\_65 | new |
| 6A\_5 | 5 | 6A | LTQ-OrbitrapO\_65 | new |
| 6A\_6 | 6 | 6A | LTQ-OrbitrapO\_65 | new |
| 6A\_7 | 7 | 6A | LTQ-OrbitrapW\_56 |  |
| 6A\_8 | 8 | 6A | LTQ-OrbitrapW\_56 |  |
| 6A\_9 | 9 | 6A | LTQ-OrbitrapW\_56 |  |
| 6B\_1 | 1 | 6B | LTQ-Orbitrap\_86 | new |
| 6B\_2 | 2 | 6B | LTQ-Orbitrap\_86 | new |
| 6B\_3 | 3 | 6B | LTQ-Orbitrap\_86 | new |
| 6B\_4 | 4 | 6B | LTQ-OrbitrapO\_65 | new |
| 6B\_5 | 5 | 6B | LTQ-OrbitrapO\_65 | new |
| 6B\_6 | 6 | 6B | LTQ-OrbitrapO\_65 | new |
| 6B\_7 | 7 | 6B | LTQ-OrbitrapW\_56 |  |
| 6B\_8 | 8 | 6B | LTQ-OrbitrapW\_56 |  |
| 6B\_9 | 9 | 6B | LTQ-OrbitrapW\_56 |  |
| 6C\_1 | 1 | 6C | LTQ-Orbitrap\_86 | new |
| 6C\_2 | 2 | 6C | LTQ-Orbitrap\_86 | new |
| 6C\_3 | 3 | 6C | LTQ-Orbitrap\_86 | new |
| 6C\_4 | 4 | 6C | LTQ-OrbitrapO\_65 | new |
| 6C\_5 | 5 | 6C | LTQ-OrbitrapO\_65 | new |
| 6C\_6 | 6 | 6C | LTQ-OrbitrapO\_65 | new |
| 6C\_7 | 7 | 6C | LTQ-OrbitrapW\_56 | new |
| 6C\_8 | 8 | 6C | LTQ-OrbitrapW\_56 | new |
| 6C\_9 | 9 | 6C | LTQ-OrbitrapW\_56 | new |

* Repeat the analysis described in chapter [5](sec-quant.html#sec-quant) using the
  extended dataset, trying to optimise true positive results and
  avoiding false positive. Think about the best experimental design
  approach, how to best process the data, visualising important steps
  along the way, to conclude with a volcano plot and a table tallying
  the number of true/false positive/negative results.