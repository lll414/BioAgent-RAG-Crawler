# Alignment-based functions

Source URL: https://jlsteenwyk.com/PhyKIT/usage/index.html
Date Scraped: 2025-12-11

---

---

### Alignment length[](#alignment-length "Link to this heading")

[![../_images/aln_len.png](../_images/aln_len.png)](../_images/aln_len.png)

Function names: alignment\_length; aln\_len; al   
Command line interface: pk\_alignment\_length; pk\_aln\_len; pk\_al

Length of an input alignment is calculated using this function.

Longer alignments are associated with strong phylogenetic signal.

Association between alignment length and phylogenetic signal
was determined by Shen et al., Genome Biology and Evolution (2016),
doi: 10.1093/gbe/evw179.

```
phykit aln_len <alignment>
```

Options:   
*<alignment>*: first argument after function name should be an alignment file

### Alignment length no gaps[](#alignment-length-no-gaps "Link to this heading")

[![../_images/aln_len_no_gaps.png](../_images/aln_len_no_gaps.png)](../_images/aln_len_no_gaps.png)

Function names: alignment\_length\_no\_gaps; aln\_len\_no\_gaps; alng   
Command line interface: pk\_alignment\_length\_no\_gaps; pk\_aln\_len\_no\_gaps; pk\_alng

Calculate alignment length excluding sites with gaps.

Longer alignments when excluding sites with gaps is
associated with strong phylogenetic signal.

PhyKIT reports three tab delimited values:
col1: number of sites without gaps
col2: total number of sites
col3: percentage of sites without gaps

Association between alignment length when excluding sites
with gaps and phylogenetic signal was determined by Shen
et al., Genome Biology and Evolution (2016),
doi: 10.1093/gbe/evw179.

```
phykit aln_len_no_gaps <alignment>
```

Options:   
*<alignment>*: first argument after function name should be an alignment file

### Alignment recoding[](#alignment-recoding "Link to this heading")

Function names: alignment\_recoding; aln\_recoding; recode   
Command line interface: pk\_alignment\_recoding; pk\_aln\_recoding; bk\_recode

Recode alignments using reduced character states.

Alignments can be recoded using established or
custom recoding schemes. Recoding schemes are
specified using the -c/–code argument. Custom
recoding schemes can be used and should be formatted
as a two column file wherein the first column is the
recoded character and the second column is the character
in the alignment.

```
phykit alignment_recoding <fasta> [-c/--code <code>]
```

Codes for which recoding scheme to use:   
**RY-nucleotide**   
R = purines (i.e., A and G)   
Y = pyrimidines (i.e., T and C)

**SandR-6**   
0 = A, P, S, and T   
1 = D, E, N, and G   
2 = Q, K, and R   
3 = M, I, V, and L   
4 = W and C   
5 = F, Y, and H

**KGB-6**   
0 = A, G, P, and S   
1 = D, E, N, Q, H, K, R, and T   
2 = M, I, and L   
3 = W   
4 = F and Y   
5 = C and V

**Dayhoff-6**   
0 = A, G, P, S, and T   
1 = D, E, N, and Q   
2 = H, K, and R   
3 = I, L, M, and V   
4 = F, W, and Y   
5 = C

**Dayhoff-9**   
0 = D, E, H, N, and Q   
1 = I, L, M, and V   
2 = F and Y   
3 = A, S, and T   
4 = K and R   
5 = G   
6 = P   
7 = C   
8 = W

**Dayhoff-12**   
0 = D, E, and Q   
1 = M, L, I, and V   
2 = F and Y   
3 = K, H, and R   
4 = G   
5 = A   
6 = P   
7 = S   
8 = T   
9 = N   
A = W   
B = C

**Dayhoff-15**   
0 = D, E, and Q   
1 = M and L   
2 = I and V   
3 = F and Y   
4 = G   
5 = A   
6 = P   
7 = S   
8 = T   
9 = N   
A = K   
B = H   
C = R   
D = W   
E = C

**Dayhoff-18**   
0 = F and Y   
1 = M and L   
2 = I   
3 = V   
4 = G   
5 = A   
6 = P   
7 = S   
8 = T   
9 = D   
A = E   
B = Q   
C = N   
D = K   
E = H   
F = R   
G = W   
H = C

Options:   
*<alignment>*: first argument after function name should be an alignment file   
*-c/--code*: argument to specify the recoding scheme to use

### Column score[](#column-score "Link to this heading")

[![../_images/column_score.png](../_images/column_score.png)](../_images/column_score.png)

Function names: column\_score; cs   
Command line interface: pk\_column\_score; pk\_cs

Calculates column score.

Column is an accuracy metric for a multiple alignment relative
to a reference alignment. It is calculated by summing the correctly
aligned columns over all columns in an alignment. Thus, values range
from 0 to 1 and higher values indicate more accurate alignments.

Column score is calculated following Thompson et al., Nucleic
Acids Research (1999), doi: 10.1093/nar/27.13.2682.

```
phykit column_score <alignment> --reference <reference_alignment>
```

Options:   
*<alignment>*: first argument after function name should be a query
fasta alignment file to be scored for accuracy   
*-r/\-\-reference*: reference alignment to compare the query alignment
to

### Compositional bias per site[](#compositional-bias-per-site "Link to this heading")

Function names: compositional\_bias\_per\_site; comp\_bias\_per\_site; cbps   
Command line interface: pk\_compositional\_bias\_per\_site; pk\_comp\_bias\_per\_site; pk\_cbps

Calculates compositional bias per site in an alignment.

Site-wise chi-squared tests are conducted in an alignment to
detect compositional biases. PhyKIT outputs four columns:   
col 1: index in alignment   
col 2: chi-squared statistic (higher values indicate greater bias)   
col 3: multi-test corrected p-value (Benjamini-Hochberg false discovery rate procedure)   
col 4: uncorrected p-value

```
phykit comp_bias_per_site <alignment>
```

Options:   
*<alignment>*: first argument after function name should be a query
fasta alignment to calculate the site-wise compositional bias of

### Create concatenation matrix[](#create-concatenation-matrix "Link to this heading")

[![../_images/create_concat_matrix.png](../_images/create_concat_matrix.png)](../_images/create_concat_matrix.png)

Function names: create\_concatenation\_matrix, create\_concat, cc   
Command line interface: pk\_create\_concatenation\_matrix, pk\_create\_concat, pk\_cc

Create a concatenated alignment file. This function is
used to help in the construction of multi-locus data
matrices.

PhyKIT will output three files:
1) A fasta file with ‘.fa’ appended to the prefix specified with the -p/\-\-prefix parameter.
2) A partition file ready for input into RAxML or IQ-tree.
3) An occupancy file that summarizes the taxon occupancy per sequence.

```
phykit create_concat -a <file> -p <string>
```

Options:   
*-a/\-\-alignment*: alignment list file. File should contain a single column list of alignment
sequence files to concatenate into a single matrix. Provide path to files relative to
working directory or provide absolute path.   
*-p/\-\-prefix*: prefix of output files

### Evolutionary Rate per Site[](#evolutionary-rate-per-site "Link to this heading")

Function names: evolutionary\_rate\_per\_site; evo\_rate\_per\_site; erps   
Command line interface: pk\_evolutionary\_rate\_per\_site; pk\_evo\_rate\_per\_site; pk\_erps

Estimate evolutionary rate per site.

Evolutionary rate per site is one minus the sum of squared frequency of different
characters at a given site. Values may range from 0 (slow evolving; no diversity
at the given site) to 1 (fast evolving; all characters appear only once).

```
phykit evo_rate_per_site <alignment>
```

Options:   
*<alignment>*: first argument after function name should be a query
fasta alignment to calculate the site-wise evolutionary rate of

### Faidx[](#faidx "Link to this heading")

[![../_images/faidx.png](../_images/faidx.png)](../_images/faidx.png)

Function names: faidx; get\_entry; ge   
Command line interface: pk\_faidx; pk\_get\_entry; pk\_ge

Extracts sequence entry from fasta file.

This function works similarly to the faidx function
in samtools, but does not requiring an indexing step.

To obtain multiple entries, input multiple entries separated
by a comma (,). For example, if you want entries
named “seq\_0” and “seq\_1”, the string “seq\_0,seq\_1”
should be associated with the -e argument.

```
phykit faidx <fasta> -e/--entry <fasta entry>
```

Options:   
*<fasta>*: first argument after function name should be a fasta file   
*-v/\-\-verbose*: entry name to be extracted from the inputted fasta file
entry

### Guanine-cytosine (GC) content[](#guanine-cytosine-gc-content "Link to this heading")

[![../_images/gc_content.png](../_images/gc_content.png)](../_images/gc_content.png)

Function names: gc\_content; gc   
Command line interface: pk\_gc\_content; pk\_gc

Calculate GC content of a fasta file.

GC content is negatively correlated with phylogenetic signal.

If there are multiple entries, use the -v/\-\-verbose option
to determine the GC content of each fasta entry separately.
Association between GC content and phylogenetic signal was
determined by Shen et al., Genome Biology and Evolution (2016),
doi: 10.1093/gbe/evw179.

```
phykit gc_content <fasta> [-v/--verbose]
```

Options:   
*<fasta>*: first argument after function name should be a fasta file   
*-v/\-\-verbose*: optional argument to print the GC content of each fasta
entry

### Pairwise identity[](#pairwise-identity "Link to this heading")

[![../_images/pairwise_identity.png](../_images/pairwise_identity.png)](../_images/pairwise_identity.png)

Function names: pairwise\_identity; pairwise\_id, pi   
Command line interface: pk\_pairwise\_identity; pk\_pairwise\_id, pk\_pi

Calculate the average pairwise identity among sequences.

Pairwise identities can be used as proxies for the evolutionary rate of sequences.

Pairwise identity is defined as the number of identical
columns (including gaps) between two aligned sequences divided
by the number of columns in the alignment. Summary statistics
are reported unless used with the verbose option in which
all pairwise identities will be reported.

An example of pairwise identities being used as a proxy
for evolutionary rate can be found here: Chen et al.
Genome Biology and Evolution (2017), doi: 10.1093/gbe/evx147.

```
phykit pairwise_identity <alignment> [-v/--verbose]
```

Options:   
*<alignment>*: first argument after function name should be an alignment file   
*-v/\-\-verbose*: optional argument to print identity per pair|br|
*-e/–exclude\_gaps*: if a site has a gap, ignore it

### Parsimony informative sites[](#parsimony-informative-sites "Link to this heading")

Function names: parsimony\_informative\_sites; pis   
Command line interface: pk\_parsimony\_informative\_sites; pk\_pis

Calculate the number and percentage of parismony
informative sites in an alignment.

The number of parsimony informative sites in an alignment
is associated with strong phylogenetic signal.

PhyKIT reports three tab delimited values:
col1: number of parsimony informative sites
col2: total number of sites
col3: percentage of parsimony informative sites

Association between the number of parsimony informative
sites and phylogenetic signal was determined by Shen
et al., Genome Biology and Evolution (2016),
doi: 10.1093/gbe/evw179 and Steenwyk et al., PLOS Biology
(2020), doi: 10.1371/journal.pbio.3001007.

```
phykit parsimony_informative_sites <alignment>
```

Options:   
*<alignment>*: first argument after function name should be an alignment file

### Protein-to-nucleotide alignment[](#protein-to-nucleotide-alignment "Link to this heading")

Function names: thread\_dna; pal2nal, p2n   
Command line interface: pk\_thread\_dna; pk\_pal2nal, pk\_p2n

Thread DNA sequence onto a protein alignment to create a
codon-based alignment.

This function requires input alignments are in fasta format.
Codon alignments are then printed to stdout. Note, paired
sequences are assumed to have the same name between the
protein and nucleotide file. The order does not matter.

To thread nucleotide sequences over a trimmed amino acid
alignment, provide PhyKIT with a log file specifying which
sites have been trimmed and which have been kept. The log
file must be formatted the same as the log files outputted
by the alignment trimming toolkit ClipKIT (see -l in ClipKIT
documentation.) Details about ClipKIT can be seen here:
<https://github.com/JLSteenwyk/ClipKIT>.

If using a ClipKIT log file, the untrimmed protein alignment
should be provided in the -p/–protein argument.

```
phykit thread_dna -p <file> -n <file> [-s]
```

Options:   
*-p/\-\-protein*: protein alignment file   
*-n/\-\-nucleotide*: nucleotide sequence file   
*-c/\-\-clipkit\_log*: clipkit outputted log file   
*-s/\-\-stop*: boolean for whether or not stop codons should be kept.
If used, stop codons will be removed.

### Relative composition variability[](#relative-composition-variability "Link to this heading")

Function names: relative\_composition\_variability; rel\_comp\_var; rcv   
Command line interface: pk\_relative\_composition\_variability; pk\_rel\_comp\_var; pk\_rcv

Calculate RCV (relative composition variability) for an alignment.

Lower RCV values are thought to be desirable because they represent
a lower composition bias in an alignment. Statistically, RCV describes
the average variability in sequence composition among taxa.

RCV is calculated following Phillips and Penny, Molecular Phylogenetics
and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

```
phykit relative_composition_variability <alignment>
```

Options:   
*<alignment>*: first argument after function name should be an alignment file

### Relative composition variability, taxon[](#relative-composition-variability-taxon "Link to this heading")

Function names: relative\_composition\_variability\_taxon; rel\_comp\_var\_taxon; rcvt   
Command line interface: pk\_relative\_composition\_variability\_taxon; pk\_rel\_comp\_var\_taxon; pk\_rcvt

Calculate RCVT (relative composition variability, taxon) for an alignment.

RCVT is the relative composition variability metric for individual taxa.
This facilitates identifying specific taxa that may have compositional
biases. Lower RCVT values are more desirable because they indicate
a lower composition bias for a given taxon in an alignment.

```
phykit relative_composition_variability_taxon <alignment>
```

Options:   
*<alignment>*: first argument after function name should be an alignment file

### Rename FASTA entries[](#rename-fasta-entries "Link to this heading")

Function names: rename\_fasta\_entries; rename\_fasta   
Command line interface: pk\_rename\_fasta\_entries; pk\_rename\_fasta

Renames fasta entries.

Renaming fasta entries will follow the scheme of a tab-delimited
file wherein the first column is the current fasta entry name and
the second column is the new fasta entry name in the resulting
output alignment. Note, the input fasta file does not need to be
an alignment file.

```
phykit rename_fasta_entries <fasta> -i/--idmap <idmap> [-o/--output <output_file>]
```

Options:   
*<alignment>*: first argument after function name should be an alignment file   
*-i/\-\-idmap*: identifier map of current FASTA names (col1) and desired FASTA names (col2)

### Sum-of-pairs score[](#sum-of-pairs-score "Link to this heading")

Function names: sum\_of\_pairs\_score; sops; sop   
Command line interface: pk\_sum\_of\_pairs\_score; pk\_sops; pk\_sop

Calculates sum-of-pairs score.

Sum-of-pairs is an accuracy metric for a multiple alignment relative
to a reference alignment. It is calculated by summing the correctly
aligned residue pairs over all pairs of sequences. Thus, values range
from 0 to 1 and higher values indicate more accurate alignments.

Column score is calculated following Thompson et al., Nucleic
Acids Research (1999), doi: 10.1093/nar/27.13.2682.

```
phykit sum_of_pairs_score <alignment> --reference <reference_alignment>
```

Options:   
*<alignment>*: first argument after function name should be a query
fasta alignment file to be scored for accuracy   
*-r/\-\-reference*: reference alignment to compare the query alignment
to

### Variable sites[](#variable-sites "Link to this heading")

Function names: variable\_sites; vs   
Command line interface: pk\_variable\_sites; pk\_vs

Calculate the number of variable sites in an alignment.

The number of variable sites in an alignment is
associated with strong phylogenetic signal.
PhyKIT reports three tab delimited values:
col1: number of variable sites
col2: total number of sites
col3: percentage of variable sites

Association between the number of variable sites and
phylogenetic signal was determined by Shen et al.,
Genome Biology and Evolution (2016),
doi: 10.1093/gbe/evw179.

```
phykit variable_sites <alignment>
```

Options:   
*<alignment>*: first argument after function name should be an alignment file