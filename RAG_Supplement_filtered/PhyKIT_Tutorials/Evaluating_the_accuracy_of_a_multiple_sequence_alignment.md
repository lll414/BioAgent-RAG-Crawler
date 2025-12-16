# Evaluating the accuracy of a multiple sequence alignment

Source URL: https://jlsteenwyk.com/PhyKIT/tutorials/index.html
Date Scraped: 2025-12-11

---

---

Evaluating the accuracy of multiple sequence alignments is an appropriate way to benchmark multiple sequence alignment strategies.
Two popular methods to assess multiple sequence alignment accuracy are sum-of-pairs score and column score, which were introduced by
Thompson et al., Nucleic Acids Research (1999), doi: 10.1093/nar/27.13.2682. Sum-of-pairs is calculated by summing the correctly
aligned residue pairs over all pairs of sequences. Column score is calculated by summing the correctly aligned columns over all
columns in an alignment. Both metrics range from 0 to 1 and higher values indicate more accurate alignments. Correctly aligned
pairs or columns require knowing some ground truth of what the correct alignment is. Thus, a reference alignment that is
perfectly (or near-perfectly) aligned is required. A reference alignment can be generated using simulations or be obtained from
publicly available databases such as BAliBASE 4 (<http://www.lbgi.fr/balibase/>). For this tutorial we will use a reference alignment
from BAliBASE.

**Download test data:
[`msa_accuracy.tar.gz`](../_downloads/1fd76132a6e0df73727f2e3628e402bb/msa_accuracy.tar.gz)**

### Step 0: Generate query alignments[](#step-0-generate-query-alignments "Link to this heading")

In the *msa\_accuracy* directory, there are two fasta files: *BBA0001\_query.faa*, an unaligned set of sequences and *BBA0001\_reference.faa*,
the reference alignment. We will align *BBA0001\_query.faa* using three different strategies implemented in Mafft, v.7.475
(<https://mafft.cbrc.jp/alignment/software/>), and evaluate the accuracy of each strategy.

To do so, please ensure Mafft is installed and then execute the following commands:

```
# first alignment
mafft --localpair BBA0001_query.faa > BBA0001_query.localpair.faa

# second alignment
mafft --genafpair BBA0001_query.faa > BBA0001_query.genafpair.faa

# third alignment
mafft --globalpair BBA0001_query.faa > BBA0001_query.globalpair.faa
```

To gain some initial insight as to whether the alignments differ, we can look at the length
of each alignment using the *aln\_len* function

```
for i in $(ls *pair.faa) ; do phykit aln_len $i ; done
1560
1464
1497
```

However, alignment length is not a measurement of accuracy. Thus, we will score each alignment
using the *sum\_of\_pairs\_score* and *column\_score* functions.

### Step 1: Score each alignment[](#step-1-score-each-alignment "Link to this heading")

For both functions, the first argument is the query alignment and the *-r/\-\-reference* argument specifies the reference alignment.
We will programmatically score each alignment using the same for loop that was used to calculate alignment length.

```
echo -e "BAliBASE_id_and_aln_strategy\tsop\tcs"
for i in $(ls *pair.faa)
do
   sop=$(phykit sum_of_pairs_score $i -r BBA0001_reference.faa)
   cs=$(phykit column_score $i -r BBA0001_reference.faa)
   echo -e "$i\t$sop\t$cs"
done
BAliBASE_id_and_aln_strategy sop     cs
BBA0001_query.genafpair.faa  0.8964  0.2943
BBA0001_query.globalpair.faa 0.9025  0.292
BBA0001_query.localpair.faa  0.8992  0.2959
```

Examination of the output reveals that the *globalpair* strategy has more correctly aligned pairs because it has a higher sum-of-pairs
score whereas the *localpair* strategy has more correctly aligned columns. Of note, the column scores are generally low, which
reflects a potential limitation of column score wherein column score is sensitive to alignment errors.

In summary, calculating sum-of-pairs score and column score can help assess the accuracy of multiple sequence alignment strategies.