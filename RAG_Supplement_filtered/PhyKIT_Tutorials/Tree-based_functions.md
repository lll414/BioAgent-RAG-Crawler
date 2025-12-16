# Tree-based functions

Source URL: https://jlsteenwyk.com/PhyKIT/usage/index.html
Date Scraped: 2025-12-11

---

---

### Bipartition support statistics[](#bipartition-support-statistics "Link to this heading")

Function names: bipartition\_support\_stats; bss   
Command line interface: pk\_bipartition\_support\_stats; pk\_bss

Calculate summary statistics for bipartition support.

High bipartition support values are thought to be desirable because
they are indicative of greater certainty in tree topology.

To obtain all bipartition support values, use the -v/\-\-verbose option.
In addition to support values for each node, the names of all terminal
branches tips are also included. Each terminal branch name is separated
with a semi-colon (;).

```
phykit bipartition_support_stats <tree> [-v/--verbose]
```

Options:   
*<alignment>*: first argument after function name should be a tree file   
*-v/\-\-verbose*: optional argument to print all bipartition support values

### Branch length multiplier[](#branch-length-multiplier "Link to this heading")

Function names: branch\_length\_multiplier; blm   
Command line interface: pk\_branch\_length\_multiplier; pk\_blm

Multiply branch lengths in a phylogeny by a given factor.

This can help modify reference trees when conducting simulations
or other analyses.

```
phykit branch_length_multiplier <tree> -f n [-o--output <output_file>]
```

Options:   
*<alignment>*: first argument after function name should be a tree file   
*-f/\-\-factor*: factor to multiply branch lengths by   
*-o/\-\-output*: optional argument to name the outputted tree file. Default
output will have the same name as the input file but with the suffix “.factor\_(n).tre”

### Collapse bipartitions[](#collapse-bipartitions "Link to this heading")

Function names: collapse\_branches, collapse, cb   
Command line interface: pk\_collapse\_branches, pk\_collapse, pk\_cb

Collapse branches on a phylogeny according to bipartition support.

Bipartitions will be collapsed if they are less than the user specified
value.

```
phykit collapse_branches <tree> -s/--support n [-o/--output <output_file>]
```

Options:   
*<alignment>*: first argument after function name should be a tree file   
*-s/\-\-support*: bipartitions with support less than this value will be
collapsed   
*-o/\-\-output*: optional argument to name the outputted tree file. Default
output will have the same name as the input file but with the suffix
“.collapsed\_(support).tre”

### Covarying evolutionary rates[](#covarying-evolutionary-rates "Link to this heading")

Function names: covarying\_evolutionary\_rates; cover   
Command line interface: pk\_covarying\_evolutionary\_rates; pk\_cover

Determine if two genes have a signature of covariation with one another.
Genes that have covarying evolutionary histories tend to have
similar functions and expression levels.

Input two phylogenies and calculate the correlation among relative
evolutionary rates between the two phylogenies. The two input trees
do not have to have the same taxa. This function will first prune both
trees to have the same tips. To transform branch lengths into relative
rates, PhyKIT uses the putative species tree’s branch lengths, which is
inputted by the user. As recommended by the original method developers,
outlier branche lengths are removed. Outlier branches have a relative
evolutionary rate greater than five.

PhyKIT reports two tab delimited values:
col1: correlation coefficient
col2: p-value

Method is empirically evaluated by Clark et al., Genome Research
(2012), doi: 10.1101/gr.132647.111. Normalization method using a
species tree follows Sato et al., Bioinformatics (2005), doi:
10.1093/bioinformatics/bti564.

```
phykit covarying_evolutionary_rates <tree_file_zero> <tree_file_one> -r/--reference <reference_tree_file> [-v/--verbose]
```

Options:   
*<tree\_file\_zero>*: first argument after function name should be an alignment file   
*<tree\_file\_one>*: first argument after function name should be an alignment file   
*-r/\-\-reference*: a tree to correct branch lengths by in the two input trees. Typically,
this is a putative species tree.   
*-v/\-\-verbose*: print out corrected branch lengths shared between tree 0 and tree 1

### Degree of violation of the molecular clock[](#degree-of-violation-of-the-molecular-clock "Link to this heading")

Function names: degree\_of\_violation\_of\_a\_molecular\_clock, dvmc   
Command line interface: pk\_degree\_of\_violation\_of\_a\_molecular\_clock, pk\_dvmc

Calculate degree of violation of a molecular clock (or DVMC) in a phylogeny.

Lower DVMC values are thought to be desirable because they are indicative
of a lower degree of violation in the molecular clock assumption.

Typically, outgroup taxa are not included in molecular clock analysis. Thus,
prior to calculating DVMC from a single gene tree, users may want to prune
outgroup taxa from the phylogeny. To prune tips from a phylogeny, see the
prune\_tree function.

Calculate DVMC in a tree following Liu et al., PNAS (2017), doi: 10.1073/pnas.1616744114.

```
phykit degree_of_violation_of_a_molecular_clock <tree>
```

Options:   
*<tree>*: input file tree name

### Evolutionary rate[](#evolutionary-rate "Link to this heading")

Function names: evolutionary\_rate, evo\_rate   
Command line interface: pk\_evolutionary\_rate, pk\_evo\_rate

Calculate a tree-based estimation of the evolutionary rate of a gene.

Evolutionary rate is the total tree length divided by the number
of terminals.

Calculate evolutionary rate following Telford et al., Proceedings
of the Royal Society B (2014).

```
phykit evolutionary_rate <tree>
```

Options:   
*<tree>*: input file tree name

### Hidden paralogy check[](#hidden-paralogy-check "Link to this heading")

Function names: hidden\_paralogy\_check, clan\_check   
Command line interface: pk\_hidden\_paralogy\_check, pk\_clan\_check

Scan tree for evidence of hidden paralogy.

This analysis can be used to identify hidden paralogy.
Specifically, this method will examine if a set of
well known monophyletic taxa are, in fact, monophyletic.
If they are not, the evolutionary history of the gene may
be subject to hidden paralogy. This analysis is typically
done with single-copy orthologous genes.

Requires a clade file, which species which monophyletic
lineages to check for. Multiple monophyletic
lineages can be specified. Each lineage should
be specified on a single line and each tip name
(or taxon name) should be separated by a space.
For example, if it is anticipated that tips
“A”, “B”, and “C” are monophyletic and “D”,
“E”, and “F” are expected to be monophyletic, the
clade file should be formatted as follows:   
“   
A B C   
D E F   
“

The output will report if the specified taxa were monophyletic
or not. The number of rows will reflect how many groups of taxa
were checked for monophyly. For example,
if there were three rows of clades in the -c file, there will be
three rows in the output
where the first row in the output corresponds to the
results of the first row in the clade file.

The concept behind this analysis follows
Siu-Ting et al., Molecular Biology and Evolution (2019),
doi: 10.1093/molbev/msz067.

```
phykit hidden_paralogy_check <tree> -c/--clade <clade_file>
```

Options:   
*-t/\-\-tree*: input file tree name
*-c/\-\-clade*: clade file detailing which monophyletic lineages should
be scanned for

### Internal branch statistics[](#internal-branch-statistics "Link to this heading")

Function names: internal\_branch\_stats; ibs   
Command line interface: pk\_internal\_branch\_stats; pk\_ibs

Calculate summary statistics for internal branch lengths in a phylogeny.

Internal branch lengths can be useful for phylogeny diagnostics.

To obtain all internal branch lengths, use the -v/\-\-verbose option.

```
phykit internal_branch_stats <tree> [-v/--verbose]
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*-v/\-\-verbose*: optional argument to print all internal branch lengths

### Internode labeler[](#internode-labeler "Link to this heading")

Function names: internode\_labeler; il   
Command line interface: pk\_internode\_labeler; pk\_il

Appends numerical identifiers to bipartitions in place of support values.
This is helpful for pointing to specific internodes in supplementary files
or otherwise.

```
phykit internode_labeler <tree> [-o/--output <file>]
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*-o/\-\-output*: optional argument to name the outputted tree file

### Last common ancestor subtree[](#last-common-ancestor-subtree "Link to this heading")

Function names: last\_common\_ancestor\_subtree; lca\_subtree   
Command line interface: pk\_last\_common\_ancestor\_subtree; pk\_lca\_subtree

Obtains subtree from a phylogeny by getting the last common ancestor
from a list of taxa.

```
phykit last_common_ancestor_subtree <file> <list_of_taxa> [-o/--output <file>]
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*<list\_of\_taxa>*: second argument after function name should be a single column
file with the list of taxa to get the last common ancestor subtree for
*-o/\-\-output*: optional argument to print all LB score values

### Long branch score[](#long-branch-score "Link to this heading")

Function names: long\_branch\_score; lb\_score; lbs   
Command line interface: pk\_long\_branch\_score; pk\_lb\_score; pk\_lbs

Calculate long branch (LB) scores in a phylogeny.

Lower LB scores are thought to be desirable because
they are indicative of taxa or trees that likely do
not have issues with long branch attraction.

LB score is the mean pairwise patristic distance of
taxon i compared to all other taxa over the average
pairwise patristic distance.

PhyKIT reports summary statistics. To obtain LB scores
for each taxa, use the -v/–verbose option.

LB scores are calculated following Struck, Evolutionary
Bioinformatics (2014), doi: 10.4137/EBO.S14239.

```
phykit long_branch_score <tree> [-v/--verbose]
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*-v/\-\-verbose*: optional argument to print all LB score values

### Monophyly check[](#monophyly-check "Link to this heading")

Function names: monophyly\_check; is\_monophyletic   
Command line interface: pk\_monophyly\_check; pk\_is\_monophyletic

This analysis can be used to determine if a set of
taxa are exclusively monophyletic. By exclusively monophyletic,
if other taxa are in the same clade, the lineage will not be
considered exclusively monophyletic.

Requires a taxa file, which species which tip names
are expected to be monophyletic. File format is a
single column file with tip names. Tip names not
present in the tree will not be considered when
examining monophyly.

The output will have six columns.
col 1: if the clade was or wasn’t monophyletic
col 2: average bipartition support value in the clade of interest
col 3: maximum bipartition support value in the clade of interest
col 4: minimum bipartition support value in the clade of interest
col 5: standard deviation of bipartition support values in the clade of interest
col 6: tip names of taxa monophyletic with the lineage of interest excluding those that are listed in the taxa\_of\_interest file

```
phykit monophyly_check <tree> <list_of_taxa>
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*<list\_of\_taxa>*: single column file with list of tip names to
examine the monophyly of

### Nearest neighbor interchange[](#nearest-neighbor-interchange "Link to this heading")

Function names: nearest\_neighbor\_interchange; nni   
Command line interface: pk\_nearest\_neighbor\_interchange; pk\_nni

Generate all nearest neighbor interchange moves for a binary
rooted tree.

By default, the output file will have the same name as the input
file but with the suffix “.nnis”

The output file will also include the original phylogeny.

```
phykit nearest_neighbor_interchange <tree> [-o/--output]
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*-o/\-\-output*: optional argument to specify output file name

### Patristic distances[](#patristic-distances "Link to this heading")

Function names: patristic\_distances; pd   
Command line interface: pk\_patristic\_distances; pk\_pd

Calculate summary statistics among patristic distances in a phylogeny.

Patristic distances are all tip-to-tip distances in a phylogeny.

To obtain all patristic distances, use the -v/–verbose option.
With the -v option, the first column will have two taxon names
separated by a ‘-’ followed by the patristic distance. Features
will be tab separated.

```
phykit patristic_distances <tree> [-v/--verbose]
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*-o/\-\-output*: optional argument to print all tip-to-tip distances

### Polytomy testing[](#polytomy-testing "Link to this heading")

Function names: polytomy\_test; polyt\_test; polyt; ptt   
Command line interface: pk\_polytomy\_test; pk\_polyt\_test; pk\_polyt; pk\_ptt

Conduct a polytomy test for three clades in a phylogeny.

Polytomy tests can be used to identify putative radiations
as well as identify well supported alternative topologies.

The polytomy testing function takes as input a file with
the three groups of taxa to test the relationships for and
a single column file with the names of the desired tree files
to use for polytomy testing. Next, the script to examine
support for the grouping of the three taxa using triplets
and gene support frequencies.

This function can account for uncertainty in gene trees -
that is, the input phylogenies can have collapsed bipartitions.

Thereafter, a chi-squared test is conducted to determine if there
is evidence to reject the null hypothesis wherein the null
hypothesis is that the three possible topologies among the three
groups are equally supported. This test is done using gene support
frequencies.

```
phykit polytomy_test -t/--trees <trees> -g/--groups <groups>
```

Options:   
*-t/\-\-trees <trees>*: single column file with the names of
phylogenies to use for polytomy testing   
*-g/\-\-groups*: a tab-delimited file with the grouping designations
to test. Lines starting with commetns are not considered. Names of
individual taxa should be separated by a semi-colon ‘;’

For example, the groups file could look like the following:

```
#label group0  group1  group2
name_of_test    tip_name_A;tip_name_B   tip_name_C  tip_name_D;tip_name_E
```

### Print tree[](#print-tree "Link to this heading")

Function names: print\_tree; print; pt   
Command line interface: pk\_print\_tree; pk\_print; pk\_pt

Print ascii tree of input phylogeny.

Phylogeny can be printed with or without branch lengths.
By default, the phylogeny will be printed with branch lengths
but branch lengths can be removed using the -r/–remove argument.

```
phykit print_tree <tree> [-r/--remove]
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*-r/\-\-remove*: optional argument to print the phylogeny without branch
lengths

### Prune tree[](#prune-tree "Link to this heading")

Function names: prune\_tree; prune   
Command line interface: pk\_prune\_tree; pk\_prune

Prune tips from a phylogeny.

Provide a single column file with the names of the tips
in the input phylogeny you would like to prune from the
tree.

```
phykit prune_tree <tree> <list_of_taxa> [-o/--output <output_file>
-k/--keep]
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*<list\_of\_taxa>*: single column file with the names of the tips to remove
from the phylogeny   
*-o/\-\-output*: name of output file for the pruned phylogeny.
Default output will have the same name as the input file but with the suffix
“.pruned”
*-k/–keep*: optional argument. If used instead of pruning taxa in <list\_of\_taxa>,
keep them
|

### Rename tree tips[](#rename-tree-tips "Link to this heading")

Function names: rename\_tree; rename\_tips   
Command line interface: pk\_rename\_tree; pk\_rename\_tips

Renames tips in a phylogeny.

Renaming tip files will follow the scheme of a tab-delimited
file wherein the first column is the current tip name and the
second column is the desired tip name in the resulting
phylogeny.

```
phykit rename_tree_tips <tree> -i/--idmap <idmap.txt> [-o/--output <output_file>]
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*-i/\-\-idmap*: identifier map of current tip names (col1) and desired
tip names (col2)   
*-o/\-\-output*: optional argument to write the renamed tree files to. Default
output will have the same name as the input file but with the suffix “.renamed”

### Robinson-Foulds distance[](#robinson-foulds-distance "Link to this heading")

Function names: robinson\_foulds\_distance; rf\_distance; rf\_dist; rf   
Command line interface: pk\_robinson\_foulds\_distance; pk\_rf\_distance; pk\_rf\_dist; pk\_rf

Calculate Robinson-Foulds (RF) distance between two trees.

Low RF distances reflect greater similarity between two phylogenies.
This function prints out two values, the plain RF value and the
normalized RF value, which are separated by a tab. Normalized RF values
are calculated by taking the plain RF value and dividing it by 2(n-3)
where n is the number of tips in the phylogeny. Prior to calculating
an RF value, PhyKIT will first determine the number of shared tips
between the two input phylogenies and prune them to a common set of
tips. Thus, users can input trees with different topologies and
infer an RF value among subtrees with shared tips.

PhyKIT will print out
col 1; the plain RF distance and
col 2: the normalized RF distance.

RF distances are calculated following Robinson & Foulds, Mathematical
Biosciences (1981), doi: 10.1016/0025-5564(81)90043-2.

```
phykit robinson_foulds_distance <tree_file_zero> <tree_file_one>
```

Options:   
*<tree\_file\_zero>*: first argument after function name should be a tree file
*<tree\_file\_one>*: second argument after function name should be a tree file

### Root tree[](#root-tree "Link to this heading")

Function names: root\_tree; root; rt   
Command line interface: pk\_root\_tree; pk\_root; pk\_rt

Roots phylogeny using user-specified taxa.

A list of taxa to root the phylogeny on should be specified using the -r
argument. The root\_taxa file should be a single-column file with taxa names.
The outputted file will have the same name as the inputted tree file but with
the suffix “.rooted”.

```
phykit root_tree <tree> -r/--root <root_taxa> [-o/--output <output_file>]
```

Options:   
*<tree>*: first argument after function name should be a tree file to root|br|
*-r/\-\-root*: single column file with taxa names to root the phylogeny on|br|
*-o/\-\-output*: optional argument to specify the name of the output file

### Spurious homolog identification[](#spurious-homolog-identification "Link to this heading")

Function names: spurious\_sequence; spurious\_seq; ss   
Command line interface: pk\_spurious\_sequence; pk\_spurious\_seq; pk\_ss

Determines potentially spurious homologs using branch lengths.

Identifies potentially spurious sequences and reports
tips in the phylogeny that could possibly be removed
from the associated multiple sequence alignment. PhyKIT
does so by identifying and reporting long terminal branches
defined as branches that are equal to or 20 times the median
length of all branches.

PhyKIT reports the following information
col1: name of tip that is a putatively spurious sequence
col2: length of branch leading to putatively spurious sequence
col3: threshold used to identify putatively spurious sequences
col4: median branch length in the phylogeny

If there are no putatively spurious sequences, “None” is reported.

Using this method to identify potentially spurious sequences
was, to my knowledge, first introduced by Shen et al., (2018)
Cell doi: 10.1016/j.cell.2018.10.023.

```
phykit spurious_seq <file> -f/\\-\\-factor
```

Options:   
*<file>*: first argument after function name should be a tree file
*-f/\-\-factor*: factor to multiply median branch length by to calculate
the threshold of long branches. (Default: 20)

### Terminal branch statistics[](#terminal-branch-statistics "Link to this heading")

Function names: terminal\_branch\_stats; tbs   
Command line interface: pk\_terminal\_branch\_stats; pk\_tbs

Calculate summary statistics for terminal branch lengths in a phylogeny.

Terminal branch lengths can be useful for phylogeny diagnostics.

To obtain all terminal branch lengths, use the -v/\-\-verbose option.

```
phykit terminal_branch_stats <tree> [-v/--verbose]
```

Options:   
*<tree>*: first argument after function name should be a tree file   
*-v/\-\-verbose*: optional argument to print all terminal branch lengths

### Tip labels[](#tip-labels "Link to this heading")

Function names: tip\_labels; tree\_labels; labels; tl   
Command line interface: pk\_tip\_labels; pk\_tree\_labels; pk\_labels; pk\_tl

Prints the tip labels (or names) a phylogeny.

```
phykit tip_labels <tree>
```

Options:   
*<tree>*: first argument after function name should be a tree file

### Tip-to-tip distance[](#tip-to-tip-distance "Link to this heading")

Function names: tip\_to\_tip\_distance; t2t\_dist; t2t   
Command line interface: pk\_tip\_to\_tip\_distance; pk\_t2t\_dist; pk\_t2t

Calculate distance between two tips (or leaves) in a phylogeny.

Distances are in substitutions per site.

```
phykit tip_to_tip_distance <tree_file> <tip_1> <tip_2>
```

Options:   
*<tree\_file>*: first argument after function name should be a tree file   
*<tip\_1>*: second argument should be the name of the first tip of interest   
*<tip\_2>*: third argument should be the name of the second tip of interest

### Tip-to-tip node distance[](#tip-to-tip-node-distance "Link to this heading")

Function names: tip\_to\_tip\_node\_distance; t2t\_node\_dist; t2t\_nd   
Command line interface: pk\_tip\_to\_tip\_node\_distance; pk\_t2t\_node\_dist; pk\_t2t\_nd

Calculate distance between two tips (or leaves) in a phylogeny.

Distance is measured by the number of nodes between one tip
and another.

```
phykit tip_to_tip_node_distance <tree_file> <tip_1> <tip_2>
```

Options:   
*<tree\_file>*: first argument after function name should be a tree file   
*<tip\_1>*: second argument should be the name of the first tip of interest   
*<tip\_2>*: third argument should be the name of the second tip of interest

### Total tree length[](#total-tree-length "Link to this heading")

Function names: total\_tree\_length; tree\_len   
Command line interface: pk\_total\_tree\_length; pk\_tree\_len

Calculate total tree length, which is a sum of all branches.

```
phykit total_tree_length <tree>
```

Options:   
*<tree>*: first argument after function name should be a tree file

### Treeness[](#treeness "Link to this heading")

Function names: treeness; tness   
Command line interface: pk\_treeness; pk\_tness

Calculate treeness statistic for a phylogeny.

Higher treeness values are thought to be desirable because they
represent a higher signal-to-noise ratio.

Treeness is the sum of internal branch lengths divided by the total
tree length. Therefore, values range from 0 to 1. Treeness can be
used as a measure of the signal-to-noise ratio in a phylogeny.

Calculate treeness (also referred to as stemminess) following
Lanyon, The Auk (1988), doi: 10.1093/auk/105.3.565 and
Phillips and Penny, Molecular Phylogenetics and Evolution
(2003), doi: 10.1016/S1055-7903(03)00057-5.

```
phykit treeness <tree>
```

Options:   
*<tree>*: first argument after function name should be a tree file