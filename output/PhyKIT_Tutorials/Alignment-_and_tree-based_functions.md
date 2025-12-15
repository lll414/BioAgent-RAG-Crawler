# Alignment- and tree-based functions

Source URL: https://jlsteenwyk.com/PhyKIT/usage/index.html
Date Scraped: 2025-12-11

---

---

### Saturation[](#saturation "Link to this heading")

Function names: saturation; sat   
Command line interface: pk\_saturation; pk\_sat

Calculate saturation for a given tree and alignment.

Saturation is defined as sequences in multiple sequence
alignments that have undergone numerous substitutions such
that the distances between taxa are underestimated.

Data with no saturation will have a value of 1. The closer
the value is to 1, the less saturated the data.

This function outputs two values (as of v1.19.9). The first
value is the saturation value and the second column is the absolute
value of saturation minus 1. Thus, lower values in the second column
are indicative of values closer to one and, thus, less saturation.

Saturation is calculated following Philippe et al., PLoS
Biology (2011), doi: 10.1371/journal.pbio.1000602.

```
phykit saturation -a <alignment> -t <tree> [-v/--verbose]
```

Options:   
*-a/\-\-alignment*: an alignment file   
*-t/\-\-tree*: a tree file   
*-e/–exclude\_gaps*: if a site has a gap, ignore it   
*-v/\-\-verbose*: print out patristic distances and uncorrected   
distances used to determine saturation

### Treeness over RCV[](#treeness-over-rcv "Link to this heading")

Function names: treeness\_over\_rcv; toverr; tor   
Command line interface: pk\_treeness\_over\_rcv; pk\_toverr; pk\_tor

Calculate treeness/RCV for a given alignment and tree.

Higher treeness/RCV values are thought to be desirable because
they harbor a high signal-to-noise ratio are least susceptible
to composition bias.

PhyKIT reports three tab delimited values:
col1: treeness/RCV
col2: treeness
col3: RCV

Calculate treeness/RCV following Phillips and Penny, Molecular
Phylogenetics and Evolution (2003), doi: 10.1016/S1055-7903(03)00057-5.

```
phykit treeness_over_rcv -a/--alignment <alignment> -t/--tree <tree>
```

Options:   
*-a/\-\-alignment*: an alignment file   
*-t/\-\-tree*: a tree file