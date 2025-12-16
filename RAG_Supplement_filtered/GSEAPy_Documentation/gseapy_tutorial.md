# 5.A Protocol to Prepare files for GSEApy

Source URL: https://gseapy.readthedocs.io/en/latest/gseapy_tutorial.html
Date Scraped: 2025-12-11

---

As a biological researcher, I like protocols.

Here is a short tutorial for you to walk you through gseapy.

For file format explanation, please see [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats)

In order to run gseapy successfully, install gseapy use pip.

```
pip install gseapy

# if you have conda
conda install -c bioconda gseapy
```

## 5.1. Use `gsea` command, or `gsea()`[](#use-gsea-command-or-gsea "Link to this heading")

Follow the steps blow.

One thing you should know is that the gseapy input files are the same as
`GSEA` desktop required. You can use these files below to run `GSEA` desktop, too.

### 5.1.1. Prepare an tabular text file of gene expression like this:[](#prepare-an-tabular-text-file-of-gene-expression-like-this "Link to this heading")

**RNA-seq,ChIP-seq, Microarry data** are all supported.

Here is to see what the structure of expression table looks like

```
import pandas as pd
df = pd.read_table('./test/gsea_data.txt')
df.head()

#or assign dataframe to the parameter 'data'
```

|  | Gene\_Symbol | ACD2 | BCD2 | CCD2 | APYD2 | BPYD2 | CPYD2 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | SCARA3 | 6.287 | 6.821 | 6.005 | 2.525 | 1.911 | 1.500 |
| 1 | POU5F1 | 11.168 | 11.983 | 10.469 | 7.795 | 7.911 | 6.174 |
| 2 | CTLA2B | 4.362 | 5.708 | 4.633 | 1.493 | 0.000 | 1.369 |
| 3 | CRYAB | 11.339 | 11.662 | 11.714 | 7.698 | 7.928 | 7.779 |
| 4 | PMP22 | 7.259 | 7.548 | 6.803 | 4.418 | 2.239 | 3.071 |

### 5.1.2. An cls file is also expected.[](#an-cls-file-is-also-expected "Link to this heading")

This file is used to specify column attributes in step 1, just like `GSEA` asked.

An example of cls file looks like below.

```
with open('gsea/edb/C1OE.cls') as cls:
    print(cls.read())

# or assign a list object to parameter 'cls' like this
# cls=['C1OE', 'C1OE', 'C1OE', 'Vector', 'Vector', 'Vector']
```

```
6 2 1
# C1OE Vector
C1OE C1OE C1OE Vector Vector Vector
```

The first line specify the total samples and phenotype numbers. Leave number 1 always be 1.

The second line specify the phenotype class(name).

The third line specify column attributes in step 1.

So you could prepare the cls file in python like this

```
groups = ['C1OE', 'C1OE', 'C1OE', 'Vector', 'Vector', 'Vector']
with open('gsea/edb/C1OE.cls', "w") as cl:
   line = f"{len(groups)} 2 1\n# C10E Vector\n"
   cl.write(line)
   cl.write(" ".join(groups) + "\n")
```

### 5.1.3. Gene\_sets file in gmt format.[](#gene-sets-file-in-gmt-format "Link to this heading")

All you need to do is to download gene set database file from `GSEA` or `Enrichr` website.

Or you could use enrichr library. In this case, just provide library name to parameter ‘gene\_sets’

If you would like to use you own gene\_sets.gmts files, build such a file use excel:

An example of gmt file looks like below:

```
with open('gsea/edb/gene_sets.gmt') as gmt:
    print(gmt.read())
```

```
ES-SPECIFIC Arid3a_used     ACTA1   CALML4  CORO1A  DHX58   DPYS    EGR1    ESRRB   GLI2    GPX2    HCK     INHBB
HDAC-UNIQUE     Arid3a_used 1700017B05RIK   8430427H17RIK   ABCA3   ANKRD44 ARL4A   BNC2    CLDN3
XEN-SPECIFIC        Arid3a_used     1110036O03RIK   A130022J15RIK   B2M     B3GALNT1        CBX4    CITED1  CLU     CTSH    CYP26A1
GATA-SPECIFIC       Arid3a_used     1200009I06RIK   5430407P10RIK   BAIAP2L1        BMP8B   CITED1  CLDN3   COBLL1  CORO1A  CRYAB   CTDSPL  DKKL1
TS-SPECIFIC Arid3a_used     5430407P10RIK   AFAP1L1 AHNAK   ANXA2   ANXA3   ANXA5   B2M     BIK     BMP8B   CAMK1D  CBX4    CLDN3   CSRP1   DKKL1   DSC2
```

## 5.2. Use `enrichr` command, or `enrichr()`[](#use-enrichr-command-or-enrichr "Link to this heading")

The only thing you need to prepare is a gene list file.

**Note**: Enrichr uses a list of Entrez gene symbols as input.

For `enrichr` , you could assign a list object

```
# assign a list object to enrichr
l = ['SCARA3', 'LOC100044683', 'CMBL', 'CLIC6', 'IL13RA1', 'TACSTD2', 'DKKL1', 'CSF1',
     'SYNPO2L', 'TINAGL1', 'PTX3', 'BGN', 'HERC1', 'EFNA1', 'CIB2', 'PMP22', 'TMEM173']

gseapy.enrichr(gene_list=l, gene_sets='KEGG_2016', outfile='test')
```

or a gene list file in txt format(one gene id per row)

```
gseapy.enrichr(gene_list='gene_list.txt',  gene_sets='KEGG_2016', outfile='test')
```

Let’s see what the txt file looks like.

```
with open('data/gene_list.txt') as genes:
    print(genes.read())
```

```
CTLA2B
SCARA3
LOC100044683
CMBL
CLIC6
IL13RA1
TACSTD2
DKKL1
CSF1
CITED1
SYNPO2L
TINAGL1
PTX3
```

Select the library you want to do enrichment analysis. To get a list of all available libraries, run

```
#s get_library_name(), it will print out all library names.
import gseapy
names = gseapy.get_library_name()
print(names)
```

```
 ['Genome_Browser_PWMs',
'TRANSFAC_and_JASPAR_PWMs',
'ChEA_2013',
'Drug_Perturbations_from_GEO_2014',
'ENCODE_TF_ChIP-seq_2014',
'BioCarta_2013',
'Reactome_2013',
'WikiPathways_2013',
'Disease_Signatures_from_GEO_up_2014',
'KEGG_2013',
'TF-LOF_Expression_from_GEO',
'TargetScan_microRNA',
'PPI_Hub_Proteins',
'GO_Molecular_Function_2015',
'GeneSigDB',
'Chromosome_Location',
'Human_Gene_Atlas',
'Mouse_Gene_Atlas',
'GO_Cellular_Component_2015',
'GO_Biological_Process_2015',
'Human_Phenotype_Ontology',
'Epigenomics_Roadmap_HM_ChIP-seq',
'KEA_2013',
'NURSA_Human_Endogenous_Complexome',
'CORUM',
'SILAC_Phosphoproteomics',
'MGI_Mammalian_Phenotype_Level_3',
'MGI_Mammalian_Phenotype_Level_4',
'Old_CMAP_up',
'Old_CMAP_down',
'OMIM_Disease',
'OMIM_Expanded',
'VirusMINT',
'MSigDB_Computational',
'MSigDB_Oncogenic_Signatures',
'Disease_Signatures_from_GEO_down_2014',
'Virus_Perturbations_from_GEO_up',
'Virus_Perturbations_from_GEO_down',
'Cancer_Cell_Line_Encyclopedia',
'NCI-60_Cancer_Cell_Lines',
'Tissue_Protein_Expression_from_ProteomicsDB',
'Tissue_Protein_Expression_from_Human_Proteome_Map',
'HMDB_Metabolites',
'Pfam_InterPro_Domains',
'GO_Biological_Process_2013',
'GO_Cellular_Component_2013',
'GO_Molecular_Function_2013',
'Allen_Brain_Atlas_up',
'ENCODE_TF_ChIP-seq_2015',
'ENCODE_Histone_Modifications_2015',
'Phosphatase_Substrates_from_DEPOD',
'Allen_Brain_Atlas_down',
'ENCODE_Histone_Modifications_2013',
'Achilles_fitness_increase',
'Achilles_fitness_decrease',
'MGI_Mammalian_Phenotype_2013',
'BioCarta_2015',
'HumanCyc_2015',
'KEGG_2015',
'NCI-Nature_2015',
'Panther_2015',
'WikiPathways_2015',
'Reactome_2015',
'ESCAPE',
'HomoloGene',
'Disease_Perturbations_from_GEO_down',
'Disease_Perturbations_from_GEO_up',
'Drug_Perturbations_from_GEO_down',
'Genes_Associated_with_NIH_Grants',
'Drug_Perturbations_from_GEO_up',
'KEA_2015',
'Single_Gene_Perturbations_from_GEO_up',
'Single_Gene_Perturbations_from_GEO_down',
'ChEA_2015',
'dbGaP',
'LINCS_L1000_Chem_Pert_up',
'LINCS_L1000_Chem_Pert_down',
'GTEx_Tissue_Sample_Gene_Expression_Profiles_down',
'GTEx_Tissue_Sample_Gene_Expression_Profiles_up',
'Ligand_Perturbations_from_GEO_down',
'Aging_Perturbations_from_GEO_down',
'Aging_Perturbations_from_GEO_up',
'Ligand_Perturbations_from_GEO_up',
'MCF7_Perturbations_from_GEO_down',
'MCF7_Perturbations_from_GEO_up',
'Microbe_Perturbations_from_GEO_down',
'Microbe_Perturbations_from_GEO_up',
'LINCS_L1000_Ligand_Perturbations_down',
'LINCS_L1000_Ligand_Perturbations_up',
'LINCS_L1000_Kinase_Perturbations_down',
'LINCS_L1000_Kinase_Perturbations_up',
'Reactome_2016',
'KEGG_2016',
'WikiPathways_2016',
'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
'Kinase_Perturbations_from_GEO_down',
'Kinase_Perturbations_from_GEO_up',
'BioCarta_2016',
'Humancyc_2016',
'NCI-Nature_2016',
'Panther_2016']
```

For more details, please track the official links: <http://amp.pharm.mssm.edu/Enrichr/>

## 5.3. Use `replot` Command, or `replot()`[](#use-replot-command-or-replot "Link to this heading")

You may also want to use `replot()` to reproduce `GSEA` desktop plots.

The only input of `replot()` is the directory of `GSEA` desktop output.

The input directory(e.g. gsea), must contained **edb** folder, gseapy need 4 data files
inside edb folder.The gsea document tree looks like this:

```
gsea
└─edb
    └─test.cls
    └─gene_sets.gmt
    └─gsea_data.rnk
    └─results.edb
```

After this, you can start to run gseapy.

```
import gseapy
gseapy.replot(indir ='gsea', outdir = 'gseapy_out')
```

If you prefer to run in command line, it’s more simple.

```
gseapy replot -i gsea -o gseapy_out
```

For advanced usage of library, see the [Developmental Guide](run.html#run).