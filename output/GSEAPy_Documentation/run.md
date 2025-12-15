# 6.Developmental Guide

Source URL: https://gseapy.readthedocs.io/en/latest/run.html
Date Scraped: 2025-12-11

---

## 6.1. Module APIs[](#module-apis "Link to this heading")

gseapy.gsea()[[source]](_modules/gseapy.html#gsea)[](#gseapy.gsea "Link to this definition")
:   Run Gene Set Enrichment Analysis.

    Parameters:
    :   * **data** – Gene expression data table, Pandas DataFrame, gct file.
        * **gene\_sets** – Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.
        * **cls** – A list or a .cls file format required for GSEA.
        * **outdir** (*str*) – Results output directory. If None, nothing will write to disk.
        * **permutation\_num** (*int*) – Number of permutations. Default: 1000.
          Minimial possible nominal p-value is about 1/nperm.
        * **permutation\_type** (*str*) – Type of permutation reshuffling,
          choose from {“phenotype”: ‘sample.labels’ , “gene\_set” : gene.labels}.
        * **min\_size** (*int*) – Minimum allowed number of genes from gene set also the data set. Default: 15.
        * **max\_size** (*int*) – Maximum allowed number of genes from gene set also the data set. Default: 500.
        * **weight** (*float*) – Refer to `algorithm.enrichment_score()`. Default:1.
        * **method** –

          The method used to calculate a correlation or ranking. Default: ‘signal\_to\_noise’.
          Others methods are:

          1. ’signal\_to\_noise’

             You must have at least three samples for each phenotype to use this metric.
             The larger the signal-to-noise ratio, the larger the differences of the means
             (scaled by the standard deviations); that is, the more distinct
             the gene expression is in each phenotype and the more the gene acts as a “class marker.”
          2. ’t\_test’

             Uses the difference of means scaled by the standard deviation and number of samples.
             Note: You must have at least three samples for each phenotype to use this metric.
             The larger the tTest ratio, the more distinct the gene expression is in each phenotype
             and the more the gene acts as a “class marker.”
          3. ’ratio\_of\_classes’ (also referred to as fold change).

             Uses the ratio of class means to calculate fold change for natural scale data.
          4. ’diff\_of\_classes’

             Uses the difference of class means to calculate fold change for nature scale data
          5. ’log2\_ratio\_of\_classes’

             Uses the log2 ratio of class means to calculate fold change for natural scale data.
             This is the recommended statistic for calculating fold change for log scale data.
        * **ascending** (*bool*) – Sorting order of rankings. Default: False.
        * **threads** (*int*) – Number of threads you are going to use. Default: 4.
        * **figsize** (*list*) – Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
        * **format** (*str*) – Matplotlib figure format. Default: ‘pdf’.
        * **graph\_num** (*int*) – Plot graphs for top sets of each phenotype.
        * **no\_plot** (*bool*) – If equals to True, no figure will be drawn. Default: False.
        * **seed** – Random seed. expect an integer. Default:None.
        * **verbose** (*bool*) – Bool, increase output verbosity, print out progress of your job, Default: False.

    Returns:
    :   Return a GSEA obj. All results store to a dictionary, obj.results,
        where contains:

        ```
        | {
        |  term: gene set name,
        |  es: enrichment score,
        |  nes: normalized enrichment score,
        |  pval:  Nominal p-value (from the null distribution of the gene set,
        |  fdr: FDR qvalue (adjusted False Discory Rate),
        |  fwerp: Family wise error rate p-values,
        |  tag %: Percent of gene set before running enrichment peak (ES),
        |  gene %: Percent of gene list before running enrichment peak (ES),
        |  lead_genes: leading edge genes (gene hits before running enrichment peak),
        |  matched genes: genes matched to the data,
        | }
        ```

gseapy.prerank()[[source]](_modules/gseapy.html#prerank)[](#gseapy.prerank "Link to this definition")
:   Run Gene Set Enrichment Analysis with pre-ranked correlation defined by user.

    Parameters:
    :   * **rnk** – pre-ranked correlation table or pandas DataFrame. Same input with `GSEA` .rnk file.
        * **gene\_sets** – Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.
        * **outdir** – results output directory. If None, nothing will write to disk.
        * **permutation\_num** (*int*) – Number of permutations. Default: 1000.
          Minimial possible nominal p-value is about 1/nperm.
        * **min\_size** (*int*) – Minimum allowed number of genes from gene set also the data set. Default: 15.
        * **max\_size** (*int*) – Maximum allowed number of genes from gene set also the data set. Defaults: 500.
        * **weight** (*str*) – Refer to `algorithm.enrichment_score()`. Default:1.
        * **ascending** (*bool*) – Sorting order of rankings. Default: False for descending. If None, do not sort the ranking.
        * **threads** (*int*) – Number of threads you are going to use. Default: 4.
        * **figsize** (*list*) – Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [6.5,6].
        * **format** (*str*) – Matplotlib figure format. Default: ‘pdf’.
        * **graph\_num** (*int*) – Plot graphs for top sets of each phenotype.
        * **no\_plot** (*bool*) – If equals to True, no figure will be drawn. Default: False.
        * **seed** – Random seed. expect an integer. Default:None.
        * **verbose** (*bool*) – Bool, increase output verbosity, print out progress of your job, Default: False.

    Returns:
    :   Return a Prerank obj. All results store to a dictionary, obj.results,
        where contains:

        ```
        | {
        |  term: gene set name,
        |  es: enrichment score,
        |  nes: normalized enrichment score,
        |  pval:  Nominal p-value (from the null distribution of the gene set,
        |  fdr: FDR qvalue (adjusted False Discory Rate),
        |  fwerp: Family wise error rate p-values,
        |  tag %: Percent of gene set before running enrichment peak (ES),
        |  gene %: Percent of gene list before running enrichment peak (ES),
        |  lead_genes: leading edge genes (gene hits before running enrichment peak),
        |  matched genes: genes matched to the data,
        | }
        ```

gseapy.ssgsea()[[source]](_modules/gseapy.html#ssgsea)[](#gseapy.ssgsea "Link to this definition")
:   Run Gene Set Enrichment Analysis with single sample GSEA tool

    Parameters:
    :   * **data** – Expression table, pd.Series, pd.DataFrame, GCT file, or .rnk file format.
        * **gene\_sets** – Enrichr Library name or .gmt gene sets file or dict of gene sets. Same input with GSEA.
        * **outdir** – Results output directory. If None, nothing will write to disk.
        * **sample\_norm\_method** (*str*) –

          Sample normalization method. Choose from {‘rank’, ‘log’, ‘log\_rank’, None}. Default: rank.
          this argument will be used for ordering genes.

          1. ’rank’: Rank your expression data, and transform by 10000\*rank\_dat/gene\_numbers
          2. ’log’ : Do not rank, but transform data by log(data + exp(1)), while data = data[data<1] =1.
          3. ’log\_rank’: Rank your expression data, and transform by log(10000\*rank\_dat/gene\_numbers+ exp(1))
          4. None or ‘custom’: Do nothing, and use your own rank value to calculate enrichment score.

    see here: <https://github.com/GSEA-MSigDB/ssGSEAProjection-gpmodule/blob/master/src/ssGSEAProjection.Library.R>, line 86

    Parameters:
    :   * **correl\_norm\_type** (*str*) –

          correlation normalization type. Choose from {‘rank’, ‘symrank’, ‘zscore’, None}. Default: rank.
          After ordering genes by sample\_norm\_method, further data transformed could be applied to get enrichment score.

          when weight == 0, sample\_norm\_method and correl\_norm\_type do not matter;
          when weight > 0, the combination of sample\_norm\_method and correl\_norm\_type
          dictate how the gene expression values in input data are transformed
          to obtain the score – use this setting with care (the transformations
          can skew scores towards +ve or -ve values)

          sample\_norm\_method will first transformed and rank original data. the data is named correl\_vector for each sample.
          then correl\_vector is transformed again by

          1. correl\_norm\_type is None or ‘rank’ : do nothing, genes are weighted by actual correl\_vector.
          2. correl\_norm\_type ==’symrank’: symmetric ranking.
          3. correl\_norm\_type ==’zscore’: standardizes the correl\_vector before using them to calculate scores.
        * **min\_size** (*int*) – Minimum allowed number of genes from gene set also the data set. Default: 15.
        * **max\_size** (*int*) – Maximum allowed number of genes from gene set also the data set. Default: 2000.
        * **permutation\_num** (*int*) – For ssGSEA, default is 0.
          However, if you try to use ssgsea method to get pval and fdr, set to an interger.
        * **weight** (*str*) – Refer to `algorithm.enrichment_score()`. Default:0.25.
        * **ascending** (*bool*) – Sorting order of rankings. Default: False.
        * **threads** (*int*) – Number of threads you are going to use. Default: 4.
        * **figsize** (*list*) – Matplotlib figsize, accept a tuple or list, e.g. [width,height]. Default: [7,6].
        * **format** (*str*) – Matplotlib figure format. Default: ‘pdf’.
        * **graph\_num** (*int*) – Plot graphs for top sets of each phenotype.
        * **no\_plot** (*bool*) – If equals to True, no figure will be drawn. Default: False.
        * **seed** – Random seed. expect an integer. Default:None.
        * **verbose** (*bool*) – Bool, increase output verbosity, print out progress of your job, Default: False.

    Returns:
    :   Return a ssGSEA obj.
        All results store to a dictionary, access enrichment score or
        normalized enrichment score by obj.res2d or obj.results.
        if permutation\_num > 0, additional results contain:

        ```
        | {
        |  term: gene set name,
        |  es: enrichment score,
        |  nes: normalized enrichment score,
        |  pval:  Nominal p-value (from the null distribution of the gene set (if permutation_num > 0),
        |  fdr: FDR qvalue (adjusted FDR) (if permutation_num > 0),
        |  fwerp: Family wise error rate p-values (if permutation_num > 0),
        |  tag %: Percent of gene set before running enrichment peak (ES),
        |  gene %: Percent of gene list before running enrichment peak (ES),
        |  lead_genes: leading edge genes (gene hits before running enrichment peak),
        |  matched genes: genes matched to the data,
        | }
        ```

gseapy.enrichr()[[source]](_modules/gseapy.html#enrichr)[](#gseapy.enrichr "Link to this definition")
:   Enrichr API.

    Parameters:
    :   * **gene\_list** – str, list, tuple, series, dataframe. Also support input txt file with one gene id per row.
          The input identifier should be the same type to gene\_sets.
        * **gene\_sets** –

          str, list, tuple of Enrichr Library name(s).
          or custom defined gene\_sets (dict, or gmt file).

          Examples:

          Input Enrichr Libraries (<https://maayanlab.cloud/Enrichr/#stats>):
          :   str: ‘KEGG\_2016’
              list: [‘KEGG\_2016’,’KEGG\_2013’]
              Use comma to separate each other, e.g. “KEGG\_2016,huMAP,GO\_Biological\_Process\_2018”

          Input custom files:
          :   dict: gene\_sets={‘A’:[‘gene1’, ‘gene2’,…],
              :   ’B’:[‘gene2’, ‘gene4’,…], …}

              gmt: “genes.gmt”

          see also the online docs:
          <https://gseapy.readthedocs.io/en/latest/gseapy_example.html#2.-Enrichr-Example>
        * **organism** –

          Enrichr supported organism. Select from (human, mouse, yeast, fly, fish, worm).
          This argument only affects the Enrichr library names you’ve chosen.
          No any affects to gmt or dict input of gene\_sets.

          see here for more details: <https://maayanlab.cloud/modEnrichr/>.
        * **outdir** – Output file directory
        * **background** –

          int, list, str.
          Background genes. This argument works only if gene\_sets has a type Dict or gmt file.
          If your input are just Enrichr library names, this argument will be ignored.

          However, this argument is not straightforward when gene\_sets is given a custom input (a gmt file or dict).

          By default, all genes listed in the gene\_sets input will be used as background.

          There are 3 ways to tune this argument:

          1. (Recommended) Input a list of background genes: [‘gene1’, ‘gene2’,…]
             The background gene list is defined by your experment. e.g. the expressed genes in your RNA-seq.
             The gene identifer in gmt/dict should be the same type to the backgound genes.
          2. Specify a number: e.g. 20000. (the number of total expressed genes).
             This works, but not recommend. It assumes that all your genes could be found in background.
             If genes exist in gmt but not included in background provided,
             they will affect the significance of the statistical test.
          3. Set a Biomart dataset name: e.g. “hsapiens\_gene\_ensembl”
             The background will be all annotated genes from the BioMart datasets you’ve choosen.
             The program will try to retrieve the background information automatically.

             Enrichr module use the code below to get the background genes:
             :   ```
                 >>> from gseapy.parser import Biomart
                 >>> bm = Biomart()
                 >>> df = bm.query(dataset=background, #  e.g. 'hsapiens_gene_ensembl'
                              attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'],
                              filename=f'~/.cache/gseapy/{background}.background.genes.txt')
                 >>> df.dropna(subset=["entrezgene_id"], inplace=True)
                 ```

             So only genes with entrezid above will be the background genes if not input specify by user.
        * **cutoff** – Show enriched terms which Adjusted P-value < cutoff.
          Only affects the output figure, not the final output file. Default: 0.05
        * **format** – Output figure format supported by matplotlib,(‘pdf’,’png’,’eps’…). Default: ‘pdf’.
        * **figsize** – Matplotlib figsize, accept a tuple or list, e.g. (width,height). Default: (6.5,6).
        * **no\_plot** (*bool*) – If equals to True, no figure will be drawn. Default: False.
        * **verbose** (*bool*) – Increase output verbosity, print out progress of your job, Default: False.

    Returns:
    :   An Enrichr object, which obj.res2d stores your last query, obj.results stores your all queries.

gseapy.enrich()[[source]](_modules/gseapy.html#enrich)[](#gseapy.enrich "Link to this definition")
:   Perform over-representation analysis (hypergeometric test).

    Parameters:
    :   * **gene\_list** – str, list, tuple, series, dataframe. Also support input txt file with one gene id per row.
          The input identifier should be the same type to gene\_sets.
        * **gene\_sets** –

          str, list, tuple of Enrichr Library name(s).
          or custom defined gene\_sets (dict, or gmt file).

          Examples:

          > dict: gene\_sets={‘A’:[‘gene1’, ‘gene2’,…],
          > :   ’B’:[‘gene2’, ‘gene4’,…], …}
          >
          > gmt: “genes.gmt”
        * **outdir** – Output file directory
        * **background** –

          None | int | list | str.
          Background genes. This argument works only if gene\_sets has a type Dict or gmt file.

          However, this argument is not straightforward when gene\_sets is given a custom input (a gmt file or dict).

          By default, all genes listed in the gene\_sets input will be used as background.

          There are 3 ways to tune this argument:

          1. (Recommended) Input a list of background genes: [‘gene1’, ‘gene2’,…]
             The background gene list is defined by your experment. e.g. the expressed genes in your RNA-seq.
             The gene identifer in gmt/dict should be the same type to the backgound genes.
          2. Specify a number: e.g. 20000. (the number of total expressed genes).
             This works, but not recommend. It assumes that all your genes could be found in background.
             If genes exist in gmt but not included in background provided,
             they will affect the significance of the statistical test.
          3. Set a Biomart dataset name: e.g. “hsapiens\_gene\_ensembl”
             The background will be all annotated genes from the BioMart datasets you’ve choosen.
             The program will try to retrieve the background information automatically.

             Enrichr module use the code below to get the background genes:
             :   ```
                 >>> from gseapy.parser import Biomart
                 >>> bm = Biomart()
                 >>> df = bm.query(dataset=background, #  e.g. 'hsapiens_gene_ensembl'
                              attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id'],
                              filename=f'~/.cache/gseapy/{background}.background.genes.txt')
                 >>> df.dropna(subset=["entrezgene_id"], inplace=True)
                 ```

             So only genes with entrezid above will be the background genes if not input specify by user.
        * **cutoff** – Show enriched terms which Adjusted P-value < cutoff.
          Only affects the output figure, not the final output file. Default: 0.05
        * **format** – Output figure format supported by matplotlib,(‘pdf’,’png’,’eps’…). Default: ‘pdf’.
        * **figsize** – Matplotlib figsize, accept a tuple or list, e.g. (width,height). Default: (6.5,6).
        * **no\_plot** (*bool*) – If equals to True, no figure will be drawn. Default: False.
        * **verbose** (*bool*) – Increase output verbosity, print out progress of your job, Default: False.

    Returns:
    :   An Enrichr object, which obj.res2d stores your last query, obj.results stores your all queries.

gseapy.replot()[[source]](_modules/gseapy.html#replot)[](#gseapy.replot "Link to this definition")
:   The main function to reproduce GSEA desktop outputs.

    Parameters:
    :   * **indir** – GSEA desktop results directory. In the sub folder, you must contain edb file folder.
        * **outdir** – Output directory.
        * **weight** (*float*) – weighted score type. choose from {0,1,1.5,2}. Default: 1.
        * **figsize** (*list*) – Matplotlib output figure figsize. Default: [6.5,6].
        * **format** (*str*) – Matplotlib output figure format. Default: ‘pdf’.
        * **min\_size** (*int*) – Min size of input genes presented in Gene Sets. Default: 3.
        * **max\_size** (*int*) – Max size of input genes presented in Gene Sets. Default: 5000.
          You are not encouraged to use min\_size, or max\_size argument in [`replot()`](#gseapy.replot "gseapy.replot") function.
          Because gmt file has already been filtered.
        * **verbose** – Bool, increase output verbosity, print out progress of your job, Default: False.

    Returns:
    :   Generate new figures with selected figure format. Default: ‘pdf’.

## 6.2. GSEA Statistics[](#module-gseapy.gsea "Link to this heading")

*class* gseapy.gsea.GSEA(*data: DataFrame | str*, *gene\_sets: List[str] | str | Dict[str, str]*, *classes: List[str] | str | Dict[str, str]*, *outdir: str | None = None*, *min\_size: int = 15*, *max\_size: int = 500*, *permutation\_num: int = 1000*, *weight: float = 1.0*, *permutation\_type: str = 'phenotype'*, *method: str = 'signal\_to\_noise'*, *ascending: bool = False*, *threads: int = 1*, *figsize: Tuple[float, float] = (6.5, 6)*, *format: str = 'pdf'*, *graph\_num: int = 20*, *no\_plot: bool = False*, *seed: int = 123*, *verbose: bool = False*)[[source]](_modules/gseapy/gsea.html#GSEA)[](#gseapy.gsea.GSEA "Link to this definition")
:   GSEA main tool

    calc\_metric(*df: DataFrame*, *method: str*, *pos: str*, *neg: str*, *classes: Dict[str, str]*, *ascending: bool*) → Tuple[List[int], Series][[source]](_modules/gseapy/gsea.html#GSEA.calc_metric)[](#gseapy.gsea.GSEA.calc_metric "Link to this definition")
    :   The main function to rank an expression table. works for 2d array.

        Parameters:
        :   * **df** – gene\_expression DataFrame.
            * **method** –

              The method used to calculate a correlation or ranking. Default: ‘log2\_ratio\_of\_classes’.
              Others methods are:

              1. ’signal\_to\_noise’ (s2n) or ‘abs\_signal\_to\_noise’ (abs\_s2n)

                 > You must have at least three samples for each phenotype.
                 > The more distinct the gene expression is in each phenotype,
                 > the more the gene acts as a “class marker”.
              2. ’t\_test’

                 > Uses the difference of means scaled by the standard deviation and number of samples.
                 > Note: You must have at least three samples for each phenotype to use this metric.
                 > The larger the t-test ratio, the more distinct the gene expression is in each phenotype
                 > and the more the gene acts as a “class marker.”
              3. ’ratio\_of\_classes’ (also referred to as fold change).

                 > Uses the ratio of class means to calculate fold change for natural scale data.
              4. ’diff\_of\_classes’

                 > Uses the difference of class means to calculate fold change for natural scale data
              5. ’log2\_ratio\_of\_classes’

                 > Uses the log2 ratio of class means to calculate fold change for natural scale data.
                 > This is the recommended statistic for calculating fold change for log scale data.
            * **pos** (*str*) – one of labels of phenotype’s names.
            * **neg** (*str*) – one of labels of phenotype’s names.
            * **classes** (*dict*) – column id to group mapping.
            * **ascending** (*bool*) – bool or list of bool. Sort ascending vs. descending.

        Returns:
        :   returns argsort values of a tuple where
            0: argsort positions (indices)
            1: pd.Series of correlation value. Gene\_name is index, and value is rankings.

        visit here for more docs: <http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html>

    load\_classes(*classes: str | List[str] | Dict[str, Any]*)[[source]](_modules/gseapy/gsea.html#GSEA.load_classes)[](#gseapy.gsea.GSEA.load_classes "Link to this definition")
    :   Parse group (classes)

    load\_data() → Tuple[DataFrame, Dict][[source]](_modules/gseapy/gsea.html#GSEA.load_data)[](#gseapy.gsea.GSEA.load_data "Link to this definition")
    :   pre-processed the data frame.new filtering methods will be implement here.

    run()[[source]](_modules/gseapy/gsea.html#GSEA.run)[](#gseapy.gsea.GSEA.run "Link to this definition")
    :   GSEA main procedure

    to\_cls(*outdir: str*)[[source]](_modules/gseapy/gsea.html#GSEA.to_cls)[](#gseapy.gsea.GSEA.to_cls "Link to this definition")
    :   Save group information to cls file

*class* gseapy.gsea.Prerank(*rnk: DataFrame | Series | str*, *gene\_sets: List[str] | str | Dict[str, str]*, *outdir: str | None = None*, *pheno\_pos='Pos'*, *pheno\_neg='Neg'*, *min\_size: int = 15*, *max\_size: int = 500*, *permutation\_num: int = 1000*, *weight: float = 1.0*, *ascending: bool | None = False*, *threads: int = 1*, *figsize: Tuple[float, float] = (6.5, 6)*, *format: str = 'pdf'*, *graph\_num: int = 20*, *no\_plot: bool = False*, *seed: int = 123*, *verbose: bool = False*)[[source]](_modules/gseapy/gsea.html#Prerank)[](#gseapy.gsea.Prerank "Link to this definition")
:   GSEA prerank tool

    load\_ranking()[[source]](_modules/gseapy/gsea.html#Prerank.load_ranking)[](#gseapy.gsea.Prerank.load_ranking "Link to this definition")
    :   parse rnk input

    run()[[source]](_modules/gseapy/gsea.html#Prerank.run)[](#gseapy.gsea.Prerank.run "Link to this definition")
    :   GSEA prerank workflow

*class* gseapy.gsea.Replot(*indir: str*, *outdir: str = 'GSEApy\_Replot'*, *weight: float = 1.0*, *min\_size: int = 3*, *max\_size: int = 1000*, *figsize: Tuple[float, float] = (6.5, 6)*, *format: str = 'pdf'*, *verbose: bool = False*)[[source]](_modules/gseapy/gsea.html#Replot)[](#gseapy.gsea.Replot "Link to this definition")
:   To reproduce GSEA desktop output results.

    gsea\_edb\_parser(*results\_path*)[[source]](_modules/gseapy/gsea.html#Replot.gsea_edb_parser)[](#gseapy.gsea.Replot.gsea_edb_parser "Link to this definition")
    :   Parse results.edb file stored under **edb** file folder.

        Parameters:
        :   **results\_path** – the path of results.edb file.

        Returns:
        :   a dict contains { enrichment\_term: [es, nes, pval, fdr, fwer, hit\_ind]}

    run()[[source]](_modules/gseapy/gsea.html#Replot.run)[](#gseapy.gsea.Replot.run "Link to this definition")
    :   main replot function

*class* gseapy.base.GMT(*mapping: Dict[str, List[str]] | None = None*, *description: str | None = None*, *source: str | None = None*, *name: str | None = 'default'*)[[source]](_modules/gseapy/base.html#GMT)[](#gseapy.base.GMT "Link to this definition")
:   A collection of gene set dictionaries with metadata.

    Attributes:
    :   \_collections: Dict[str, Dict[str, Any]] - Stores gene set collections
        :   key: collection name
            value: {

            > ‘genes’: Dict[str, List[str]] - Gene set mappings
            > ‘description’: str - Collection description
            > ‘source’: str - Source of the gene sets

            }

    add(*mapping: Dict[str, List[str]]*, *description: str | None = None*, *source: str | None = None*, *name: str | None = 'default'*)[[source]](_modules/gseapy/base.html#GMT.add)[](#gseapy.base.GMT.add "Link to this definition")
    :   Add a gene set collection with metadata.

        Args:
        :   mapping: Gene set dictionary to add
            description: Description of the gene sets
            source: Source of the gene sets
            name: Name for this collection

    filter(*min\_size: int | None = None*, *max\_size: int | None = None*, *gene\_list: List[str] | None = None*, *collections: List[str] | None = None*) → [GMT](#gseapy.base.GMT "gseapy.base.GMT")[[source]](_modules/gseapy/base.html#GMT.filter)[](#gseapy.base.GMT.filter "Link to this definition")
    :   Filter gene sets based on size and gene membership.

        Args:
        :   min\_size: Minimum number of genes in a set
            max\_size: Maximum number of genes in a set
            gene\_list: Only keep genes present in this list
            collections: Only keep these named collections

        Returns:
        :   A new filtered GMT object

    get(*name: str = 'default'*) → Dict[str, List[str]][[source]](_modules/gseapy/base.html#GMT.get)[](#gseapy.base.GMT.get "Link to this definition")
    :   Get gene sets by collection name.

    get\_metadata(*name: str = 'default'*) → Dict[str, Any][[source]](_modules/gseapy/base.html#GMT.get_metadata)[](#gseapy.base.GMT.get_metadata "Link to this definition")
    :   Get metadata for a collection.

    items()[[source]](_modules/gseapy/base.html#GMT.items)[](#gseapy.base.GMT.items "Link to this definition")
    :   Iterate over (name, gene\_sets) pairs.

    *classmethod* read(*paths: str*, *source: str | None = None*) → [GMT](#gseapy.base.GMT "gseapy.base.GMT")[[source]](_modules/gseapy/base.html#GMT.read)[](#gseapy.base.GMT.read "Link to this definition")
    :   Read GMT files into a collection.

        Args:
        :   paths: Comma-separated list of GMT file paths
            source: Source annotation for the files

    write(*ofname: str*)[[source]](_modules/gseapy/base.html#GMT.write)[](#gseapy.base.GMT.write "Link to this definition")
    :   Write GMT file to disk.

*class* gseapy.base.GSEAbase(*outdir: str | None = None*, *gene\_sets: List[str] | str | Dict[str, str] = 'KEGG\_2016'*, *module: str = 'base'*, *threads: int = 1*, *enrichr\_url: str = 'http://maayanlab.cloud'*, *verbose: bool = False*)[[source]](_modules/gseapy/base.html#GSEAbase)[](#gseapy.base.GSEAbase "Link to this definition")
:   base class of GSEA.

    check\_uppercase(*gene\_list: List[str | int]*) → bool[[source]](_modules/gseapy/base.html#GSEAbase.check_uppercase)[](#gseapy.base.GSEAbase.check_uppercase "Link to this definition")
    :   Check whether a list of gene names are mostly in uppercase.

        ### 6. Parameters[](#parameters "Link to this heading")

        gene\_listlist, int
        :   A list of gene names or Entrez IDs

        ### 6. Returns[](#returns "Link to this heading")

        bool
        :   Whether the list of gene names are mostly in uppercase

    enrichment\_score(*gene\_list: Iterable[str]*, *correl\_vector: Iterable[float]*, *gene\_set: Dict[str, List[str]]*, *weight: float = 1.0*, *nperm: int = 1000*, *seed: int = 123*, *single: bool = False*, *scale: bool = False*)[[source]](_modules/gseapy/base.html#GSEAbase.enrichment_score)[](#gseapy.base.GSEAbase.enrichment_score "Link to this definition")
    :   This is the most important function of GSEApy. It has the same algorithm with GSEA and ssGSEA.

        Parameters:
        :   * **gene\_list** – The ordered gene list gene\_name\_list, rank\_metric.index.values
            * **gene\_set** – gene\_sets in gmt file, please use gmt\_parser to get gene\_set.
            * **weight** – It’s the same with gsea’s weighted\_score method. Weighting by the correlation
              is a very reasonable choice that allows significant gene sets with less than perfect coherence.
              options: 0(classic),1,1.5,2. default:1. if one is interested in penalizing sets for lack of
              coherence or to discover sets with any type of nonrandom distribution of tags, a value p < 1
              might be appropriate. On the other hand, if one uses sets with large number of genes and only
              a small subset of those is expected to be coherent, then one could consider using p > 1.
              Our recommendation is to use p = 1 and use other settings only if you are very experienced
              with the method and its behavior.
            * **correl\_vector** – A vector with the correlations (e.g. signal to noise scores) corresponding to the genes in
              the gene list. Or rankings, rank\_metric.values
            * **nperm** – Only use this parameter when computing esnull for statistical testing. Set the esnull value
              equal to the permutation number.
            * **seed** – Random state for initializing gene list shuffling. Default: seed=None

        Returns:

        ES: Enrichment score (real number between -1 and +1)

        ESNULL: Enrichment score calculated from random permutations.

        Hits\_Indices: Index of a gene in gene\_list, if gene is included in gene\_set.

        RES: Numerical vector containing the running enrichment score for all locations in the gene list .

    get\_libraries() → List[str][[source]](_modules/gseapy/base.html#GSEAbase.get_libraries)[](#gseapy.base.GSEAbase.get_libraries "Link to this definition")
    :   return active enrichr library name.Offical API

    load\_gmt(*gene\_list: Iterable[str]*, *gmt: List[str] | str | Dict[str, str]*) → Dict[str, List[str]][[source]](_modules/gseapy/base.html#GSEAbase.load_gmt)[](#gseapy.base.GSEAbase.load_gmt "Link to this definition")
    :   load gene set dict

    load\_gmt\_only(*gmt: List[str] | str | Dict[str, str]*) → Dict[str, List[str]][[source]](_modules/gseapy/base.html#GSEAbase.load_gmt_only)[](#gseapy.base.GSEAbase.load_gmt_only "Link to this definition")
    :   parse gene\_sets.
        gmt: List, Dict, Strings

        However,this function will merge different gene sets into one big dict to
        save computation time for later.

    make\_unique(*rank\_metric: DataFrame*, *col\_idx: int*) → DataFrame[[source]](_modules/gseapy/base.html#GSEAbase.make_unique)[](#gseapy.base.GSEAbase.make_unique "Link to this definition")
    :   make gene id column unique by adding a digit, similar to R’s make.unique

    parse\_gmt(*gmt: str*) → Dict[str, List[str]][[source]](_modules/gseapy/base.html#GSEAbase.parse_gmt)[](#gseapy.base.GSEAbase.parse_gmt "Link to this definition")
    :   gmt parser when input is a string

    plot(*terms: str | List[str]*, *colors: str | List[str] | None = None*, *legend\_kws: Dict[str, Any] | None = None*, *figsize: Tuple[float, float] = (4, 5)*, *show\_ranking: bool = True*, *ofname: str | None = None*)[[source]](_modules/gseapy/base.html#GSEAbase.plot)[](#gseapy.base.GSEAbase.plot "Link to this definition")
    :   terms: str, list. terms/pathways to show
        colors: str, list. list of colors for each term/pathway
        legend\_kws: kwargs to pass to ax.legend. e.g. loc, bbox\_to\_achor.
        ofname: savefig

    prepare\_outdir()[[source]](_modules/gseapy/base.html#GSEAbase.prepare_outdir)[](#gseapy.base.GSEAbase.prepare_outdir "Link to this definition")
    :   create temp directory.

    *property* results[](#gseapy.base.GSEAbase.results "Link to this definition")
    :   compatible to old style

    to\_df(*gsea\_summary: List[Dict]*, *gmt: Dict[str, List[str]]*, *rank\_metric: Series | DataFrame*, *indices: List | None = None*)[[source]](_modules/gseapy/base.html#GSEAbase.to_df)[](#gseapy.base.GSEAbase.to_df "Link to this definition")
    :   Convernt GSEASummary to DataFrame

        rank\_metric: if a Series, then it must be sorted in descending order already
        :   if a DataFrame, indices must not None.

        indices: Only works for DataFrame input. Stores the indices of sorted array

## 6.3. Over-representation Statistics[](#module-gseapy.stats "Link to this heading")

gseapy.stats.calc\_pvalues(*query*, *gene\_sets*, *background=20000*, *\*\*kwargs*)[[source]](_modules/gseapy/stats.html#calc_pvalues)[](#gseapy.stats.calc_pvalues "Link to this definition")
:   calculate pvalues for all categories in the graph

    Parameters:
    :   * **query** (*set*) – set of identifiers for which the p value is calculated
        * **gene\_sets** (*dict*) – gmt file dict after background was set
        * **background** (*set*) – total number of genes in your annotated database.

    Returns:
    :   pvalues
        x: overlapped gene number
        n: length of gene\_set which belongs to each terms
        hits: overlapped gene names.

    ### 6. For 2\*2 contingency table:[](#for-2-2-contingency-table "Link to this heading")

    > in query | not in query | row total

    => in gene\_set | a | b | a+b
    => not in gene\_set | c | d | c+d

    > column total | a+b+c+d = anno database

    Then, in R
    :   x=a the number of white balls drawn without replacement
        :   from an urn which contains both black and white balls.

        m=a+b the number of white balls in the urn
        n=c+d the number of black balls in the urn
        k=a+c the number of balls drawn from the urn

    In Scipy:
    for args in scipy.hypergeom.sf(k, M, n, N, loc=0):

    > M: the total number of objects,
    > n: the total number of Type I objects.
    > k: the random variate represents the number of Type I objects in N drawn
    >
    > > without replacement from the total population.

    Therefore, these two functions are the same when using parameters from 2\*2 table:
    R: > phyper(x-1, m, n, k, lower.tail=FALSE)
    Scipy: >>> hypergeom.sf(x-1, m+n, m, k)

    For Odds ratio in Enrichr (see <https://maayanlab.cloud/Enrichr/help#background&q=4>)

    > oddsRatio = (1.0 \* x \* d) / Math.max(1.0 \* b \* c, 1)

    where:

    > x are the overlapping genes,
    > b (m-x) are the genes in the annotated set - overlapping genes,
    > c (k-x) are the genes in the input set - overlapping genes,
    > d (bg-m-k+x) are the 20,000 genes (or total genes in the background) - genes in the annotated set - genes in the input set + overlapping genes

gseapy.stats.fdrcorrection(*pvals*, *alpha=0.05*)[[source]](_modules/gseapy/stats.html#fdrcorrection)[](#gseapy.stats.fdrcorrection "Link to this definition")
:   benjamini hocheberg fdr correction. inspired by statsmodels

gseapy.stats.multiple\_testing\_correction(*ps*, *alpha=0.05*, *method='benjamini-hochberg'*, *\*\*kwargs*)[[source]](_modules/gseapy/stats.html#multiple_testing_correction)[](#gseapy.stats.multiple_testing_correction "Link to this definition")
:   correct pvalues for multiple testing and add corrected q value

    Parameters:
    :   * **ps** – list of pvalues
        * **alpha** – significance level default : 0.05
        * **method** – multiple testing correction method [bonferroni|benjamini-hochberg]

    Returns (q, rej):
    :   two lists of q-values and rejected nodes

## 6.4. Enrichr API[](#module-gseapy.enrichr "Link to this heading")

*class* gseapy.enrichr.Enrichr(*gene\_list: Iterable[str]*, *gene\_sets: List[str] | str | Dict[str, str]*, *organism: str = 'human'*, *outdir: str | None = 'Enrichr'*, *background: List[str] | int | str = 'hsapiens\_gene\_ensembl'*, *cutoff: float = 0.05*, *format: str = 'pdf'*, *figsize: Tuple[float, float] = (6.5, 6)*, *top\_term: int = 10*, *no\_plot: bool = False*, *verbose: bool = False*)[[source]](_modules/gseapy/enrichr.html#Enrichr)[](#gseapy.enrichr.Enrichr "Link to this definition")
:   Enrichr API

    check\_genes(*gene\_list: List[str]*, *usr\_list\_id: str*)[[source]](_modules/gseapy/enrichr.html#Enrichr.check_genes)[](#gseapy.enrichr.Enrichr.check_genes "Link to this definition")
    :   Compare the genes sent and received to get successfully recognized genes

    check\_uppercase(*gene\_list: List[str]*)[[source]](_modules/gseapy/enrichr.html#Enrichr.check_uppercase)[](#gseapy.enrichr.Enrichr.check_uppercase "Link to this definition")
    :   Check whether a list of gene names are mostly in uppercase.

        ### 6. Parameters[](#id1 "Link to this heading")

        gene\_listlist
        :   A list of gene names

        ### 6. Returns[](#id2 "Link to this heading")

        bool
        :   Whether the list of gene names are mostly in uppercase

    enrich(*gmt: Dict[str, List[str]]*)[[source]](_modules/gseapy/enrichr.html#Enrichr.enrich)[](#gseapy.enrichr.Enrichr.enrich "Link to this definition")
    :   use local mode

        p = p-value computed using the Fisher exact test (Hypergeometric test)
        z = z-score (Odds Ratio)
        combine score = - log(p)·z

        see here: <http://amp.pharm.mssm.edu/Enrichr/help#background&q=4>

        columns contain:

        > Term Overlap P-value Odds Ratio Combinde Score Adjusted\_P-value Genes

    filter\_gmt(*gmt*, *background*)[[source]](_modules/gseapy/enrichr.html#Enrichr.filter_gmt)[](#gseapy.enrichr.Enrichr.filter_gmt "Link to this definition")
    :   the gmt values should be filtered only for genes that exist in background
        this substantially affect the significance of the test, the hypergeometric distribution.

        Parameters:
        :   * **gmt** – a dict of gene sets.
            * **background** – list, set, or tuple. A list of custom backgound genes.

    get\_background() → Set[str][[source]](_modules/gseapy/enrichr.html#Enrichr.get_background)[](#gseapy.enrichr.Enrichr.get_background "Link to this definition")
    :   get background gene

    get\_libraries() → List[str][[source]](_modules/gseapy/enrichr.html#Enrichr.get_libraries)[](#gseapy.enrichr.Enrichr.get_libraries "Link to this definition")
    :   return active enrichr library name. Official API

    get\_results(*gene\_list: List[str]*) → Tuple[AnyStr, DataFrame][[source]](_modules/gseapy/enrichr.html#Enrichr.get_results)[](#gseapy.enrichr.Enrichr.get_results "Link to this definition")
    :   Enrichr API

    parse\_background(*gmt: Dict[str, List[str]] | None = None*)[[source]](_modules/gseapy/enrichr.html#Enrichr.parse_background)[](#gseapy.enrichr.Enrichr.parse_background "Link to this definition")
    :   set background genes

    parse\_genelists() → str[[source]](_modules/gseapy/enrichr.html#Enrichr.parse_genelists)[](#gseapy.enrichr.Enrichr.parse_genelists "Link to this definition")
    :   parse gene list

    parse\_genesets(*gene\_sets=None*)[[source]](_modules/gseapy/enrichr.html#Enrichr.parse_genesets)[](#gseapy.enrichr.Enrichr.parse_genesets "Link to this definition")
    :   parse gene\_sets input file type

    prepare\_outdir()[[source]](_modules/gseapy/enrichr.html#Enrichr.prepare_outdir)[](#gseapy.enrichr.Enrichr.prepare_outdir "Link to this definition")
    :   create temp directory.

    run()[[source]](_modules/gseapy/enrichr.html#Enrichr.run)[](#gseapy.enrichr.Enrichr.run "Link to this definition")
    :   run enrichr for one sample gene list but multi-libraries

    send\_genes(*payload*, *url*) → Dict[[source]](_modules/gseapy/enrichr.html#Enrichr.send_genes)[](#gseapy.enrichr.Enrichr.send_genes "Link to this definition")
    :   send gene list to enrichr server

    set\_organism()[[source]](_modules/gseapy/enrichr.html#Enrichr.set_organism)[](#gseapy.enrichr.Enrichr.set_organism "Link to this definition")
    :   Select Enrichr organism from below:

        Human & Mouse, H. sapiens & M. musculus
        Fly, D. melanogaster
        Yeast, S. cerevisiae
        Worm, C. elegans
        Fish, D. rerio

## 6.5. BioMart API[](#module-gseapy.biomart "Link to this heading")

*class* gseapy.biomart.Biomart(*host: str = 'www.ensembl.org'*, *verbose: bool = False*)[[source]](_modules/gseapy/biomart.html#Biomart)[](#gseapy.biomart.Biomart "Link to this definition")
:   query from BioMart

    add\_filter(*name: str*, *value: Iterable[str]*)[[source]](_modules/gseapy/biomart.html#Biomart.add_filter)[](#gseapy.biomart.Biomart.add_filter "Link to this definition")
    :   key: filter names
        value: Iterable[str]

    get\_attributes(*dataset: str = 'hsapiens\_gene\_ensembl'*)[[source]](_modules/gseapy/biomart.html#Biomart.get_attributes)[](#gseapy.biomart.Biomart.get_attributes "Link to this definition")
    :   Get available attritbutes from dataset you’ve selected

    get\_datasets(*mart: str = 'ENSEMBL\_MART\_ENSEMBL'*)[[source]](_modules/gseapy/biomart.html#Biomart.get_datasets)[](#gseapy.biomart.Biomart.get_datasets "Link to this definition")
    :   Get available datasets from mart you’ve selected

    get\_filters(*dataset: str = 'hsapiens\_gene\_ensembl'*)[[source]](_modules/gseapy/biomart.html#Biomart.get_filters)[](#gseapy.biomart.Biomart.get_filters "Link to this definition")
    :   Get available filters from dataset you’ve selected

    get\_marts()[[source]](_modules/gseapy/biomart.html#Biomart.get_marts)[](#gseapy.biomart.Biomart.get_marts "Link to this definition")
    :   Get available marts and their names.

    query(*dataset: str = 'hsapiens\_gene\_ensembl'*, *attributes: List[str] | None = []*, *filters: Dict[str, Iterable[str]] | None = {}*, *filename: str | None = None*)[[source]](_modules/gseapy/biomart.html#Biomart.query)[](#gseapy.biomart.Biomart.query "Link to this definition")
    :   mapping ids using BioMart.

        Parameters:
        :   * **dataset** – str, default: ‘hsapiens\_gene\_ensembl’
            * **attributes** – str, list, tuple
            * **filters** – dict, {‘filter name’: list(filter value)}
            * **host** – www.ensembl.org, asia.ensembl.org, useast.ensembl.org

        Returns:
        :   a dataframe contains all attributes you selected.

        Example:

        ```
        >>> queries = {'ensembl_gene_id': ['ENSG00000125285','ENSG00000182968'] } # need to be a python dict
        >>> results = bm.query(dataset='hsapiens_gene_ensembl',
                               attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'go_id'],
                               filters=queries)
        ```

    query\_simple(*dataset: str = 'hsapiens\_gene\_ensembl'*, *attributes: List[str] = []*, *filters: Dict[str, Iterable[str]] = {}*, *filename: str | None = None*)[[source]](_modules/gseapy/biomart.html#Biomart.query_simple)[](#gseapy.biomart.Biomart.query_simple "Link to this definition")
    :   This function is a simple version of BioMart REST API.
        same parameter to query().

        However, you could get cross page of mapping. such as Mouse 2 human gene names

        **Note**: it will take a couple of minutes to get the results.
        A xml template for querying biomart. (see <https://gist.github.com/keithshep/7776579>)

        Example::
        :   ```
            >>> from gseapy import Biomart
            >>> bm = Biomart()
            >>> results = bm.query_simple(dataset='mmusculus_gene_ensembl',
                                          attributes=['ensembl_gene_id',
                                                      'external_gene_name',
                                                      'hsapiens_homolog_associated_gene_name',
                                                      'hsapiens_homolog_ensembl_gene'])
            ```

## 6.6. Parser[](#module-gseapy.parser "Link to this heading")

gseapy.parser.download\_library(*name: str*, *organism: str = 'human'*, *filename: str | None = None*) → Dict[str, List[str]][[source]](_modules/gseapy/parser.html#download_library)[](#gseapy.parser.download_library "Link to this definition")
:   download enrichr libraries.

    Parameters:
    :   * **name** (*str*) – the enrichr library name. see gseapy.get\_library\_name().
        * **organism** (*str*) – Select one from { ‘Human’, ‘Mouse’, ‘Yeast’, ‘Fly’, ‘Fish’, ‘Worm’ }
        * **filename** (*str*) – the file name to save if not None.

    Return dict:
    :   gene\_sets of the enrichr library from selected organism

gseapy.parser.get\_library(*name: str*, *organism: str = 'Human'*, *min\_size: int = 0*, *max\_size: int = 2000*, *save: str | None = None*, *gene\_list: List[str] | None = None*) → Dict[str, List[str]][[source]](_modules/gseapy/parser.html#get_library)[](#gseapy.parser.get_library "Link to this definition")
:   Parse gene\_sets.gmt(gene set database) file or download from enrichr server.

    Parameters:
    :   * **name** (*str*) – the gene\_sets.gmt file or an enrichr library name.
          checkout full enrichr library name here: <https://maayanlab.cloud/Enrichr/#libraries>
        * **organism** (*str*) – choose one from { ‘Human’, ‘Mouse’, ‘Yeast’, ‘Fly’, ‘Fish’, ‘Worm’ }.
          This arugment has not effect if input is a .gmt file.
        * **min\_size** – Minimum allowed number of genes for each gene set. Default: 0.
        * **max\_size** – Maximum allowed number of genes for each gene set. Default: 2000.
        * **save** (*str*) – the path to save the filtered gene set database.
        * **gene\_list** – if input a gene list, min and max overlapped genes between gene set and gene\_list are kept.

    Return dict:
    :   Return a filtered gene set database dictionary.

    Note: **DO NOT** filter gene sets, when use `replot()`. Because `GSEA` Desktop have already
    done this for you.

gseapy.parser.get\_library\_name(*organism: str = 'Human'*) → List[str][[source]](_modules/gseapy/parser.html#get_library_name)[](#gseapy.parser.get_library_name "Link to this definition")
:   return enrichr active enrichr library name.
    see also: <https://maayanlab.cloud/modEnrichr/>

    Parameters:
    :   **organism** (*str*) – Select one from { ‘Human’, ‘Mouse’, ‘Yeast’, ‘Fly’, ‘Fish’, ‘Worm’ }

    Returns:
    :   a list of enrichr libraries from selected database

gseapy.parser.gsea\_cls\_parser(*cls: str*) → Tuple[str][[source]](_modules/gseapy/parser.html#gsea_cls_parser)[](#gseapy.parser.gsea_cls_parser "Link to this definition")
:   Extract class(phenotype) name from .cls file.

    Parameters:
    :   **cls** – the a class list instance or .cls file which is identical to GSEA input .

    Returns:
    :   phenotype name and a list of class vector.

gseapy.parser.gsea\_edb\_parser(*results\_path: str*) → Dict[str, List[str]][[source]](_modules/gseapy/parser.html#gsea_edb_parser)[](#gseapy.parser.gsea_edb_parser "Link to this definition")
:   Parse results.edb file stored under **edb** file folder.

    Parameters:
    :   **results\_path** – the path of results.edb file.

    Returns:
    :   a dict contains { enrichment\_term: [es, nes, pval, fdr, fwer, hit\_ind]}

gseapy.parser.read\_gmt(*path: str*) → Dict[str, List[str]][[source]](_modules/gseapy/parser.html#read_gmt)[](#gseapy.parser.read_gmt "Link to this definition")
:   Read GMT file

    Parameters:
    :   **path** (*str*) – the path to a gmt file.

    Returns:
    :   a dict object

## 6.7. Visualization[](#module-gseapy.plot "Link to this heading")

*class* gseapy.plot.MidpointNormalize(*vmin=None*, *vmax=None*, *vcenter=None*, *clip=False*)[[source]](_modules/gseapy/plot.html#MidpointNormalize)[](#gseapy.plot.MidpointNormalize "Link to this definition")
:   inverse(*value*)[[source]](_modules/gseapy/plot.html#MidpointNormalize.inverse)[](#gseapy.plot.MidpointNormalize.inverse "Link to this definition")
    :   Maps the normalized value (i.e., index in the colormap) back to image
        data value.

        ### 6. Parameters[](#id3 "Link to this heading")

        value
        :   Normalized value.

gseapy.plot.barplot(*df: DataFrame*, *column: str = 'Adjusted P-value'*, *group: str | None = None*, *title: str = ''*, *cutoff: float = 0.05*, *top\_term: int = 10*, *ax: Axes | None = None*, *figsize: Tuple[float, float] = (4, 6)*, *color: str | List[str] | Dict[str, str] = 'salmon'*, *ofname: str | None = None*, *\*\*kwargs*)[[source]](_modules/gseapy/plot.html#barplot)[](#gseapy.plot.barplot "Link to this definition")
:   Visualize GSEApy Results.
    When multiple datasets exist in the input dataframe, the group argument is your friend.

    Parameters:
    :   * **df** – GSEApy DataFrame results.
        * **column** – column name in df to map the x-axis data. Default: Adjusted P-value
        * **group** – group by the variable in df that will produce bars with different colors.
        * **title** – figure title.
        * **cutoff** – terms with column value < cut-off are shown. Work only for
          (“Adjusted P-value”, “P-value”, “NOM p-val”, “FDR q-val”)
        * **top\_term** – number of top enriched terms grouped by hue are shown.
        * **ax** – Matplotlib axes. If None, create a new figure.
        * **figsize** – tuple, matplotlib figsize. only used when ax is None.
        * **color** – color or list or dict of matplotlib.colors. Must be reconigzed by matplotlib.
          if dict input, dict keys must be found in the group
        * **ofname** – output file name. If None, don’t save figure

    Returns:
    :   matplotlib.Axes. return None if given ofname.
        Only terms with column <= cut-off are plotted.

gseapy.plot.dotplot(*df: DataFrame*, *column: str = 'Adjusted P-value'*, *x: str | None = None*, *y: str = 'Term'*, *x\_order: List[str] | bool = False*, *y\_order: List[str] | bool = False*, *title: str = ''*, *cutoff: float = 0.05*, *top\_term: int = 10*, *size: float = 5*, *ax: Axes | None = None*, *figsize: Tuple[float, float] = (4, 6)*, *cmap: str = 'viridis\_r'*, *ofname: str | None = None*, *xticklabels\_rot: float | None = None*, *yticklabels\_rot: float | None = None*, *marker: str = 'o'*, *show\_ring: bool = False*, *\*\*kwargs*)[[source]](_modules/gseapy/plot.html#dotplot)[](#gseapy.plot.dotplot "Link to this definition")
:   Visualize GSEApy Results with categorical scatterplot
    When multiple datasets exist in the input dataframe, the x argument is your friend.

    Parameters:
    :   * **df** – GSEApy DataFrame results.
        * **column** – column name in df that map the dot colors. Default: Adjusted P-value.
        * **x** – Categorical variable in df that map the x-axis data. Default: None.
        * **y** – Categorical variable in df that map the y-axis data. Default: Term.
        * **x\_order** – bool, array-like list. Default: False.
          If True, peformed hierarchical\_clustering on X-axis.
          or input a array-like list of x categorical levels.
        * **x\_order** – bool, array-like list. Default: False.
          If True, peformed hierarchical\_clustering on Y-axis.
          or input a array-like list of y categorical levels.
        * **title** – Figure title.
        * **cutoff** – Terms with column value < cut-off are shown. Work only for
          (“Adjusted P-value”, “P-value”, “NOM p-val”, “FDR q-val”)
        * **top\_term** – Number of enriched terms to show (based on values in the column (colormap)).
        * **size** – float, scale the dot size to get proper visualization.
        * **ax** – Matplotlib axes.
        * **figsize** – tuple, matplotlib figure size, only used when ax is None.
        * **cmap** – Matplotlib colormap for mapping the column semantic.
        * **ofname** – Output file name. If None, don’t save figure
        * **marker** – The matplotlib.markers. See <https://matplotlib.org/stable/api/markers_api.html>
        * **bool** (*show\_ring*) – Whether to draw outer ring.

    Returns:
    :   matplotlib.Axes if ofname is None.
        Only terms with column <= cut-off are plotted.

gseapy.plot.enrichment\_map(*df: DataFrame*, *column: str = 'Adjusted P-value'*, *cutoff: float = 0.05*, *top\_term: int = 10*, *\*\*kwargs*) → Tuple[DataFrame, DataFrame][[source]](_modules/gseapy/plot.html#enrichment_map)[](#gseapy.plot.enrichment_map "Link to this definition")
:   Visualize GSEApy Results.
    Node size corresponds to the percentage of gene overlap in a certain term of interest.
    Colour of the node corresponds to the significance of the enriched terms.
    Edge size corresponds to the number of genes that overlap between the two connected nodes.
    Gray edges correspond to both nodes when it is the only colour edge.
    When there are two different edge colours, red corresponds to positve nodes and blue corresponds to negative nodes.

    Parameters:
    :   * **df** – GSEApy DataFrame results.
        * **column** – column name in df to map the node colors. Default: Adjusted P-value or FDR q-val.
          choose from (“Adjusted P-value”, “P-value”, “FDR q-val”, “NOM p-val”).
        * **group** – group by the variable in df that will produce bars with different colors.
        * **title** – figure title.
        * **cutoff** – nodes with column value < cut-off are shown. Work only for
          (“Adjusted P-value”, “P-value”, “NOM p-val”, “FDR q-val”)
        * **top\_term** – number of top enriched terms are selected as nodes.

    Returns:
    :   tuple of dataframe (nodes, edges)

gseapy.plot.gseaplot(*term: str*, *hits: Sequence[int]*, *nes: float*, *pval: float*, *fdr: float*, *RES: Sequence[float]*, *rank\_metric: Sequence[float] | None = None*, *pheno\_pos: str = ''*, *pheno\_neg: str = ''*, *color: str = '#88C544'*, *figsize: Tuple[float, float] = (6, 5.5)*, *cmap: str = 'seismic'*, *ofname: str | None = None*, *\*\*kwargs*) → List[Axes] | None[[source]](_modules/gseapy/plot.html#gseaplot)[](#gseapy.plot.gseaplot "Link to this definition")
:   This is the main function for generating the gsea plot.

    Parameters:
    :   * **term** – gene\_set name
        * **hits** – hits indices of rank\_metric.index presented in gene set S.
        * **nes** – Normalized enrichment scores.
        * **pval** – nominal p-value.
        * **fdr** – false discovery rate.
        * **RES** – running enrichment scores.
        * **rank\_metric** – pd.Series for rankings, rank\_metric.values.
        * **pheno\_pos** – phenotype label, positive correlated.
        * **pheno\_neg** – phenotype label, negative correlated.
        * **color** – color for RES and hits.
        * **figsize** – matplotlib figsize.
        * **ofname** – output file name. If None, don’t save figure

    return matplotlib.Figure.

gseapy.plot.gseaplot2(*terms: List[str]*, *hits: List[Sequence[int]]*, *RESs: List[Sequence[float]]*, *rank\_metric: Sequence[float] | None = None*, *colors: str | List[str] | None = None*, *figsize: Tuple[float, float] = (6, 4)*, *legend\_kws: Dict[str, Any] | None = None*, *ofname: str | None = None*, *\*\*kwargs*) → List[Axes] | None[[source]](_modules/gseapy/plot.html#gseaplot2)[](#gseapy.plot.gseaplot2 "Link to this definition")
:   Trace plot for combining multiple terms/pathways into one plot
    :param terms: list of terms to show in trace plot
    :param hits: list of hits indices correspond to each term.
    :param RESs: list of running enrichment scores correspond to each term.
    :param rank\_metric: Optional, rankings.
    :param figsize: matplotlib figsize.
    :legend\_kws: Optional, contol the location of lengends
    :param ofname: output file name. If None, don’t save figure

    return matplotlib.Figure.

gseapy.plot.heatmap(*df: DataFrame*, *z\_score: int | None = None*, *title: str = ''*, *figsize: Tuple[float, float] = (5, 5)*, *cmap: str | None = None*, *xticklabels: bool = True*, *yticklabels: bool = True*, *ofname: str | None = None*, *ax: Axes | None = None*, *\*\*kwargs*)[[source]](_modules/gseapy/plot.html#heatmap)[](#gseapy.plot.heatmap "Link to this definition")
:   Visualize the dataframe.

    Parameters:
    :   * **df** – DataFrame from expression table.
        * **z\_score** – 0, 1, or None. z\_score axis{0, 1}. If None, not scale.
        * **title** – figure title.
        * **figsize** – heatmap figsize.
        * **cmap** – matplotlib colormap. e.g. “RdBu\_r”.
        * **xticklabels** – bool, whether to show xticklabels.
        * **xticklabels** – bool, whether to show xticklabels.
        * **ofname** – output file name. If None, don’t save figure.
        * **ax** – matplotlib axes. Default: None.

    Returns:
    :   ax if ofname is None.

gseapy.plot.ringplot(*df: DataFrame*, *column: str = 'Adjusted P-value'*, *x: str | None = None*, *title: str = ''*, *cutoff: float = 0.05*, *top\_term: int = 10*, *size: float = 5*, *figsize: Tuple[float, float] = (4, 6)*, *cmap: str = 'viridis\_r'*, *ofname: str | None = None*, *xticklabels\_rot: float | None = None*, *yticklabels\_rot: float | None = None*, *marker='o'*, *show\_ring: bool = True*, *\*\*kwargs*)[[source]](_modules/gseapy/plot.html#ringplot)[](#gseapy.plot.ringplot "Link to this definition")
:   ringplot is deprecated, use dotplot instead

    Parameters:
    :   * **df** – GSEApy DataFrame results.
        * **x** – Group by the variable in df that will produce categorical scatterplot.
        * **column** – column name in df to map the dot colors. Default: Adjusted P-value
        * **title** – figure title
        * **cutoff** – terms with column value < cut-off are shown. Work only for
          (“Adjusted P-value”, “P-value”, “NOM p-val”, “FDR q-val”)
        * **top\_term** – number of enriched terms to show.
        * **size** – float, scale the dot size to get proper visualization.
        * **figsize** – tuple, matplotlib figure size.
        * **cmap** – matplotlib colormap for mapping the column semantic.
        * **ofname** – output file name. If None, don’t save figure
        * **marker** – the matplotlib.markers. See <https://matplotlib.org/stable/api/markers_api.html>
        * **bool** (*show\_ring*) – whether to show outer ring.

    Returns:
    :   matplotlib.Axes. return None if given ofname.
        Only terms with column <= cut-off are plotted.

gseapy.plot.zscore(*data2d: DataFrame*, *axis: int | None = 0*)[[source]](_modules/gseapy/plot.html#zscore)[](#gseapy.plot.zscore "Link to this definition")
:   Standardize the mean and variance of the data axis Parameters.

    Parameters:
    :   * **data2d** – DataFrame to normalize.
        * **axis** – int, Which axis to normalize across. If 0, normalize across rows,
          if 1, normalize across columns. If None, don’t change data

    Returns:
    :   Normalized DataFrame. Normalized data with a mean of 0 and variance of 1
        across the specified axis.

## 6.8. Scientific Journal and Sci- themed Color Palettes[](#module-gseapy.scipalette "Link to this heading")

## 6.9. Utils[](#utils "Link to this heading")