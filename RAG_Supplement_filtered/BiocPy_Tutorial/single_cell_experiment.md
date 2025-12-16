Source URL: https://biocpy.github.io/tutorial/chapters/experiments/single_cell_experiment.html
Date Scraped: 2025-12-15

---

This package provides container class to represent single-cell experimental data as 2-dimensional matrices. In these matrices, the rows typically denote features or genomic regions of interest, while columns represent cells. In addition, a `SingleCellExperiment` (SCE) object may contain low-dimensionality embeddings, alternative experiments performed on same sample or set of cells.

Important

The design of `SingleCellExperiment` class and its derivates adheres to the R/Bioconductor specification, where rows correspond to features, and columns represent cells.

Note

These classes follow a functional paradigm for accessing or setting properties, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.

To get started, install the package from [PyPI](https://pypi.org/project/singlecellexperiment/)

```
pip install singlecellexperiment
```

The `SingleCellExperiment` extends `RangeSummarizedExperiment` and contains additional attributes:

* `reduced_dims`: Slot for low-dimensionality embeddings for each cell.
* `alternative_experiments`: Manages multi-modal experiments performed on the same sample or set of cells.
* `row_pairs` or `column_pairs`: Stores relationships between features or cells.

Note

In contrast to R, matrices in Python are unnamed and do not contain row or column names. Hence, these matrices cannot be directly used as values in assays or alternative experiments. We strictly enforce type checks in these cases. To relax these restrictions for alternative experiments, set `type_check_alternative_experiments` to `False`.

Important

If you are using the `alternative_experiments` slot, the number of cells must match the parent experiment. Otherwise, the expectation is that the cells do not share the same sample or annotations and cannot be set in alternative experiments!

Before we construct a `SingleCellExperiment` object, lets generate information about rows, columns and a mock experimental data from single-cell rna-seq experiments:

Show the code

```
import pandas as pd
import numpy as np
from scipy import sparse as sp
from biocframe import BiocFrame
from genomicranges import GenomicRanges
from random import random 

nrows = 200
ncols = 6
counts = sp.rand(nrows, ncols, density=0.2, format="csr")
row_data = BiocFrame(
    {
        "seqnames": [
            "chr1",
            "chr2",
            "chr2",
            "chr2",
            "chr1",
            "chr1",
            "chr3",
            "chr3",
            "chr3",
            "chr3",
        ]
        * 20,
        "starts": range(100, 300),
        "ends": range(110, 310),
        "strand": ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"] * 20,
        "score": range(0, 200),
        "GC": [random() for _ in range(10)] * 20,
    }
)

col_data = pd.DataFrame(
    {
        "celltype": ["cluster1", "cluster2"] * 3,
    }
)

counts
```

```
<200x6 sparse matrix of type '<class 'numpy.float64'>'
    with 240 stored elements in Compressed Sparse Row format>
```

Now lets create the `SingleCellExperiment` instance:

```
from singlecellexperiment import SingleCellExperiment

sce = SingleCellExperiment(
    assays={"counts": counts}, row_data=row_data, column_data=col_data,
    reduced_dims = {"random_embeds": np.random.rand(ncols, 4)}
)

print(sce)
```

```
class: SingleCellExperiment
dimensions: (200, 6)
assays(1): ['counts']
row_data columns(6): ['seqnames', 'starts', 'ends', 'strand', 'score', 'GC']
row_names(0):  
column_data columns(1): ['celltype']
column_names(6): ['0', '1', '2', '3', '4', '5']
main_experiment_name:  
reduced_dims(1): ['random_embeds']
alternative_experiments(0): []
row_pairs(0): []
column_pairs(0): []
metadata(0):
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

Tip

You can also use delayed or file-backed arrays for representing experimental data, check out [this section](https://biocpy.github.io/tutorial/chapters/experiments/summarized_experiment.html#delayed-or-file-backed-arrays) from summarized experiment.

We provide convenient methods for loading an `AnnData` or `h5ad` file into `SingleCellExperiment` objects.

For example, lets create an `AnnData` object,

```
import anndata as ad
from scipy import sparse as sp

counts = sp.csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
adata = ad.AnnData(counts)
adata
```

```
AnnData object with n_obs × n_vars = 100 × 2000
```

Converting `AnnData` as `SingleCellExperiment` is straightforward:

```
sce_adata = SingleCellExperiment.from_anndata(adata)
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

and vice-verse. All assays from SCE are represented in the `layers` slot of the `AnnData` object:

```
adata2 = sce_adata.to_anndata()
print(adata2)
```

```
(AnnData object with n_obs × n_vars = 100 × 2000
    layers: 'X', None)
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/anndata/_core/aligned_df.py:67: ImplicitModificationWarning: Transforming to str index.
  warnings.warn("Transforming to str index.", ImplicitModificationWarning)
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/anndata/_core/aligned_df.py:67: ImplicitModificationWarning: Transforming to str index.
  warnings.warn("Transforming to str index.", ImplicitModificationWarning)
```

Similarly, one can load a h5ad file:

```
from singlecellexperiment import read_h5ad
sce_h5 = read_h5ad("../../assets/data/adata.h5ad")
print(sce_h5)
```

```
class: SingleCellExperiment
dimensions: (20, 30)
assays(3): ['array', 'sparse', 'X']
row_data columns(5): ['var_cat', 'cat_ordered', 'int64', 'float64', 'uint8']
row_names(20): ['gene0', 'gene1', 'gene2', ..., 'gene17', 'gene18', 'gene19']
column_data columns(5): ['obs_cat', 'cat_ordered', 'int64', 'float64', 'uint8']
column_names(30): ['cell0', 'cell1', 'cell2', ..., 'cell27', 'cell28', 'cell29']
main_experiment_name:  
reduced_dims(0): []
alternative_experiments(0): []
row_pairs(0): []
column_pairs(0): []
metadata(2): O_recarray nested
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

In addition, we also provide convenient methods to load a [10X Genomics HDF5 Feature-Barcode Matrix Format](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices) file.

```
from singlecellexperiment import read_tenx_h5
sce_h5 = read_tenx_h5("../../assets/data/tenx.sub.h5")
print(sce_h5)
```

```
class: SingleCellExperiment
dimensions: (1000, 3005)
assays(1): ['counts']
row_data columns(2): ['id', 'name']
row_names(0):  
column_data columns(0): []
column_names(0):  
main_experiment_name:  
reduced_dims(0): []
alternative_experiments(0): []
row_pairs(0): []
column_pairs(0): []
metadata(0):
```

Note

Methods are also available to read a 10x matrix market directory using the `read_tenx_mtx` function.

Getters are available to access various attributes using either the property notation or functional style.

Since `SingleCellExperiment` extends `RangedSummarizedExperiment`, all getters and setters from the base class are accessible here; more details [here](https://biocpy.github.io/tutorial/chapters/experiments/summarized_experiment.html).

```
# access assay names
print("reduced dim names (as property): ", sce.reduced_dim_names)
print("reduced dim names (functional style): ", sce.get_reduced_dim_names())

# access row data
print(sce.row_data)
```

```
reduced dim names (as property):  ['random_embeds']
reduced dim names (functional style):  ['random_embeds']
BiocFrame with 200 rows and 6 columns
      seqnames  starts    ends strand   score                  GC
        <list> <range> <range> <list> <range>              <list>
  [0]     chr1     100     110      -       0  0.9428265564541277
  [1]     chr2     101     111      +       1 0.22619486204847783
  [2]     chr2     102     112      +       2  0.3756405606579122
           ...     ...     ...    ...     ...                 ...
[197]     chr3     297     307      +     197 0.16951649789798318
[198]     chr3     298     308      -     198 0.39887354687383936
[199]     chr3     299     309      -     199 0.11610477514345241
```

One can access an reduced dimension by index or name:

```
sce.reduced_dim(0) # same as se.reduced_dim("random_embeds")
```

```
array([[0.44447259, 0.44771424, 0.65664489, 0.72822602],
       [0.3439318 , 0.83002796, 0.047597  , 0.58570195],
       [0.57008475, 0.39113037, 0.99168239, 0.38135755],
       [0.17634127, 0.56165685, 0.75319488, 0.32916635],
       [0.6977452 , 0.92124643, 0.93556774, 0.71834461],
       [0.79196687, 0.16772082, 0.45144653, 0.41957007]])
```

You can subset experimental data by using the subset (`[]`) operator. This operation accepts different slice input types, such as a boolean vector, a `slice` object, a list of indices, or names (if available) to subset.

In our previous example, we didn’t include row or column names. Let’s create another `SingleCellExperiment` object that includes names.

```
subset_sce = sce[0:10, 0:3]
print(subset_sce)
```

```
class: SingleCellExperiment
dimensions: (10, 3)
assays(1): ['counts']
row_data columns(6): ['seqnames', 'starts', 'ends', 'strand', 'score', 'GC']
row_names(0):  
column_data columns(1): ['celltype']
column_names(3): ['0', '1', '2']
main_experiment_name:  
reduced_dims(1): ['random_embeds']
alternative_experiments(0): []
row_pairs(0): []
column_pairs(0): []
metadata(0):
```

`SingleCellExperiment` implements methods for the `combine` generic from [**BiocUtils**](https://github.com/BiocPy/biocutils).

These methods enable the merging or combining of multiple `SingleCellExperiment` objects, allowing users to aggregate data from different experiments or conditions. Note: `row_pairs` and `column_pairs` are not ignored as part of this operation.

To demonstrate, let’s create multiple `SingleCellExperiment` objects (read more about this in [combine section from `SummarizedExperiment`](https://biocpy.github.io/tutorial/chapters/experiments/summarized_experiment.html#combining-experiments)).

Show the code

```
ncols = 10
nrows = 100
sce1 = SingleCellExperiment(
    assays={"counts": np.random.poisson(lam=10, size=(nrows, ncols))},
    row_data=BiocFrame({"A": [1] * nrows}),
    column_data=BiocFrame({"A": [1] * ncols}),
)

sce2 = SingleCellExperiment(
    assays={
        "counts": np.random.poisson(lam=10, size=(nrows, ncols)),
        # "normalized": np.random.normal(size=(nrows, ncols)),
    },
    row_data=BiocFrame({"A": [3] * nrows}),
    column_data=BiocFrame({"A": [3] * ncols}),
)

rowdata1 = pd.DataFrame(
    {
        "seqnames": ["chr_5", "chr_3", "chr_2"],
        "start": [500, 300, 200],
        "end": [510, 310, 210],
    },
    index=["HER2", "BRCA1", "TPFK"],
)
coldata1 = pd.DataFrame(
    {
        "sample": ["SAM_1", "SAM_2", "SAM_3"],
        "disease": ["True", "True", "True"],
        "doublet_score": [0.15, 0.62, 0.18],
    },
    index=["cell_1", "cell_2", "cell_3"],
)
sce_alts1 = SingleCellExperiment(
    assays={
        "counts": np.random.poisson(lam=5, size=(3, 3)),
        "lognorm": np.random.lognormal(size=(3, 3)),
    },
    row_data=rowdata1,
    column_data=coldata1,
    row_names=["HER2", "BRCA1", "TPFK"],
    column_names=["cell_1", "cell_2", "cell_3"],
    metadata={"seq_type": "paired"},
    reduced_dims={"PCA": np.random.poisson(lam=10, size=(3, 5))},
    alternative_experiments={
        "modality1": SingleCellExperiment(
            assays={"counts2": np.random.poisson(lam=10, size=(3, 3))},
        )
    },
)

rowdata2 = pd.DataFrame(
    {
        "seqnames": ["chr_5", "chr_3", "chr_2"],
        "start": [500, 300, 200],
        "end": [510, 310, 210],
    },
    index=["HER2", "BRCA1", "TPFK"],
)
coldata2 = pd.DataFrame(
    {
        "sample": ["SAM_4", "SAM_5", "SAM_6"],
        "disease": ["True", "False", "True"],
        "doublet_score": [0.05, 0.23, 0.54],
    },
    index=["cell_4", "cell_5", "cell_6"],
)
sce_alts2 = SingleCellExperiment(
    assays={
        "counts": np.random.poisson(lam=5, size=(3, 3)),
        # "lognorm": np.random.lognormal(size=(3, 3)),
    },
    row_data=rowdata2,
    column_data=coldata2,
    metadata={"seq_platform": "Illumina NovaSeq 6000"},
    reduced_dims={"PCA": np.random.poisson(lam=5, size=(3, 5))},
    alternative_experiments={
        "modality1": SingleCellExperiment(
            assays={"counts2": np.random.poisson(lam=5, size=(3, 3))},
        )
    },
)
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

The `combine_rows` or `combine_columns` operations, expect all experiments to contain the same assay names. To combine experiments by row:

```
from biocutils import relaxed_combine_columns, combine_columns, combine_rows, relaxed_combine_rows
sce_combined = combine_rows(sce2, sce1)
print(sce_combined)
```

```
class: SingleCellExperiment
dimensions: (200, 10)
assays(1): ['counts']
row_data columns(1): ['A']
row_names(0):  
column_data columns(1): ['A']
column_names(0):  
main_experiment_name:  
reduced_dims(0): []
alternative_experiments(0): []
row_pairs(0): []
column_pairs(0): []
metadata(0):
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/singlecellexperiment/SingleCellExperiment.py:1089: UserWarning: 'row_pairs' and 'column_pairs' are currently ignored during this operation.
  warn(
```

Similarly to combine by column:

```
sce_combined = combine_columns(sce2, sce1)
print(sce_combined)
```

```
class: SingleCellExperiment
dimensions: (100, 20)
assays(1): ['counts']
row_data columns(1): ['A']
row_names(0):  
column_data columns(1): ['A']
column_names(0):  
main_experiment_name:  
reduced_dims(0): []
alternative_experiments(0): []
row_pairs(0): []
column_pairs(0): []
metadata(0):
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/singlecellexperiment/SingleCellExperiment.py:1131: UserWarning: 'row_pairs' and 'column_pairs' are currently ignored during this operation.
  warn(
```

Note

You can use `relaxed_combine_columns` or `relaxed_combined_rows` when there’s mismatch in the number of features or samples. Missing rows or columns in any object are filled in with appropriate placeholder values before combining, e.g. missing assay’s are replaced with a masked numpy array.

```
# sce_alts1 contains an additional assay not present in sce_alts2
sce_relaxed_combine = relaxed_combine_columns(sce_alts1, sce_alts2)
print(sce_relaxed_combine)
```

```
class: SingleCellExperiment
dimensions: (3, 6)
assays(2): ['counts', 'lognorm']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(3): ['HER2', 'BRCA1', 'TPFK']
column_data columns(3): ['sample', 'disease', 'doublet_score']
column_names(6): ['cell_1', 'cell_2', 'cell_3', 'cell_4', 'cell_5', 'cell_6']
main_experiment_name:  
reduced_dims(1): ['PCA']
alternative_experiments(1): ['modality1']
row_pairs(0): []
column_pairs(0): []
metadata(1): seq_type
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/singlecellexperiment/SingleCellExperiment.py:1241: UserWarning: 'row_pairs' and 'column_pairs' are currently ignored during this operation.
  warn("'row_pairs' and 'column_pairs' are currently ignored during this operation.")
```

The package also provides methods to convert a `SingleCellExperiment` object into a `MuData` representation:

```
mdata = sce.to_mudata()
mdata
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/anndata/_core/aligned_df.py:67: ImplicitModificationWarning: Transforming to str index.
  warnings.warn("Transforming to str index.", ImplicitModificationWarning)
```

```
MuData object with n_obs × n_vars = 6 × 200
  var:  'seqnames', 'starts', 'ends', 'strand', 'score', 'GC'
  1 modality
    Unknown Modality:   6 x 200
      obs:  'celltype'
      var:  'seqnames', 'starts', 'ends', 'strand', 'score', 'GC'
      obsm: 'random_embeds'
      layers:   'counts'
```

or coerce to an `AnnData` object:

```
adata, alts = sce.to_anndata()
print("main experiment: ", adata)
print("alternative experiments: ", alts)
```

```
main experiment:  AnnData object with n_obs × n_vars = 6 × 200
    obs: 'celltype'
    var: 'seqnames', 'starts', 'ends', 'strand', 'score', 'GC'
    obsm: 'random_embeds'
    layers: 'counts'
alternative experiments:  None
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/anndata/_core/aligned_df.py:67: ImplicitModificationWarning: Transforming to str index.
  warnings.warn("Transforming to str index.", ImplicitModificationWarning)
```

---

* Check out reference [documentation](https://biocpy.github.io/SingleCellExperiment/index.html) for more details.
* R/Bioconductor’s [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) package.