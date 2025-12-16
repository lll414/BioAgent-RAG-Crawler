Source URL: https://biocpy.github.io/tutorial/chapters/experiments/summarized_experiment.html
Date Scraped: 2025-12-15

---

This package provides containers to represent genomic experimental data as 2-dimensional matrices. In these matrices, the rows typically denote features or genomic regions of interest, while columns represent samples or cells.

The package currently includes representations for both `SummarizedExperiment` and `RangedSummarizedExperiment`. A distinction lies in the fact that the rows of a `RangedSummarizedExperiment` object are expected to be `GenomicRanges` (tutorial [here](https://biocpy.github.io/tutorial/chapters/representations/genomic_ranges.html)), representing genomic regions of interest.

Important

The design of `SummarizedExperiment` class and its derivates adheres to the R/Bioconductor specification, where rows correspond to features, and columns represent samples or cells.

Note

These classes follow a functional paradigm for accessing or setting properties, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.

To get started, install the package from [PyPI](https://pypi.org/project/summarizedexperiment/)

```
pip install summarizedexperiment
```

A `SummarizedExperiment` contains three key attributes,

* `assays`: A dictionary of matrices with assay names as keys, e.g. counts, logcounts etc.
* `row_data`: Feature information e.g. genes, transcripts, exons, etc.
* `column_data`: Sample information about the columns of the matrices.

Important

Both `row_data` and `column_data` are expected to be [BiocFrame](https://biocpy.github.io/tutorial/chapters/representations/biocframe.html) objects and will be coerced to a `BiocFrame` for consistent downstream operations.

In addition, these classes can optionally accept `row_names` and `column_names`. Since `row_data` and `column_data` may also contain names, the following rules are used in the implementation:

* On **construction**, if `row_names` or `column_names` are not provided, these are automatically inferred from `row_data` and `column_data` objects.
* On **accessors** of these objects, the `row_names` in `row_data` and `column_data` are replaced by the equivalents from the SE level.
* On **setters** for these attributes, especially with the functional style (`set_row_data` and `set_column_data` methods), additional options are available to replace the names in the SE object.

Caution

These rules help avoid unexpected mdifications in names, when either `row_data` or `column_data` objects are modified.

To construct a `SummarizedExperiment`, we’ll first generate a matrix of read counts, representing the read counts from a series of RNA-seq experiments. Following that, we’ll create a `BiocFrame` object to denote feature information and a table for column annotations. This table may include the names for the columns and any other values we wish to represent.

Show the code

```
from random import random
import pandas as pd
import numpy as np
from biocframe import BiocFrame

nrows = 200
ncols = 6
counts = np.random.rand(nrows, ncols)
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
        "treatment": ["ChIP", "Input"] * 3,
    }
)
```

Note

The inputs `row_data` and `column_data` are expected to be `BiocFrame` objects and will be coerced to a `BiocFrame` if a pandas `DataFrame` is supplied.

Now, we can construct a `SummarizedExperiment` from this information.

```
from summarizedexperiment import SummarizedExperiment

se = SummarizedExperiment(
    assays={"counts": counts}, row_data=row_data, column_data=col_data
)

print(se)
```

```
class: SummarizedExperiment
dimensions: (200, 6)
assays(1): ['counts']
row_data columns(6): ['seqnames', 'starts', 'ends', 'strand', 'score', 'GC']
row_names(0):  
column_data columns(1): ['treatment']
column_names(6): ['0', '1', '2', '3', '4', '5']
metadata(0):
```

Similarly, we can use the same information to construct a `RangeSummarizedExperiment`. We convert feature information into a `GenomicRanges` object and provide this as `row_ranges`:

```
from genomicranges import GenomicRanges
from summarizedexperiment import RangedSummarizedExperiment

gr = GenomicRanges.from_pandas(row_data.to_pandas())

rse = RangedSummarizedExperiment(
    assays={"counts": counts}, row_data=row_data, row_ranges=gr, column_data=col_data
)
print(rse)
```

```
class: RangedSummarizedExperiment
dimensions: (200, 6)
assays(1): ['counts']
row_data columns(6): ['seqnames', 'starts', 'ends', 'strand', 'score', 'GC']
row_names(0):  
column_data columns(1): ['treatment']
column_names(6): ['0', '1', '2', '3', '4', '5']
metadata(0):
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

The general idea is that `DelayedArray`’s are a drop-in replacement for NumPy arrays, at least for [BiocPy](https://github.com/BiocPy) applications. Learn more about [delayed arrays here](https://biocpy.github.io/tutorial/chapters/representations/delayed_arrays.html).

For example, we can use the `DelayedArray` inside a `SummarizedExperiment`:

```
import numpy
import delayedarray

# create a delayed array, can also be file-backed
x = numpy.random.rand(100, 20)
d = delayedarray.wrap(x)

# operate over delayed arrays
filtered = d[1:100:2,1:8]
total = filtered.sum(axis=0)
normalized = filtered / total
transformed = numpy.log1p(normalized)

import summarizedexperiment as SE
se_delayed = SE.SummarizedExperiment({ "counts": filtered, "lognorm": transformed })
print(se_delayed)
```

```
class: SummarizedExperiment
dimensions: (50, 7)
assays(2): ['counts', 'lognorm']
row_data columns(0): []
row_names(0):  
column_data columns(0): []
column_names(0):  
metadata(0):
```

Converting a `SummarizedExperiment` to an `AnnData` representation is straightforward:

```
adata = se.to_anndata()
print(adata)
```

```
AnnData object with n_obs × n_vars = 6 × 200
    obs: 'treatment'
    var: 'seqnames', 'starts', 'ends', 'strand', 'score', 'GC'
    layers: 'counts'
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/anndata/_core/aligned_df.py:67: ImplicitModificationWarning: Transforming to str index.
  warnings.warn("Transforming to str index.", ImplicitModificationWarning)
```

Tip

To convert an `AnnData` object to a BiocPy representation, utilize the `from_anndata` method in the [SingleCellExperiment](https://biocpy.github.io/tutorial/chapters/experiments/single_cell_experiment.html) class. This minimizes the loss of information when converting between these two representations.

Getters are available to access various attributes using either the property notation or functional style.

```
# access assay names
print("assay names (as property): ", se.assay_names)
print("assay names (functional style): ", se.get_assay_names())

# access row data
print(se.row_data)
```

```
assay names (as property):  ['counts']
assay names (functional style):  ['counts']
BiocFrame with 200 rows and 6 columns
      seqnames  starts    ends strand   score                  GC
        <list> <range> <range> <list> <range>              <list>
  [0]     chr1     100     110      -       0  0.5937744055105046
  [1]     chr2     101     111      +       1  0.9968901673418803
  [2]     chr2     102     112      +       2 0.17011805868533125
           ...     ...     ...    ...     ...                 ...
[197]     chr3     297     307      +     197  0.1530596336715876
[198]     chr3     298     308      -     198 0.45833737361602656
[199]     chr3     299     309      -     199  0.9847343304818247
```

One can access an assay by index or name:

```
se.assay(0) # same as se.assay("counts")
```

```
array([[0.99010486, 0.3455029 , 0.02584926, 0.43959487, 0.73111188,
        0.07249016],
       [0.76088651, 0.7316889 , 0.05156813, 0.28567382, 0.61107981,
        0.54028849],
       [0.25542983, 0.30257744, 0.92346279, 0.93999678, 0.13602946,
        0.4609789 ],
       ...,
       [0.73177911, 0.03745444, 0.53115654, 0.02262766, 0.04611963,
        0.08600043],
       [0.82160941, 0.91571671, 0.32419894, 0.74134188, 0.08032582,
        0.98589182],
       [0.66215336, 0.22253834, 0.77108608, 0.32142434, 0.41430143,
        0.84835408]])
```

Important

All property-based setters are `in_place` operations, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.

```
modified_column_data = se.column_data.set_column("score", range(10,16))
modified_se = se.set_column_data(modified_column_data)
print(modified_se)
```

```
class: SummarizedExperiment
dimensions: (200, 6)
assays(1): ['counts']
row_data columns(6): ['seqnames', 'starts', 'ends', 'strand', 'score', 'GC']
row_names(0):  
column_data columns(2): ['treatment', 'score']
column_names(6): ['0', '1', '2', '3', '4', '5']
metadata(0):
```

Now, lets check the `column_data` on the original object.

```
print(se.column_data)
```

```
BiocFrame with 6 rows and 1 column
  treatment
     <list>
0      ChIP
1     Input
2      ChIP
3     Input
4      ChIP
5     Input
```

You can subset experimental data by using the subset (`[]`) operator. This operation accepts different slice input types, such as a boolean vector, a `slice` object, a list of indices, or names (if available) to subset.

In our previous example, we didn’t include row or column names. Let’s create another `SummarizedExperiment` object that includes names.

Show the code

```
row_data = BiocFrame({
    "seqnames": ["chr_5", "chr_3", "chr_2"],
    "start": [100, 200, 300],
    "end": [110, 210, 310],
})

col_data = BiocFrame({
    "sample": ["SAM_1", "SAM_3", "SAM_3"],
    "disease": ["True", "True", "True"],
    },
    row_names=["cell_1", "cell_2", "cell_3"],
)

se_with_names = SummarizedExperiment(
    assays={
        "counts": np.random.poisson(lam=5, size=(3, 3)),
        "lognorm": np.random.lognormal(size=(3, 3)),
    },
    row_data=row_data,
    column_data=col_data,
    row_names=["HER2", "BRCA1", "TPFK"],
    column_names=["cell_1", "cell_2", "cell_3"],
)

print(se_with_names)
```

```
class: SummarizedExperiment
dimensions: (3, 3)
assays(2): ['counts', 'lognorm']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(3): ['HER2', 'BRCA1', 'TPFK']
column_data columns(2): ['sample', 'disease']
column_names(3): ['cell_1', 'cell_2', 'cell_3']
metadata(0):
```

A straightforward slice operation:

```
subset_se = se_with_names[0:10, 0:3]
print(subset_se)
```

```
class: SummarizedExperiment
dimensions: (3, 3)
assays(2): ['counts', 'lognorm']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(3): ['HER2', 'BRCA1', 'TPFK']
column_data columns(2): ['sample', 'disease']
column_names(3): ['cell_1', 'cell_2', 'cell_3']
metadata(0):
```

Either one or both of the slices can contain names. These names are mapped to `row_names` and `column_names` of the `SummarizedExperiment` object.

```
subset_se = se_with_names[:2, ["cell_1", "cell_3"]]
print(subset_se)
```

```
class: SummarizedExperiment
dimensions: (2, 2)
assays(2): ['counts', 'lognorm']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(2): ['HER2', 'BRCA1']
column_data columns(2): ['sample', 'disease']
column_names(2): ['cell_1', 'cell_3']
metadata(0):
```

An `Exception` is raised if a names does not exist.

Similarly, you can also slice by a boolean array.

Important

Note that the boolean vectors should contain the same number of features for the row slice and the same number of samples for the column slice.

```
subset_se_with_bools = se_with_names[[True, True, False], [True, False, True]]
print(subset_se_with_bools)
```

```
class: SummarizedExperiment
dimensions: (2, 2)
assays(2): ['counts', 'lognorm']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(2): ['HER2', 'BRCA1']
column_data columns(2): ['sample', 'disease']
column_names(2): ['cell_1', 'cell_3']
metadata(0):
```

This is a feature not a bug :), you can specify an empty list to completely remove all rows or samples.

Warning

An empty array (`[]`) is not the same as an empty slice (`:`). This helps us avoid unintented operations.

```
subset = se_with_names[:2, []]
print(subset)
```

```
class: SummarizedExperiment
dimensions: (2, 0)
assays(2): ['counts', 'lognorm']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(2): ['HER2', 'BRCA1']
column_data columns(2): ['sample', 'disease']
column_names(0): []
metadata(0):
```

Additionally, since `RangeSummarizedExperiment` contains `row_ranges`, this allows us to perform a number of range-based operations that are possible on a `GenomicRanges` object.

For example, to subset `RangeSummarizedExperiment` with a **query** set of regions:

```
from iranges import IRanges
query = GenomicRanges(seqnames=["chr2"], ranges=IRanges([4], [6]), strand=["+"])

result = rse.subset_by_overlaps(query)
print(result)
```

```
class: RangedSummarizedExperiment
dimensions: (0, 6)
assays(1): ['counts']
row_data columns(6): ['seqnames', 'starts', 'ends', 'strand', 'score', 'GC']
row_names(0):  
column_data columns(1): ['treatment']
column_names(6): ['0', '1', '2', '3', '4', '5']
metadata(0):
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

Additionally, RSE supports many other interval based operations. Checkout the [documentation](https://biocpy.github.io/SummarizedExperiment/api/modules.html) for more details.

`SummarizedExperiment` implements methods for the `combine` generics from [**BiocUtils**](https://github.com/BiocPy/biocutils).

These methods enable the merging or combining of multiple `SummarizedExperiment` objects, allowing users to aggregate data from different experiments or conditions. To demonstrate, let’s create multiple `SummarizedExperiment` objects.

Show the code

```
rowData1 = pd.DataFrame(
    {
        "seqnames": ["chr_5", "chr_3", "chr_2"],
        "start": [10293804, 12098948, 20984392],
        "end": [28937947, 3872839, 329837492]
    },
    index=["HER2", "BRCA1", "TPFK"],
)
colData1 = pd.DataFrame(
    {
        "sample": ["SAM_1", "SAM_3", "SAM_3"],
        "disease": ["True", "True", "True"],
    },
    index=["cell_1", "cell_2", "cell_3"],
)
se1 = SummarizedExperiment(
    assays={
        "counts": np.random.poisson(lam=5, size=(3, 3)),
        "lognorm": np.random.lognormal(size=(3, 3))
    },
    row_data=rowData1,
    column_data=colData1,
    metadata={"seq_type": "paired"},
)

rowData2 = pd.DataFrame(
    {
        "seqnames": ["chr_5", "chr_3", "chr_2"],
        "start": [10293804, 12098948, 20984392],
        "end": [28937947, 3872839, 329837492]
    },
    index=["HER2", "BRCA1", "TPFK"],
)
colData2 = pd.DataFrame(
    {
        "sample": ["SAM_4", "SAM_5", "SAM_6"],
        "disease": ["True", "False", "True"],
    },
    index=["cell_4", "cell_5", "cell_6"],
)
se2 = SummarizedExperiment(
    assays={
        "counts": np.random.poisson(lam=5, size=(3, 3)),
        "lognorm": np.random.lognormal(size=(3, 3))
    },
    row_data=rowData2,
    column_data=colData2,
    metadata={"seq_platform": "Illumina NovaSeq 6000"},
)

rowData3 = pd.DataFrame(
    {
        "seqnames": ["chr_7", "chr_1", "chr_Y"],
        "start": [1084390, 1874937, 243879798],
        "end": [243895239, 358908298, 390820395]
    },
    index=["MYC", "BRCA2", "TPFK"],
)
colData3 = pd.DataFrame(
    {
        "sample": ["SAM_7", "SAM_8", "SAM_9"],
        "disease": ["True", "False", "False"],
        "doublet_score": [.15, .62, .18]
    },
    index=["cell_7", "cell_8", "cell_9"],
)
se3 = SummarizedExperiment(
    assays={
        "counts": np.random.poisson(lam=5, size=(3, 3)),
        "lognorm": np.random.lognormal(size=(3, 3)),
        "beta": np.random.beta(a=1, b=1, size=(3, 3))
    },
    row_data=rowData3,
    column_data=colData3,
    metadata={"seq_platform": "Illumina NovaSeq 6000"},
)

print(se1)
print(se2)
print(se3)
```

```
class: SummarizedExperiment
dimensions: (3, 3)
assays(2): ['counts', 'lognorm']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(3): ['HER2', 'BRCA1', 'TPFK']
column_data columns(2): ['sample', 'disease']
column_names(3): ['cell_1', 'cell_2', 'cell_3']
metadata(1): seq_type

class: SummarizedExperiment
dimensions: (3, 3)
assays(2): ['counts', 'lognorm']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(3): ['HER2', 'BRCA1', 'TPFK']
column_data columns(2): ['sample', 'disease']
column_names(3): ['cell_4', 'cell_5', 'cell_6']
metadata(1): seq_platform

class: SummarizedExperiment
dimensions: (3, 3)
assays(3): ['counts', 'lognorm', 'beta']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(3): ['MYC', 'BRCA2', 'TPFK']
column_data columns(3): ['sample', 'disease', 'doublet_score']
column_names(3): ['cell_7', 'cell_8', 'cell_9']
metadata(1): seq_platform
```

Important

The `combine_rows` or `combine_columns` operations, expect all experiments to contain the same assay names.

To combine experiments by row:

```
from biocutils import relaxed_combine_columns, combine_columns, combine_rows, relaxed_combine_rows
se_combined = combine_rows(se2, se1)
print(se_combined)
```

```
class: SummarizedExperiment
dimensions: (6, 3)
assays(2): ['lognorm', 'counts']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(6): ['HER2', 'BRCA1', 'TPFK', 'HER2', 'BRCA1', 'TPFK']
column_data columns(2): ['sample', 'disease']
column_names(3): ['cell_4', 'cell_5', 'cell_6']
metadata(1): seq_platform
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/summarizedexperiment/BaseSE.py:96: UserWarning: 'row_data' does not contain unique 'row_names'.
  warn("'row_data' does not contain unique 'row_names'.", UserWarning)
```

Similarly to combine by column:

```
se_combined = combine_columns(se2, se1)
print(se_combined)
```

```
class: SummarizedExperiment
dimensions: (3, 6)
assays(2): ['lognorm', 'counts']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(3): ['HER2', 'BRCA1', 'TPFK']
column_data columns(2): ['sample', 'disease']
column_names(6): ['cell_4', 'cell_5', 'cell_6', 'cell_1', 'cell_2', 'cell_3']
metadata(1): seq_platform
```

Important

You can use `relaxed_combine_columns` or `relaxed_combined_rows` when there’s mismatch in the number of features or samples. Missing rows or columns in any object are filled in with appropriate placeholder values before combining, e.g. missing assay’s are replaced with a masked numpy array.

```
# se3 contains an additional assay not present in se1
se_relaxed_combine = relaxed_combine_columns(se3, se1)
print(se_relaxed_combine)
```

```
class: SummarizedExperiment
dimensions: (3, 6)
assays(3): ['lognorm', 'beta', 'counts']
row_data columns(3): ['seqnames', 'start', 'end']
row_names(3): ['MYC', 'BRCA2', 'TPFK']
column_data columns(3): ['sample', 'disease', 'doublet_score']
column_names(6): ['cell_7', 'cell_8', 'cell_9', 'cell_1', 'cell_2', 'cell_3']
metadata(1): seq_platform
```

Both these classes can also contain no experimental data, and they tend to be useful when integrated into more extensive data structures but do not contain any data themselves.

To create an empty `SummarizedExperiment`:

```
empty_se = SummarizedExperiment()
print(empty_se)
```

```
class: SummarizedExperiment
dimensions: (0, 0)
assays(0): []
row_data columns(0): []
row_names(0):  
column_data columns(0): []
column_names(0):  
metadata(0):
```

Similarly an empty `RangeSummarizedExperiment`:

```
empty_rse = RangedSummarizedExperiment()
print(empty_rse)
```

```
class: RangedSummarizedExperiment
dimensions: (0, 0)
assays(0): []
row_data columns(0): []
row_names(0):  
column_data columns(0): []
column_names(0):  
metadata(0):
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

---

* Check out reference [documentation](https://biocpy.github.io/SummarizedExperiment/index.html) for more details.
* R/Bioconductor’s [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) package.