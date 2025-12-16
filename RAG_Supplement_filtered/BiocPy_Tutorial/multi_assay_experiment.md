Source URL: https://biocpy.github.io/tutorial/chapters/experiments/multi_assay_experiment.html
Date Scraped: 2025-12-15

---

`MultiAssayExperiment` (MAE) simplifies the management of multiple experimental assays conducted on a shared set of specimens.

Note

These classes follow a functional paradigm for accessing or setting properties, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.

To get started, install the package from [PyPI](https://pypi.org/project/multiassayexperiment/)

```
pip install multiassayexperiment
```

An MAE contains three main entities,

* **Primary information** (`column_data`): Bio-specimen/sample information. The `column_data` may provide information about patients, cell lines, or other biological units. Each row in this table represents an independent biological unit. It must contain an `index` that maps to the ‘primary’ in `sample_map`.
* **Experiments** (`experiments`): Genomic data from each experiment. either a `SingleCellExperiment`, `SummarizedExperiment`, `RangedSummarizedExperiment` or any class that extends a `SummarizedExperiment`.
* **Sample Map** (`sample_map`): Map biological units from `column_data` to the list of `experiments`. Must contain columns,

  + **assay** provides the names of the different experiments performed on the biological units. All experiment names from experiments must be present in this column.
  + **primary** contains the sample name. All names in this column must match with row labels from col\_data.
  + **colname** is the mapping of samples/cells within each experiment back to its biosample information in col\_data.

  Each sample in `column_data` may map to one or more columns per assay.

Let’s start by first creating few experiments:

Show the code

```
from random import random

import numpy as np
from biocframe import BiocFrame
from genomicranges import GenomicRanges
from iranges import IRanges

nrows = 200
ncols = 6
counts = np.random.rand(nrows, ncols)
gr = GenomicRanges(
    seqnames=[
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
        ] * 20,
    ranges=IRanges(range(100, 300), range(110, 310)),
    strand = ["-", "+", "+", "*", "*", "+", "+", "+", "-", "-"] * 20,
    mcols=BiocFrame({
        "score": range(0, 200),
        "GC": [random() for _ in range(10)] * 20,
    })
)

col_data_sce = BiocFrame({"treatment": ["ChIP", "Input"] * 3},
    row_names=[f"sce_{i}" for i in range(6)],
)

col_data_se = BiocFrame({"treatment": ["ChIP", "Input"] * 3},
    row_names=[f"se_{i}" for i in range(6)],
)
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

More importantly, we need to provide `sample_map` information:

```
sample_map = BiocFrame({
    "assay": ["sce", "se"] * 6,
    "primary": ["sample1", "sample2"] * 6,
    "colname": ["sce_0", "se_0", "sce_1", "se_1", "sce_2", "se_2", "sce_3", "se_3", "sce_4", "se_4", "sce_5", "se_5"]
})

sample_data = BiocFrame({"samples": ["sample1", "sample2"]}, row_names= ["sample1", "sample2"])

print(sample_map)
```

```
BiocFrame with 12 rows and 3 columns
      assay primary colname
     <list>  <list>  <list>
 [0]    sce sample1   sce_0
 [1]     se sample2    se_0
 [2]    sce sample1   sce_1
        ...     ...     ...
 [9]     se sample2    se_4
[10]    sce sample1   sce_5
[11]     se sample2    se_5
```

Finally, we can create an `MultiAssayExperiment` object:

```
from multiassayexperiment import MultiAssayExperiment
from singlecellexperiment import SingleCellExperiment
from summarizedexperiment import SummarizedExperiment

tsce = SingleCellExperiment(
    assays={"counts": counts}, row_data=gr.to_pandas(), column_data=col_data_sce
)

tse2 = SummarizedExperiment(
    assays={"counts": counts.copy()},
    row_data=gr.to_pandas().copy(),
    column_data=col_data_se.copy(),
)

mae = MultiAssayExperiment(
    experiments={"sce": tsce, "se": tse2},
    column_data=sample_data,
    sample_map=sample_map,
    metadata={"could be": "anything"},
)

print(mae)
```

```
class: MultiAssayExperiment containing 2 experiments
[0] sce: SingleCellExperiment with 200 rows and 6 columns 
[1] se: SummarizedExperiment with 200 rows and 6 columns 
column_data columns(1): ['samples']
sample_map columns(3): ['assay', 'primary', 'colname']
metadata(1): could be
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

If both `column_data` and `sample_map` are `None`, the constructor naively creates sample mapping, with each `experiment` considered to be a independent `sample`. We add a sample to `column_data` in this pattern - `unknown_sample_{experiment_name}`.

All cells from the each experiment are considered to be from the same sample and is reflected in `sample_map`.

Important

***This is not a recommended approach, but if you don’t have sample mapping, then it doesn’t matter***.

```
mae = MultiAssayExperiment(
    experiments={"sce": tsce, "se": tse2},
    metadata={"could be": "anything"},
)

print(mae)
```

```
class: MultiAssayExperiment containing 2 experiments
[0] sce: SingleCellExperiment with 200 rows and 6 columns 
[1] se: SummarizedExperiment with 200 rows and 6 columns 
column_data columns(1): ['samples']
sample_map columns(3): ['assay', 'primary', 'colname']
metadata(1): could be
```

We provide convenient methods to easily convert a `MuData` object into an `MultiAssayExperiment`.

Let’s create a mudata object:

Show the code

```
import numpy as np
from anndata import AnnData

np.random.seed(1)

n, d, k = 1000, 100, 10

z = np.random.normal(loc=np.arange(k), scale=np.arange(k) * 2, size=(n, k))
w = np.random.normal(size=(d, k))
y = np.dot(z, w.T)

adata = AnnData(y)
adata.obs_names = [f"obs_{i+1}" for i in range(n)]
adata.var_names = [f"var_{j+1}" for j in range(d)]

d2 = 50
w2 = np.random.normal(size=(d2, k))
y2 = np.dot(z, w2.T)

adata2 = AnnData(y2)
adata2.obs_names = [f"obs_{i+1}" for i in range(n)]
adata2.var_names = [f"var2_{j+1}" for j in range(d2)]

from mudata import MuData
mdata = MuData({"rna": adata, "spatial": adata2})

print(mdata)
```

```
MuData object with n_obs × n_vars = 1000 × 150
  2 modalities
    rna:    1000 x 100
    spatial:    1000 x 50
```

Lets convert this object to an MAE:

```
from multiassayexperiment import MultiAssayExperiment

mae_obj = MultiAssayExperiment.from_mudata(input=mdata)
print(mae_obj)
```

```
class: MultiAssayExperiment containing 2 experiments
[0] rna: SingleCellExperiment with 100 rows and 1000 columns 
[1] spatial: SingleCellExperiment with 50 rows and 1000 columns 
column_data columns(1): ['samples']
sample_map columns(3): ['assay', 'primary', 'colname']
metadata(0):
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

Getters are available to access various attributes using either the property notation or functional style.

```
# access assays
print("experiment names (as property): ", mae.experiment_names)
print("experiment names (functional style): ", mae.get_experiment_names())

# access sample data
print(mae.column_data)
```

```
experiment names (as property):  ['sce', 'se']
experiment names (functional style):  ['sce', 'se']
BiocFrame with 2 rows and 1 column
                              samples
                               <list>
unknown_sample_sce unknown_sample_sce
 unknown_sample_se  unknown_sample_se
```

Check out the c[lass documentation](https://biocpy.github.io/MultiAssayExperiment/api/multiassayexperiment.html#multiassayexperiment.MultiAssayExperiment.MultiAssayExperiment) for the full list of accessors and setters.

A helper method is available to easily access row or column names across all experiments. This method returns a dictionary with experiment names as keys and the corresponding values, which can be either the row or column names depending on the function:

```
from rich import print as pprint
pprint("row names:", mae.get_row_names())
pprint("column names:", mae.get_column_names())
```

```
row names:
{
    'sce': Names(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17',
'18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', 
'37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', 
'56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', 
'75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', 
'94', '95', '96', '97', '98', '99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', 
'111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '124', '125', '126', 
'127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '140', '141', '142', 
'143', '144', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157', '158', 
'159', '160', '161', '162', '163', '164', '165', '166', '167', '168', '169', '170', '171', '172', '173', '174', 
'175', '176', '177', '178', '179', '180', '181', '182', '183', '184', '185', '186', '187', '188', '189', '190', 
'191', '192', '193', '194', '195', '196', '197', '198', '199']),
    'se': Names(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', 
'18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', 
'37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', 
'56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', 
'75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', 
'94', '95', '96', '97', '98', '99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', 
'111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '124', '125', '126', 
'127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '140', '141', '142', 
'143', '144', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157', '158', 
'159', '160', '161', '162', '163', '164', '165', '166', '167', '168', '169', '170', '171', '172', '173', '174', 
'175', '176', '177', '178', '179', '180', '181', '182', '183', '184', '185', '186', '187', '188', '189', '190', 
'191', '192', '193', '194', '195', '196', '197', '198', '199'])
}
```

```
column names:
{
    'sce': Names(['sce_0', 'sce_1', 'sce_2', 'sce_3', 'sce_4', 'sce_5']),
    'se': Names(['se_0', 'se_1', 'se_2', 'se_3', 'se_4', 'se_5'])
}
```

One can access an experiment by name:

```
print(mae.experiment("se"))
```

```
class: SummarizedExperiment
dimensions: (200, 6)
assays(1): ['counts']
row_data columns(7): ['starts', 'widths', 'ends', 'seqnames', 'strand', 'score', 'GC']
row_names(200): ['0', '1', '2', ..., '197', '198', '199']
column_data columns(1): ['treatment']
column_names(6): ['se_0', 'se_1', 'se_2', 'se_3', 'se_4', 'se_5']
metadata(0):
```

Additionally you may access an experiment with the sample information included in the column data of the experiment:

Note

This creates a copy of the experiment.

```
expt_with_sample_info = mae.experiment("se", with_sample_data=True)
print(expt_with_sample_info)
```

```
class: SummarizedExperiment
dimensions: (200, 6)
assays(1): ['counts']
row_data columns(7): ['starts', 'widths', 'ends', 'seqnames', 'strand', 'score', 'GC']
row_names(200): ['0', '1', '2', ..., '197', '198', '199']
column_data columns(5): ['assay', 'primary', 'colname', 'treatment', 'samples']
column_names(6): ['se_0', 'se_1', 'se_2', 'se_3', 'se_4', 'se_5']
metadata(0):
```

Note

For consistency with the R MAE’s interface, we also provide `get_with_column_data` method, that performs the same operation.

Important

All property-based setters are `in_place` operations, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.

```
modified_column_data = mae.column_data.set_column("score", range(len(mae.column_data)))
modified_mae = mae.set_column_data(modified_column_data)
print(modified_mae)
```

```
class: MultiAssayExperiment containing 2 experiments
[0] sce: SingleCellExperiment with 200 rows and 6 columns 
[1] se: SummarizedExperiment with 200 rows and 6 columns 
column_data columns(2): ['samples', 'score']
sample_map columns(3): ['assay', 'primary', 'colname']
metadata(1): could be
```

Now, lets check the `column_data` on the original object.

```
print(mae.column_data)
```

```
BiocFrame with 2 rows and 1 column
                              samples
                               <list>
unknown_sample_sce unknown_sample_sce
 unknown_sample_se  unknown_sample_se
```

You can subset `MultiAssayExperiment` by using the subset (`[]`) operator. This operation accepts different slice input types, such as a boolean vector, a `slice` object, a list of indices, or names (if available) to subset.

`MultiAssayExperiment` allows subsetting by three dimensions: `rows`, `columns`, and `experiments`. ***`sample_map` is automatically filtered during this operation***.

```
subset_mae = mae[1:5, 0:4]
print(subset_mae)
```

```
class: MultiAssayExperiment containing 2 experiments
[0] sce: SingleCellExperiment with 4 rows and 6 columns 
[1] se: SummarizedExperiment with 4 rows and 6 columns 
column_data columns(1): ['samples']
sample_map columns(3): ['assay', 'primary', 'colname']
metadata(1): could be
```

The following creates a subset based on the experiments dimension:

```
subset_mae = mae[1:5, 0:1, ["se"]]
print(subset_mae)
```

```
class: MultiAssayExperiment containing 1 experiments
[0] se: SummarizedExperiment with 4 rows and 0 columns 
column_data columns(1): ['samples']
sample_map columns(3): ['assay', 'primary', 'colname']
metadata(1): could be
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/multiassayexperiment/MultiAssayExperiment.py:63: UserWarning: 'primary' from `sample_map` & `column_data` mismatch.
  warn("'primary' from `sample_map` & `column_data` mismatch.", UserWarning)
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/multiassayexperiment/MultiAssayExperiment.py:74: UserWarning: 'experiments' contains names not represented in 'sample_map' or vice-versa.
  warn(
```

Note

If you’re wondering about why the experiment “se” has 0 columns, it’s important to note that our MAE implementation does not remove columns from an experiment solely because none of the columns map to the samples of interest. This approach aims to prevent unexpected outcomes in complex subset operations.

The `MultiAssayExperiment` class also provides a few methods for sample management.

The `complete_cases` function is designed to identify samples that contain data across all experiments. It produces a boolean vector with the same length as the number of samples in `column_data`. Each element in the vector is `True` if the sample is present in all experiments, or `False` otherwise.

```
print(mae.complete_cases())
```

```
[False, False]
```

You can use this boolean vector to select samples with complete data across all assays or experiments.

```
subset_mae = mae[:, mae.complete_cases(),]
print(subset_mae)
```

```
class: MultiAssayExperiment containing 2 experiments
[0] sce: SingleCellExperiment with 200 rows and 0 columns 
[1] se: SummarizedExperiment with 200 rows and 0 columns 
column_data columns(1): ['samples']
sample_map columns(3): ['assay', 'primary', 'colname']
metadata(1): could be
```

This method identifies ‘samples’ with replicates within each experiment. The result is a dictionary where experiment names serve as keys, and the corresponding values indicate whether the sample is replicated within each experiment.

```
from rich import print as pprint # mainly for pretty printing
pprint(mae.replicated())
```

```
{
    'sce': {
        'unknown_sample_sce': [True, True, True, True, True, True],
        'unknown_sample_se': [False, False, False, False, False, False]
    },
    'se': {
        'unknown_sample_sce': [False, False, False, False, False, False],
        'unknown_sample_se': [True, True, True, True, True, True]
    }
}
```

The `intersect_rows` finds common `row_names` across all experiments and returns a `MultiAssayExperiment` with those rows.

```
common_rows_mae = mae.intersect_rows()
print(common_rows_mae)
```

```
class: MultiAssayExperiment containing 2 experiments
[0] sce: SingleCellExperiment with 200 rows and 6 columns 
[1] se: SummarizedExperiment with 200 rows and 6 columns 
column_data columns(1): ['samples']
sample_map columns(3): ['assay', 'primary', 'colname']
metadata(1): could be
```

If you are only interested in finding common `row_names` across all experiments:

```
common_rows = mae.find_common_row_names()
print(common_rows)
```

```
{'82', '96', '162', '126', '184', '93', '106', '99', '183', '163', '38', '61', '157', '12', '83', '7', '121', '0', '141', '24', '154', '117', '56', '46', '186', '50', '66', '91', '151', '180', '68', '74', '88', '105', '112', '94', '129', '8', '49', '130', '108', '103', '17', '181', '39', '169', '13', '148', '142', '75', '124', '14', '47', '30', '34', '36', '118', '127', '10', '59', '51', '45', '60', '41', '119', '9', '58', '85', '104', '125', '1', '29', '76', '97', '164', '177', '63', '43', '152', '173', '134', '69', '62', '135', '33', '133', '28', '147', '87', '171', '194', '187', '155', '67', '19', '70', '57', '172', '71', '102', '116', '72', '26', '22', '114', '178', '188', '199', '139', '168', '196', '113', '128', '42', '23', '149', '156', '122', '165', '159', '84', '179', '101', '131', '89', '185', '143', '31', '197', '90', '190', '40', '107', '15', '109', '18', '136', '48', '80', '193', '92', '115', '4', '52', '120', '5', '55', '145', '98', '78', '79', '20', '110', '11', '137', '146', '160', '198', '25', '81', '191', '192', '2', '138', '170', '140', '132', '167', '37', '53', '111', '100', '21', '174', '195', '35', '176', '32', '3', '123', '77', '182', '65', '150', '161', '73', '153', '144', '175', '6', '16', '27', '54', '64', '44', '95', '158', '189', '166', '86'}
```

While the necessity of an empty `MultiAssayExperiment` might not be apparent, for the sake of consistency with the rest of the tutorials:

```
mae = MultiAssayExperiment(experiments={})
print(mae)
```

```
class: MultiAssayExperiment containing 0 experiments
column_data columns(1): ['samples']
sample_map columns(3): ['assay', 'primary', 'colname']
metadata(0):
```

---

* Check out reference [documentation](https://biocpy.github.io/MultiAssayExperiment/index.html) for more details.
* R/Bioconductor’s [MultiAssayExperiment](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html) package.