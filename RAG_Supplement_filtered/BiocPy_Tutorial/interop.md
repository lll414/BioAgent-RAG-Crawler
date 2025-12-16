Source URL: https://biocpy.github.io/tutorial/chapters/interop.html
Date Scraped: 2025-12-15

---

The [rds2py](https://github.com/BiocPy/rds2py) package serves as a Python interface to the [rds2cpp](https://github.com/LTLA/rds2cpp) library, enabling direct reading of RDS files within Python. This eliminates the need for additional data conversion tools or intermediate formats, streamlining the transition between Python and R for seamless analysis.

One notable feature is the use of memory views (excluding strings) to access the same memory from C++ in Python, facilitated through Cython. This approach is particularly advantageous for handling large datasets, as it avoids unnecessary duplication of data.

What sets `rds2py` apart from other similar parsers is its capability to read **S4** classes. This unique feature allows the parsing of Bioconductor data types directly from R into Python.

To get started, install the package from [PyPI](https://pypi.org/project/rds2py/)

```
pip install rds2py
```

Reading an RDS file in Python involves a two-step process. First, we parse the serialized **RDS** into a readable Python object, typically a dictionary. This object contains both the ***data*** and relevant ***metadata*** about the structure and internal representation of the R object. Subsequently, we use one of the available functions to convert this object into a Python representation.

Before we begin, let’s create a test dataset from R. In this example, we’ll download the “zeisel brain” dataset from the [scRNAseq package](https://bioconductor.org/packages/release/data/experiment/html/scRNAseq.html). For tutorial purposes, we’ll filter the dataset to the first 1000 rows and save it as an RDS file.

```
library(scRNAseq)
sce <- ZeiselBrainData()
sub <- sce[1:1000,]
saveRDS(sub, "../assets/data/zeisel-brain-subset.rds")
```

Now, we can read the RDS file in Python using the `read_rds` function, which parses the file contents and returns a dictionary of the R object.

```
from rds2py import read_rds

r_object = read_rds("../assets/data/zeisel-brain-subset.rds")

from rich import print as pprint
pprint(r_object) # hiding the response
```

Note

The output of the above code block is hidden to maintain the cleanliness and visual appeal of this document :)

Once we have a realized structure, we can convert this object into useful Python representations. It contains **two** keys: - **data**: If atomic entities, contains the NumPy view of the memory space. - **attributes**: Additional properties available for the object.

The package provides functions to convert these R objects into useful Python representations.

```
from rds2py import as_summarized_experiment

# to convert an robject to SCE
sce = as_summarized_experiment(r_object)

print(sce)
```

```
class: SingleCellExperiment
dimensions: (1000, 3005)
assays(1): ['counts']
row_data columns(1): ['featureType']
row_names(1000): ['0', '1', '2', ..., '997', '998', '999']
column_data columns(10): ['tissue', 'group #', 'total mRNA mol', 'well', 'sex', 'age', 'diameter', 'cell_id', 'level1class', 'level2class']
column_names(3005): ['1772071015_C02', '1772071017_G12', '1772071017_A05', ..., '1772063068_D01', '1772066098_A12', '1772058148_F03']
main_experiment_name:  
reduced_dims(0): []
alternative_experiments(2): ['ERCC', 'repeat']
row_pairs(0): []
column_pairs(0): []
metadata(0):
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/genomicranges/SeqInfo.py:348: UserWarning: 'seqnames' is deprecated, use 'get_seqnames' instead
  warn("'seqnames' is deprecated, use 'get_seqnames' instead", UserWarning)
```

And that’s it! It’s as straightforward as that. The `as_summarized_experiment` function serves as an example of how to convert complex R structures into Python representations. Similarly, the package offers parsers for atomic vectors, lists, sparse/dense matrices, data frames, and most R data structures.

You can continue to convert this object into `AnnData` representation and perform analysis. For more details on `SingleCellExperiment`, refer to the documentation [here](https://biocpy.github.io/tutorial/chapters/experiments/single_cell_experiment.html).

```
sce.to_anndata()
```

```
(AnnData object with n_obs × n_vars = 3005 × 1000
     obs: 'tissue', 'group #', 'total mRNA mol', 'well', 'sex', 'age', 'diameter', 'cell_id', 'level1class', 'level2class'
     var: 'featureType'
     layers: 'counts',
 None)
```

Well, that’s all. Dive in, explore, and create more base representations to encapsulate complex R structures. If you wish to add more representations, we are more than happy to accept contributions.