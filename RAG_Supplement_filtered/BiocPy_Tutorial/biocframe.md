Source URL: https://biocpy.github.io/tutorial/chapters/representations/biocframe.html
Date Scraped: 2025-12-15

---

`BiocFrame` class is a Bioconductor-friendly data frame class. Its primary advantage lies in not making assumptions about the types of the columns - as long as an object has a length (`__len__`) and supports slicing methods (`__getitem__`), it can be used inside a `BiocFrame`.

This flexibility allows us to accept arbitrarily complex objects as columns, which is often the case in Bioconductor objects. Also check out Bioconductor’s [**S4Vectors**](https://bioconductor.org/packages/S4Vectors) package, which implements the `DFrame` class on which `BiocFrame` was based.

Note

These classes follow a functional paradigm for accessing or setting properties, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.

To get started, install the package from [PyPI](https://pypi.org/project/biocframe/)

```
pip install biocframe
```

One of the core principles guiding the implementation of the `BiocFrame` class is “***what you put is what you get***”. Unlike Pandas `DataFrame`, `BiocFrame` makes no assumptions about the types of the columns provided as input. Some key differences to highlight the advantages of using `BiocFrame` are especially in terms of modifications to column types and handling nested dataframes.

As an example, Pandas `DataFrame` modifies the types of the input data. These assumptions may cause issues when interoperating between R and Python.

```
import pandas as pd
import numpy as np
from array import array

df = pd.DataFrame({
    "numpy_vec": np.zeros(10),
    "list_vec": [1]* 10,
    "native_array_vec": array('d', [3.14] * 10) # less used but native python arrays
})

print("type of numpy_vector column:", type(df["numpy_vec"]), df["numpy_vec"].dtype)
print("type of list_vector column:", type(df["list_vec"]), df["list_vec"].dtype)
print("type of native_array_vector column:", type(df["native_array_vec"]), df["native_array_vec"].dtype)

print(df)
```

```
type of numpy_vector column: <class 'pandas.core.series.Series'> float64
type of list_vector column: <class 'pandas.core.series.Series'> int64
type of native_array_vector column: <class 'pandas.core.series.Series'> float64
   numpy_vec  list_vec  native_array_vec
0        0.0         1              3.14
1        0.0         1              3.14
2        0.0         1              3.14
3        0.0         1              3.14
4        0.0         1              3.14
5        0.0         1              3.14
6        0.0         1              3.14
7        0.0         1              3.14
8        0.0         1              3.14
9        0.0         1              3.14
```

With `BiocFrame`, no assumptions are made, and the input data is not cast into (un)expected types:

```
from biocframe import BiocFrame
import numpy as np
from array import array

bframe_types = BiocFrame({
    "numpy_vec": np.zeros(10),
    "list_vec": [1]* 10,
    "native_array_vec": array('d', [3.14] * 10)
})

print("type of numpy_vector column:", type(bframe_types["numpy_vec"]))
print("type of list_vector column:", type(bframe_types["list_vec"]))
print("type of native_array_vector column:", type(bframe_types["native_array_vec"]))

print(bframe_types)
```

```
type of numpy_vector column: <class 'numpy.ndarray'>
type of list_vector column: <class 'list'>
type of native_array_vector column: <class 'array.array'>
BiocFrame with 10 rows and 3 columns
             numpy_vec list_vec native_array_vec
    <ndarray[float64]>   <list>          <array>
[0]                0.0        1             3.14
[1]                0.0        1             3.14
[2]                0.0        1             3.14
[3]                0.0        1             3.14
[4]                0.0        1             3.14
[5]                0.0        1             3.14
[6]                0.0        1             3.14
[7]                0.0        1             3.14
[8]                0.0        1             3.14
[9]                0.0        1             3.14
```

Note

This behavior remains consistent when extracting, slicing, combining, or performing any supported operation on `BiocFrame` objects.

Pandas `DataFrame` does not support nested structures; therefore, running the snippet below will result in an error:

```
df = pd.DataFrame({
    "ensembl": ["ENS00001", "ENS00002", "ENS00002"],
    "symbol": ["MAP1A", "BIN1", "ESR1"],
    "ranges": pd.DataFrame({
        "chr": ["chr1", "chr2", "chr3"],
        "start": [1000, 1100, 5000],
        "end": [1100, 4000, 5500]
    }),
})
print(df)
```

However, it is handled seamlessly with `BiocFrame`:

```
bframe_nested = BiocFrame({
    "ensembl": ["ENS00001", "ENS00002", "ENS00002"],
    "symbol": ["MAP1A", "BIN1", "ESR1"],
    "ranges": BiocFrame({
        "chr": ["chr1", "chr2", "chr3"],
        "start": [1000, 1100, 5000],
        "end": [1100, 4000, 5500]
    }),
})

print(bframe_nested)
```

```
BiocFrame with 3 rows and 3 columns
     ensembl symbol         ranges
      <list> <list>    <BiocFrame>
[0] ENS00001  MAP1A chr1:1000:1100
[1] ENS00002   BIN1 chr2:1100:4000
[2] ENS00002   ESR1 chr3:5000:5500
```

Note

This behavior remains consistent when extracting, slicing, combining, or performing any other supported operations on `BiocFrame` objects.

Creating a `BiocFrame` object is straightforward; just provide the `data` as a dictionary.

```
from biocframe import BiocFrame

obj = {
    "ensembl": ["ENS00001", "ENS00002", "ENS00003"],
    "symbol": ["MAP1A", "BIN1", "ESR1"],
}
bframe = BiocFrame(obj)
print(bframe)
```

```
BiocFrame with 3 rows and 2 columns
     ensembl symbol
      <list> <list>
[0] ENS00001  MAP1A
[1] ENS00002   BIN1
[2] ENS00003   ESR1
```

Tip

You can specify complex objects as columns, as long as they have some “length” equal to the number of rows. For example, we can embed a `BiocFrame` within another `BiocFrame`.

```
obj = {
    "ensembl": ["ENS00001", "ENS00002", "ENS00002"],
    "symbol": ["MAP1A", "BIN1", "ESR1"],
    "ranges": BiocFrame({
        "chr": ["chr1", "chr2", "chr3"],
        "start": [1000, 1100, 5000],
        "end": [1100, 4000, 5500]
    }),
}

bframe2 = BiocFrame(obj, row_names=["row1", "row2", "row3"])
print(bframe2)
```

```
BiocFrame with 3 rows and 3 columns
      ensembl symbol         ranges
       <list> <list>    <BiocFrame>
row1 ENS00001  MAP1A chr1:1000:1100
row2 ENS00002   BIN1 chr2:1100:4000
row3 ENS00002   ESR1 chr3:5000:5500
```

The `row_names` parameter is analogous to index in the pandas world and should not contain missing strings. Additionally, you may provide:

* `column_data`: A `BiocFrame`object containing metadata about the columns. This must have the same number of rows as the numbers of columns.
* `metadata`: Additional metadata about the object, usually a dictionary.
* `column_names`: If different from the keys in the `data`. If not provided, this is automatically extracted from the keys in the `data`.

`BiocFrame` is intended for accurate representation of Bioconductor objects for interoperability with R, many users may prefer working with **pandas** `DataFrame` objects for their actual analyses. This conversion is easily achieved:

```
from biocframe import BiocFrame
bframe3 = BiocFrame(
    {
        "foo": ["A", "B", "C", "D", "E"],
        "bar": [True, False, True, False, True]
    }
)

df = bframe3.to_pandas()
print(type(df))
print(df)
```

```
<class 'pandas.core.frame.DataFrame'>
  foo    bar
0   A   True
1   B  False
2   C   True
3   D  False
4   E   True
```

Converting back to a `BiocFrame` is similarly straightforward:

```
out = BiocFrame.from_pandas(df)
print(out)
```

```
BiocFrame with 5 rows and 2 columns
     foo    bar
  <list> <list>
0      A   True
1      B  False
2      C   True
3      D  False
4      E   True
```

BiocPy classes follow a functional paradigm for accessing or setting properties, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.

Properties can be directly accessed from the object:

```
print("shape:", bframe.shape)
print("column names (functional style):", bframe.get_column_names())
print("column names (as property):", bframe.column_names) # same as above
```

```
shape: (3, 2)
column names (functional style): ['ensembl', 'symbol']
column names (as property): ['ensembl', 'symbol']
```

We can fetch individual columns:

```
print("functional style:", bframe.get_column("ensembl"))
print("w/ accessor", bframe["ensembl"])
```

```
functional style: ['ENS00001', 'ENS00002', 'ENS00003']
w/ accessor ['ENS00001', 'ENS00002', 'ENS00003']
```

And we can get individual rows as a dictionary:

```
bframe.get_row(2)
```

```
{'ensembl': 'ENS00003', 'symbol': 'ESR1'}
```

To retrieve a subset of the data in the `BiocFrame`, we use the subset (`[]`) operator. This operator accepts different subsetting arguments, such as a boolean vector, a `slice` object, a sequence of indices, or row/column names.

```
sliced_with_bools = bframe[1:2, [True, False, False]]
print("Subset using booleans: \n", sliced_with_bools)

sliced_with_names = bframe[[0,2], ["symbol", "ensembl"]]
print("\nSubset using column names: \n", sliced_with_names)

# Short-hand to get a single column:
print("\nShort-hand to get a single column: \n", bframe["ensembl"])
```

```
Subset using booleans: 
 BiocFrame with 1 row and 1 column
     ensembl
      <list>
[0] ENS00002

Subset using column names: 
 BiocFrame with 2 rows and 2 columns
    symbol  ensembl
    <list>   <list>
[0]  MAP1A ENS00001
[1]   ESR1 ENS00003

Short-hand to get a single column: 
 ['ENS00001', 'ENS00002', 'ENS00003']
```

For setting properties, we encourage a **functional style** of programming to avoid mutating the object directly. This helps prevent inadvertent modifications of `BiocFrame` instances within larger data structures.

```
modified = bframe.set_column_names(["column1", "column2"])
print(modified)
```

```
BiocFrame with 3 rows and 2 columns
     column1 column2
      <list>  <list>
[0] ENS00001   MAP1A
[1] ENS00002    BIN1
[2] ENS00003    ESR1
```

Now let’s check the column names of the original object,

```
# Original is unchanged:
print(bframe.get_column_names())
```

```
['ensembl', 'symbol']
```

To add new columns, or replace existing ones:

```
modified = bframe.set_column("symbol", ["A", "B", "C"])
print(modified)

modified = bframe.set_column("new_col_name", range(2, 5))
print(modified)
```

```
BiocFrame with 3 rows and 2 columns
     ensembl symbol
      <list> <list>
[0] ENS00001      A
[1] ENS00002      B
[2] ENS00003      C
BiocFrame with 3 rows and 3 columns
     ensembl symbol new_col_name
      <list> <list>      <range>
[0] ENS00001  MAP1A            2
[1] ENS00002   BIN1            3
[2] ENS00003   ESR1            4
```

Change the row or column names:

```
modified = bframe.\
    set_column_names(["FOO", "BAR"]).\
    set_row_names(['alpha', 'bravo', 'charlie'])
print(modified)
```

```
BiocFrame with 3 rows and 2 columns
             FOO    BAR
          <list> <list>
  alpha ENS00001  MAP1A
  bravo ENS00002   BIN1
charlie ENS00003   ESR1
```

Tip

The functional style allows you to chain multiple operations.

We also support Bioconductor’s metadata concepts, either along the columns or for the entire object:

```
modified = bframe.\
    set_metadata({ "author": "Jayaram Kancherla" }).\
    set_column_data(BiocFrame({"column_source": ["Ensembl", "HGNC" ]}))
print(modified)
```

```
BiocFrame with 3 rows and 2 columns
     ensembl symbol
      <list> <list>
[0] ENS00001  MAP1A
[1] ENS00002   BIN1
[2] ENS00003   ESR1
------
column_data(1): column_source
metadata(1): author
```

Properties can also be set by direct assignment for in-place modification. We prefer not to do it this way as it can silently mutate `BiocFrame` instances inside other data structures. Nonetheless:

```
testframe = BiocFrame({ "A": [1,2,3], "B": [4,5,6] })
testframe.column_names = ["column1", "column2" ]
print(testframe)
```

```
BiocFrame with 3 rows and 2 columns
    column1 column2
     <list>  <list>
[0]       1       4
[1]       2       5
[2]       3       6
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/biocframe/BiocFrame.py:468: UserWarning: Setting property 'column_names' is an in-place operation, use 'set_column_names' instead
  warn(
```

Caution

Warnings are raised when properties are directly mutated. These assignments are the same as calling the corresponding `set_*()` methods with `in_place = True`. It is best to do this only if the `BiocFrame` object is not being used anywhere else; otherwise, it is safer to just create a (shallow) copy via the default `in_place = False`.

Similarly, we could set or replace columns directly:

```
testframe["column2"] = ["A", "B", "C"]
testframe[1:3, ["column1","column2"]] = BiocFrame({"x":[4, 5], "y":["E", "F"]})
```

```
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/biocframe/BiocFrame.py:833: UserWarning: This method performs an in-place operation, use 'set_column' instead
  warn(
/opt/hostedtoolcache/Python/3.9.18/x64/lib/python3.9/site-packages/biocframe/BiocFrame.py:827: UserWarning: This method performs an in-place operation, use 'set_slice' instead
  warn(
```

You can iterate over the rows of a `BiocFrame` object. `name` is `None` if the object does not contain any `row_names`. To iterate over the first two rows:

```
for name, row in bframe[:2,]:
    print(name, row)
```

```
None {'ensembl': 'ENS00001', 'symbol': 'MAP1A'}
None {'ensembl': 'ENS00002', 'symbol': 'BIN1'}
```

`BiocFrame` implements methods for the various `combine` generics from [**BiocUtils**](https://github.com/BiocPy/biocutils). For example, to combine by row:

```
import biocutils

bframe1 = BiocFrame({
    "odd": [1, 3, 5, 7, 9],
    "even": [0, 2, 4, 6, 8],
})

bframe2 = BiocFrame({
    "odd": [11, 33, 55, 77, 99],
    "even": [0, 22, 44, 66, 88],
})

combined = biocutils.combine_rows(bframe1, bframe2)
print(combined)
```

```
BiocFrame with 10 rows and 2 columns
       odd   even
    <list> <list>
[0]      1      0
[1]      3      2
[2]      5      4
[3]      7      6
[4]      9      8
[5]     11      0
[6]     33     22
[7]     55     44
[8]     77     66
[9]     99     88
```

Similarly, to combine by column:

```
bframe3 = BiocFrame({
    "foo": ["A", "B", "C", "D", "E"],
    "bar": [True, False, True, False, True]
})

combined = biocutils.combine_columns(bframe1, bframe3)
print(combined)
```

```
BiocFrame with 5 rows and 4 columns
       odd   even    foo    bar
    <list> <list> <list> <list>
[0]      1      0      A   True
[1]      3      2      B  False
[2]      5      4      C   True
[3]      7      6      D  False
[4]      9      8      E   True
```

By default, the combine methods assume that the number and identity of columns (for `combine_rows()`) or rows (for `combine_columns()`) are the same across objects. In situations where this is not the case, such as having different columns across objects, we can use `relaxed_combine_rows()` instead:

```
from biocframe import relaxed_combine_rows

modified2 = bframe2.set_column("foo", ["A", "B", "C", "D", "E"])

combined = biocutils.relaxed_combine_rows(bframe1, modified2)
print(combined)
```

```
BiocFrame with 10 rows and 3 columns
       odd   even    foo
    <list> <list> <list>
[0]      1      0   None
[1]      3      2   None
[2]      5      4   None
[3]      7      6   None
[4]      9      8   None
[5]     11      0      A
[6]     33     22      B
[7]     55     44      C
[8]     77     66      D
[9]     99     88      E
```

Similarly, if the rows are different, we can use `BiocFrame`’s `merge` function. This function uses the *row\_names* as the index to perform this operation; you can specify an alternative set of keys through the `by` parameter.

```
from biocframe import merge

modified1 = bframe1.set_row_names(["A", "B", "C", "D", "E"])
modified3 = bframe3.set_row_names(["C", "D", "E", "F", "G"])

combined = merge([modified1, modified3], by=None, join="outer")
print(combined)
```

```
BiocFrame with 7 rows and 4 columns
     odd   even    foo    bar
  <list> <list> <list> <list>
A      1      0   None   None
B      3      2   None   None
C      5      4      A   True
D      7      6      B  False
E      9      8      C   True
F   None   None      D  False
G   None   None      E   True
```

We can create empty `BiocFrame` objects that only specify the number of rows. This is beneficial in scenarios where `BiocFrame` objects are incorporated into larger data structures but do not contain any data themselves.

```
empty = BiocFrame(number_of_rows=100)
print(empty)
```

```
BiocFrame with 100 rows and 0 columns
```

Most operations detailed in this document can be performed on an empty `BiocFrame` object.

```
print("Column names:", empty.column_names)

subset_empty = empty[1:10,:]
print("\nSubsetting an empty BiocFrame: \n", subset_empty)
```

```
Column names: []

Subsetting an empty BiocFrame: 
 BiocFrame with 9 rows and 0 columns
```

---

* Explore more details [the reference documentation](https://biocpy.github.io/BiocFrame/).
* Additionally, take a look at Bioconductor’s [**S4Vectors**](https://bioconductor.org/packages/S4Vectors) package, which implements the `DFrame` class upon which `BiocFrame` was built.