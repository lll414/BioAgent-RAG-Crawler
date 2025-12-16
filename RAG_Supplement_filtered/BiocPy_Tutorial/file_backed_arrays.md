Source URL: https://biocpy.github.io/tutorial/chapters/representations/file_backed_arrays.html
Date Scraped: 2025-12-15

---

This is the Python equivalent of Bioconductor’s [HDF5Array](https://bioconductor.org/packages/HDF5Array) package, providing a representation of HDF5-backed arrays within the [delayedarray framework](https://biocpy.github.io/tutorial/chapters/representations/delayed_arrays.html). The idea is to allow users to store, manipulate and operate on large datasets without loading them into memory, in a manner that is trivially compatible with other data structures in the [BiocPy ecosystem](https://github.com/BiocPy/).

Note

These classes follow a functional paradigm for accessing or setting properties, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.

This package is published to [PyPI](https://pypi.org/project/delayedarray/) and can be installed via the usual methods:

```
pip install hdf5array tiledbarray
```

Let’s mock up a dense array:

```
import numpy
data = numpy.random.rand(40, 50)

import h5py
with h5py.File("whee.h5", "w") as handle:
    handle.create_dataset("yay", data=data)
```

We can now represent it as a `Hdf5DenseArray`:

```
import hdf5array
arr = hdf5array.Hdf5DenseArray("whee.h5", "yay", native_order=True)
arr
```

```
<40 x 50> Hdf5DenseArray object of type 'float64'
[[0.6137166 , 0.527271  , 0.93094047, ..., 0.89824058, 0.3751469 ,
  0.93230924],
 [0.97032952, 0.17491716, 0.9653401 , ..., 0.16868539, 0.86533484,
  0.86432493],
 [0.25271991, 0.21226941, 0.63550637, ..., 0.68211344, 0.25971667,
  0.51211434],
 ...,
 [0.98859864, 0.74363204, 0.47433832, ..., 0.94123875, 0.45731422,
  0.72858729],
 [0.72920257, 0.41178416, 0.79736108, ..., 0.19227502, 0.49228222,
  0.09636427],
 [0.56070914, 0.39983991, 0.48079399, ..., 0.76084754, 0.43936791,
  0.91062151]]
```

This is just a subclass of a `DelayedArray` and can be used anywhere in the BiocPy framework. Parts of the NumPy API are also supported - for example, we could apply a variety of [delayed operations](https://biocpy.github.io/tutorial/chapters/representations/delayed_arrays.html):

```
scaling = numpy.random.rand(100)
transformed = numpy.log1p(arr / scaling)
transformed
```

```
<40 x 50> DelayedArray object of type 'float64'
[[8.48910063, 0.45282309, 0.86869162, ..., 0.64237837, 0.36916216,
  0.69562291],
 [8.9471275 , 0.17395586, 0.88991513, ..., 0.15632231, 0.70802249,
  0.65838789],
 [7.60214298, 0.20748297, 0.66507891, ..., 0.52129451, 0.26936297,
  0.43955721],
 ...,
 [8.9657778 , 0.59209249, 0.53361019, ..., 0.66481297, 0.43458535,
  0.57962064],
 [8.66148636, 0.36969894, 0.78172189, ..., 0.17635798, 0.46117954,
  0.09882514],
 [8.39878913, 0.36069396, 0.5392223 , ..., 0.56712022, 0.420657  ,
  0.68389445]]
```

Note

Check out the [documentation](https://biocpy.github.io/hdf5array/) for more details.

We support a variety of compressed sparse formats where the non-zero elements are held inside three separate datasets - usually `data`, `indices` and `indptr`, based on the [10X Genomics sparse HDF5 format](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices). To demonstrate, let’s mock up some sparse data using **scipy**:

```
import scipy.sparse
mock = scipy.sparse.random(1000, 200, 0.1).tocsc()

with h5py.File("sparse_whee.h5", "w") as handle:
    handle.create_dataset("sparse_blah/data", data=mock.data, compression="gzip")
    handle.create_dataset("sparse_blah/indices", data=mock.indices, compression="gzip")
    handle.create_dataset("sparse_blah/indptr", data=mock.indptr, compression="gzip")
```

We can then create a sparse HDF5-backed matrix.

Note

Note that there is some variation in this HDF5 compressed sparse format, notably where the dimensions are stored and whether it is column/row-major.

The constructor will not do any auto-detection so we need to provide this information explicitly:

```
import hdf5array
arr = hdf5array.Hdf5CompressedSparseMatrix(
    "sparse_whee.h5", 
    "sparse_blah", 
    shape=(100, 200), 
    by_column=True
)
arr
```

```
<100 x 200> sparse Hdf5CompressedSparseMatrix object of type 'float64'
[[0.        , 0.        , 0.        , ..., 0.        , 0.        ,
  0.        ],
 [0.        , 0.        , 0.        , ..., 0.        , 0.        ,
  0.        ],
 [0.        , 0.72135589, 0.        , ..., 0.64795697, 0.        ,
  0.        ],
 ...,
 [0.        , 0.        , 0.        , ..., 0.        , 0.86939557,
  0.95820628],
 [0.57983664, 0.        , 0.        , ..., 0.        , 0.        ,
  0.        ],
 [0.        , 0.        , 0.        , ..., 0.        , 0.80669413,
  0.        ]]
```

Let’s mock up a dense array:

```
import numpy
import tiledb

data = numpy.random.rand(40, 50)
tiledb.from_numpy("dense.tiledb", data)
```

```
DenseArray(uri='dense.tiledb', mode=r, ndim=2)
```

We can now represent it as a `TileDbArray`:

```
import tiledbarray
arr = tiledbarray.TileDbArray("dense.tiledb", attribute_name="")
```

This is just a subclass of a `DelayedArray` and can be used anywhere in the BiocPy framework. Parts of the NumPy API are also supported - for example, we could apply a variety of delayed operations:

```
scaling = numpy.random.rand(100)
transformed = numpy.log1p(arr / scaling)
transformed
```

```
<40 x 50> DelayedArray object of type 'float64'
[[0.38570168, 1.34329305, 0.37061962, ..., 0.18149232, 0.48812983,
  3.36119847],
 [0.01007671, 1.1368437 , 0.9035832 , ..., 0.16488582, 0.60526886,
  3.20586654],
 [0.48338776, 0.71156784, 0.72885705, ..., 0.69442619, 0.48340426,
  3.71077405],
 ...,
 [0.31601366, 0.35049416, 0.58996299, ..., 0.51546406, 0.14903065,
  3.70188055],
 [0.14035231, 0.85252485, 0.67058939, ..., 0.67844625, 0.02391537,
  3.32386021],
 [0.69106365, 0.11115126, 0.23616254, ..., 0.69022239, 0.08169617,
  3.77476317]]
```

Check out the [documentation](https://biocpy.github.io/tiledbarray/) for more details.

We can perform similar operations on a sparse matrix as well. Lets mock a sparse matrix and store it as a tiledb file.

```
dir_path = "sparse_array.tiledb"
dom = tiledb.Domain(
     tiledb.Dim(name="rows", domain=(0, 4), tile=4, dtype=numpy.int32),
     tiledb.Dim(name="cols", domain=(0, 4), tile=4, dtype=numpy.int32),
)
schema = tiledb.ArraySchema(
     domain=dom, sparse=True, attrs=[tiledb.Attr(name="", dtype=numpy.int32)]
)
tiledb.SparseArray.create(f"{dir_path}", schema)

tdb = tiledb.SparseArray(f"{dir_path}", mode="w")
i, j = [1, 2, 2], [1, 4, 3]
data = numpy.array(([1, 2, 3]))
tdb[i, j] = data
```

We can now represent this as a `TileDbArray`:

```
import tiledbarray
arr = tiledbarray.TileDbArray(dir_path, attribute_name="")

slices = (slice(0, 2), [2, 3])

import delayedarray
subset = delayedarray.extract_sparse_array(arr, (*slices,))
subset
```

```
<2 x 2> SparseNdarray object of type 'int32'
[[0, 0],
 [0, 0]]
```

---

* Just like [delayedarrays](https://biocpy.github.io/tutorial/chapters/representations/delayed_arrays.html#interoperability-with-biocpy-packages), you can use file-backed matrices to represent experimental data in [summarized experiment](https://biocpy.github.io/tutorial/chapters/experiments/summarized_experiment.html#delayed-arrays) and its derivates, reducing the *in-memory* footprint to load large datasets.
* Check out the [documentation](https://biocpy.github.io/hdf5array/) for more details.