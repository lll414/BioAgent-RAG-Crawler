Source URL: https://biocpy.github.io/tutorial/chapters/representations/atomics.html
Date Scraped: 2025-12-15

---

The [BiocUtils](https://github.com/BiocPy/BiocUtils) package offers essential utilities designed for universal use across all packages, with a focus on emulating convenient features of base R. In particular, this package addresses challenges associated with Python lists, which lack type specificity, leading to the need for inference when dealing with lists containing booleans, numbers, floats or strings.

To begin using the package, you can install it from [PyPI](https://pypi.org/project/biocutils/)

```
pip install biocutils
```

The package provides several atomic lists that are coerced into appropriate types. These include `BooleanList`, `FloatList`, `NamedList`, `IntegerList`, and `StringList`.

Let’s explore `BooleanList`, which resembles a regular Python list but coercing anything added to it into a `boolean`. Additionally, `None` values are accepted and treated as missing booleans.

This list may also be named (see `NamedList`), which provides dictionary-like functionality.

```
from biocutils import BooleanList, NamedList

x = BooleanList([ True, False, False, True ])
print(x)
```

```
[True, False, False, True]
```

Similarly, one can create atomic lists for other types, such as `FloatList`:

```
from biocutils import FloatList

x = FloatList([ 1.1, 2, 3, 4 ])
print(x)
```

```
[1.1, 2.0, 3.0, 4.0]
```

Accessing these vectors is similar to any other list:

```
print("2nd element:", x[2])

print("reassign value")
x[1] = 50
print("x: ", x)
```

```
2nd element: 3.0
```

```
reassign value
```

```
x:  [1.1, 50.0, 3.0, 4.0]
```

To convert objects back to Python lists:

```
print(list(x))
```

```
[1.1, 50.0, 3.0, 4.0]
```

The `Factor` class is analogous to R’s factor. It comprises a vector of integer `codes`, each corresponding to an index within a list of unique strings (`levels`). The purpose is to encode a list of strings as integers for streamlined numerical analysis.

The most straightforward way to create a `Factor` is from an existing list of strings:

```
from biocutils import Factor

f = Factor.from_sequence(["A", "B", "A", "B", "E"])
print(f)
```

```
Factor of length 5 with 3 levels
values: A, B, A, B, E
levels: A, B, E
ordered: False
```

Alternatively, if you already have a list of `codes` and associated `levels`:

```
f = Factor([0, 1, 2, 0, 2, 4], levels=["A", "B", "C", "D", "E"])
print(f)
```

```
Factor of length 6 with 5 levels
values: A, B, C, A, C, E
levels: A, B, C, D, E
ordered: False
```

To convert a `Factor` back to a Python list:

```
print(list(f))
```

```
['A', 'B', 'C', 'A', 'C', 'E']
```

The Biocutils package introduces a `subset` generic function designed to handle n-dimensional objects, where n > 1 (i.e., objects with a shape property of length greater than 1). When applied, the function first verifies the dimensionality of the input objects. If they are n-dimensional, it invokes `subset_rows()` to perform the subsetting along the first dimension. On the other hand, if the objects are deemed vector-like, the function utilizes `subset_sequence()` for the subsetting operation.

```
from biocutils import subset

x = [1, 2, 3, 4, 5]
print(subset(x, [0, 2, 4]))
```

```
[1, 3, 5]
```

The `combine` generic function in Biocutils is designed to accommodate objects of varying dimensions. It begins by examining the dimensionality of the input objects: if they are n-dimensional for n > 1 (i.e., possessing a shape property of length greater than 1), the function utilizes `combine_rows()` to merge them along the first dimension. Conversely, if the objects exhibit a vector-like structure, the function employs `combine_sequences()` for the combination process.

```
import numpy as np
from biocutils import combine

x = [1, 2, 3]
y = [0.1, 0.2]
xd = np.array(x)

combine(xd, y)
```

```
array([1. , 2. , 3. , 0.1, 0.2])
```

Note

The `combine` generic, usually returns an object that is same type as the first argument.

```
import biocutils
biocutils.match(["A", "C", "E"], ["A", "B", "C", "D", "E"])
```

```
array([0, 2, 4], dtype=int8)
```

```
import biocutils
biocutils.intersect(["A", "B", "C", "D"], ["D", "A", "E"])
```

```
['A', 'D']
```

```
import biocutils
biocutils.union(["A", "B", "C", "D"], ["D", "A", "E"])
```

```
['A', 'B', 'C', 'D', 'E']
```

Checks if all elements of a list or tuple are of the same type.

```
import biocutils
import numpy as np

x = [np.random.rand(3), np.random.rand(3, 2)]
biocutils.is_list_of_type(x, np.ndarray)
```

```
True
```

---

Refer to the [documentation](https://biocpy.github.io/BiocUtils/) for comprehensive details on the functionality the package offers.