Source URL: https://biocpy.github.io/tutorial/chapters/extras/iranges.html
Date Scraped: 2025-12-15

---

Python implementation of the [**IRanges**](https://bioconductor.org/packages/IRanges) Bioconductor package.

An `IRanges` holds a **start** position and a **width**, and is typically used to represent coordinates along a genomic sequence. The interpretation of the **start** position depends on the application; for sequences, the **start** is usually a 1-based position, but other use cases may allow zero or even negative values, e.g., circular genomes.

`IRanges` uses [nested containment lists](https://github.com/pyranges/ncls) under the hood to perform fast overlap and search based operations.

Note

These classes follow a functional paradigm for accessing or setting properties, with further details discussed in [functional paradigm](https://biocpy.github.io/tutorial/chapters/philosophy.html#functional-discipline) section.

To get started, install the package from [PyPI](https://pypi.org/project/IRanges/)

```
pip install iranges
```

Note

The descriptions for some of these methods come from the [Bioconductor documentation](https://bioconductor.org/packages/release/bioc/html/IRanges.html).

An `IRanges` holds a **start** position and a **width**, and is most typically used to represent coordinates along some genomic sequence. The interpretation of the start position depends on the application; for sequences, the start is usually a 1-based position, but other use cases may allow zero or even negative values (e.g. circular genomes).

```
from iranges import IRanges

starts = [-2, 6, 9, -4, 1, 0, -6, 10]
widths = [5, 0, 6, 1, 4, 3, 2, 3]
ir = IRanges(starts, widths)

print(ir)
```

```
IRanges object with 8 ranges and 0 metadata columns
               start              end            width
    <ndarray[int32]> <ndarray[int32]> <ndarray[int32]>
[0]               -2                3                5
[1]                6                6                0
[2]                9               15                6
[3]               -4               -3                1
[4]                1                5                4
[5]                0                3                3
[6]               -6               -4                2
[7]               10               13                3
```

Properties can be accessed directly from the object:

```
print("Number of intervals:", len(ir))

print("start positions:", ir.get_start())
print("width of each interval:", ir.get_width())
print("end positions:", ir.get_end())
```

```
Number of intervals: 8
start positions: [-2  6  9 -4  1  0 -6 10]
width of each interval: [5 0 6 1 4 3 2 3]
end positions: [ 3  6 15 -3  5  3 -4 13]
```

Tip

Just like BiocFrame, these classes offer both functional-style and property-based getters and setters.

```
print("start positions:", ir.start)
print("width of each interval:", ir.width)
print("end positions:", ir.end)
```

```
start positions: [-2  6  9 -4  1  0 -6 10]
width of each interval: [5 0 6 1 4 3 2 3]
end positions: [ 3  6 15 -3  5  3 -4 13]
```

`reduce` method reduces the intervals to an `IRanges` where the intervals are:

* not empty
* not overlapping
* ordered from left to right
* not even adjacent (i.e. there must be a non empty gap between 2 consecutive ranges).

```
reduced = ir.reduce()
print(reduced)
```

```
IRanges object with 4 ranges and 0 metadata columns
               start              end            width
    <ndarray[int32]> <ndarray[int32]> <ndarray[int32]>
[0]               -6               -3                3
[1]               -2                5                7
[2]                6                6                0
[3]                9               15                6
```

`IRanges` uses [nested containment lists](https://github.com/pyranges/ncls) under the hood to perform fast overlap and search-based operations.

```
subject = IRanges([2, 2, 10], [1, 2, 3])
query = IRanges([1, 4, 9], [5, 4, 2])

overlap = subject.find_overlaps(query)
print(overlap)
```

```
[[1, 0], [], [2]]
```

The `nearest`, `precede` or `follow` methods finds the nearest overlapping range along the specified direction.

```
query = IRanges([1, 3, 9], [2, 5, 2])
subject = IRanges([3, 5, 12], [1, 2, 1])

nearest = subject.nearest(query, select="all")
print(nearest)
```

```
[[0], [0, 1], [2]]
```

Note

These methods typically return a list of indices from `subject` for each interval in `query`.

The `coverage` method counts the number of overlaps for each position.

```
cov = subject.coverage()
print(cov)
```

```
[0 0 1 0 1 1 0 0 0 0 0 1]
```

`shift` adjusts the start positions by their **shift**.

```
shifted = ir.shift(shift=10)
print(shifted)
```

```
IRanges object with 8 ranges and 0 metadata columns
               start              end            width
    <ndarray[int32]> <ndarray[int32]> <ndarray[int32]>
[0]                8               13                5
[1]               16               16                0
[2]               19               25                6
[3]                6                7                1
[4]               11               15                4
[5]               10               13                3
[6]                4                6                2
[7]               20               23                3
```

Other range transformation methods include `narrow`, `resize`, `flank`, `reflect` and `restrict`. For example `narrow` supports the adjustment of `start`, `end` and `width` values, which should be relative to each range.

```
narrowed = ir.narrow(start=4, width=2)
print(narrowed)
```

```
IRanges object with 8 ranges and 0 metadata columns
               start              end            width
    <ndarray[int64]> <ndarray[int64]> <ndarray[int64]>
[0]                1                3                2
[1]                9               11                2
[2]               12               14                2
[3]               -1                1                2
[4]                4                6                2
[5]                3                5                2
[6]               -3               -1                2
[7]               13               15                2
```

Well as the name says, computes disjoint intervals.

```
disjoint = ir.disjoin()
print(disjoint)
```

```
IRanges object with 9 ranges and 0 metadata columns
               start              end            width
    <ndarray[int32]> <ndarray[int32]> <ndarray[int32]>
[0]               -6               -4                2
[1]               -4               -3                1
[2]               -2                0                2
[3]                0                1                1
[4]                1                3                2
[5]                3                5                2
[6]                9               10                1
[7]               10               13                3
[8]               13               15                2
```

`reflect` reverses each range within a set of common reference bounds.

```
starts = [2, 5, 1]
widths = [2, 3, 3]
x = IRanges(starts, widths)
bounds = IRanges([0, 5, 3], [11, 2, 7])

res = x.reflect(bounds=bounds)
print(res)
```

```
IRanges object with 3 ranges and 0 metadata columns
               start              end            width
    <ndarray[int32]> <ndarray[int32]> <ndarray[int32]>
[0]                7                9                2
[1]                4                7                3
[2]                9               12                3
```

`flank` returns ranges of a specified width that flank, to the left (default) or right, each input range. One use case of this is forming promoter regions for a set of genes.

```
starts = [2, 5, 1]
widths = [2, 3, 3]
x = IRanges(starts, widths)

res = x.flank(2, start=False)
print(res)
```

```
IRanges object with 3 ranges and 0 metadata columns
               start                end              width
    <ndarray[int32]> <ndarray[float64]> <ndarray[float64]>
[0]                4                6.0                2.0
[1]                8               10.0                2.0
[2]                4                6.0                2.0
```

`IRanges` supports most interval set operations. For example, to compute `gaps`:

```
gaps = ir.gaps()
print(gaps)
```

```
IRanges object with 2 ranges and 0 metadata columns
               start              end            width
    <ndarray[int32]> <ndarray[int32]> <ndarray[int32]>
[0]               -3               -2                1
[1]                5                9                4
```

Or Perform interval set operations, e..g `union`, `intersection`, `disjoin`:

```
x = IRanges([1, 5, -2, 0, 14], [10, 5, 6, 12, 4])
y = IRanges([14, 0, -5, 6, 18], [7, 3, 8, 3, 3])

intersection = x.intersect(y)
print(intersection)
```

```
IRanges object with 3 ranges and 0 metadata columns
               start              end            width
    <ndarray[int32]> <ndarray[int32]> <ndarray[int32]>
[0]               -2                3                5
[1]                6                9                3
[2]               14               18                4
```

---

* [IRanges reference](https://biocpy.github.io/IRanges/api/iranges.html#iranges-package)
* [Bioc/IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html)