Source URL: https://biocpy.github.io/tutorial/chapters/representations/
Date Scraped: 2025-12-15

---

BiocPy implements essential foundational data structures that serve as the building blocks for extensive and complex representations. This includes the [BiocFrame package](https://biocpy.github.io/tutorial/chapters/representations/biocframe.html) providing a Bioconductor-like data frame class, and the [GenomicRanges package](https://biocpy.github.io/tutorial/chapters/representations/genomic_ranges.html) to aid in representing genomic regions and facilitating analysis. The [BiocUtils](https://biocpy.github.io/tutorial/chapters/representations/atomics.html) package provides many atomic data type classes, defines generics and efficiently manages most common operations across these packages.

The [biocpy](https://github.com/BiocPy/BiocPy) package serves as a convenient wrapper that installs all the core packages within the ecosystem.

```
pip install biocpy
```

Alternatively, you can install specific packages as required. For example:

```
pip install biocframe # <package-name>
```

To update packages, use the following command:

```
pip install -U biocframe # or <package-name>
```