Source URL: https://bioconductor.org/books/release/BiocBookDemo/pages/biocbook-vs-rebook.html
Date Scraped: 2025-12-15

---

{`rebook`} is another book rendering package currently used by Bioconductor to build book packages such as `OSCA.*` and `SingleRBook`.

## 4.1 Differences with {`rebook`}-based books

`OSCA` books provided by Bioconductor are rendered upon building by the *Bioconductor Build System* (BBS). They rely on {`rebook`} to orchestrate the rendering. Briefly, `OSCA` books have their book pages in `/inst/book/`, while the `vignettes/` folder contains (1) a `Makefile` and (2) a dummy `stub.Rmd` vignette (required to trigger vignette rendering and thus `make`).

When the BBS triggers `OSCA.intro` package building, upon vignette rendering, the `Makefile` triggers the following commands:

```
R
```

```
work.dir <- rebook::bookCache('OSCA.intro')
handle <- rebook::preCompileBook('../inst/book', work.dir=work.dir, desc='../DESCRIPTION')
old.dir <- setwd(work.dir)
bookdown::render_book('index.Rmd')
setwd(old.dir)
rebook::postCompileBook(work.dir=work.dir, final.dir='../inst/doc/book', handle=handle)
```

The resulting book, pre-compiled by {`rebook`} and assembled by {`bookdown`}, is eventually served by Bioconductor from the `/inst/doc/book/` folder.

{`BiocBook`}-based packages follow a strategy similar to that of `OSCA` books: they provide a `Makefile` in the `vignettes/` folder to trigger book rendering when building the package. However, the executed command does not rely on {`rebook`} and {`bookdown`}, but on a `render` command from the `quarto` software.

```
shell
```

```
quarto render ../inst/
mv ../inst/docs ../inst/doc/book
```

The resulting book, fully compiled by native `quarto`, is also located in `/inst/doc/book/` once the package is built (and in `doc/book` in the library directory once installed).

This implies that `quarto` (>= 1.3) has to be installed in the system building the package!

## 4.2 `BiocBook` features missing from `rebook`

* A `BiocBook` can be readily initiated using `BiocBook::init()`;
* It relies on modern `.qmd` files supported by `Quarto`;
* It can work as a standalone Github-hosted package, without necessarily having to be submitted to/built by Bioconductor. The book should be rendered exactly the same way through Github or by the BBS;
* It supports versioning of the online book served by `gh-pages` **through the author Github account**;
* It distributes versioned Dockerfiles **through the author Github account**;
* `BiocBook`-based packages can actually provide fully-fledged functions in `R/`, manual pages and vignettes. They can be installed exactly the same way than other software packages.

## 4.3 `rebook` features missing from `BiocBook`

* Smart reuse of objects generated in one book in another book;

Tip

This can still be achieved **within** a book by saving a data object as an `.rds` file and loading it in a subsequent chapter.

```
isthisworking <- readRDS('isthisworking.rds')
isthisworking
##  [1] "yes"
```

* Native support of cross-references across books.