Source URL: https://bioconductor.org/books/release/BiocBookDemo/
Date Scraped: 2025-12-15

---

![](assets/cover.png "Write Quarto books with Bioconductor")

**Package:** BiocBookDemo  
**Authors:** Jacques Serizay [aut, cre]  
**Compiled:** 2025-11-03  
**Package version:** 1.8.0  
**R version:** **R version 4.5.1 Patched (2025-08-23 r88802)**  
**BioC version:** **3.22**  
**License:** MIT + file LICENSE

# What are `BiocBook`s?

`BiocBook`s are **package-based, versioned online books** with a **supporting `Docker` image** for each book version.

A `BiocBook` can be created by authors (e.g. `R` developers, but also scientists, teachers, communicators, …) who wish to:

1. *Write*: compile a **body of biological and/or bioinformatics knowledge**;
2. *Containerize*: provide **Docker images** to reproduce the examples illustrated in the compendium;
3. *Publish*: deploy an **online book** to disseminate the compendium;
4. *Version*: **automatically** generate specific online book versions and Docker images for specific [Bioconductor releases](https://contributions.bioconductor.org/use-devel.html).

Tip

A {`BiocBook`}-based package hosted on **GitHub** with a branch named `RELEASE_X_Y` provides:

* A **Docker image**: hosted on [ghcr.io](https://github.com/features/packages);
* An **online book** (a.k.a website): hosted on the GitHub repository `gh-pages` branch;

Both are built against the specific Bioconductor release `X.Y`.

A {`BiocBook`}-based package submitted to **Bioconductor** also lead to the online book being independently built by the **Bioconductor Build System (BBS)** and deployed to `https://bioconductor.org/books/<bioc_version>/<pkg>/`

![](pages/images/workflow.jpg)

# What is this {`BiocBook`} package?

The [{BiocBook} **package**](https://github.com/js2264/BiocBook) offers a streamlined approach to creating `BiocBook`s, with several important benefits:

* The author creates a {`BiocBook`}-based package without leaving R;
* The author writes book chapters in `pages/*.qmd` files using enhanced markdown;
* The author can submit its {`BiocBook`}-based package to Bioconductor.

The containerization and publishing of the new {`BiocBook`}-based package is automated:

* A Github Action workflow **generates different Docker images** for different Bioconductor releases, with the packages used in the book pre-installed;
* A Github Action workflow **publishes different book versions** for different Bioconductor releases.

# Main features of `BiocBook`s

## Fully compatible with the *Bioconductor Build System*

When a {`BiocBook`}-based package is accepted into Bioconductor, it is automatically integrated into the *Bionconductor Build System* (BBS).

This means that it is getting built using `R CMD build --keep-empty-dirs --no-resave-data .`. This triggers the rendering of the book contained in `/inst/`. Book packages built by the BBS are then automatically deployed and are eventually available at `https://bioconductor.org/books/<bioc_version>/<pkg>/`.

## Automated versioning of Docker images

A separate Docker image is built for each branch (named `devel` or `RELEASE_X_Y`) of a {`BiocBook`}-based **Github repository**.

Each Docker image provides pre-installed `R` packages:

* Bioconductor release `X.Y`;
* Specific book dependencies from Bioconductor release `X.Y` (listed in `DESCRIPTION`);
* The book package itself

The Docker images also include a `micromamba`-based environment, named `BiocBook`, in which all the softwares listed in `requirements.yml` are installed.

For example, `Docker` images built from the {`BiocBookDemo`} package repository are available here:

👉 [ghcr.io/js2264/biocbookdemo](https://ghcr.io/js2264/biocbookdemo) 🐳

Get started now 🎉

You can get access to all the packages used in this book in < 1 minute, using this command in a terminal:

```
bash
```

```
docker run -it ghcr.io/js2264/biocbookdemo:devel R
```

## Automated versioning of the online book

Regardless of whether the book package is submitted to Bioconductor, a Github Actions workflow publishes individual online books for each branch (named `devel` or `RELEASE_X_Y`) of a `BiocBook`-based **Github repository**.

For example, the online book version matching the `devel` version of the {`BiocBook`} package is available from:

👉 <http://js2264.github.io/BiocBookDemo/devel/> 📘

## RStudio Server

An RStudio Server instance based on a specific Bioconductor `<version>` (`devel` or `RELEASE_X_Y`) can be initiated from the corresponding `Docker` image as follows:

```
bash
```

```
docker run \
    --volume <local_folder>:<destination_folder> \
    -e PASSWORD=OHCA \
    -p 8787:8787 \
    ghcr.io/<github_user>/<biocbook_repo>:<version>
```

The initiated RStudio Server instance will be available at [https://localhost:8787](https://localhost:8787/).

Further instructions regarding Bioconductor-based Docker images are available [here](https://bioconductor.org/help/docker/).

# Acknowledgments

This works was inspired by and closely follows the strategy used in coordination by the Bioconductor core team and Aaron Lun to submit book-containing packages (from the `OSCA` series as well as `SingleR` and `csaw` books).

* *Orchestrating Single-Cell Analysis with Bioconductor* ([n.d.](https://bioconductor.org/books/release/BiocBookDemo/#ref-OSCA))
* *Assigning cell types with SingleR* ([2023](https://bioconductor.org/books/release/BiocBookDemo/#ref-SingleR))
* *The csaw Book* ([2023](https://bioconductor.org/books/release/BiocBookDemo/#ref-csaw))

This package was also inspired by the `*down` package series, including:

* Xie ([2014](https://bioconductor.org/books/release/BiocBookDemo/#ref-knitr))
* Xie ([2016](https://bioconductor.org/books/release/BiocBookDemo/#ref-bookdown))
* Wickham, Hesselberth, and Salmon ([2022](https://bioconductor.org/books/release/BiocBookDemo/#ref-pkgdown))

# References