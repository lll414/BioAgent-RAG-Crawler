---
source_url: https://r4ds.hadley.nz/import.html
download_date: 2025-11-20
---

# Import

In this part of the book, you’ll learn how to import a wider range of data into R, as well as how to get it into a form useful for analysis. Sometimes this is just a matter of calling a function from the appropriate data import package. But in more complex cases it might require both tidying and transformation in order to get to the tidy rectangle that you’d prefer to work with.

![Our data science model with import highlighted in blue. ](diagrams/data-science/import.png)

Figure 1: Data import is the beginning of the data science process; without data you can’t do data science!

In this part of the book you’ll learn how to access data stored in the following ways:

* In [20  Spreadsheets](spreadsheets.html), you’ll learn how to import data from Excel spreadsheets and Google Sheets.
* In [21  Databases](databases.html), you’ll learn about getting data out of a database and into R (and you’ll also learn a little about how to get data out of R and into a database).
* In [22  Arrow](arrow.html), you’ll learn about Arrow, a powerful tool for working with out-of-memory data, particularly when it’s stored in the parquet format.
* In [23  Hierarchical data](rectangling.html), you’ll learn how to work with hierarchical data, including the deeply nested lists produced by data stored in the JSON format.
* In [24  Web scraping](webscraping.html), you’ll learn web “scraping”, the art and science of extracting data from web pages.

There are two important tidyverse packages that we don’t discuss here: haven and xml2. If you’re working with data from SPSS, Stata, and SAS files, check out the **haven** package, <https://haven.tidyverse.org>. If you’re working with XML data, check out the **xml2** package, <https://xml2.r-lib.org>. Otherwise, you’ll need to do some research to figure which package you’ll need to use; google is your friend here 😃.
