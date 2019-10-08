
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis build
status](https://travis-ci.org/heike/extracat.svg?branch=master)](https://travis-ci.org/heike/extracat)

# extracat

`extracat` is an R package for working with categorical variables.

## Installation

The package `extracat` has been archived from CRAN, but you can install
it from this github repository:

``` r
devtools::install_github("heike/extracat")
```

## Simple use case

Several data sets are provided within the package, the command `visna`
provides a good first overview of the where missing values are in a data
set:

``` r
library(extracat)
visna(GeneEx)
```

![](man/figures/README-unnamed-chunk-3-1.png)<!-- -->

Using additional sorting methods, we can work out patterns of missing
values:

``` r
# reorder both rows and columns:
visna(GeneEx, sort = "b", sort.method="optile")
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->
