# epimedtools

[![Build Status](https://api.travis-ci.org/fchuffar/epimedtools.png?branch=master)](https://travis-ci.org/fchuffar/epimedtools)
[![Coverage Status](https://img.shields.io/codecov/c/github/fchuffar/epimedtools/master.svg)](https://codecov.io/github/fchuffar/epimedtools?branch=master)


The __epimedtools__ package provides set of useful statistical functions. This package rely on 'GEOquery' and 'affy' packages, R Reference Classes, caching features and memoisation to simplify, homogenize and speed up multi-omic analysis.

## Data Struture

The __epimedtools__ package aims to manage a `study` object composed of the three data structure: the `data` matrix and the two `exp_grp` and `platform` dataframe.

Take care: 

  - study$data is a matrix of numerics
  - `colnames(study$data)` contain `rownames(study$exp_grp)`
  - `rownames(study$data)` contains `rownames(study$platform)`

```R
> study
#' $data
#'       ctrl1 ctrl2 ctrl3 ctrl4 ctrl5 ctrl6 case1 case2 case3 case4 case5 case6
#' prb1  1.801 3.036 1.856 2.717 2.374 3.755 0.000 3.656 0.725 5.043 3.265 3.539
#' prb2  4.369 3.190 3.025 3.591 2.363 2.120 3.503 3.685 2.948 2.751 3.656 2.525
#' prb3  2.650 3.198 2.660 1.827 4.225 3.284 0.677 1.382 3.741 1.568 3.287 2.930
#' prb4  3.546 3.589 1.904 1.925 2.968 2.349 2.608 3.467 2.098 3.541 2.854 1.796
#' prb5  1.721 2.419 3.216 1.936 3.419 2.881 1.959 2.040 2.325 2.973 4.222 2.918
#' prb6  3.963 4.579 4.468 3.106 4.373 3.699 2.575 2.107 3.588 3.329 3.626 0.257
#' prb7  3.305 2.019 4.151 3.240 3.758 3.150 5.418 1.672 3.073 3.425 2.434 0.734
#' prb8  3.657 1.829 5.768 2.217 3.377 2.493 4.815 2.900 4.521 4.502 2.144 3.973
#' prb9  2.158 3.413 3.564 2.754 2.634 2.556 2.868 1.514 2.705 4.202 2.030 1.197
#' prb10 3.390 2.958 2.204 4.523 2.152 4.561 1.629 0.891 2.519 2.783 1.541 2.375
#'
#' $exp_grp
#'          sex age      tabac treatment histo
#' ctrl1 Female  23     Smoker      0 ug  lung
#' ctrl2   Male  28     Smoker      0 ug  lung
#' ctrl3 Female  26     Smoker      0 ug  lung
#' ctrl4 Female  30 Non Smoker      0 ug  lung
#' ctrl5   Male  22 Non Smoker      0 ug  lung
#' ctrl6   Male  23 Non Smoker      0 ug  lung
#' case1   Male  25     Smoker     15 ug  lung
#' case2 Female  27     Smoker     15 ug  lung
#' case3   Male  29     Smoker     15 ug  lung
#' case4 Female  23 Non Smoker     15 ug  lung
#' case5   Male  23 Non Smoker     15 ug  lung
#' case6   Male  25 Non Smoker     15 ug  lung
#'
#' $platform
#'         GB_ACC      Species Gene Symbol ENTREZ_GENE_ID
#' prb1   1U48705 Homo sapiens        DDR1            780
#' prb2   2M87338 Homo sapiens        RFC2           5982
#' prb3   3X51757 Homo sapiens       HSPA6           3310
#' prb4   4X69699 Homo sapiens        PAX8           7849
#' prb5   5L36861 Homo sapiens      GUCA1A           2978
#' prb6   6L13852 Homo sapiens        UBA7           7318
#' prb7   7X55005 Homo sapiens        THRA           7067
#' prb8   8X79510 Homo sapiens      PTPN21          11099
#' prb9   9M21121 Homo sapiens        CCL5           6352
#' prb10  1J02843 Homo sapiens      CYP2E1           1571

```

## Treatment and Analysis Features

The __epimedtools__ package provides not only a set of treatment features like data fetching or normalization, but also genome and genes anlysis tools.

## Installation

To get the current development version from github, you need first to install following packages from bioconductor:

  * ``Biobase`
  * ``affy``
  * ``GEOquery``

Then, install ``epimedtool``:

```R
install.packages("devtools")
devtools::install_github("fchuffar/epimedtools")
```


## Vignettes

To browse available vignettes:

```R
browseVignettes(package = 'epimedtools')
```

