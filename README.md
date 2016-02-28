# epimedtools

[![Build Status](https://api.travis-ci.org/fchuffar/epimedtools.png?branch=master)](https://travis-ci.org/fchuffar/epimedtools)
[![Coverage Status](https://img.shields.io/codecov/c/github/fchuffar/epimedtools/master.svg)](https://codecov.io/github/fchuffar/epimedtools?branch=master)


The __epimedtools__ package provides set of useful tools encapsulating some 'GEOquery' functions. This package rely on RC, caching features and memoisation to simplify, homogenize and speed up multi-omic analysis.

## Installation

To get the current development version from github:

```R
install.packages("devtools")
devtools::install_github("fchuffar/epimedtools")
```


## Vignettes

To browse available vignettes:

```R
browseVignettes(package = 'epimedtools')
```

## Data Struture

The `epimedtools` package manage a `study` object composed of the three data structure: `data`, `exp_grp` and `platform`.

```R
> study
$data
      ctrl1 ctrl2 ctrl3 ctrl4 ctrl5 ctrl6 case1 case2 case3 case4 case5 case6
prb1  3.773 0.016 3.040 1.356 3.196 3.457 2.539 1.092 2.244 0.000 4.056 2.595
prb2  3.194 3.813 3.519 0.821 3.535 1.454 1.265 4.435 2.541 1.275 2.270 2.336
prb3  3.074 3.098 2.000 1.657 1.874 3.567 3.876 0.313 3.974 5.980 4.072 2.730
prb4  3.743 4.026 2.285 1.631 2.184 2.930 2.975 3.069 2.605 3.533 3.415 1.984
prb5  2.226 3.780 3.081 2.025 2.465 2.610 3.132 3.242 2.255 3.837 3.413 2.373
prb6  3.302 3.323 1.569 3.059 2.212 2.877 3.114 1.632 3.817 2.955 1.943 1.481
prb7  1.857 2.852 2.293 3.708 2.975 5.154 2.461 3.699 2.249 2.662 0.778 3.696
prb8  2.412 3.330 2.970 2.807 2.071 3.148 2.364 2.163 4.072 0.000 3.122 2.198
prb9  2.831 2.611 4.469 3.053 2.237 4.002 1.983 2.724 2.380 2.797 3.132 2.017
prb10 3.244 3.322 3.490 0.788 1.955 3.345 0.280 4.069 3.259 2.175 1.150 3.355

$exp_grp
      sex age tabac treatment histo   orig
ctrl1   F  28   Smk      0 ug  lung France
ctrl2   F  21   Smk      0 ug  lung France
ctrl3   F  20   Smk      0 ug  lung France
ctrl4   F  23 N-Smk      0 ug  lung France
ctrl5   M  21 N-Smk      0 ug  lung France
ctrl6   F  20 N-Smk      0 ug  lung France
case1   M  21   Smk     15 ug  lung France
case2   M  26   Smk     15 ug  lung France
case3   F  21   Smk     15 ug  lung France
case4   F  28 N-Smk     15 ug  lung France
case5   F  21 N-Smk     15 ug  lung France
case6   M  22 N-Smk     15 ug  lung France

$platform
      gene_name GOID
prb1        ...  ...
prb2        ...  ...
prb3        ...  ...
prb4        ...  ...
prb5        ...  ...
prb6        ...  ...
prb7        ...  ...
prb8        ...  ...
prb9        ...  ...
prb10       ...  ...

```