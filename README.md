# epimedtools

[![Build Status](https://api.travis-ci.org/fchuffar/epimedtools.png?branch=master)](https://travis-ci.org/fchuffar/epimedtools)
[![Coverage Status](https://img.shields.io/codecov/c/github/fchuffar/epimedtools/master.svg)](https://codecov.io/github/fchuffar/epimedtools?branch=master)


The __epimedtools__ package provides set of useful tools encapsulating some 'GEOquery' functions. This package rely on RC, caching features and memoisation to simplify, homogenize and speed up multi-omic analysis.

## Installation

To get the current development version from github:

```R
install.packages("devtools")
devtools::install_github("fchuffar/epimedtools", build_vignettes=TRUE)
```


## Vignettes

To list available vignettes:

```R
vignette(package="epimedtools")
```

To open a vigentte:
```R
vignette(topic="GSE26471_local",package="epimedtools")
```

