## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(epimedtools)

## ------------------------------------------------------------------------
study = create_study()

## ------------------------------------------------------------------------
series_matrix_filename = system.file(
  "extdata/GSE26471", 
  "GSE26471_series_matrix.txt.gz", 
  package = "epimeddata"
)
print(series_matrix_filename)
study$series_matrix_filename = series_matrix_filename

## ------------------------------------------------------------------------
study$get_gset()

## ------------------------------------------------------------------------
dim(study$get_data())

## ------------------------------------------------------------------------
dim(study$get_exp_grp())

## ---- fig.width=7, fig.height=5------------------------------------------
plot(density(study$get_data()[,1]))

## ------------------------------------------------------------------------
study$save("/tmp/study.rds")

