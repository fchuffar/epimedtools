## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(epimedtools)

## ------------------------------------------------------------------------
study = create_study()

## ------------------------------------------------------------------------
study$platform_name = "GPL570"

## ------------------------------------------------------------------------
exp_grp_filename = system.file(
  "extdata/GSE26471", 
  "GSE26471_exp_grp.tab.gz", 
  package = "epimedtools"
)
print(exp_grp_filename)
study$exp_grp = read.table(file=gzfile(exp_grp_filename), stringsAsFactors=FALSE)

## ------------------------------------------------------------------------
data_filename = system.file(
  "extdata/GSE26471", 
  "GSE26471_data.tab.gz", 
  package = "epimedtools"
)
print(data_filename)
study$data = as.matrix(read.table(file=gzfile(data_filename)))

## ------------------------------------------------------------------------
dim(study$get_data())

## ------------------------------------------------------------------------
dim(study$get_exp_grp())

## ---- fig.width=7, fig.height=5------------------------------------------
plot(density(study$get_data()[,1]))

## ------------------------------------------------------------------------
study$save("/tmp/study.rds")

