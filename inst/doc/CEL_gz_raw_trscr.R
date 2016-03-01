## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(epimedtools)

## ------------------------------------------------------------------------
study = create_study()

## ------------------------------------------------------------------------
kc_cel_filedir = system.file(
  "extdata/trscr_raw_kc", 
  package = "epimeddata"
)
print(kc_cel_filedir)
ctrl_cel_filedir = system.file(
  "extdata/trscr_raw_ctrl", 
  package = "epimeddata"
)
print(ctrl_cel_filedir)
study$cel_filedirs = c(kc_cel_filedir, ctrl_cel_filedir)
print(study$cel_filedirs)

## ------------------------------------------------------------------------
exp_grp_filename = system.file(
  "extdata/trscr_raw_kc", 
  "expgrp_kc.csv.gz", 
  package = "epimeddata"
)
exp_grp = read.table(file=gzfile(exp_grp_filename), stringsAsFactors=FALSE, sep=";", header=TRUE)
rownames(exp_grp) = paste(exp_grp$sample, ".CEL.gz", sep="")
study$exp_grp = exp_grp

## ------------------------------------------------------------------------
head(study$get_data())

## ------------------------------------------------------------------------
study$get_exp_grp()

## ------------------------------------------------------------------------
study$platform_name = "GPL570"

## ---- fig.width=7, fig.height=5------------------------------------------
plot(density(study$get_data()[,1]))

## ---- fig.width=7, fig.height=5------------------------------------------
ctrl_sample_names = rownames(study$exp_grp)[study$get_exp_grp()$orig == "trscr_raw_ctrl"]
plot(density(study$get_ratio(ctrl_sample_names)[,1]))

## ------------------------------------------------------------------------
study$save("/tmp/study.rds")

