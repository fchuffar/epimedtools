## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(epimedtools)

## ------------------------------------------------------------------------
study = create_study()

## ------------------------------------------------------------------------
kc_cel_filedir = "../inst/extdata/trscr_raw_kc"
ctrl_cel_filedir = "../inst/extdata/trscr_raw_ctrl"
study$cel_filedirs = c(kc_cel_filedir, ctrl_cel_filedir)
print(study$cel_filedirs)

## ------------------------------------------------------------------------
head(study$get_data())

## ------------------------------------------------------------------------
print(study$get_exp_grp())

## ------------------------------------------------------------------------
# selecting intersting probes
probe_names = rownames(study$get_data())[1:10]
# subset dataset
ctrl_exp_grp_key = "orig"
case_exp_grp_key = "orig"
ctrl_factor_name = "trscr_raw_ctrl"
case_factor_name = "trscr_raw_kc"

mw_pval = study$do_mw_test(probe_names, ctrl_exp_grp_key, case_exp_grp_key, ctrl_factor_name, case_factor_name)


