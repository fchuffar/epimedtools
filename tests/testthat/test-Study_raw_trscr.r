context("Study_raw_trscr")

test_that("Study_raw_trscr$data could be compute from .CEL.gz", {
  study = create_study(Study_RC_name="Study_raw_trscr")
  kc_cel_filedir = "../../inst/extdata/trscr_raw_kc"
  ctrl_cel_filedir = "../../inst/extdata/trscr_raw_ctrl"
  study$cel_filedirs = c(kc_cel_filedir, ctrl_cel_filedir)
  library(affy)
  data = study$get_data()
  exp_grp = study$get_exp_grp()
  ratio = study$get_ratio()
  expect_equal(dim(data), c(54675,6))
  expect_equal(dim(ratio), c(54675,6))
  expect_equal(dim(exp_grp), c(6,1))
})

test_that("exp_grp could be merge", {
  #exp_grp1
  exp_grp1 = read.table(file=gzfile("../../inst/extdata/trscr_raw_kc/expgrp_kc.csv.gz"), stringsAsFactors=FALSE, sep=";", header=TRUE)
  rownames(exp_grp1) = paste(exp_grp1$sample, ".CEL.gz", sep="")
  #exp_grp2
  kc_cel_filedir = "../../inst/extdata/trscr_raw_kc"
  ctrl_cel_filedir = "../../inst/extdata/trscr_raw_ctrl"
  study = create_study(Study_RC_name="Study_raw_trscr")
  study$exp_grp = exp_grp1
  study$cel_filedirs = c(kc_cel_filedir, ctrl_cel_filedir)
  data = study$get_data()
  library(affy)
  exp_grp2 = study$get_exp_grp()
  expect_equal(dim(exp_grp2), c(6,16))
  #exp_grp3
  exp_grp_filename = system.file(
    "extdata/GSE26471",
    "GSE26471_exp_grp.tab.gz",
    package = "epimedtools"
  )
  exp_grp3 = read.table(file=gzfile(exp_grp_filename), stringsAsFactors=FALSE)
  exp_grp3$orig="GSE26471"
  exp_grp3 = exp_grp3[36:37]
  fused_exp_grp2 = fuse_exp_grp(exp_grp2, exp_grp3, by=c("row.names", "orig"))
  expect_equal(dim(fused_exp_grp2), c(7,17))
})




