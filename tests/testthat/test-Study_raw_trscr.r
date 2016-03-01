context("Study_raw_trscr")

test_that(".CEL.gz files could be retrieve from geo and exp_grp properly fused.", {
  # Case_staudy from GEO
  case_study_dir = "/tmp"
  study_case = create_study()
  study_case$gse = "GSE49506"
  cel_files = study_case$get_cel_files(dest_dir=case_study_dir)
  case_exp_grp = study_case$get_exp_grp()
  case_cel_dir = paste(case_study_dir, "/", study_case$gse, sep="")
  # Go!
  study = create_study()
  study$exp_grp = case_exp_grp
  # ctrl_cel_filedir = "../../inst/extdata/trscr_raw_ctrl"
  ctrl_cel_filedir = system.file(
    "extdata/trscr_raw_ctrl", 
    package = "epimeddata"
  )
  study$cel_filedirs = c(case_cel_dir, ctrl_cel_filedir)
  library(affy)
  data = study$get_data()
  exp_grp = study$get_exp_grp()
  expect_equal(dim(data), c(54675,8))
  expect_equal(dim(exp_grp), c(8,37))
})

test_that("Study_raw_trscr$data could be compute from .CEL.gz", {
  study = create_study()
  kc_cel_filedir = system.file(
    "extdata/trscr_raw_kc", 
    package = "epimeddata"
  )
  ctrl_cel_filedir = system.file(
    "extdata/trscr_raw_ctrl", 
    package = "epimeddata"
  )
  study$cel_filedirs = c(kc_cel_filedir, ctrl_cel_filedir)
  library(affy)
  data = study$get_data()
  expect_equal(dim(data), c(54675,6))
  exp_grp = study$get_exp_grp()
  expect_equal(dim(exp_grp), c(6,1))
  expect_equal(sum(rownames(exp_grp) %in% colnames(data)), 6)
  ratio = study$get_ratio()
  expect_equal(dim(ratio), c(54675,6))
  probe_names=rownames(data)[1:10]
  exp_grp_key="orig"
  ctrl_name="trscr_raw_ctrl"
  nb_perm=10
  m2s = study$do_m2s_analysis(probe_names=probe_names, exp_grp_key=exp_grp_key, ctrl_name=ctrl_name, nb_perm=nb_perm)
  expect_equal(dim(m2s), c(10,7))
})

test_that("that hooks could be activated.", {
  study = create_study()
  kc_cel_filedir = system.file(
    "extdata/trscr_raw_kc", 
    package = "epimeddata"
  )
  ctrl_cel_filedir = system.file(
    "extdata/trscr_raw_ctrl", 
    package = "epimeddata"
  )
  study$cel_filedirs = c(kc_cel_filedir, ctrl_cel_filedir)
  expect_error(study$get_data(CACHE=FALSE, hook = "stop", "It's just a hook test!"), regexp="It's just a hook test!")
})

test_that("exp_grp could be merge", {
  #exp_grp1  
  exp_grp1_filename = system.file(
    "extdata/trscr_raw_kc", 
    "expgrp_kc.csv.gz", 
    package = "epimeddata"
  )
  exp_grp1 = read.table(file=gzfile(exp_grp1_filename), stringsAsFactors=FALSE, sep=";", header=TRUE)
  rownames(exp_grp1) = paste(exp_grp1$sample, ".CEL.gz", sep="")
  #exp_grp2
  kc_cel_filedir = system.file(
    "extdata/trscr_raw_kc", 
    package = "epimeddata"
  )
  ctrl_cel_filedir = system.file(
    "extdata/trscr_raw_ctrl", 
    package = "epimeddata"
  )
  study = create_study()
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
    package = "epimeddata"
  )
  exp_grp3 = read.table(file=gzfile(exp_grp_filename), stringsAsFactors=FALSE)
  exp_grp3$orig="GSE26471"
  exp_grp3 = exp_grp3[36:37]
  fused_exp_grp2 = fuse_exp_grp(exp_grp2, exp_grp3, by=c("row.names", "orig"))
  expect_equal(dim(fused_exp_grp2), c(7,17))
})




