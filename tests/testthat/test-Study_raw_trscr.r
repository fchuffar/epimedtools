context("Study_raw_trscr")

test_that("Study_loc objects could be plot", {
  study = create_study(Study_RC_name="Study_raw_trscr")
  kc_cel_filedir = system.file(
    "extdata/trscr_raw_kc", 
    package = "epimedtools"
  )
  ctrl_cel_filedir = system.file(
    "extdata/trscr_raw_ctrl", 
    package = "epimedtools"
  )
  study$cel_filedirs = c(kc_cel_filedir, ctrl_cel_filedir)
  data = study$get_data(celfile.path="")
  exp_grp = study$get_exp_grp()
  expect_equal(dim(data), c(54675,6))
  expect_equal(dim(exp_grp), c(6,1))
})


