context("Study_loc")

test_that("Platform name could be retrieved", {
  study = create_study()
  study$platform_name = "GPL570"
  platform_name = study$get_platform_name()
  expect_equal(platform_name, "GPL570")
})

test_that("Experimental grouping could be retrieved", {
  study = create_study()
  exp_grp_filename = system.file(
    "extdata/GSE26471", 
    "GSE26471_exp_grp.tab.gz", 
    package = "epimeddata"
  )
  study$exp_grp = read.table(file=gzfile(exp_grp_filename), stringsAsFactors=FALSE)
  exp_grp = study$get_exp_grp()
  expect_equal(dim(exp_grp), c(1,36))
})

test_that("Data could be retrieved", {
  study = create_study()
  data_filename = system.file(
    "extdata/GSE26471", 
    "GSE26471_data.tab.gz", 
    package = "epimeddata"
  )
  study$data = as.matrix(read.table(file=gzfile(data_filename)))
  data = study$get_data()
  expect_equal(dim(data), c(54675,1))
})

# test_that("Platform description could be retrieved", {
#   study = create_study()
#   study$platform_filename = "../../inst/extdata/GSE26471/GPL570.soft.gz"
#   expect_equal(dim(study$get_platform()), c(54675,16))
# })

test_that("Study_loc objects could be plot", {
  study = create_study()
  study$platform_name = "GPL570"
  exp_grp_filename = system.file(
    "extdata/GSE26471", 
    "GSE26471_exp_grp.tab.gz", 
    package = "epimeddata"
  )
  study$exp_grp = read.table(file=gzfile(exp_grp_filename), stringsAsFactors=FALSE)
  data_filename = system.file(
    "extdata/GSE26471", 
    "GSE26471_data.tab.gz", 
    package = "epimeddata"
  )
  study$data = as.matrix(read.table(file=gzfile(data_filename)))
  expect_silent(study$plot_qc())
})


