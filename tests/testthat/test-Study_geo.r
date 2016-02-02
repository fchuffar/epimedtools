context("Study_geo")

test_that("Data are Laoded using GSE number", {
  study = create_study()
  study$gse = "GSE26471"
  extdata_dir = system.file("extdata", package = "epimedtools")
  print(extdata_dir)
  gset = study$get_gset(dest_dir=extdata_dir)
  expect_equal(dim(gset), c(Features=54675, Samples=1))
})

test_that("Data are loaded using series_matrix_filename", {
  study = create_study()
  series_matrix_filename = system.file(
    "extdata/GSE26471", 
    "GSE26471_series_matrix.txt.gz", 
    package = "epimedtools"
  )
  study$series_matrix_filename = series_matrix_filename
  gset = study$get_gset()
  expect_equal(dim(gset), c(Features=54675, Samples=1))
})

test_that("Error is return when no GSE number is given", {
  study = create_study()
  expect_error(study$get_gset())
})

test_that("Error is return when a wrong series_matrix_filename is given", {
  study = create_study()
  expect_error(study$get_gset())
  study$series_matrix_filename = "path/to/a/wrong/series_matrix.txt.gz"
  expect_error(study$get_gset())
})

test_that("Platform name could be retrieved", {
  study = create_study()
  study$gse = "GSE26471"
  extdata_dir = system.file("extdata", package = "epimedtools")
  platform_name = study$get_platform_name(dest_dir=extdata_dir)
  expect_equal(platform_name, "GPL570")
})

test_that("Experimental grouping could be retrieved", {
  study = create_study()
  study$gse = "GSE26471"
  extdata_dir = system.file("extdata", package = "epimedtools")
  exp_grp = study$get_exp_grp(dest_dir=extdata_dir)
  expect_equal(dim(exp_grp), c(1,36))
})

test_that("Data could be retrieved", {
  study = create_study()
  study$gse = "GSE26471"
  extdata_dir = system.file("extdata", package = "epimedtools")
  data = study$get_data(dest_dir=extdata_dir)
  expect_equal(dim(data), c(54675,1))
})


