context("Study_abstract features")

test_that("Data are Laoded using GSE number", {
  study = create_study()
  study$gse = "GSE26471"
  extdata_dir = system.file("extdata", package = "epimedtools")
  ratio = study$get_ratio(dest_dir=extdata_dir)
  expect_equal(dim(ratio), c(54675,1))
})


