context("Study RC")

test_that("Caching features work", {
  tmp_cache_filename = "tmp_cached_study.Rds"
  if (file.exists(tmp_cache_filename)) {
    file.remove(tmp_cache_filename)
  }
  s1 = create_study(tmp_cache_filename)
  s1$series_matrix_filename = "../../inst/extdata/GSE26471/GSE26471_series_matrix.txt.gz"
  s1$gse = "GSE26471"
  s1$platform_filename = "foo/bar/baz"
  s1$cache_it()
  s2 = create_study(tmp_cache_filename)
  for (f in names(s2$getRefClass()$fields())) {
    expect_equal(s2[[f]], s1[[f]])    
  }
  tmp2_cache_filename = "tmp2_cached_study.Rds"
  if (file.exists(tmp2_cache_filename)) {
    file.remove(tmp2_cache_filename)
  }
  s2$cache_it(tmp2_cache_filename)
  s3 = create_study(tmp2_cache_filename)
  for (f in names(s3$getRefClass()$fields())) {
    expect_equal(s2[[f]], s3[[f]])    
  }  
  expect_true(s2$cache_filename == s3$cache_filename)
  expect_false(s1$cache_filename == s3$cache_filename)
  file.remove(c(tmp2_cache_filename, tmp_cache_filename))
})

test_that("Data are Laoded using GSE number", {
  study = create_study()
  study$gse = "GSE26471"
  gset = study$get_gset(dest_dir="../../inst/extdata")
  expect_equal(dim(gset), c(Features=54675, Samples=1))
})

test_that("Data are loaded using series_matrix_filename", {
  study = create_study()
  study$series_matrix_filename = "../../inst/extdata/GSE26471/GSE26471_series_matrix.txt.gz"
  gset = study$get_gset()
  expect_equal(dim(gset), c(Features=54675, Samples=1))
})

test_that("Platform name could be retrieved", {
  study = create_study()
  study$gse = "GSE26471"
  platform_name = study$get_platform_name(dest_dir="../../inst/extdata")
  expect_equal(platform_name, "GPL570")
})

test_that("Experimental grouping could be retrieved", {
  study = create_study()
  study$gse = "GSE26471"
  exp_grp = study$get_exp_grp(dest_dir="../../inst/extdata")
  expect_equal(dim(exp_grp), c(1,36))
})

test_that("Data could be retrieved", {
  study = create_study()
  study$gse = "GSE26471"
  data = study$get_data(dest_dir="../../inst/extdata")
  expect_equal(dim(data), c(54675,1))
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


