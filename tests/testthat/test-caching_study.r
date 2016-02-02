context("Caching Study RC")

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