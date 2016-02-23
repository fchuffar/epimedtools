context("utils")

test_that("monitored_apply and monitored_sapply work as sapply and sapply", {
  mat = matrix(rnorm(50), 10)
  func = function(v) {return(c(mean(v), sd(v)))}
  foo = monitored_apply(mat, 1, func)
  bar = apply(mat, 1, func)
  expect_equivalent(foo, bar)
  func2 = function(i) {return(c(i*i, i+i))}
  foo = monitored_sapply(1:10, func2)
  bar = sapply(1:10, func2)
  expect_equivalent(foo, bar)
})

test_that("simplify_sample_names simplify sample names", {
  sample_names = c("kc_GSM748053", "kc_GSM748054", "kc_GSM748055.CEL.gz", "ctrl_GSM80562.cel.gz",  "ctrl_GSM80576.CEL.gz",  "ctrl_GSM80577.CEL.gz") 
  sample_names = simplify_sample_names(sample_names)
  expect_equal(sum(duplicated(sample_names)), 0)
})

