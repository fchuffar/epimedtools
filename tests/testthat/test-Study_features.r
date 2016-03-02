context("Study_features")

test_that("do_gm2sd_analysis works", {
  study = create_study()
  s = get_fake_study()
  study$data = s$data
  study$exp_grp = s$exp_grp

  probe_names = rownames(study$get_data())
  # subset dataset
  ctrl_exp_grp_key = "tabac"
  case_exp_grp_key = "tabac"
  ctrl_factor_name = "Non Smoker"
  case_factor_name = "Smoker"
  # case_factor_name = "Placenta"

  gm2sd = study$do_gm2sd_analysis(probe_names, ctrl_exp_grp_key, case_exp_grp_key, ctrl_factor_name, case_factor_name, ctrl_thres_func=m2sd, case_value_func=mean, comp_func=get("<"), nb_perm=100, MONITORED=FALSE)
  expect_equal(dim(gm2sd), c(10,3))

  gm2sd0 = study$do_gm2sd_analysis(probe_names, ctrl_exp_grp_key, case_exp_grp_key, ctrl_factor_name, case_factor_name, ctrl_thres_func=m2sd, case_value_func=mean, comp_func=get("<"), nb_perm=0, MONITORED=FALSE)
  expect_equal(dim(gm2sd0), c(10,3))
  expect_equal(gm2sd0[1, 2], 1)
})

test_that("do_mw_test works", {
  study = create_study()
  s = get_fake_study()
  study$data = s$data
  study$exp_grp = s$exp_grp

  probe_names = rownames(study$get_data())
  # subset dataset
  ctrl_exp_grp_key = "tabac"
  case_exp_grp_key = "tabac"
  ctrl_factor_name = "Non Smoker"
  case_factor_name = "Smoker"
  # case_factor_name = "Placenta"

  mw = study$do_mw_test(probe_names, ctrl_exp_grp_key, case_exp_grp_key, ctrl_factor_name, case_factor_name)
  expect_equal(dim(mw), c(10,2))
})

test_that("do_pca works", {
  study = create_study()
  s = get_fake_study()
  study$data = s$data
  study$exp_grp = s$exp_grp

  probe_names = rownames(study$get_data())
  # subset dataset
  ctrl_exp_grp_key = "tabac"
  case_exp_grp_key = "tabac"
  ctrl_factor_name = "Non Smoker"
  case_factor_name = "Smoker"
  # case_factor_name = "Placenta"

  pca_res = study$do_pca()
  expect_silent(study$plot_pca_eig(pca_res))
  expect_silent(study$plot_pca_pc(pca_res, col=1))
})
