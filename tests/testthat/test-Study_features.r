context("Study_features")

test_that("compute_survival_table plot_survival_panel perform_anova_gen work", {
  study = get_fake_study() 

  suffix = "tabac"
  exp_grp_key = "tabac"
  sample_names = rownames(study$exp_grp[!is.na(study$exp_grp$ss),])
  probe_names = rownames(study$data)
  USE_CACHE=FALSE
  PLOT_SCURVE=FALSE

  survival_res = compute_survival_table(probe_names, sample_names, exp_grp_key, study, suffix, USE_CACHE=USE_CACHE)

  anova_res = perform_anova_gen(study$exp_grp, val~tabac, study$data)

  mw_res = study$do_mw_test(probe_names, exp_grp_key, ctrl_fctr="Non_Smoker")

  probe_name = probe_names[1]
  anova_mw_res = cbind(anova_res[probe_names,], mw_res[probe_names,])
  gene_pf_colname="gene_name"
  plot_survival_panel(probe_name, sample_names, exp_grp_key, study, anova_mw_res=anova_mw_res, gene_pf_colname=gene_pf_colname)

  case_fctr = "Smoker"
  mw_pval="mw_pval_adj"
  main = ""
  fc_thres = c(-1.5, -1.2, 1.2, 1.5)
  # plot_hm(exp_grp_key, case_fctr, anova_mw_res, study, main=main, mw_pval=mw_pval, fc_thres=fc_thres)
  res_key_lr = "logratio_tabac_Smoker"
  res_key_pval = "adj_pval_tabac"
  plot_volcano(res_key_lr, res_key_pval, anova_mw_res)   

  expect_equal(dim(survival_res), c(nrow(study$data),5))
  expect_equal(dim(anova_res), c(nrow(study$data),10))
})



test_that("do_gm2sd_analysis works", {
  study = get_fake_study()

  probe_names = rownames(study$get_data())
  # subset dataset
  ctrl_key = "tabac"
  case_key = "tabac"
  ctrl_fctr = "Non_Smoker"
  case_fctr = "Smoker"
  # case_fctr = "Placenta"

  gm2sd = study$do_gm2sd_analysis(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr, ctrl_thres_func=m2sd, case_value_func=mean, comp_func=get("<"), nb_perm=100, MONITORED=FALSE)
  # print(gm2sd)
  expect_equal(dim(gm2sd), c(10,3))

  gm2sd0 = study$do_gm2sd_analysis(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr, ctrl_thres_func=m2sd, case_value_func=mean, comp_func=get("<"), nb_perm=0, MONITORED=FALSE)
  expect_equal(dim(gm2sd0), c(10,3))
  expect_equal(gm2sd0[1, 2], 1)
})

test_that("do_mw_test works", {
  study = get_fake_study()

  probe_names = rownames(study$get_data())
  # subset dataset
  ctrl_key = "tabac"
  case_key = "tabac"
  ctrl_fctr = "Non_Smoker"
  case_fctr = "Smoker"

  mw = study$do_mw_test(probe_names, ctrl_key, case_key, ctrl_fctr, case_fctr)
  # print(mw)
  expect_equal(dim(mw), c(10,4))
})


test_that("simplified_XXX_names do the job.", {
  colnames_mw_res = c("histo_other.histo_neuro.mean_fc.res", "histo_other.histo_neuro.mw_pval.res")
  simplified_colnames_mw_res = simplify_column_names(colnames_mw_res)
  expect_equal(simplified_colnames_mw_res, c("histo_other_histo_neuro_mean_fc_res", "histo_other_histo_neuro_mw_pval_res"))
  # factor_names = c("histo_other_histo_neuro_mean_fc_res", "histo_other_histo_neuro_mw_pval_res")
  simplified_factnames_mw_res = simplify_factor_names(simplified_colnames_mw_res, "_")
  expect_equal(simplified_factnames_mw_res, c("mean_fc", "mw_pval"))
  expect_equal(simplify_factor_names(c("a_bb_c", "a_b_c"), "_"), c("bb", "b"))
  expect_equal(simplify_factor_names(c("a_bb_c", "a_b_c_c"), "_"), c("bb", "b_c"))
  expect_equal(simplify_factor_names(c("bb", "b"), "_"), c("bb", "b"))
  expect_equal(simplify_factor_names(c("bb", "b_c"), "_"), c("bb", "b_c"))
  factor_names = c("aa", "aaa")
  simplify_factor_names(c("aa", "aaa"))
})

test_that("do_mw_test works", {
  study = get_fake_study()

  probe_names = rownames(study$get_data())
  # subset dataset
  ctrl_key = "treatment"
  ctrl_fctr = "0ug"
  case_key = "tabac"
  study$exp_grp[study$exp_grp[[ctrl_key]] == ctrl_fctr, case_key] = NA
  # print(study$exp_grp)

  mw = study$do_mw_test(probe_names, ctrl_key=ctrl_key, case_key=case_key, ctrl_fctr=ctrl_fctr)
  # mw = study$do_mw_test(probe_names, ctrl_key=ctrl_key, case_key=case_key)
  # print(mw)
  expect_equal(dim(mw), c(10,8))
  expect_gt(mw[1,1], 0)
})


test_that("do_pca works", {
  study = get_fake_study()

  probe_names = rownames(study$get_data())
  # subset dataset
  ctrl_key = "tabac"
  case_key = "tabac"
  ctrl_fctr = "Non_Smoker"
  case_fctr = "Smoker"
  # case_fctr = "Placenta"

  pca_res = study$do_pca()
  expect_silent(study$plot_pca_eig(pca_res))
  expect_silent(study$plot_pca_pc(pca_res, col=1))
})
