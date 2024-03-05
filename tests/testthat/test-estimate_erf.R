test_that("estimate erf works as expected", {

  set.seed(7312)
  m_d <- generate_syn_data(sample_size = 400)

  m_xgboost <- function(nthread = 1,
                        ntrees = 35,
                        shrinkage = 0.3,
                        max_depth = 5,
                        ...) {SuperLearner::SL.xgboost(
                          nthread = nthread,
                          ntrees = ntrees,
                          shrinkage=shrinkage,
                          max_depth=max_depth,
                          ...)}

  assign("m_xgboost", m_xgboost, envir = .GlobalEnv)

  data_with_gps <- estimate_gps(.data = m_d,
                                .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                                sl_lib = c("m_xgboost"),
                                gps_density = "kernel")


  cw_object_matching <- compute_counter_weight(gps_obj = data_with_gps,
                                               ci_appr = "matching",
                                               bin_seq = NULL,
                                               nthread = 1,
                                               delta_n = 0.1,
                                               dist_measure = "l1",
                                               scale = 1)


  cw_object_weighting <- compute_counter_weight(gps_obj = data_with_gps,
                                                ci_appr = "weighting",
                                                bin_seq = NULL,
                                                nthread = 1,
                                                delta_n = 0.1,
                                                dist_measure = "l1",
                                                scale = 0.5)

  pseudo_pop_weighting <- generate_pseudo_pop(.data = m_d,
                                              cw_obj = cw_object_weighting,
                                              covariate_col_names = c("cf1", "cf2", "cf3",
                                                                      "cf4", "cf5", "cf6"),
                                              covar_bl_trs = 0.1,
                                              covar_bl_trs_type = "maximal",
                                              covar_bl_method = "absolute")


  pseudo_pop_matching <- generate_pseudo_pop(.data = m_d,
                                             cw_obj = cw_object_matching,
                                             covariate_col_names = c("cf1", "cf2", "cf3",
                                                                     "cf4", "cf5", "cf6"),
                                             covar_bl_trs = 0.1,
                                             covar_bl_trs_type = "maximal",
                                             covar_bl_method = "absolute")



  # parametric -----------------------------------------------------------------
  erf_obj_parametric_matching <- estimate_erf(
    .data = pseudo_pop_matching$.data,
    .formula = Y ~ w,
    weights_col_name = "counter_weight",
    model_type = "parametric",
    w_vals = seq(2,20,0.5),
    .family = "gaussian")


  expect_true(".data_original" %in% names(erf_obj_parametric_matching))
  expect_true(".data_prediction" %in% names(erf_obj_parametric_matching))
  expect_true("params" %in% names(erf_obj_parametric_matching))
  expect_equal(erf_obj_parametric_matching$params$model_type, "parametric")
  expect_equal(nrow(erf_obj_parametric_matching$.data_original), 400L)
  expect_equal(nrow(erf_obj_parametric_matching$.data_prediction), 37L)
  expect_equal(length(erf_obj_parametric_matching$.data_original), 3L)
  expect_equal(length(erf_obj_parametric_matching$.data_prediction), 2L)


  erf_obj_parametric_weighting <- estimate_erf(
    .data = pseudo_pop_weighting$.data,
    .formula = Y ~ w,
    weights_col_name = "counter_weight",
    model_type = "parametric",
    w_vals = seq(2,20,0.5),
    .family = "gaussian")

  expect_true(".data_original" %in% names(erf_obj_parametric_weighting))
  expect_true(".data_prediction" %in% names(erf_obj_parametric_weighting))
  expect_true("params" %in% names(erf_obj_parametric_weighting))
  expect_equal(erf_obj_parametric_weighting$params$model_type, "parametric")
  expect_equal(nrow(erf_obj_parametric_weighting$.data_original), 400L)
  expect_equal(nrow(erf_obj_parametric_weighting$.data_prediction), 37L)
  expect_equal(length(erf_obj_parametric_weighting$.data_original), 3L)
  expect_equal(length(erf_obj_parametric_weighting$.data_prediction), 2L)


  # semi-parametric -----------------------------------------------------------------
  erf_obj_semiparametric_matching <- estimate_erf(
    .data = pseudo_pop_matching$.data,
    .formula = Y ~ w,
    weights_col_name = "counter_weight",
    model_type = "semiparametric",
    w_vals = seq(2,20,0.5),
    .family = "gaussian")


  expect_true(".data_original" %in% names(erf_obj_semiparametric_matching))
  expect_true(".data_prediction" %in% names(erf_obj_semiparametric_matching))
  expect_true("params" %in% names(erf_obj_semiparametric_matching))
  expect_equal(erf_obj_semiparametric_matching$params$model_type, "semiparametric")
  expect_equal(nrow(erf_obj_semiparametric_matching$.data_original), 400L)
  expect_equal(nrow(erf_obj_semiparametric_matching$.data_prediction), 37L)
  expect_equal(length(erf_obj_semiparametric_matching$.data_original), 3L)
  expect_equal(length(erf_obj_semiparametric_matching$.data_prediction), 2L)


  erf_obj_semiparametric_weighting <- estimate_erf(
    .data = pseudo_pop_weighting$.data,
    .formula = Y ~ w,
    weights_col_name = "counter_weight",
    model_type = "semiparametric",
    w_vals = seq(2,20,0.5),
    .family = "gaussian")

  expect_true(".data_original" %in% names(erf_obj_semiparametric_weighting))
  expect_true(".data_prediction" %in% names(erf_obj_semiparametric_weighting))
  expect_true("params" %in% names(erf_obj_semiparametric_weighting))
  expect_equal(erf_obj_semiparametric_weighting$params$model_type, "semiparametric")
  expect_equal(nrow(erf_obj_semiparametric_weighting$.data_original), 400L)
  expect_equal(nrow(erf_obj_semiparametric_weighting$.data_prediction), 37L)
  expect_equal(length(erf_obj_semiparametric_weighting$.data_original), 3L)
  expect_equal(length(erf_obj_semiparametric_weighting$.data_prediction), 2L)

  # non-parametric -----------------------------------------------------------------
  erf_obj_nonparametric <- estimate_erf(.data = pseudo_pop_weighting$.data,
                                        .formula = Y ~ w,
                                        weights_col_name = "counter_weight",
                                        model_type = "nonparametric",
                                        w_vals = seq(2,20,0.5),
                                        bw_seq = seq(0.2,2,0.2),
                                        kernel_appr = "kernsmooth")


  expect_true(".data_original" %in% names(erf_obj_nonparametric))
  expect_true(".data_prediction" %in% names(erf_obj_nonparametric))
  expect_true("params" %in% names(erf_obj_nonparametric))
  expect_equal(erf_obj_nonparametric$params$model_type, "nonparametric")
  expect_equal(nrow(erf_obj_nonparametric$.data_original), 400L)
  expect_equal(nrow(erf_obj_nonparametric$.data_prediction), 37L)
  expect_equal(length(erf_obj_nonparametric$.data_original), 3L)
  expect_equal(length(erf_obj_nonparametric$.data_prediction), 2L)

})
