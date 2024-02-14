test_that("check_kolmogorov_smirnov works as expected.", {
  skip_on_cran()
  data.table::setDTthreads(1)
  set.seed(8422)
  n <- 500
  s_data <- generate_syn_data(sample_size = n,
                              outcome_sd = 10,
                              gps_spec = 1,
                              cova_spec = 1)

  s_data$id <- seq_along(1:nrow(s_data))
  m_xgboost <- function(nthread = 4,
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

  data_with_gps_1 <- estimate_gps(
    .data = s_data,
    .formula = w ~ I(cf1^2) + cf2 + I(cf3^2) + cf4 + cf5 + cf6,
    sl_lib = c("m_xgboost"),
    gps_density = "normal")

  cw_object_matching <- compute_counter_weight(gps_obj = data_with_gps_1,
                                               ci_appr = "matching",
                                               bin_seq = NULL,
                                               nthread = 1,
                                               delta_n = 0.1,
                                               dist_measure = "l1",
                                               scale = 0.5)

  ps_pop1 <- generate_pseudo_pop(.data = s_data,
                                 cw_obj = cw_object_matching,
                                 covariate_col_names = c("cf1", "cf2",
                                                         "cf3", "cf4",
                                                         "cf5", "cf6"),
                                 covar_bl_trs = 0.1,
                                 covar_bl_trs_type = "maximal",
                                 covar_bl_method = "absolute")

  output <- CausalGPS:::check_kolmogorov_smirnov(w = ps_pop1$.data[, c("w")],
                                                 c = ps_pop1$.data[, ps_pop1$params$covariate_col_names],
                                                 counter = ps_pop1$.data[, c("counter_weight")],
                                                 ci_appr="matching")

  expect_equal(length(output), 2L)
  expect_equal(length(output$ks_stat), 7L)
  expect_equal(length(output$stat_vals), 3L)
  expect_equal(output$ks_stat[["w"]], 0.041004785, tolerance = 0.000001)
  expect_equal(output$ks_stat[["cf1"]], 0.008191388, tolerance = 0.000001)
  expect_equal(output$stat_vals[["maximal_val"]], 0.04100478, tolerance = 0.000001)
})
