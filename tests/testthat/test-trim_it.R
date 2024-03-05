test_that("trim_it works as expected", {

  set.seed(967)
  m_d <- generate_syn_data(sample_size = 1000)

  trimmed_data <- trim_it(data_obj = m_d,
                          trim_quantiles = c(0.05, 0.95),
                          variable = "w")

  expect_equal(length(trimmed_data), 9L)
  expect_equal(nrow(trimmed_data), 900L)

  m_xgboost <- function(nthread = 1,
                        ntrees = 100,
                        shrinkage = 0.3,
                        max_depth = 6,
                        minobspernode = 1,
                        verbose = 0,
                        ...) {SuperLearner::SL.xgboost(
                          nthread = nthread,
                          ntrees = ntrees,
                          shrinkage=shrinkage,
                          max_depth=max_depth,
                          mibobspernode=minobspernode,
                          verbose=verbose,
                          ...)}


  assign("m_xgboost", m_xgboost, envir = .GlobalEnv)

  data_with_gps_1 <- estimate_gps(
    .data = m_d,
    .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
    sl_lib = c("m_xgboost"))

  trimmed_gps_obj <- trim_it(data_with_gps_1,
                             trim_quantiles = c(0.15, 0.90),
                             variable = "gps")

  expect_equal(nrow(trimmed_gps_obj$.data), 750L)
})
