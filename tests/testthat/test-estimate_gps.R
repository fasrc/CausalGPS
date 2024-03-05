test_that("estimate_gps works as expected.", {

  #skip_on_cran()
  set.seed(895)
  m_d <- generate_syn_data(sample_size = 100)
  m_d$id <- seq_along(1:nrow(m_d))

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

  expect_equal(length(data_with_gps_1$.data), 6L)
  expect_equal(nrow(data_with_gps_1$.data), 100)
  expect_true(any(colnames(data_with_gps_1$.data) %in% "id"))
  expect_true(".formula" %in% names(data_with_gps_1$params))

  data_with_gps_2 <- estimate_gps(
                       .data = m_d,
                       .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                       sl_lib = c("m_xgboost")
  )

  required_elements <- c( ".data",
                          "params")
  expect_equal(length(intersect(c(attributes(data_with_gps_2)$names),
                                required_elements)), 2L)

  # Missing values
  set.seed(1789)
  m_d_2 <- generate_syn_data(sample_size = 100)
  m_d_2$w[20] <- NA
  m_d_2$id <- seq_along(1:nrow(m_d_2))
  # Missing value in target
  # Error because SL does not support missing data.
  expect_error(estimate_gps(
                 .data = m_d_2,
                 .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                 sl_lib = c("m_xgboost")))

})


test_that("estimate_gps input data should include id column.", {
  skip_on_cran()
  set.seed(895)
  m_d <- generate_syn_data(sample_size = 1000)
  m_d$id <- NULL

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
  expect_error(estimate_gps(.data = m_d,
                            .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                            sl_lib = c("m_xgboost")),
               regexp = "data should include id column.")
})


test_that("estimate_gps residuals works as expected.", {
  skip_on_cran()
  set.seed(895)

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
  m_d <- generate_syn_data(sample_size = 1000)
  data_with_gps <- estimate_gps(.data = m_d,
                                .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                                sl_lib = c("m_xgboost"),
                                gps_density = "normal")
  std_pred <- data_with_gps$.data$e_gps_std_pred
  expect_equal(std_pred[2], std_pred[30], tolerance = 0.0000001)
  expect_equal(std_pred[200], std_pred[650], tolerance = 0.0000001)


  set.seed(274)
  m_d <- generate_syn_data(sample_size = 1000)
  data_with_gps <- estimate_gps(.data = m_d,
                                .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                                sl_lib = c("m_xgboost"),
                                gps_density = "kernel")
  std_pred <- data_with_gps$.data$e_gps_std_pred
  expect_true(std_pred[2] != std_pred[30])
  expect_true(std_pred[200] != std_pred[650])
})

