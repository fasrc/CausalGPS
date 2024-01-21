test_that("estimate_gps works as expected.", {

  #skip_on_cran()
  set.seed(895)
  data.table::setDTthreads(1)
  m_d <- generate_syn_data(sample_size = 100)
  m_d$id <- seq_along(1:nrow(m_d))

  data_with_gps_1 <- estimate_gps(
                       data = m_d,
                       formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                       sl_lib = c("m_xgboost"))

  expect_equal(length(data_with_gps_1$dataset), 6)
  expect_equal(nrow(data_with_gps_1$dataset), 100)
  expect_equal(data_with_gps_1$dataset$gps[2], 20.991916, tolerance = 0.00001)
  expect_true(any(colnames(data_with_gps_1$dataset) %in% "id"))
  expect_true("formula" %in% names(data_with_gps_1))

  data_with_gps_2 <- estimate_gps(
                       data = m_d,
                       formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                       sl_lib = c("m_xgboost")
  )

  required_elements <- c( "dataset",
                          "gps_mx",
                          "w_mx" )
  expect_equal(length(intersect(c(attributes(data_with_gps_2)$names),
                                required_elements)), 3L)
  expect_equal(data_with_gps_2$dataset$e_gps_pred[58,], 19.07269287,
               tolerance = 0.00001)

  # Missing values
  set.seed(1789)
  m_d_2 <- generate_syn_data(sample_size = 100)
  m_d_2$w[20] <- NA
  m_d_2$id <- seq_along(1:nrow(m_d_2))
  # Missing value in target
  # Error because SL does not support missing data.
  expect_error(estimate_gps(
                 data = m_d_2,
                 formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                 sl_lib = c("m_xgboost")))

})


test_that("estimate_gps input data should include id column.", {
  skip_on_cran()
  set.seed(895)
  m_d <- generate_syn_data(sample_size = 1000)
  m_d$id <- NULL


  expect_error(estimate_gps(data = m_d,
                            formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                            sl_lib = c("m_xgboost")),
               regexp = "data should include id column.")
})


test_that("estimate_gps residuals works as expected.", {
  skip_on_cran()
  set.seed(895)
  m_d <- generate_syn_data(sample_size = 1000)
  data_with_gps <- estimate_gps(data = m_d,
                                formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                                sl_lib = c("m_xgboost"),
                                params = list(xgb_max_depth = c(3,4,5),
                                              xgb_rounds = c(10,20,30,40)),
                                gps_density = "normal")
  std_pred <- data_with_gps$dataset$e_gps_std_pred
  expect_equal(std_pred[2], std_pred[30], tolerance = 0.0000001)
  expect_equal(std_pred[200], std_pred[650], tolerance = 0.0000001)


  set.seed(274)
  m_d <- generate_syn_data(sample_size = 1000)
  data_with_gps <- estimate_gps(data = m_d,
                                formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                                sl_lib = c("m_xgboost"),
                                params = list(xgb_max_depth = c(3,4,5),
                                              xgb_rounds = c(10,20,30,40)),
                                gps_density = "kernel")
  std_pred <- data_with_gps$dataset$e_gps_std_pred
  expect_true(std_pred[2] != std_pred[30])
  expect_true(std_pred[200] != std_pred[650])
})

