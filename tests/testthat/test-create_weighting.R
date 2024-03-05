test_that("create_weighting works as expected.", {

  skip_on_cran()
  set.seed(481)
  data.table::setDTthreads(1)
  mydata <- generate_syn_data(sample_size = 100)

  mydata$id <- seq_along(1:nrow(mydata))

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
    .data = mydata,
    .formula = w ~ I(cf1^2) + cf2 + I(cf3^2) + cf4 + cf5 + cf6,
    sl_lib = c("m_xgboost"),
    gps_density = "normal")



  result <- create_weighting(dataset = data_with_gps_1$.data,
                             exposure_col_name = "w")


  expect_equal(length(result), 7)
  expect_equal(nrow(result), 100)
})
