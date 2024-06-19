test_that("matching_l1 functions as expected.", {

  skip_on_cran()
  #data.table::setDTthreads(1)
  set.seed(721)
  s_data <- generate_syn_data(sample_size=200,
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
    .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
    sl_lib = c("m_xgboost"),
    gps_density = "normal")


  m_d <- data_with_gps_1
  exposure_col_name = "w"
  gps_mx <- compute_min_max(m_d$.data[["gps"]])
  w_mx <- compute_min_max(m_d$.data[[exposure_col_name]])

  val <- matching_fn(w = 10,
                     dataset = m_d$.data,
                     exposure_col_name = exposure_col_name,
                     e_gps_pred = as.vector(m_d$.data$e_gps_pred),
                     e_gps_std_pred = as.vector(m_d$.data$e_gps_std_pred),
                     w_resid = as.vector(m_d$.data$w_resid),
                     gps_mx = gps_mx,
                     w_mx = w_mx,
                     dist_measure = "l1",
                     gps_density = "kernel",
                     delta_n = 1,
                     scale = 0.5,
                     nthread = 1)



   expect_equal(nrow(val), 11)
   expect_equal(length(val), 2)
})
