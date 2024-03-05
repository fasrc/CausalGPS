test_that("estimate_npmetric_erf works as expected", {
  # skip_on_cran()
  # set.seed(347)
  # n <- 1000
  #
  # mydata <- generate_syn_data(sample_size=n)
  # year <- sample(x=c("2001","2002","2003","2004","2005"),
  #                size = n, replace = TRUE)
  # region <- sample(x=c("North", "South", "East", "West"),
  #                  size = n, replace = TRUE)
  #
  # mydata$year <- as.factor(year)
  # mydata$region <- as.factor(region)
  # mydata$cf5 <- as.factor(mydata$cf5)
  #
  # # mydata$id <- seq_along(1:nrow(mydata))
  #
  # m_xgboost <- function(nthread = 4,
  #                       ntrees = 35,
  #                       shrinkage = 0.3,
  #                       max_depth = 5,
  #                       ...) {SuperLearner::SL.xgboost(
  #                         nthread = nthread,
  #                         ntrees = ntrees,
  #                         shrinkage=shrinkage,
  #                         max_depth=max_depth,
  #                         ...)}
  #
  # assign("m_xgboost", m_xgboost, envir = .GlobalEnv)
  #
  # data_with_gps_1 <- estimate_gps(
  #   .data = mydata,
  #   .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6 + year + region,
  #   sl_lib = c("m_xgboost"),
  #   gps_density = "normal")
  #
  # cw_object_matching <- compute_counter_weight(gps_obj = data_with_gps_1,
  #                                              ci_appr = "matching",
  #                                              bin_seq = NULL,
  #                                              nthread = 1,
  #                                              delta_n = 0.1,
  #                                              dist_measure = "l1",
  #                                              scale = 0.5)
  #
  # pseudo_pop <- generate_pseudo_pop(.data = mydata,
  #                                cw_obj = cw_object_matching,
  #                                covariate_col_names = c("cf1", "cf2",
  #                                                        "cf3", "cf4",
  #                                                        "cf5", "cf6",
  #                                                        "year", "region"),
  #                                covar_bl_trs = 0.1,
  #                                covar_bl_trs_type = "maximal",
  #                                covar_bl_method = "absolute")
  #
  # # pseudo_pop <- generate_pseudo_pop(m_d[, c("id", "w")],
  # #                                   m_d[, c("id", "cf1", "cf2", "cf3",
  # #                                           "cf4", "cf5", "cf6")],
  # #                                   ci_appr = "matching",
  # #                                   gps_density = "kernel",
  # #                                   exposure_trim_qtls = c(0.01,0.99),
  # #                                   sl_lib = c("SL.xgboost"),
  # #                                   covar_bl_method = "absolute",
  # #                                   covar_bl_trs = 0.1,
  # #                                   covar_bl_trs_type = "mean",
  # #                                   max_attempt = 1,
  # #                                   dist_measure = "l1",
  # #                                   delta_n = 1,
  # #                                   scale = 1)
  #
  # min_w <- min(pseudo_pop$pseudo_pop$w)
  # max_w <- max(pseudo_pop$pseudo_pop$w)
  #
  # data <- pseudo_pop$pseudo_pop
  #
  # res <- estimate_npmetric_erf(data$Y,
  #                              data$w,
  #                              data$counter_weight,
  #                              bw_seq=seq(0.2,0.2,0.2),
  #                              w_vals=seq(min_w,max_w,0.5),
  #                              nthread = 1)
  #
  # expect_equal(class(res),"gpsm_erf")
  # expect_equal(length(res$params$bw_seq), 1)
  # expect_equal(length(res$params$w_vals), length(res$erf))
  # #expect_equal(res$risk_val[1], 1305125, tolerance = 0.00001)
})


test_that("estimate_npmetric_erf works as expected (with earth)", {

  # skip_if_not_installed("earth")
  # skip_on_cran()
  #
  # set.seed(347)
  # m_d <- generate_syn_data(sample_size = 400)
  # m_d$id <- seq_along(1:nrow(m_d))
  #
  # pseudo_pop <- generate_pseudo_pop(m_d[, c("id", "w")],
  #                                   m_d[, c("id", "cf1", "cf2", "cf3", "cf4",
  #                                           "cf5", "cf6")],
  #                                   ci_appr = "matching",
  #                                   pred_model = "sl",
  #                                   gps_density = "kernel",
  #                                   exposure_trim_qtls = c(0.01,0.99),
  #                                   sl_lib = c("SL.xgboost",
  #                                              "SL.earth",
  #                                              "SL.gam"),
  #                                   covar_bl_method = "absolute",
  #                                   covar_bl_trs = 0.1,
  #                                   covar_bl_trs_type = "mean",
  #                                   max_attempt = 1,
  #                                   dist_measure = "l1",
  #                                   delta_n = 1,
  #                                   scale = 0.5)
  #
  # min_w <- min(pseudo_pop$pseudo_pop$w)
  # max_w <- max(pseudo_pop$pseudo_pop$w)
  #
  # data <- merge(m_d[, c("id", "Y")], pseudo_pop$pseudo_pop, by="id")
  #
  # res <- estimate_npmetric_erf(data$Y,
  #                              data$w,
  #                              data$counter_weight,
  #                              bw_seq=seq(0.2,0.4,0.2),
  #                              w_vals=seq(min_w,max_w,0.5),
  #                              nthread = 1)
  #
  # expect_equal(class(res),"gpsm_erf")
  # expect_equal(length(res$params$bw_seq), 2)
  # expect_equal(length(res$params$w_vals), length(res$erf))
  # #expect_equal(res$risk_val[1], 1305125, tolerance = 0.00001)

})

