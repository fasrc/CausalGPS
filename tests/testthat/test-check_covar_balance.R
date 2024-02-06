test_that("Covariate balance check works as expected", {
skip_on_cran()
data.table::setDTthreads(1)
set.seed(532)
s_data <- generate_syn_data(sample_size = 200,
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



confounders <- paste0("cf", seq(1,6))
val1 <- check_covar_balance(w = ps_pop1$.data[, c("w")],
                            c = ps_pop1$.data[, confounders],
                            ci_appr = "matching",
                            counter_weight = ps_pop1$.data[, c("counter_weight")],
                            covar_bl_method="absolute",
                            covar_bl_trs=0.3,
                            covar_bl_trs_type="mean")







expect_true(val1$pass)
val2 <- check_covar_balance(w = ps_pop1$.data[, c("w")],
                            c = ps_pop1$.data[, confounders],
                            ci_appr = "matching",
                            counter_weight = ps_pop1$.data[, c("counter_weight")],
                            covar_bl_method="absolute",
                            covar_bl_trs=0.1,
                            covar_bl_trs_type="mean")

expect_false(val2$pass)
})
#
# test_that("Covariate balance check works as expected - part 2", {
# skip_on_cran()
# data.table::setDTthreads(1)
# set.seed(987)
# s_data <- generate_syn_data(sample_size = 500,
#                             outcome_sd = 10,
#                             gps_spec = 1,
#                             cova_spec = 1)
#
# year <- c(rep(c("2001"), each=100),
#           rep(c("2002"), each=100),
#           rep(c("2003"), each=100),
#           rep(c("2004"), each=100),
#           rep(c("2005"), each=100))
#
# region <- rep(c(rep("North",each=25),
#                 rep("South",each=25),
#                 rep("East", each=25),
#                 rep("West", each=25)), each=5)
#
# s_data$year <- as.factor(year)
# s_data$region <- as.factor(region)
#
# weight_test <- generate_pseudo_pop(
#   s_data[, c("id", "w")],
#   s_data[, c("id", "cf1", "cf2", "cf3",
#              "cf4","cf5","cf6", "year", "region")],
#   ci_appr = "weighting",
#   gps_density = "normal",
#   use_cov_transform = TRUE,
#   transformers = list("pow2", "pow3"),
#   sl_lib = c("SL.xgboost"),
#   params = list(xgb_nrounds = 50,
#                 xgb_max_depth = 6,
#                 xgb_eta = 0.3,
#                 xgb_min_child_weight = 1),
#   nthread = 1,
#   covar_bl_method = "absolute",
#   covar_bl_trs = 0.1,
#   covar_bl_trs_type = "mean",
#   exposure_trim_qtls = c(0.01, 0.99),
#   gps_trim_qtls = c(0.0, 1.0),
#   max_attempt = 1,
#   delta_n = 1,
#   scale = 0.5)
#
#   w_1 <- weight_test$pseudo_pop[, c("w")]
#   c_1 <- weight_test$pseudo_pop[, c("cf1", "cf2", "cf3",
#                                     "cf4", "cf5", "cf6",
#                                     "year", "region")]
#   cw <- weight_test$pseudo_pop[, c("counter_weight")]
#
#   val3 <- check_covar_balance(w = w_1,
#                               c = c_1,
#                               counter_weight = cw,
#                               ci_appr = "weighting",
#                               nthread = 1,
#                               covar_bl_method = "absolute",
#                               covar_bl_trs = 0.1,
#                               covar_bl_trs_type = "mean")
#
#   expect_false(val3$pass)
#
#   val4 <- check_covar_balance(w = w_1,
#                               c = c_1,
#                               counter_weight = cw,
#                               ci_appr = "weighting",
#                               nthread = 1,
#                               covar_bl_method="absolute",
#                               covar_bl_trs=0.5,
#                               covar_bl_trs_type="mean")
#   expect_true(val4$pass)
# })
