set.seed(712)
m_d <- generate_syn_data(sample_size = 1000)

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
                                             scale = 0.5)




plot(cw_object_matching, subset_ids = c(100, 300))


cw_object_weighting <- compute_counter_weight(gps_obj = data_with_gps,
                                              ci_appr = "weighting",
                                              bin_seq = NULL,
                                              nthread = 1,
                                              delta_n = 0.1,
                                              dist_measure = "l1",
                                              scale = 0.5)


plot(cw_object_weighting)


pseudo_pop_weighting <- generate_pseudo_pop(.data = m_d,
                                            cw_obj = cw_object_weighting,
                                            covariate_col_names = c("cf1", "cf2", "cf3",
                                                                    "cf4", "cf5", "cf6"),
                                            covar_bl_trs = 0.1,
                                            covar_bl_trs_type = "maximal",
                                            covar_bl_method = "absolute")

plot(pseudo_pop_weighting, include_details = TRUE)

pseudo_pop_matching <- generate_pseudo_pop(.data = m_d,
                                           cw_obj = cw_object_matching,
                                           covariate_col_names = c("cf1", "cf2", "cf3",
                                                                   "cf4", "cf5", "cf6"),
                                           covar_bl_trs = 0.1,
                                           covar_bl_trs_type = "maximal",
                                           covar_bl_method = "absolute")

plot(pseudo_pop_matching, include_details = TRUE)


# Exposure response function


# parametric
erf_obj_parametric <- estimate_erf(.data = pseudo_pop_matching$.data,
                        .formula = Y ~ w,
                        weights_col_name = "counter_weight",
                        model_type = "parametric",
                        w_vals = seq(2,20,0.5),
                        .family = "gaussian")



plot(erf_obj_parametric)


# semiparametric
erf_obj_semiparametric <- estimate_erf(.data = pseudo_pop_weighting$.data,
                                       .formula = Y ~ w,
                                       weights_col_name = "counter_weight",
                                       model_type = "semiparametric",
                                       w_vals = seq(2,20,0.5),
                                       .family = "gaussian")



plot(erf_obj_semiparametric)


# non-parametric
erf_obj_nonparametric <- estimate_erf(.data = pseudo_pop_weighting$.data,
                                       .formula = Y ~ w,
                                       weights_col_name = "counter_weight",
                                       model_type = "nonparametric",
                                       w_vals = seq(2,20,0.5),
                                       bw_seq = seq(0.2,2,0.2),
                                       kernel_appr = "kernsmooth")



plot(erf_obj_nonparametric)
