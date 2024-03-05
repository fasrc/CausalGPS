set.seed(967)
m_d <- generate_syn_data(sample_size = 10000)

trimmed_data <- trim_it(data_obj = m_d,
                        trim_quantiles = c(0.05, 0.95),
                        variable = "w")

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

data_with_gps <- estimate_gps(.data = trimmed_data,
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




plot(cw_object_matching)


cw_object_weighting <- compute_counter_weight(gps_obj = data_with_gps,
                                              ci_appr = "weighting",
                                              bin_seq = NULL,
                                              nthread = 1,
                                              delta_n = 0.1,
                                              dist_measure = "l1",
                                              scale = 0.5)




plot(cw_object_weighting)


pseudo_pop_weighting <- generate_pseudo_pop(.data = trimmed_data,
                                  cw_obj = cw_object_weighting,
                                  covariate_col_names = c("cf1", "cf2", "cf3",
                                                          "cf4", "cf5", "cf6"),
                                  covar_bl_trs = 0.1,
                                  covar_bl_trs_type = "maximal",
                                  covar_bl_method = "absolute")

plot(pseudo_pop_weighting, include_details = TRUE)

pseudo_pop_matching <- generate_pseudo_pop(.data = trimmed_data,
                                            cw_obj = cw_object_matching,
                                            covariate_col_names = c("cf1", "cf2", "cf3",
                                                                    "cf4", "cf5", "cf6"),
                                            covar_bl_trs = 0.1,
                                            covar_bl_trs_type = "maximal",
                                            covar_bl_method = "absolute")

plot(pseudo_pop_matching, include_details = TRUE)
