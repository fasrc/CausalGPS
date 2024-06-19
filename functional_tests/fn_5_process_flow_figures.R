set.seed(513)
m_d <- generate_syn_data(sample_size = 5000)
m_d <- trim_it(m_d, c(0.01,0.99), variable = "w")

m_xgboost <- function(nthread = 6,
                      ntrees = 30,
                      shrinkage = 0.19,
                      max_depth = 5,
                      ...) {SuperLearner::SL.xgboost(
                        nthread = nthread,
                        ntrees = ntrees,
                        shrinkage=shrinkage,
                        max_depth=max_depth,
                        ...)}

gps_obj <- estimate_gps(.data = m_d,
                        .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                        sl_lib = c("m_xgboost"),
                        gps_density = "normal")


pdf("process_flow_fig_1_gps.pdf", width = 8, height = 6)
plot(gps_obj)
dev.off()

cw_object_matching <- compute_counter_weight(gps_obj = gps_obj,
                                             ci_appr = "matching",
                                             bin_seq = NULL,
                                             nthread = 6,
                                             delta_n = 0.1,
                                             dist_measure = "l1",
                                             scale = 0.5)

pdf("process_flow_fig_2_cw.pdf", width = 8, height = 6)
plot(cw_object_matching, subset_ids = c(200, 300), every_n = 10)
dev.off()

pseudo_pop_matching <- generate_pseudo_pop(.data = m_d,
                                           cw_obj = cw_object_matching,
                                           covariate_col_names = c("cf1", "cf2", "cf3",
                                                                   "cf4", "cf5", "cf6"),
                                           covar_bl_trs = 0.1,
                                           covar_bl_trs_type = "maximal",
                                           covar_bl_method = "absolute")

pdf("process_flow_fig_3_pspop.pdf", width = 8, height = 6)
plot(pseudo_pop_matching)
dev.off()

# Exposure response function

# non-parametric
erf_obj_nonparametric <- estimate_erf(.data = pseudo_pop_matching$.data,
                                       .formula = Y ~ w,
                                       weights_col_name = "counter_weight",
                                       model_type = "nonparametric",
                                       w_vals = seq(2,20,0.5),
                                       bw_seq = seq(0.2,2,0.2),
                                       kernel_appr = "kernsmooth",
                                       nthread = 6)


pdf("process_flow_fig_4_erf_np.pdf", width = 8, height = 6)
plot(erf_obj_nonparametric)
dev.off()
