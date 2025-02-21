# (OPTIONAL) UPDATE TO MOST RECENT VERSION OF CAUSALGPS
library("devtools")
devtools::install_github("NSAPH-Software/CausalGPS", ref = "paper_examples")

library(CausalGPS)
# Confirm version installed
#> packageVersion("CausalGPS")
#[1] ‘0.5.0.9000’

# GENERTE THE DATA
m_d <- generate_syn_data(sample_size = 10000)
m_d$id <- seq_along(1:nrow(m_d))

# DEFINE THE SMOOTHER
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

# ESTIMATE THE GPS
data_with_gps_1 <- estimate_gps(
  .data = m_d,
  .formula = w ~ I(cf1^2) + cf2 + I(cf3^2) + cf4 + cf5 + cf6,
  sl_lib = c("m_xgboost"),
  gps_density = "normal")

# COMPUTE THE COUNTER WEIGHTS
cw_object_matching <- compute_counter_weight(gps_obj = data_with_gps_1,
  ci_appr = "matching",
  bin_seq = NULL,
  nthread = 1,
  delta_n = 0.1,
  dist_measure = "l1",
  scale = 1)

# (OPTIONAL) TRIM DATA TO REMOVE TAILED DISTRIBUTION WITH POTENTIAL EXTREME VALUES
trimmed_m_d <- trim_it(m_d, c(0.05, 0.95), "w")

# GENERATE PSEUDO POPULATION
# CHECK THE COVARIATE BALANCE
pseudo_pop <- generate_pseudo_pop(.data = trimmed_m_d,
  cw_obj = cw_object_matching,
  covariate_col_names = c("cf1", "cf2",
    "cf3", "cf4",
    "cf5", "cf6"),
  covar_bl_trs = 0.1,
  covar_bl_trs_type = "maximal", #or covar_bl_trs_type = "mean"
  covar_bl_method = "absolute")

# ESTIMATE THE ERF
# PARAMETRIC
erf_obj_parametric_1 <- estimate_erf(.data=pseudo_pop$.data, 
                                   .formula= as.formula("Y ~ w"), 
                                   weights_col_name="counter_weight", 
                                   model_type="parametric",
                                   w_vals = seq(2,20,0.5), 
                                   .family="gaussian")

plot(erf_obj_parametric_1)

erf_obj_parametric_2 <- estimate_erf(.data=pseudo_pop$.data, 
                                   .formula= as.formula("Y ~ w + I(w^2)"), 
                                   weights_col_name="counter_weight", 
                                   model_type="parametric",
                                   w_vals = seq(2,20,0.5), 
                                   .family="gaussian")

plot(erf_obj_parametric_2)

erf_obj_parametric_3 <- estimate_erf(.data=pseudo_pop$.data, 
                                   .formula= as.formula("Y ~ w + I(w^2) + I(w^3)"), 
                                   weights_col_name="counter_weight", 
                                   model_type="parametric",
                                   w_vals = seq(2,20,0.5), 
                                   .family="gaussian")

plot(erf_obj_parametric_3)

erf_obj_parametric_4 <- estimate_erf(.data=pseudo_pop$.data, 
                                   .formula= as.formula("Y ~ w + I(w^2) + I(w^3) + exp(w)"), 
                                   weights_col_name="counter_weight", 
                                   model_type="parametric",
                                   w_vals = seq(2,20,0.5), 
                                   .family="gaussian")

plot(erf_obj_parametric_4)

erf_obj_parametric_5 <- estimate_erf(.data=pseudo_pop$.data, 
                                   .formula= as.formula("Y ~ w + I(w^2) + I(w^3) + exp(w) + log(w)"), 
                                   weights_col_name="counter_weight", 
                                   model_type="parametric",
                                   w_vals = seq(2,20,0.5), 
                                   .family="gaussian")

plot(erf_obj_parametric_5)

# SEMI-PARAMETRIC
erf_obj_semiparametric_1 <- estimate_erf(.data=pseudo_pop$.data, 
                                       .formula= as.formula("Y ~ w"), 
                                       weights_col_name="counter_weight", 
                                       model_type="semiparametric",
                                       w_vals = seq(2,20,0.5), 
                                       .family="gaussian")

plot(erf_obj_semiparametric_1)

erf_obj_semiparametric_2 <- estimate_erf(.data=pseudo_pop$.data, 
                                       .formula= as.formula("Y ~ s(w, 3)"), 
                                       weights_col_name="counter_weight", 
                                       model_type="semiparametric",
                                       w_vals = seq(2,20,0.5), 
                                       .family="gaussian")

plot(erf_obj_semiparametric_2)


erf_obj_semiparametric_3 <- estimate_erf(.data=pseudo_pop$.data, 
                                   .formula= as.formula("Y ~ s(w, 100)"), 
                                   weights_col_name="counter_weight", 
                                   model_type="semiparametric",
                                   w_vals = seq(2,20,0.5), 
                                   .family="gaussian")

plot(erf_obj_semiparametric_3)

# NON-PARAMETRIC
erf_obj_nonparametric <- estimate_erf(.data=pseudo_pop$.data, 
                                      .formula= Y ~ w, 
                                      weights_col_name="counter_weight", 
                                      model_type = "nonparametric",
                                      w_vals = seq(2,20,0.5),
                                      bw_seq = 20,
                                      kernel_appr = "kernsmooth")
plot(erf_obj_nonparametric)
