##      author:  Naeem Khoshnevis
##      created: September 2024
##      purpose: Reproducing examples in the paper.


# Load libraries
library(ggplot2)
library(CausalGPS)


# Load cw object and data
load("pseudo_pop_matching_object.RData")


erf_obj_parametric_matching <- estimate_erf(
  .data = pseudo_pop_matching_object$.data,
  .formula = education ~ w,
  weights_col_name = "counter_weight",
  w_vals = seq(2,20,0.5),
  model_type = "parametric",
  .family = "gaussian")

pdf("figure_paper_8_erf_obj_parametric_matching.pdf")
plot(erf_obj_parametric_matching)
dev.off()


erf_obj_semiparametric_matching <- estimate_erf(
  .data = pseudo_pop_matching_object$.data,
  .formula = education ~ w,
  weights_col_name = "counter_weight",
  w_vals = seq(2,20,0.5),
  model_type = "semiparametric",
  .family = "gaussian")

pdf("figure_paper_8_erf_obj_semiparametric_matching.pdf")
plot(erf_obj_semiparametric_matching)
dev.off()
