##      author:  Naeem Khoshnevis
##      created: March 2024
##      purpose: Reproducing examples in the paper.


# Load libraries
library(ggplot2)
library(CausalGPS)


# Load cw object and data
load("cw_matching_object.RData")
load("study_data.RData")

confounders   <- c("mean_bmi", "smoke_rate",
                   "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "popdensity",
                   "pct_owner_occ", "summer_tmmx", "winter_tmmx",
                   "summer_rmax", "winter_rmax", "year")


pseudo_pop_weighting_object <- generate_pseudo_pop(
                                           .data = data,
                                           cw_obj = cw_matching_object,
                                           covariate_col_names = confounders,
                                           covar_bl_trs = 0.1,
                                           covar_bl_trs_type = "maximal",
                                           covar_bl_method = "absolute")


# save(pseudo_pop_matching_object, file = "pseudo_pop_matching_object.RData")
#
# pdf("figure_paper_6_pseudo_pop_matching_object.pdf")
# plot(pseudo_pop_matching_object)
# dev.off()
