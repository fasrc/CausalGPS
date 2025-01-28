##      author:  Naeem Khoshnevis
##      created: March 2024 (Updated: September 2024)
##      purpose: Reproducing examples in the paper.


# Load libraries
library(ggplot2)
library(CausalGPS)


# Load gps object
load("data_with_gps_normal.RData")

cw_weighting_object <- compute_counter_weight(gps_obj = data_with_gps_normal,
                                              ci_appr = "weighting",
                                              bin_seq = NULL,
                                              nthread = 6,
                                              delta_n = 0.1,
                                              dist_measure = "l1",
                                              scale = 1)


save(cw_weighting_object, file = "cw_weighting_object.RData")


