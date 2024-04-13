##      author:  Naeem Khoshnevis
##      created: March 2024
##      purpose: Reproducing examples in the paper.


# Load libraries
library(ggplot2)
library(CausalGPS)


# Load gps object
load("data_with_gps_normal.RData")


cw_matching_object <- compute_counter_weight(gps_obj = data_with_gps_normal,
                                             ci_appr = "matching",
                                             bin_seq = NULL,
                                             nthread = 6,
                                             delta_n = 0.1,
                                             dist_measure = "l1",
                                             scale = 0.5)


save(cw_matching_object, file = "cw_matching_object.RData")
