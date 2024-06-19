set.seed(967)
m_d <- generate_syn_data(sample_size = 500)

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


# estimate_erf as a wrapper function.

print(paste0("Size of original data: ", nrow(m_d)))

trimmed_data <- trim_it(data_obj = m_d,
                        trim_quantiles = c(0.05, 0.95),
                        variable = "w")

print(paste0("Size of data after trimming for exposure: ", nrow(trimmed_data)))

data_with_gps_1 <- estimate_gps(.data = trimmed_data,
                                .formula = w ~ I(cf1^2) + cf2 + I(cf3^2) + cf4 + cf5 + cf6,
                                sl_lib = c("m_xgboost"),
                                gps_density = "normal")

print(paste0("Size of GPS data based on exposure trime: ",
             nrow(data_with_gps_1$.data)))

summary(data_with_gps_1)
plot(data_with_gps_1)

trimmed_gps_data <- trim_it(data_with_gps_1, c(0.1, 0.9), variable = "gps")


summary(trimmed_gps_data)
plot(trimmed_gps_data)

# Estimating GPS with kernel approach

data_with_gps_kernel <- estimate_gps(.data = trimmed_data,
                                .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
                                sl_lib = c("m_xgboost"),
                                gps_density = "kernel")


plot(data_with_gps_kernel)
