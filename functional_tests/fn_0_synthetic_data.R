library(CausalGPS)


mydata <- generate_syn_data(sample_size = 1000,
                            outcome_sd = 10,
                            gps_spec = 1,
                            cova_spec = 1,
                            vectorized_y = FALSE)


head(mydata)


# More investigation about synthetic data
