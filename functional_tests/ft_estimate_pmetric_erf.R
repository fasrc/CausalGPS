set.seed(422)
n <- 2000
mydata <- generate_syn_data(sample_size=n)
year <- sample(x=c("2001","2002","2003","2004","2005"),size = n, replace = TRUE)
region <- sample(x=c("North", "South", "East", "West"),size = n, replace = TRUE)
mydata$year <- as.factor(year)
mydata$region <- as.factor(region)
mydata$cf5 <- as.factor(mydata$cf5)
set_logger(logger_level = "INFO")


id <- seq_along(1:nrow(mydata))

Y <- data.frame(id = id, out = mydata$Y)
w <- data.frame(id = id, w = mydata$w)
c <- data.frame(id = id, mydata[c("cf1","cf2","cf3","cf4","cf5","cf6","year","region")])

pseudo_pop <- generate_pseudo_pop(Y,
                                  w,
                                  c,
                                  ci_appr = "matching",
                                  sl_lib = c("m_xgboost"),
                                  params = list(xgb_nrounds=c(10,20,30),
                                                xgb_eta=c(0.1,0.2,0.3)),
                                  nthread = 1,
                                  covar_bl_method = "absolute",
                                  covar_bl_trs = 0.1,
                                  covar_bl_trs_type = "mean",
                                  max_attempt = 1,
                                  dist_measure = "l1",
                                  delta_n = 1,
                                  scale = 0.5)

outcome_m <- estimate_pmetric_erf (formula = out ~ w,
                                   family = gaussian,
                                   data = pseudo_pop$pseudo_pop)
