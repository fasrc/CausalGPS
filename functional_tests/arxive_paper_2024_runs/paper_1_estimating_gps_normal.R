##      author:  Naeem Khoshnevis
##      created: September 2023 (Updated: September 2024)
##      purpose: Reproducing examples in the paper.


# Load libraries
library(ggplot2)
library(CausalGPS)
library(data.table)

# Load data --------------------------------------------------------------------

# Study data:
# Multifactorial Zip Code-Year Dataset: Socio-Economic, Demographic, and
# Environmental Variables in the Contiguous United States (2000-2016)
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/5XBJBM


data_file <- "zip_data.RData"

load(data_file)
data.table::setDF(zip_data)
data <- zip_data

# Add id to the data
data$id <- 1:nrow(data)
data$w <- data$pm25
data$pm25 <- NULL

save(data, file = "study_data.RData")

# Estimate GPS -----------------------------------------------------------------

## Super learner wrapper
m_xgboost <- function(nthread = 6,
                      ntrees = 50,
                      shrinkage = 0.3,
                      max_depth = 6,
                      minobspernode = 1,
                      verbose = 1,
                      ...) {SuperLearner::SL.xgboost(
                        nthread = nthread,
                        ntrees = ntrees,
                        shrinkage=shrinkage,
                        max_depth=max_depth,
                        mibobspernode=minobspernode,
                        verbose=verbose,
                        ...)}

exposure <- "w"
confounders   <- c("mean_bmi", "smoke_rate",
                   "hispanic", "pct_blk", "medhouseholdincome",
                   "medianhousevalue", "poverty", "popdensity",
                   "pct_owner_occ", "summer_tmmx", "winter_tmmx",
                   "summer_rmax", "winter_rmax", "year")

formula_str <- paste(exposure, " ~ ", paste(confounders, collapse = " + "))


data_with_gps_normal <- estimate_gps(.data = data,
                                     .formula = as.formula(formula_str),
                                     gps_density = "normal",
                                     sl_lib = c("m_xgboost")
)


pdf("figure_paper_1_estimating_gps_normal.pdf")
plot(data_with_gps_normal)
dev.off()

save(data_with_gps_normal, file = "data_with_gps_normal.RData")
