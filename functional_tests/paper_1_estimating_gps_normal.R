##      author:  Naeem Khoshnevis
##      created: September 2023
##      purpose: Reproducing examples in the paper.


# Load libraries
library(ggplot2)
library(CausalGPS)
library(data.table)

# Load data --------------------------------------------------------------------
data_file <- "zip_data.RData"
if (!file.exists(data_file)) {
  stop(paste0("Download the study data file from the following link:\n",
              "https://drive.google.com/file/d/",
              "1QFdbVU8Qir1gWf96c5h_ZhT-aPjhHpqn/view?usp=share_link"))
} else {
  load(data_file)
}

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
