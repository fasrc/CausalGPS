#' @title
#' Generate pseudo population
#'
#' @description
#' Generates pseudo population data set based on user-defined causal inference
#' approach. The function uses an adaptive approach to satisfies covariate
#' balance requirements. The function terminates either by satisfying covariate
#' balance or completing the requested number of iteration, whichever comes
#' first.
#'
#' @param .data A data.frame of observation data with `id` column.
#' @param cw_obj An S3 object of counter_weight.
#' @param covariate_col_names A list of covariate columns.
#' @param covar_bl_trs Covariate balance threshold
#' @param covar_bl_trs_type Type of the covariance balance threshold.
#' @param covar_bl_method Covariate balance method.
#'
#' @return
#' Returns a pseudo population (gpsm_pspop) object that is generated
#' or augmented based on the selected causal inference approach (ci_appr). The
#' object includes the following objects:
#' - params
#'   - ci_appr
#'   - params
#' - pseudo_pop
#' - adjusted_corr_results
#' - original_corr_results
#' - best_gps_used_params
#' - effect size of generated pseudo population
#'
#' @export
#' @examples
#' \donttest{
#'
#' set.seed(967)
#'
#' m_d <- generate_syn_data(sample_size = 200)
#' m_d$id <- seq_along(1:nrow(m_d))
#'
#' m_xgboost <- function(nthread = 4,
#'                       ntrees = 35,
#'                       shrinkage = 0.3,
#'                       max_depth = 5,
#'                       ...) {SuperLearner::SL.xgboost(
#'                         nthread = nthread,
#'                         ntrees = ntrees,
#'                         shrinkage=shrinkage,
#'                         max_depth=max_depth,
#'                         ...)}
#'
#' data_with_gps_1 <- estimate_gps(
#'   .data = m_d,
#'   .formula = w ~ I(cf1^2) + cf2 + I(cf3^2) + cf4 + cf5 + cf6,
#'   sl_lib = c("m_xgboost"),
#'   gps_density = "normal")
#'
#' cw_object_matching <- compute_counter_weight(gps_obj = data_with_gps_1,
#'                                              ci_appr = "matching",
#'                                              bin_seq = NULL,
#'                                              nthread = 1,
#'                                              delta_n = 0.1,
#'                                              dist_measure = "l1",
#'                                              scale = 0.5)
#'
#' pseudo_pop <- generate_pseudo_pop(.data = m_d,
#'                                   cw_obj = cw_object_matching,
#'                                   covariate_col_names = c("cf1", "cf2",
#'                                                           "cf3", "cf4",
#'                                                           "cf5", "cf6"),
#'                                   covar_bl_trs = 0.1,
#'                                   covar_bl_trs_type = "maximal",
#'                                   covar_bl_method = "absolute")
#'
#'}
generate_pseudo_pop <- function(.data,
                                cw_obj,
                                covariate_col_names,
                                covar_bl_trs = 0.1,
                                covar_bl_trs_type = "maximal",
                                covar_bl_method = "absolute"){

  # Log system info
  log_system_info()

  # timing the function
  st_time_gpp <- proc.time()

  # function call
  fcall <- match.call()

  # Generate output set ------------------------------------
  counter <- 0

  # collect exposure and covariate columns
  exposure_col <- cw_obj$params$exposure_col
  covariate_cols <- covariate_col_names

  # join data based on id
  merged_data <- merge(.data, cw_obj$.data, by="id")

  # Check covariate balance for unweighted/unmatched data, but trimmed if any
  original_corr_obj <- check_covar_balance(
    w = merged_data[, c(exposure_col)],
    c = merged_data[, c(covariate_cols)],
    counter_weight = NULL,
    ci_appr = cw_obj$params$ci_appr,
    covar_bl_method = covar_bl_method,
    covar_bl_trs = covar_bl_trs,
    covar_bl_trs_type = covar_bl_trs_type)

  covar_bl_t <- paste0(covar_bl_trs_type, "_absolute_corr")
  message(paste0(covar_bl_trs_type, " absolute correlation (original): ",
                getElement(original_corr_obj$corr_results, covar_bl_t),
                "| Covariate balance threshold: ", covar_bl_trs))

  # Check covariate balance for weighted/matched data, and trimmed if any
  adjusted_corr_obj <- check_covar_balance(
    w = merged_data[, c(exposure_col)],
    c = merged_data[, c(covariate_cols)],
    counter_weight = merged_data$counter_weight,
    ci_appr = cw_obj$params$ci_appr,
    covar_bl_method = covar_bl_method,
    covar_bl_trs = covar_bl_trs,
    covar_bl_trs_type = covar_bl_trs_type)

  message(paste0(covar_bl_trs_type, " absolute correlation (adjusted): ",
                 getElement(adjusted_corr_obj$corr_results, covar_bl_t),
                 "| Covariate balance threshold: ", covar_bl_trs))

  # check Kolmogorov-Smirnov statistics
  ks_stats <- check_kolmogorov_smirnov(w = merged_data[, c(exposure_col)],
                                       c = merged_data[, covariate_cols],
                                       counter_weight = merged_data[,
                                                          c("counter_weight")],
                                       ci_appr = cw_obj$params$ci_appr)


  # compute effective sample size
  ess_recommended <- length(merged_data[, c(exposure_col)]) / 10
  ess <- ((sum(merged_data$counter_weight) ^ 2) /
            sum(merged_data$counter_weight ^ 2))
  if (ess < ess_recommended) {
    logger::log_warn("Effective sample size is less than recommended.",
                     "Current: {ess}, recommended min value:",
                     " {ess_recommended}.")
  }

  result <- list()
  class(result) <- "cgps_pspop"

  result$params$ci_appr <- cw_obj$params$ci_appr

  result$.data <- merged_data
  result$params$adjusted_corr_results <- adjusted_corr_obj$corr_results
  result$params$original_corr_results <- original_corr_obj$corr_results
  result$params$ks_stats <- ks_stats
  result$params$fcall <- fcall
  result$params$passed_covar_test <- adjusted_corr_obj$pass
  result$params$ci_appr <- cw_obj$params$ci_appr
  result$params$covariate_col_names <- unlist(covariate_cols)
  result$params$ess <- ess
  result$params$ess_recommended <- ess_recommended
  result$params$covar_bl_trs <- covar_bl_trs
  result$params$covar_bl_trs_type <- covar_bl_trs_type
  result$params$covar_bl_method <- covar_bl_method
  result$params$cw_obj_params <- cw_obj$params

  end_time_gpp <- proc.time()

  invisible(result)
  }

