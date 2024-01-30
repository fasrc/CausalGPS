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
#' @param w A data.frame comprised of two columns: one contains the observed
#' exposure variable, and the other is labeled as 'id'. The column for the
#' outcome variable can be assigned any name as per your requirements.
#' @param c A data.frame of includes observed covariate variables. It should
#' also consist of a column named 'id'.
#' @param ci_appr The causal inference approach. Possible values are:
#'   - "matching": Matching by GPS
#'   - "weighting": Weighting by GPS
#' @param gps_density Model type which is used for estimating GPS value,
#' including `normal` (default) and `kernel`.
#' @param use_cov_transform If TRUE, the function uses transformer to meet the
#'  covariate balance.
#' @param transformers A list of transformers. Each transformer should be a
#' unary function. You can pass name of customized function in the quotes.
#' Available transformers:
#'   - pow2: to the power of 2
#'   - pow3: to the power of 3
#' @param bin_seq Sequence of w (treatment) to generate pseudo population. If
#' NULL is passed the default value will be used, which is
#' `seq(min(w)+delta_n/2,max(w), by=delta_n)`.
#' @param exposure_trim_qtls A numerical vector of two. Represents the trim quantile
#' level for exposure values. Both numbers should be in the range of \[0,1] and
#' in increasing order (default: c(0.01, 0.99)).
#' @param gps_trim_qtls A numerical vector of two. Represents the trim quantile
#' level for the gps values. Both numbers should be in the range of \[0,1] and
#' in increasing order (default: c(0.0, 1.0)).
#' @param params Includes list of params that is used internally. Unrelated
#'  parameters will be ignored.
#' @param sl_lib A vector of prediction algorithms.
#' @param nthread An integer value that represents the number of threads to be
#' used by internal packages.
#' @param include_original_data If TRUE, includes the original data in the
#' outcome.
#' @param ...  Additional arguments passed to different models.
#' @details
#' ## Additional parameters
#' ### Causal Inference Approach (ci_appr)
#' - if ci_appr = 'matching':
#'   - *dist_measure*: Matching function. Available options:
#'     - l1: Manhattan distance matching
#'   - *delta_n*: caliper parameter.
#'   - *scale*: a specified scale parameter to control the relative weight that
#'  is attributed to the distance measures of the exposure versus the GPS.
#'   - *covar_bl_method*: covariate balance method. Available options:
#'      - 'absolute'
#'   - *covar_bl_trs*: covariate balance threshold
#'   - *covar_bl_trs_type*: covariate balance type (mean, median, maximal)
#'   - *max_attempt*: maximum number of attempt to satisfy covariate balance.
#'   - See [create_matching()] for more details about the parameters and default
#'   values.
#' - if ci_appr = 'weighting':
#'   - *covar_bl_method*: Covariate balance method.
#'   - *covar_bl_trs*: Covariate balance threshold
#'   - *max_attempt*: Maximum number of attempt to satisfy covariate balance.
#'
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
#' m_d <- generate_syn_data(sample_size = 100)
#' pseuoo_pop <- generate_pseudo_pop(m_d[, c("id", "w")],
#'                                   m_d[, c("id", "cf1","cf2","cf3","cf4","cf5","cf6")],
#'                                   ci_appr = "matching",
#'                                   gps_density = "normal",
#'                                   bin_seq = NULL,
#'                                   expos_trim_qlts = c(0.01,0.99),
#'                                   gps_trim_qlts = c(0.01,0.99),
#'                                   use_cov_transform = FALSE,
#'                                   transformers = list(),
#'                                   params = list(xgb_nrounds=c(10,20,30),
#'                                                 xgb_eta=c(0.1,0.2,0.3)),
#'                                   sl_lib = c("m_xgboost"),
#'                                   nthread = 1,
#'                                   covar_bl_method = "absolute",
#'                                   covar_bl_trs = 0.1,
#'                                   covar_bl_trs_type= "mean",
#'                                   max_attempt = 1,
#'                                   dist_measure = "l1",
#'                                   delta_n = 1,
#'                                   scale = 0.5)
#'}
generate_pseudo_pop <- function(.data,
                                cw_obj,
                                covariates,
                                covar_bl_trs = 0.1,
                                covar_bl_trs_type = "maximal",
                                covar_bl_method = "absolute"){

  # Log system info
  log_system_info()

  # timing the function
  st_time_gpp <- proc.time()

  # function call
  fcall <- match.call()

  # Check arguments ----------------------------------------
  # check_args(ci_appr, use_cov_transform, transformers,
  #            gps_density, exposure_trim_qtls, ...)

  # Generate output set ------------------------------------
  counter <- 0

  ## collect additional arguments
  # dot_args <- list(...)
  # arg_names <- names(dot_args)
  #
  # for (i in arg_names){
  #   assign(i, unlist(dot_args[i], use.names = FALSE))
  # }

  # collect exposure and covariate columns
  exposure_col <- cw_obj$params$exposure_col
  covariate_cols <- covariates

  # join data based on id
  merged_data <- merge(.data, cw_obj$.data, by="id")

  # Check covariate balance for unweighted/unmatched data, but trimmed if any
  original_corr_obj <- check_covar_balance(
    w = merged_data[, c(exposure_col)],
    c = merged_data[, c(covariate_cols)],
    counter_weight = NULL,
    ci_appr = cw_obj$params$ci_appr,
    nthread = nthread,
    covar_bl_method = covar_bl_method,
    covar_bl_trs = covar_bl_trs,
    covar_bl_trs_type = covar_bl_trs_type)

  # Check covariate balance for weighted/matched data, and trimmed if any
  adjusted_corr_obj <- check_covar_balance(
    w = merged_data[, c(exposure_col)],
    c = merged_data[, c(covariate_cols)],
    counter_weight = merged_data$counter_weight,
    ci_appr = cw_obj$params$ci_appr,
    nthread = nthread,
    covar_bl_method = covar_bl_method,
    covar_bl_trs = covar_bl_trs,
    covar_bl_trs_type = covar_bl_trs_type)

  # check Kolmogorov-Smirnov statistics
  ks_stats <- check_kolmogorov_smirnov(w = merged_data[, c(exposure_col)],
                                       c = merged_data[, covariate_cols],
                                       counter_weight = merged_data[,
                                                          c("counter_weight")],
                                       ci_appr = cw_obj$params$ci_appr,
                                       nthread = nthread)


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
  class(result) <- "gpsm_pspop"

  result$params$ci_appr <- cw_obj$params$ci_appr
  #result$params$params <- params
  # for (item in arg_names){
  #   result$params[[item]] <- get(item)
  # }


  # if (include_original_data){
  #   result$original_data <- original_data
  # }

  # result$original_data_size <- nrow(original_data)

  result$pseudo_pop <- merged_data
  result$adjusted_corr_results <- adjusted_corr_obj$corr_results
  result$original_corr_results <- original_corr_obj$corr_results
  result$ks_stats <- ks_stats
  result$fcall <- fcall
  result$passed_covar_test <- adjusted_corr_obj$pass
  result$ci_appr <- cw_obj$ci_appr
  result$covariate_cols_name <- unlist(covariate_cols)
  result$ess <- ess
  result$ess_recommended <- ess_recommended

  end_time_gpp <- proc.time()

  # logger::log_debug("Wall clock time to run generate_pseudo_pop:",
  #                   " {(end_time_gpp -   st_time_gpp)[[3]]} seconds.")
  # logger::log_debug("Covariate balance condition has been met (TRUE/FALSE):",
  #                   " {adjusted_corr_obj$pass}, (iteration:",
  #                   " {counter} / {max_attempt})")
  invisible(result)
  }

