#' @title
#' Compile pseudo population
#'
#' @description
#' Compiles pseudo population based on the original population and estimated GPS
#' value.
#'
#' @param data_obj A S3 object including the following:
#'   - Original data set + GPS values
#'   - e_gps_pred
#'   - e_gps_std_pred
#'   - w_resid
#'   - gps_mx (min and max of gps)
#'   - w_mx (min and max of w).
#' @param ci_appr Causal inference approach.
#' @param gps_density Model type which is used for estimating GPS value,
#' including `normal` and `kernel`.
#' @param exposure_col_name Exposure data column name.
#' @param nthread An integer value that represents the number of threads to be
#' used by internal packages.
#' @param ... Additional parameters.
#'
#' @details
#' For matching approach, use an extra parameter, `bin_seq`, which is sequence
#' of w (treatment) to generate pseudo population. If `NULL` is passed the
#' default value will be used, which is
#' `seq(min(w)+delta_n/2,max(w), by=delta_n)`.
#'
#'
#' @export
#'
#' @return
#' `compile_pseudo_pop` returns the pseudo population data that is compiled based
#' on the selected causal inference approach.
#'
#' @examples
#' \donttest{
#' set.seed(112)
#' m_d <- generate_syn_data(sample_size = 100)
#'
#' m_xgboost <- function(nthread = 1,
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
#' data_with_gps <- estimate_gps(.data = m_d,
#'                               .formula = w ~ cf1 + cf2 + cf3 +
#'                                              cf4 + cf5 + cf6,
#'                               gps_density = "normal",
#'                               sl_lib = c("m_xgboost")
#'                              )
#'
#'
#' pd <- compile_pseudo_pop(data_obj = data_with_gps,
#'                          ci_appr = "matching",
#'                          gps_density = "normal",
#'                          bin_seq = NULL,
#'                          exposure_col_name = c("w"),
#'                          nthread = 1,
#'                          dist_measure = "l1",
#'                          covar_bl_method = 'absolute',
#'                          covar_bl_trs = 0.1,
#'                          covar_bl_trs_type= "mean",
#'                          delta_n = 0.5,
#'                          scale = 1)
#'}
compile_pseudo_pop <- function(data_obj,
                               ci_appr,
                               gps_density,
                               exposure_col_name,
                               nthread,
                               ...) {

  dist_measure <- delta_n <- bin_seq <- NULL


  if (!(is.object(data_obj) && !isS4(data_obj))) {
    stop("data_obj should be a S3 object.")
  }

  if (!(is.element(".data", attributes(data_obj)$names))) {
    stop("data_obj should have the required .data field.")
  }

  logger::log_info("Starting compiling pseudo population ",
                    " (original data size: {nrow(data_obj$.data)}) ... ")

  ## collect additional arguments
  dot_args <- list(...)
  arg_names <- names(dot_args)

  for (i in arg_names){
    assign(i, unlist(dot_args[i], use.names = FALSE))
  }


  auxilary_columns <- c("e_gps_pred", "e_gps_std_pred", "w_resid")
  if (ci_appr == 'matching'){
      matched_set <- create_matching(.data = data_obj$.data,
                                     exposure_col_name = exposure_col_name,
                                     matching_fn = matching_fn,
                                     dist_measure = dist_measure,
                                     gps_density = gps_density,
                                     delta_n = delta_n,
                                     scale = scale,
                                     bin_seq = bin_seq,
                                     nthread = nthread)
      logger::log_info("Finished compiling pseudo population ",
                      " (Pseudo population data size: {nrow(matched_set)})")
      matched_set[, (auxilary_columns) := NULL]
      return(matched_set)

  } else if (ci_appr == 'weighting'){
    weighted_set <- create_weighting(data_obj$.data,
                                     exposure_col_name)
    logger::log_info("Finished compiling pseudo population ",
                     " (Pseudo population data size: {nrow(weighted_set)})")
    weighted_set[, (auxilary_columns) := NULL]
    return(weighted_set)

  } else {

  stop(paste('The code should not get here.',
             'Something is wrong with checking arguments.'))
  }
}
