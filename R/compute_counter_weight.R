#' @title
#' Compute counter or weight of data samples
#'
#' @description
#' Computes counter (for matching approach) or weight (for weighting) approach.
#'
#' @param gps_obj A gps object that is generated with `estimate_gps` function.
#' If it is provided, the number of iteration will forced to 1 (Default: NULL).
#' @param ci_appr The causal inference approach. Possible values are:
#'   - "matching": Matching by GPS
#'   - "weighting": Weighting by GPS
#' @param nthread An integer value that represents the number of threads to be
#' used by internal packages.
#' @param ...  Additional arguments passed to different models.
#' @details
#' ## Additional parameters
#' ### Causal Inference Approach (ci_appr)
#' - if ci_appr = 'matching':
#'   - *bin_seq*: A sequence of w (treatment) to generate pseudo population.
#'   If `NULL` is passed the default value will be used, which is
#'   `seq(min(w)+delta_n/2,max(w), by=delta_n)`.
#'   - *dist_measure*: Matching function. Available options:
#'     - l1: Manhattan distance matching
#'   - *delta_n*: caliper parameter.
#'   - *scale*: a specified scale parameter to control the relative weight that
#'  is attributed to the distance measures of the exposure versus the GPS.
#'
#' @return
#' Returns a counter_weight (cgps_cw) object that includes `.data` and `params`
#' attributes.
#' - `.data`: includes `id` and `counter_weight` columns. In case of `matching`
#' the `counter_weight` column is integer values, which represent how many times
#' the provided observational data was mached during the matching process. In
#' case of `weighting` the column is double values.
#'
#' - `params`: Include related parameters that is used for the process.
#'
#' @export
#' @examples
#' \donttest{
#' m_d <- generate_syn_data(sample_size = 100)
#' gps_obj <- estimate_gps(.data = m_d,
#'                         .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
#'                         gps_density = "normal",
#'                         sl_lib = c("SL.xgboost"))
#'
#' cw_object <- compute_counter_weight(gps_obj = gps_obj,
#'                                     ci_appr = "matching",
#'                                     bin_seq = NULL,
#'                                     nthread = 1,
#'                                     delta_n = 0.1,
#'                                     dist_measure = "l1",
#'                                     scale = 0.5)
#'}
compute_counter_weight <- function(gps_obj,
                                   ci_appr,
                                   nthread = 1,
                                   ...){

  # Passing packaging check() ------------------------------
  delta_n <- NULL
  dist_measure <- NULL
  scale <- NULL
  bin_seq <- NULL

  # Log system info
  log_system_info()

  # timing the function
  st_time_gpp <- proc.time()

  # function call
  fcall <- match.call()

  ## collect additional arguments
  dot_args <- list(...)
  arg_names <- names(dot_args)

  for (i in arg_names){
    assign(i, unlist(dot_args[i], use.names = FALSE))
  }

  if (!is.null(gps_obj)){
    if (!inherits(gps_obj, "cgps_gps")){
      stop("Provided gps_obj is not an standard gps object.")
    }
  }

  exposure_col <- all.vars(gps_obj$params$.formula)[1]
  .data <- gps_obj$.data

  zero_initialize <- rep(0, nrow(.data))
  gps_obj$.data$counter_weight <- zero_initialize


  if (ci_appr == "matching"){
    if (is.null(delta_n)){
      stop("delta_n input param is not provided for matching approach.")
    }
    if (is.null(scale)){
      stop("scale input param is not provided for matching approach.")
    }
    if (is.null(dist_measure)){
      stop("dist_measure input param is not provided for matching approach.")
    }
    if (!is.numeric(delta_n)){
      stop("delta_n should be a numeric value.")
    }
    if (!is.numeric(scale) || scale < 0 || scale > 1){
      stop("scale should be a numerical value in [0, 1] range.")
    }
  }

  if (is.null(bin_seq) && ci_appr == "matching"){
    min_w <- min(.data[[exposure_col]])
    max_w <- max(.data[[exposure_col]])
    start_val <- min_w + delta_n/2
    end_val <- max_w
    if ((start_val < end_val && delta_n < 0) ||
        (start_val > end_val && delta_n > 0)) {
      stop(paste("Inconsistent values for sequencing.",
                 " start val: ", start_val,
                 " end val: ", end_val,
                 " delta_n/2: ", delta_n / 2,
                 "\n delta_n should be less than: ", (max_w - min_w) / 2 ))
    }
  }

  counter_weighted_data <- compile_pseudo_pop(
                                   data_obj = gps_obj,
                                   ci_appr = ci_appr,
                                   gps_density = gps_obj$params$gps_density,
                                   exposure_col_name = exposure_col,
                                   nthread = nthread,
                                   ...)

  counter_weighted_data <- data.frame(
    id = counter_weighted_data$id,
    counter_weight = counter_weighted_data$counter_weight)

  result <- list()
  class(result) <- "cgps_cw"
  result$.data <- counter_weighted_data
  result$params$ci_appr <- ci_appr
  result$params$exposure_col <- exposure_col

  if (ci_appr == "matching"){
    result$params$delta_n <- delta_n
    result$params$scale <- scale
    result$params$dist_measure <- dist_measure
  }

  return(result)
}





