#' @title
#' Compute counter our weight of data samples
#'
#' @description
#' TBD
#'
#' @param gps_obj A gps object that is generated with `estimate_gps` function.
#' If it is provided, the number of iteration will forced to 1 (Default: NULL).
#' @param ci_appr The causal inference approach. Possible values are:
#'   - "matching": Matching by GPS
#'   - "weighting": Weighting by GPS
#' @param bin_seq Sequence of w (treatment) to generate pseudo population. If
#' NULL is passed the default value will be used, which is
#' `seq(min(w)+delta_n/2,max(w), by=delta_n)`.
#' @param exposure_trim_qtls A numerical vector of two. Represents the trim quantile
#' level for exposure values. Both numbers should be in the range of \[0,1] and
#' in increasing order (default: c(0.01, 0.99)).
#' @param gps_trim_qtls A numerical vector of two. Represents the trim quantile
#' level for the gps values. Both numbers should be in the range of \[0,1] and
#' in increasing order (default: c(0.0, 1.0)).
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
#' - if ci_appr = 'weighting':
#'   - *covar_bl_method*: Covariate balance method.
#'   - *covar_bl_trs*: Covariate balance threshold
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
compute_counter_weight <- function(gps_obj = NULL,
                                   ci_appr,
                                   bin_seq = NULL,
                                   nthread = 1,
                                   include_original_data = FALSE,
                                   ...){

  # Passing packaging check() ------------------------------
  # max_attempt <- NULL
  # covar_bl_trs <- NULL
  # covar_bl_trs_type <- NULL
  # delta_n <- NULL

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
  # counter <- 0

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
                                   bin_seq = bin_seq,
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


  return(result)
}


#' @title
#' Trim a data frame or an S3 object
#'
#' @description
#' Trims a data frame or an S3 object's `.data` attributs.
#'
#'
#' @param data_obj A data frame or an S3 object containing the data to be
#' trimmed. For a data frame, the function operates directly on it. For an S3
#' object, the function expects a `.data` attribute containing the data.
#' @param trim_quantiles A numeric vector of length 2 specifying the lower and
#' upper quantiles used for trimming the data.
#' @param variable The name of the variable in the data on which the trimming is
#'  to be applied.
#'
#' @return
#' Returns a trimmed data frame or an S3 object with the $.data attribute
#' trimmed, depending on the input type.
#'
#' @export
#'
#' @examples
#'
#' # Example usage with a data frame
#' df <- data.frame(id = 1:10, value = rnorm(100))
#' trimmed_df <- trim_it(df, c(0.1, 0.9), "value")
#'
#' # Example usage with an S3 object
#' data_obj <- list()
#' class(data_obj) <- "myobject"
#' data_obj$.data <- df
#' trimmed_data_obj <- trim_it(data_obj, c(0.1, 0.9), "value")
#'
trim_it <- function(data_obj, trim_quantiles, variable){

  type_flag <- NULL
  if (class(data_obj) == "data.frame"){
    data <- data_obj
    type_flag <- "data.frame"
  } else if (is.object(data_obj) && !isS4(data_obj)){
    data <- data_obj$.data
    type_flag <- "s3_obj"
  } else {
    stop(paste0("The data_obj with type: ", class(data_obj),
                " is not supported." ))
  }

  if ((trim_quantiles[1] < 0 || trim_quantiles[1] > 1) ||
      (trim_quantiles[2] < 0 || trim_quantiles[2] > 1) ||
      (trim_quantiles[1] > trim_quantiles[2])) {
    stop(paste("trim_quntiles should be in the [0,1] range,",
               " and the first quantile should be less than the second one."))
  }

  id_exist <- any(colnames(data) %in% "id")
  if (!id_exist) stop("data should include id column.")

  id_exist <- any(colnames(data) %in% variable)
  if (!id_exist) stop(paste0("data should include the ", variable, " column."))

  # get trim quantiles and trim data
  q1 <- stats::quantile(data[[variable]], trim_quantiles[1])
  q2 <- stats::quantile(data[[variable]], trim_quantiles[2])

  data <- data[stats::complete.cases(data), ]
  data <- data[data[[variable]] <= q2  & data[[variable]] >= q1, ]

  if (type_flag == "data.frame"){
    data_obj <- data
  } else if (type_flag == "s3_obj") {
    data_obj$.data <- data
  }

  return(data_obj)
}


