#' @title
#' Estimate generalized propensity score (GPS) values
#'
#' @description
#' Estimates GPS value for each observation using normal or kernel
#' approaches.
#'
#' @param .data A data frame of observed continuous exposure variable and
#' observed covariates variable. Also includes `id` column for future
#' references.
#' @param .formula A formula specifying the relationship between the exposure
#' variable and the covariates. For example, w ~ I(cf1^2) + cf2.
#' @param gps_density Model type which is used for estimating GPS value,
#' including `normal` (default) and `kernel`.
#' @param sl_lib A vector of prediction algorithms to be used by the
#' SuperLearner packageg.
#' @param ...  Additional arguments passed to the model.
#'
#' @return
#' The function returns a S3 object. Including the following:
#'   - `.data `: `id`, `exposure_var`, `gps`, `e_gps_pred`, `e_gps_std_pred`,
#'   `w_resid`
#'   - `params`: Including the following fields:
#'     - gps_mx (min and max of gps)
#'     - w_mx (min and max of w).
#'     - .formula
#'     - gps_density
#'     - sl_lib
#'     - fcall (function call)
#'
#'
#' @export
#'
#' @examples
#' \donttest{
#' m_d <- generate_syn_data(sample_size = 100)
#' data_with_gps <- estimate_gps(.data = m_d,
#'                               .formula = w ~ cf1 + cf2 + cf3 + cf4 + cf5 + cf6,
#'                               gps_density = "normal",
#'                               sl_lib = c("SL.xgboost")
#'                              )
#'}
estimate_gps <- function(.data,
                         .formula,
                         gps_density = "normal",
                         sl_lib = c("SL.xgboost"),
                         ...) {

  start_time <- proc.time()

  id_exist <- any(colnames(.data) %in% "id")
  if (!id_exist) stop(".data should include id column.")

  dot_args <- list(...)
  arg_names <- names(dot_args)

  for (i in arg_names){
    assign(i, unlist(dot_args[i], use.names = FALSE))
  }

  # function call
  fcall <- match.call()

  # Check if data has missing value(s) -----------------------------------------
  if (sum(is.na(.data)) > 0){
    logger::log_warn(
      "data data.frame has {sum(is.na(.data))} missing values.")
  }


  response_var = all.vars(.formula)[1]
  model_data <- stats::model.matrix(object = .formula,
                                    data = .data)

  response_data = .data[[response_var]]

  if (gps_density == "normal"){
    e_gps <- train_it(target = response_data,
                      input = model_data,
                      sl_lib_internal = sl_lib,
                      ...)

    e_gps_pred <- e_gps$SL.predict
    e_gps_std_pred <- stats::sd(response_data - e_gps_pred)
    w_resid <- compute_resid(response_data,
                             e_gps_pred,
                             e_gps_std_pred)
    gps <- stats::dnorm(response_data,
                        mean = e_gps_pred,
                        sd = e_gps_std_pred)

  } else if (gps_density == "kernel"){

    e_gps <- train_it(target = response_data,
                      input = model_data,
                      sl_lib_internal = sl_lib, ...)
    e_gps_pred <- e_gps$SL.predict
    e_gps_std <- train_it(target = abs(response_data - e_gps_pred),
                          input = model_data,
                          sl_lib_internal = sl_lib, ...)
    e_gps_std_pred <- e_gps_std$SL.predict
    w_resid <- compute_resid(response_data,
                             e_gps_pred,e_gps_std_pred)
    gps <- compute_density(w_resid, w_resid)

  } else {

    logger::log_error("Code should nevet get here. Doublecheck check_arguments.")
    stop(paste("Invalide gps_density: ", gps_density,
               ". Use normal or kernel."))
  }

  w_mx <- compute_min_max(response_data)
  gps_mx <- compute_min_max(gps)

  # create new data.frame for output
  dataset <- stats::setNames(data.frame(.data$id,
                                        response_data,
                                        gps),
                      c("id", response_var, "gps"))
  dataset$e_gps_pred <- e_gps_pred
  if (length(e_gps_std_pred) == 1){
    e_gps_std_pred <- rep(e_gps_std_pred, nrow(dataset))
  }
  dataset$e_gps_std_pred <- e_gps_std_pred
  dataset$w_resid <- w_resid

  # Logging for debugging purposes
  logger::log_debug("Min Max of treatment: {paste(w_mx, collapse = ', ')}")
  logger::log_debug("Min Max of gps: {paste(gps_mx, collapse = ', ')}")
  logger::log_debug("Weights for the select libraries in predicting e_gps:",
          " {paste(names(e_gps$coef), collapse = ', ')}",
          " {paste(e_gps$coef, collapse = ', ')}",
          " | Overal Risk: {sum(e_gps$coef * e_gps$cvRisk)/length(e_gps$coef)}")
  logger::log_debug("Wall clock time to estimate e_gps:",
                    " {e_gps$times$everything[3]} seconds.")
  if (gps_density == "kernel"){
    logger::log_debug("Weights for the select libraries in predicting residuals:",
            " {paste(names(e_gps_std$coef), collapse = ', ')}",
            " {paste(e_gps_std$coef, collapse = ', ')} | Overal risk:",
            " {sum(e_gps_std$coef * e_gps_std$cvRisk)/length(e_gps_std$coef)}")
    logger::log_debug("Wall clock time to estimate residuals:",
                      " {e_gps_std$times$everything[3]} seconds.")
  }

  end_time <- proc.time()

  logger::log_debug("Wall clock time to run estimate_gps function: ",
                    " {(end_time - start_time)[[3]]} seconds.")


  result <- list()
  class(result) <- "cgps_gps"
  result$.data <- dataset
  result$params$gps_mx <- gps_mx
  result$params$w_mx <- w_mx
  result$params$.formula <- .formula
  result$params$gps_density <- gps_density
  result$params$sl_lib <- sl_lib
  result$params$fcall <- fcall

  invisible(result)
}
