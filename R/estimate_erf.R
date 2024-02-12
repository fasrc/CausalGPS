#' @title
#' Estimate Exposure Response Function
#'
#' @description
#' A short description...
#'
#'
#' @param .data TBD
#' @param .formula TBD
#' @param weights_col_name TBD
#' @param model_type tBD
#' @param ... TBD
#'
#' @return
#' TBD
#' @export
#'
#' @examples
#' TBD
estimate_erf <- function(.data,
                         .formula,
                         weights_col_name,
                         model_type,
                         w_vals,
                         ...) {

  # NULL and Defaults
  .family <- NULL

  # Collect additional arguments -----------------------------------------------
  dot_args <- list(...)
  arg_names <- names(dot_args)

  for (i in arg_names){
    assign(i, unlist(dot_args[i], use.names = FALSE))
  }


  # Check input parameters -----------------------------------------------------
  model_type <- tolower(model_type)

  if (!(model_type %in% c("parametric", "semiparametric", "nonparametric"))){
    stop(paste0("The provided model type: ", model_type, ", is not supported. ",
                "Acceptable models: parametric, semiparametric, ",
                "and nonparametric."))
  }

  # Estimate exposure response -------------------------------------------------

  result <- list()
  class(result) <- "cgps_erf"


  if (model_type == "parametric") {

    if (any(.data[[weights_col_name]] < 0)){
      stop("Negative weights are not accepted.\n")
    }

    if (sum(.data[[weights_col_name]]) == 0) {
      .data[[weights_col_name]] <- .data[[weights_col_name]] + 1
      logger::log_debug("Giving equal weight for all samples.")
    }

    if (is.null(.family)){
      stop(paste0("Please provide `family` in the additional argument."))
    }

    gnm_model <- do.call(gnm::gnm, c(
                           list("formula" = .formula,
                           "family" = .family,
                           "data" = .data,
                           "weights" = .data[[weights_col_name]]),
                           ...))

    if (is.null(gnm_model)) {
      stop("gnm model is null. Did not converge.")
    }


    formula_string <- deparse(gnm_model$formula)
    parts <- strsplit(formula_string, "~")[[1]]
    outcome <- trimws(parts[1])
    predictor <- trimws(parts[2])


    x <- gnm_model$x[,2]
    names(x) <- NULL

    y_original <- gnm_model$y
    names(y_original) <- NULL


    w_pred <- data.frame(w = w_vals)
    names(w_pred) <- predictor

    y_pred <- stats::predict(gnm_model, w_pred)
    names(y_pred) <- NULL

    # normalized weight
    min_weight <- min(.data[[weights_col_name]])
    max_weight <- max(.data[[weights_col_name]])
    normalized_weight <- (.data[[weights_col_name]] - min_weight)/(max_weight - min_weight)

    result_data_original <- data.frame(x = x,
                                       y_original = y_original,
                                       normalized_weight = normalized_weight)

    result_data_prediction <- data.frame(w_vals = w_vals,
                                         y_pred = y_pred)

    result$params$gnm_model <- gnm_model

  } else if (model_type == "semiparametric") {

    if (any(.data[[weights_col_name]] < 0)){
      stop("Negative weights are not accepted.\n")
    }

    if (sum(.data[[weights_col_name]]) == 0) {
      .data[[weights_col_name]] <- .data[[weights_col_name]] + 1
      logger::log_debug("Giving equal weight for all samples.")
    }

    if (is.null(.family)){
      stop(paste0("Please provide `family` in the additional argument."))
    }

    gam_model <- do.call(gam::gam, c(
      list("formula" = .formula,
           "family" = .family,
           "data" = .data,
           "weights" = .data[[weights_col_name]]),
      ...))

    if (is.null(gam_model)) {
      stop("gnm model is null. Did not converge.")
    }

    formula_string <- deparse(gam_model$formula)
    parts <- strsplit(formula_string, "~")[[1]]
    predictor <- trimws(parts[2])

    x <- gam_model$data[[predictor]]

    y_original <- gam_model$y
    names(y_original) <- NULL

    w_pred <- data.frame(w = w_vals)
    names(w_pred) <- predictor
    y_pred <- predict(gam_model, w_pred)
    names(y_pred) <- NULL

    # normalized weight
    min_weight <- min(.data[[weights_col_name]])
    max_weight <- max(.data[[weights_col_name]])
    normalized_weight <- (.data[[weights_col_name]] - min_weight)/(max_weight - min_weight)

    result_data_original <- data.frame(x = x,
                                       y_original = y_original,
                                       normalized_weight = normalized_weight)

    result_data_prediction <- data.frame(w_vals = w_vals,
                                         y_pred = y_pred)

    result$params$gam_model <- gam_model




  } else if (model_type == "nonparametric") {

    formula_string <- deparse(.formula)
    parts <- strsplit(formula_string, "~")[[1]]
    outcome <- trimws(parts[1])
    predictor <- trimws(parts[2])

    erf_np <- estimate_npmetric_erf(m_Y = .data[[outcome]],
                                    m_w = .data[[predictor]],
                                    counter_weight = .data[[weights_col_name]],
                                    bw_seq=seq(0.2,2,0.2),
                                    w_vals = seq(2,20,0.5),
                                    nthread = 1,
                                    kernel_appr = "locpol")



    formula_string <- deparse(as.formula(.formula))
    parts <- strsplit(formula_string, "~")[[1]]
    predictor <- trimws(parts[2])

    x <- erf_np$params$m_w
    y_original <- erf_np$params$m_Y
    y_pred <- erf_np$erf


    # normalized weight
    min_weight <- min(.data[[weights_col_name]])
    max_weight <- max(.data[[weights_col_name]])
    normalized_weight <- (.data[[weights_col_name]] - min_weight)/(max_weight - min_weight)

    result_data_original <- data.frame(x = x,
                                       y_original = y_original,
                                       normalized_weight = normalized_weight)

    result_data_prediction <- data.frame(w_vals = w_vals,
                                         y_pred = y_pred)

    print(erf_np)

  } else {
    stop("The code should never get here. Double check.")
  }



  result$.data_original <- result_data_original
  result$.data_prediction <- result_data_prediction
  result$params$model_type <- model_type

  return(result)

}
