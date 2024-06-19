#' @title
#' Estimate Exposure Response Function
#'
#' @description
#' Estimates the exposure-response function (ERF) for a matched and weighted
#' dataset using parametric, semiparametric, and nonparametric models.
#'
#' @param .data A data frame containing an observed continuous exposure variable, weights,
#' and an observed outcome variable. Includes an `id` column for future
#' reference.
#' @param .formula A formula specifying the relationship between the exposure
#' variable and the outcome variable. For example, Y ~ w.
#' @param weights_col_name A string representing the weight or counter column
#' name in `.data`.
#' @param model_type A string representing the model type based on preliminary
#' assumptions, including `parametric`, `semiparametric`, and `nonparametric`
#' models.
#' @param w_vals A numeric vector of values at which you want to calculate the
#' ERF.
#' @param ... Additional arguments passed to the model.
#'
#' @return
#' Returns an S3 object containing the following data and parameters:
#'   - .data_original <- result_data_original
#'   - .data_prediction <- result_data_prediction
#'   - params
#'
#' @export
#'
estimate_erf <- function(.data,
                         .formula,
                         weights_col_name,
                         model_type,
                         w_vals,
                         ...) {

  # NULL and Defaults
  .family <- kernel_appr <- bw_seq <- NULL
  nthread <- 1

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

  if (!is.data.frame(.data)) stop("`.data` must be a data frame.")
  if (!weights_col_name %in% names(.data)) {
    stop(paste0(weights_col_name, " does not exist in data" ))
  }

  # Estimate exposure response -------------------------------------------------

  result <- list()
  class(result) <- "cgps_erf"

  if (any(.data[[weights_col_name]] < 0)){
    stop("Negative weights are not accepted.\n")
  }

  if (sum(.data[[weights_col_name]]) == 0) {
    .data[[weights_col_name]] <- .data[[weights_col_name]] + 1
    logger::log_debug("Giving equal weight for all samples.")
  }

  # For better visualization of weights for each data sample.
  min_weight <- min(.data[[weights_col_name]])
  max_weight <- max(.data[[weights_col_name]])
  normalized_weight <- (.data[[weights_col_name]] - min_weight)/(max_weight - min_weight)


  if (model_type == "parametric") {

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

    result_data_original <- data.frame(x = x,
                                       y_original = y_original,
                                       normalized_weight = normalized_weight)

    result_data_prediction <- data.frame(w_vals = w_vals,
                                         y_pred = y_pred)

    result$params$gnm_model <- gnm_model

  } else if (model_type == "semiparametric") {

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
    y_pred <- stats::predict(gam_model, w_pred)
    names(y_pred) <- NULL

    result_data_original <- data.frame(x = x,
                                       y_original = y_original,
                                       normalized_weight = normalized_weight)

    result_data_prediction <- data.frame(w_vals = w_vals,
                                         y_pred = y_pred)

    result$params$gam_model <- gam_model




  } else if (model_type == "nonparametric") {

    if (is.null(bw_seq)){
      stop("Please provide `bw_seq` in the additional argument.")
    }

    if (is.null(kernel_appr)){
      stop("Please provide `kernel_appar` in the additional argument.")
    }

    formula_string <- deparse(.formula)
    parts <- strsplit(formula_string, "~")[[1]]
    outcome <- trimws(parts[1])
    predictor <- trimws(parts[2])

    erf_np <- estimate_npmetric_erf(m_Y = .data[[outcome]],
                                    m_w = .data[[predictor]],
                                    counter_weight = .data[[weights_col_name]],
                                    bw_seq=bw_seq,
                                    w_vals = w_vals,
                                    nthread = nthread,
                                    kernel_appr = kernel_appr)

    formula_string <- deparse(stats::as.formula(.formula))
    parts <- strsplit(formula_string, "~")[[1]]
    predictor <- trimws(parts[2])

    x <- erf_np$params$m_w
    y_original <- erf_np$params$m_Y
    y_pred <- erf_np$erf

    result_data_original <- data.frame(x = x,
                                       y_original = y_original,
                                       normalized_weight = normalized_weight)

    result_data_prediction <- data.frame(w_vals = w_vals,
                                         y_pred = y_pred)

  } else {
    stop("The code should never get here. Double check.")
  }

  result$.data_original <- result_data_original
  result$.data_prediction <- result_data_prediction
  result$params$model_type <- model_type

  return(result)
}


fit_model <- function(data, formula, weights, model_type, ...) {
  if (model_type == "parametric") {
    model <- do.call(gnm::gnm, c(list(formula = formula, data = data, weights = weights), ...))
  } else if (model_type == "semiparametric") {
    model <- do.call(gam::gam, c(list(formula = formula, data = data, weights = weights), ...))
  } else if (model_type == "nonparametric") {
    # Assuming a function for nonparametric fitting exists
    model <- do.call(estimate_npmetric_erf, c(list(data = data, weights = weights), ...))
  } else {
    stop("Unsupported model type")
  }
  return(model)
}

