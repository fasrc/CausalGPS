#' @title
#' Estimate Parametric Exposure Response Function
#'
#' @description
#' Estimate a constant effect size for matched and weighted data set using
#' parametric models
#'
#' @param formula a vector of outcome variable in matched set.
#' @param family a description of the error distribution (see ?gnm)
#' @param data dataset that formula is build upon (Note that there should be a
#' `counter_weight` column in this data.)
#' @param ... Additional parameters for further fine tuning the gnm model.
#'
#' @details
#' This method uses generalized nonlinear model (gnm) from gnm package.
#'
#' @return
#' returns an object of class gnm
#'
#' @keywords internal
#'
estimate_pmetric_erf <- function(formula, family, data, ...) {


  ## collect additional arguments
  dot_args <- list(...)
  named_args <- stats::setNames(dot_args, names(dot_args))

  if (any(data$counter_weight < 0)){
    stop("Negative weights are not allowed.")
  }

  if (sum(data$counter_weight) == 0) {
    data$counter_weight <- data$counter_weight + 1
    logger::log_debug("Giving equal weight for all samples.")
  }

  gnm_model <- do.call(gnm::gnm, c(list("formula" = formula,
                                        "family" = family,
                                        "data" = data,
                                        "weights" = data$counter_weight),
                                   named_args))

  if (is.null(gnm_model)) {
    stop("gnm model is null. Did not converge.")
  }

  return(gnm_model)
}
