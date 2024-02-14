#' @title
#' Estimate semi-exposure-response function (semi-ERF).
#'
#' @description
#' Estimates the smoothed exposure-response function using a generalized
#' additive model with splines.
#'
#' @param formula a vector of outcome variable in matched set.
#' @param family a description of the error distribution (see ?gam).
#' @param data dataset that formula is build upon Note that there should be a
#' `counter_weight` column in this data.).
#' @param ... Additional parameters for further fine tuning the gam model.
#'
#' @details
#' This approach uses Generalized Additive Model (gam) using mgcv package.
#'
#' @return
#' returns an object of class gam
#'
#' @keywords internal
#'
estimate_semipmetric_erf <- function(formula, family, data, ...) {


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

  gam_model <- do.call(gam::gam, c(list("formula" = formula,
                                        "family" = family,
                                        "data" = data,
                                        "weights" = data$counter_weight),
                                   named_args))

  if (is.null(gam_model)) {
    stop("gam model is null. Did not converge.")
  }

  return(gam_model)
}
