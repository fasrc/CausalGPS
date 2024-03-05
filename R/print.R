#' @title
#' Extend print function for cgps_erf object
#'
#' @param x A cgps_erf object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' No return value. This function is called for side effects.
#'
#' @export
#'
print.cgps_erf <- function(x, ...) {

  x <- unclass(x)

  cat(" CausalGPS exposure rate function object\n")
  cat(" function call: \n")
  cat("      ***       \n")
  print(x$fcall, ...)
  cat("      ***       \n")
  cat(" Output data can be accessed at $erf \n")
  cat(" Look at summary for more details.\n")
}



#' @title
#' print summary of cgps_erf object

#'
#' @param object A cgps_erf object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' Returns summary of data
#' @export
summary.cgps_erf <- function(object, ...) {

  # cat_list <- function(input) {
  #   cat(paste("   size: ", length(input),
  #             ", class: ", class(input),
  #             ", missing value(s): ", sum(is.na(input)),
  #             sep = ""))
  #   if (is.numeric(input)) {
  #     cat(paste("\n   min: ", sprintf("%.3f", min(input, na.rm = TRUE)),
  #               "\n   max: ", sprintf("%.3f", max(input, na.rm = TRUE)),
  #               "\n   mean: ", sprintf("%.3f", mean(input, na.rm = TRUE)),
  #               sep = ""))
  #   }
  # }
  #
  # object <- unclass(object)
  # cat("Input data: \n")
  # for (item in names(object$params)){
  #   cat(paste(" ", item, "\n"))
  #   cat_list(object$params[[item]])
  #   cat("\n")
  # }
  # cat("\nOutput data: \n")
  # cat(paste("  erf\n"))
  # cat_list(object$erf)
  # cat("\n")
}


#' @title
#' Extend print function for cgps_pspop object
#'
#' @param x A cgps_pspop object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' No return value. This function is called for side effects.
#'
#' @export
#'
print.cgps_pspop <- function(x, ...) {

  x <- unclass(x)

  cat(" CausalGPS pseudo population object\n")
  cat(" function call: \n")
  cat("      ***       \n")
  print(x$fcall, ...)
  cat("      ***       \n")
  cat(" Output data can be accessed at $pseudo_pop \n")
  cat(" Look at summary for more details. \n")
}

#' @title
#' print summary of cgps_pspop object

#'
#' @param object A cgps_pspop object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' Returns summary of data
#' @export
summary.cgps_pspop <- function(object, ...) {

  cat("--- CausalGPS pseudo population object summary --- \n")
  cat(paste("Pseudo population met the covariate balance requirement: ",
            object$passed_covar_test, "\n"))
  cat(paste("Absolute correlation of the original data: \n",
            "  mean:    ", sprintf("%.3f",
                                   object$original_corr_results$mean_absolute_corr),
            "\n",
            "  median:  ", sprintf("%.3f",
                                   object$original_corr_results$median_absolute_corr),
            "\n",
            "  maximal: ", sprintf("%.3f",
                                   object$original_corr_results$maximal_absolute_corr),
            "\n"
  ))
  cat(paste("\n", names(object$original_corr_results$absolute_corr), ":",
            sprintf("%.3f",object$original_corr_results$absolute_corr)))
  cat(paste("\n\n Absolute correlation of the pseudo population: \n",
            "  mean:    ", sprintf("%.3f",
                                   object$adjusted_corr_results$mean_absolute_corr),
            "\n",
            "  median:  ", sprintf("%.3f",
                                   object$adjusted_corr_results$median_absolute_corr),
            "\n",
            "  maximal: ", sprintf("%.3f",
                                   object$adjusted_corr_results$maximal_absolute_corr),
            "\n"
  ))
  cat(paste("\n", names(object$adjusted_corr_results$absolute_corr), ":",
            sprintf("%.3f",object$adjusted_corr_results$absolute_corr)))
  cat(paste("\n\n Hyperparameters used for the select population:"))
  cat(paste("\n", names(object$best_gps_used_params), ":",
            object$best_gps_used_params))
  cat("\n\n")
  cat(paste("Number of data samples: ", nrow(object$pseudo_pop), "\n"))
  cat(paste("Number of iterations: ", object$counter, "\n"))
  cat("Effective sample size: \n")
  cat(paste("  Achieved: ", object$ess, "\n"))
  cat(paste("  Min recommended: ", object$ess_recommended, "\n"))
  cat("Kolmogorov-Smirnov (KS) statistics:")
  if (is.null(object$ks_stats)){
    cat("\n  Not computed. \n")
  } else {
    cat(paste("\n", " ", names(object$ks_stats$ks_stat), ":",
              sprintf("%.3f", object$ks_stats$ks_stat)))
    cat(paste("\n summary: \n",
              "  mean:    ",
              sprintf("%.3f", object$ks_stats$stat_vals[["mean_val"]]), "\n",
              "  median:  ",
              sprintf("%.3f", object$ks_stats$stat_vals[["median_val"]]), "\n",
              "  maximal: ",
              sprintf("%.3f", object$ks_stats$stat_vals[["maximal_val"]]), "\n"
    ))
  }
  cat("--- *** --- \n")
}



#' @title
#' Extend print function for cgps_gps object
#'
#' @param x A cgps_gps object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' No return value. This function is called for side effects.
#'
#' @export
#'
print.cgps_gps <- function(x, ...) {

  x <- unclass(x)

  cat(" CausalGPS gps object\n")
  cat(" function call: \n")
  cat("      ***       \n")
  print(x$params$fcall, ...)
  cat("      ***       \n")
  cat(" Output data can be accessed at $.data \n")
  cat(" Look at summary for more details. \n")
}


#' @title
#' print summary of cgps_gps object

#'
#' @param object A cgps_gps object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' Returns summary of data
#' @export
summary.cgps_gps <- function(object, ...) {

  cat("--- CausalGPS gps object summary --- \n")
  cat("      ***       \n")
  data_size <- nrow(object$.data)
  cat(paste0("Number of data samples: ", data_size, "\n"))
  cat(paste0("Used formula: ", deparse(object$params$.formula), "\n"))
  cat("--- *** --- \n")

}



#' @title
#' Extend print function for cgps_cw object
#'
#' @param x A cgps_cw object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' No return value. This function is called for side effects.
#'
#' @export
#'
print.cgps_cw <- function(x, ...) {

  x <- unclass(x)

  cat(" TBD\n")
}


#' @title
#' print summary of cgps_cw object

#'
#' @param object A cgps_cw object.
#' @param ... Additional arguments passed to customize the results.
#'
#' @return
#' Returns summary of data
#' @export
summary.cgps_cw <- function(object, ...) {

  object <- unclass(object)
  cat("TBD \n")
}




