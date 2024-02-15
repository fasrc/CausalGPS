#' @title
#' A helper function for cgps_pspop object
#'
#' @description
#' A helper function to plot cgps_pspop object using ggplot2 package.
#'
#' @param object A cgps_pspop object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @return
#' Returns a ggplot object.
#'
#'
#' @keywords internal
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#'
autoplot.cgps_pspop <- function(object, ...){

  gg_labs <- gg_title <- include_details <- NULL

  if (object$params$ci_appr == "matching"){
    default_gg_labs <- list(x = "Absolute Correlation", y = "Covariates")
  } else if (object$params$ci_appr == "weighting"){
    default_gg_labs <- list(x = "Absolute Weighted Correlation", y="Covariates")
  }

  default_gg_title <- "Covariate Balance Test"

  ## collect additional arguments
  dot_args <- list(...)
  arg_names <- names(dot_args)

  for (i in arg_names){
    assign(i,unlist(dot_args[i], use.names = FALSE))
  }

  # create a data.frame from two data.
  balance <- data.frame(original = object$params$original_corr_results$absolute_corr,
                        adjusted = object$params$adjusted_corr_results$absolute_corr)

  # Convert row names in to attribute.
  balance$covar_label <- row.names(balance)

  # sort data.frame based on original data correlation values
  balance <- balance[order(balance$original), ]
  covar_label <- balance$covar_label
  row.names(balance) <- NULL

  n_cov <- length(balance$original)
  m_balance <- melt(as.data.table(balance), measure.vars = c("original",
                                                             "adjusted"))
  m_balance$covariates <- rep(seq(1, n_cov, 1), 2)

  g <- ggplot2::ggplot(data = m_balance,
                       ggplot2::aes(x=.data$value,
                                    y=.data$covariates,
                                    color=.data$variable,
                                    shape=.data$variable)) + # Add shape aesthetic
    ggplot2::geom_point() +
    ggplot2::geom_path() +
    ggplot2::scale_color_manual(values = c("original" = "blue", "adjusted" = "orange")) + # Define colors manually
    ggplot2::scale_shape_manual(values=c(16, 15)) +  # Differentiate shapes: 16 is a solid circle, 24 is a triangle
    ggplot2::scale_y_discrete(limit = factor(1:n_cov),labels = covar_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_vline(xintercept = object$params$covar_bl_trs) +
    ggplot2::geom_vline(xintercept = object$params$adjusted_corr_results$mean_absolute_corr,
                        linetype="dotdash") +
    ggplot2::labs(x = default_gg_labs$x, y = default_gg_labs$y) +
    ggplot2::ggtitle(default_gg_title)

  if (!is.null(gg_labs)) {
    g <- g + ggplot2::labs(x=gg_labs[1], y=gg_labs[2])
  }

  if (!is.null(gg_title)) {
    g <- g + ggplot2::ggtitle(gg_title)
  }

  if (!is.null(include_details)){
    if (!is.logical(include_details)){
      stop(paste("include_details should be logical value.",
                 "Current: ", include_details))
      }
  } else {
    include_details <- FALSE
  }


  # legend position
  if (include_details) {
    g <- g + ggplot2::theme(legend.position = c(0.8, 0.1))
  }

  # Object details
  if (include_details){
    object_sum <- "Object Summary "
    object_sum <- paste0(object_sum, "\n", "--------------")
    object_sum <- paste0(object_sum, "\n", "Approach: ",
                         object$params$ci_appr)
    if (object$params$ci_appr == "matching"){
      object_sum <- paste0(object_sum, "\n", " scale: ",
                           round(object$params$cw_obj_params$scale, 3))
      object_sum <- paste0(object_sum, "\n", " delta: ",
                           round(object$params$cw_obj_params$delta, 3))
      object_sum <- paste0(object_sum, "\n", " dist measure: ",
                           object$params$cw_obj_params$dist_measure)
    }
    object_sum <- paste0(object_sum, "\n", "--------------")
    object_sum <- paste0(object_sum, "\n", "Covar balance method: ",
                         object$params$covar_bl_method)
    object_sum <- paste0(object_sum, "\n", " Threshhold type: ",
                         object$params$covar_bl_trs_type)
    object_sum <- paste0(object_sum, "\n", " Threshhold value: ",
                         round(object$params$covar_bl_trs, 3))
    object_sum <- paste0(object_sum, "\n", " Passed covariate test: ",
                         object$params$passed_covar_test)
    object_sum <- paste0(object_sum, "\n", " Maximal abs. cov.: ",
                         round(object$params$adjusted_corr_results$maximal_absolute_corr, 3))
    object_sum <- paste0(object_sum, "\n", " Median abs. cov.: ",
                         round(object$params$adjusted_corr_results$median_absolute_corr, 3))
    object_sum <- paste0(object_sum, "\n", " Mean abs. cov.: ",
                         round(object$params$adjusted_corr_results$mean_absolute_corr, 3))
    object_sum <- paste0(object_sum, "\n", "--------------")
    object_sum <- paste0(object_sum, "\n", "Kolmogorov Smirnov Test")
    object_sum <- paste0(object_sum, "\n", " ess: ", round(object$params$ess,3))
    object_sum <- paste0(object_sum, "\n", " ess (min recomennded):", round(object$params$ess_recommended,3))


    text_grob <- cowplot::ggdraw() + cowplot::draw_label(object_sum,
                                                         fontface = 'plain',
                                                         hjust = 0,
                                                         vjust = 1,
                                                         x = 0.01,
                                                         y = 0.95,
                                                         size = 7,
                                                         lineheight = 1.2)
  }

  if (include_details){
    g <- cowplot::plot_grid(
      g, text_grob,
      ncol = 2,
      rel_widths = c(2, 1)
    )
  }

  return(g)
}


#' @title
#' Extend generic plot functions for cgps_pspop class
#'
#' @description
#' A wrapper function to extend generic plot functions for cgps_pspop class.
#'
#' @param x  A cgps_pspop object.
#' @param ... Additional arguments passed to customize the plot.
#' @details
#' ## Additional parameters
#' - *include_details*: If set to TRUE, the plot will include run details (Default = FALSE).
#'
#' @return
#' Returns a ggplot2 object, invisibly. This function is called for side effects.
#'
#' @export
#'
plot.cgps_pspop <- function(x, ...) {
  g <- ggplot2::autoplot(x, ...)
  print(g)
  invisible(g)
}


#' @title
#' A helper function for cgps_gps object
#'
#' @description
#' A helper function to plot cgps_gps object using ggplot2 package.
#'
#' @param object A cgps_gps object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @return
#' Returns a ggplot object.
#'
#'
#' @keywords internal
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#'
autoplot.cgps_gps <- function(object, ...){

  ## collect additional arguments
  dot_args <- list(...)
  arg_names <- names(dot_args)

  for (i in arg_names) {
    assign(i, unlist(dot_args[i], use.names = FALSE))
  }


  # create a density plot
  dataset <- object$.data

  g <- ggplot2::ggplot(data = dataset,
                       ggplot2::aes(x = .data$gps)) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_density(fill = "blue", alpha = 0.5) +
    ggplot2::ggtitle("Density Plot of GPS Values")
  g <- g + ggplot2::labs(x = "GPS value", y = "Density")

  return(g)

}


#' @title
#' Extend generic plot functions for cgps_gps class
#'
#' @description
#' A wrapper function to extend generic plot functions for cgps_gps class.
#'
#' @param x  A cgps_gps object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @return
#' Returns a ggplot2 object, invisibly. This function is called for side effects.
#'
#' @export
#'
plot.cgps_gps <- function(x, ...) {
  g <- ggplot2::autoplot(x, ...)
  print(g)
  invisible(g)
}


#' @title
#' A helper function for cgps_cw object
#'
#' @description
#' A helper function to plot cgps_cw object using ggplot2 package.
#'
#' @param object A cgps_cw object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @return
#' Returns a ggplot object.
#'
#'
#' @keywords internal
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#'
autoplot.cgps_cw <- function(object, ...){

  id <- counter_weight <- NULL

  dataset <- object$.data

  # Default values
  every_n <- 10
  subset_ids <- NULL


  ## collect additional arguments
  dot_args <- list(...)
  arg_names <- names(dot_args)

  for (i in arg_names) {
    assign(i, unlist(dot_args[i], use.names = FALSE))
  }


  min_id <- min(dataset$id)
  max_id <- max(dataset$id)

  if (is.null(subset_ids)){
    subset_data <- dataset
    used_min_id <- min_id
    used_max_id <- max_id
  } else {
    if (!all(is.integer(as.integer(subset_ids))) || (length(subset_ids) != 2)) {
      stop("subset_ids should be a numerical vector of length 2.")
    }

    if ((subset_ids[1] > subset_ids[2])
        || (subset_ids[1] > max_id)
        || (subset_ids[2] < min_id)){
      stop(paste0("Provided subset_ids: [", subset_ids[1], ",", subset_ids[2],
                  "] is out of valid range of ",
                  " [", min_id, "," , max_id, "]"))
    }

    if (subset_ids[1] < min_id){
      warning(paste0("The minimum id value in the data is: ", min_id))
      used_min_id <- min_id
    } else {
      used_min_id <- subset_ids[1]
    }
    if (subset_ids[2] > max_id){
      warning(paste0("The maximum id value in the data is: ", max_id))
      used_max_id <- max_id
    } else {
      used_max_id <- subset_ids[2]
    }
    subset_data <- dataset[dataset$id >= used_min_id & dataset$id <= used_max_id,]
  }

  if (object$params$ci_appr == "matching") {
    hist_density_plot <- ggplot2::ggplot(subset_data,
                                         ggplot2::aes(x = counter_weight)) +
      ggplot2::geom_histogram(bandwidth = 1,
                              fill = "blue",
                              color = "black") +
      ggplot2::labs(x = "Count", y = "Frequency") +
      ggplot2::theme_minimal()
    y_label <- "Count"
  } else if (object$params$ci_appr == "weighting") {
    hist_density_plot <- ggplot2::ggplot(subset_data,
                                         ggplot2::aes(x = counter_weight)) +
      ggplot2::geom_density(fill = "blue", color = "black") +
      ggplot2::labs(x = "Weight  (displayed in log10 scale)",
                    y = "Density") +
      ggplot2::theme_minimal() + ggplot2::scale_x_log10()
    y_label <- "Weight"
  } else {
    stop("The cgsp_cw object is not generated properly.")
  }

  plot_title <- paste0("Displaying the ", y_label)

  g <- ggplot2::ggplot(data = subset_data, ggplot2::aes(x = as.factor(id),
                                                        y = counter_weight)) +
    ggplot2::geom_point(shape = "|", color = "blue", size = 2) +
    ggplot2::geom_segment(ggplot2::aes(xend = as.factor(id), yend = 0),
                          color = "red", linewidth = 0.2) +
    ggplot2::labs(x = "ID", y = y_label, title = plot_title) +
    ggplot2::theme_minimal() +
    ggplot2::coord_flip()

  if (every_n > 1){
    g <- g + ggplot2::scale_x_discrete(breaks = unique(as.factor(subset_data$id))[c(rep(FALSE, every_n-1), TRUE)])  # Show ID labels every 20 data points
  }

  g <- cowplot::plot_grid(
    g, hist_density_plot,
    ncol = 2,
    rel_widths = c(1, 1)
  )


  return(g)
}


#' @title
#' Extend generic plot functions for cgps_cw class
#'
#' @description
#' A wrapper function to extend generic plot functions for cgps_cw class.
#'
#' @param x  A cgps_cw object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @details
#' Additional parameters:
#' - *every_n*: Puts label to ID at every n interval (default = 10)
#' - *subset_id*: A vector of range of ids to be included in the plot
#' (default = NULL)
#'
#' @return
#' Returns a ggplot2 object, invisibly. This function is called for side effects.
#'
#' @export
#'
plot.cgps_cw <- function(x, ...) {
  g <- ggplot2::autoplot(x, ...)
  print(g)
  invisible(g)
}




#' @title
#' A helper function for cgps_erf object
#'
#' @description
#' A helper function to plot cgps_erf object using ggplot2 package.
#'
#' @param object A cgps_erf object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @return
#' Returns a ggplot object.
#'
#'
#' @keywords internal
#' @importFrom ggplot2 autoplot
#' @importFrom rlang .data
#'
autoplot.cgps_erf <- function(object, ...){


  x <- y_original <- normalized_weight <- w_vals <- y_pred <- NULL


  df1 <- object$.data_original
  df2 <- object$.data_prediction

  ## collect additional arguments
  dot_args <- list(...)
  arg_names <- names(dot_args)

  for (i in arg_names) {
    assign(i, unlist(dot_args[i], use.names = FALSE))
  }

  g <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df1, ggplot2::aes(x,
                                                 y_original,
                                                 size=normalized_weight),
                        color="blue", alpha=0.1) +
    ggplot2::geom_line(data = df2, ggplot2::aes(w_vals, y_pred), color = "orange") +
    ggplot2::geom_point(data = df2, ggplot2::aes(w_vals, y_pred), color = "orange") +
    ggplot2::labs(x = "Exposure",
                  y = "Outcome") +
    ggplot2::ggtitle(paste0("Exposure Response Curve for ",
                            object$params$model_type, " model")) +
    ggplot2::theme_bw()

  return(g)
}


#' @title
#' Extend generic plot functions for cgps_cw class
#'
#' @description
#' A wrapper function to extend generic plot functions for cgps_cw class.
#'
#' @param x  A cgps_erf object.
#' @param ... Additional arguments passed to customize the plot.
#'
#' @details
#' TBD
#'
#' @return
#' Returns a ggplot2 object, invisibly. This function is called for side effects.
#'
#' @export
#'
plot.cgps_erf <- function(x, ...) {
  g <- ggplot2::autoplot(x, ...)
  print(g)
  invisible(g)
}

