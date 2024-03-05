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
  if (is.data.frame(data_obj)){
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
