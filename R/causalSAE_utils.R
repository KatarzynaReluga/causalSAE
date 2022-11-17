#
# Format matrix in a suitable way
#

format_data_matrix <- function(data, select_row = 1:nrow(data),
                               select_col = NULL,
                               name_col = "X"){

  if (is.null(select_col)){
    data_new <- data.matrix(data[select_row, ])
  } else {
    data_new <- data.matrix(data[select_row, grep(select_col, colnames(data))])
  }
  if (ncol(data_new) == 1) {
    colnames(data_new) <- paste(name_col)
  } else {
    colnames(data_new) <- paste(name_col, 1:ncol(data_new), sep = "")
  }
  data_new
}

#
# Format data
#

format_data <- function(model_formula,
                        data_sample,
                        data_out_of_sample = NULL) {

  # Get the (predictor) variables
  vars <- attr(terms(model_formula), which = "term.labels")

  # Get the response
  response <- as.character(attr(terms(model_formula), which = "variables")[[2]])

  # Get the predictors without group variable
  predictor <- paste(vars[-length(vars)], collapse = " + ")

  X = as.matrix(data_sample[, vars[ - length(vars)],  drop = F])
  Y = data_sample[, response]
  clusters_sample  = as.numeric(data_sample$group)

  if (is.null(data_out_of_sample)) {
    X_newdata = NULL
    clusters_out_of_sample = NULL
  } else {
    X_newdata = as.matrix(data_out_of_sample[, vars[ - length(vars)],  drop = F])
    clusters_out_of_sample = as.numeric(data_out_of_sample$group)
  }


  formatted_data <- list(X = X,
                         X_newdata = X_newdata,
                         Y = Y,
                         clusters_sample = clusters_sample,
                         clusters_out_of_sample = clusters_out_of_sample)

  return(formatted_data)
}
