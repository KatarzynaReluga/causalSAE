
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
