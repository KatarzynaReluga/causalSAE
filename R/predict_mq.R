#' Making predictions for mquantfit 
#'
#' @param object A \code{gab} object
#' @param newdata New data (default is the training data for \code{object})
#' @param regression_type Type of regression "continuous", "binary", "nb", "poisson"
#' @param ... Additional parameters.
#' 
#' @return a matrix, each row corresponds to one iteration.
#'
#' @importFrom stats predict
#' @importFrom dplyr select	%>%
#'
#' @export
#'
predict.mquantfit <- function(object, newdata = object$data, 
                              regression_type = c("continuous", "binary", "nb", "poisson"),
                              ...) {
  
  regression_type <- match.arg(regression_type)
  
  newdata <- 
    newdata %>% select(colnames(object$coefficients)[-c(1:2)])
  
  list_coef <- list()
  for (i in 1:dim(object$coefficients)[1]) {
    list_coef[[i]] <- (object$coefficients)[i,]
  }

  predict_one_q <- function(coef_mq, newdata, regression_type) {
    if (regression_type == "continuous") {
      predict_y <- coef_mq[2] + as.matrix(newdata) %*%t(coef_mq)[-c(1:2)]
    } else if (regression_type == "poisson") {
      predict_y <- exp(coef_mq[2] + as.matrix(newdata) %*%t(coef_mq)[-c(1:2)])
    } else if (regression_type == "binary"){
      exp_y <- exp(coef_mq[2] + as.matrix(newdata) %*%t(coef_mq)[-c(1:2)])
      predict_y <- exp_y * (1 + exp_y)^{-1}
    }
    predict_y
  }

  predict_outcome <- lapply(list_coef, predict_one_q, newdata, regression_type)
  names(predict_outcome) <- paste0("q_", (object$coefficients)[,1])
  predict_outcome
    
}