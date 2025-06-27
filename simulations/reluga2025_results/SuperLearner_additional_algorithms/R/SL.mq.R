#' @title Wrapper for mq
#' @description Wrapper for MQ regression via mquantreg().
#'
#'
#' @param Y Outcome variable
#' @param X Training dataframe
#' @param newX Test dataframe
#' @param family Gaussian or binomial
#' @param obsWeights Observation-level weights
#' @param model Whether to save model.matrix of data in fit object. Set to FALSE
#' to save memory.
#' @param ... Any remaining arguments, not used.
#'
#'
#' @export
#' 
#' @importFrom mquantreg mquantreg
#' 
SL.mq <- function(Y, X, newX, family, obsWeights, model = FALSE, ...) {
  
  # X must be a dataframe, not a matrix.
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  
  vars <- names(X[ ,!names(X) %in% c("group")])
  predictors <- paste(vars, collapse = " + ")
  model_formula <- paste("Y", " ~ ", predictors)
  
  df <- data.frame(Y, X)
  

  fit = mquantreg(formula = model_formula,
                    data = df,
                    q  = 0.5,
                    method = "continuous")
  
  # newX must be a dataframe, not a matrix.
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  
  pred <- predict(fit,
                 newdata = newX,
                 regression_type = "continuous")
  
  fit <- list(object = fit, family = family)
  class(fit) <- "SL.mq"
  out <- list(pred = pred, fit = fit)
  return(out)
}
#' 
#' @title Prediction for SL.mq
#' @description Prediction for SL.mq
#'
#' @param object SL.mq object
#' @param newdata Dataframe to generate predictions
#' @param ... Unused additional arguments
#'
#' @export
#' 
#' 
predict.SL.mq <- function(object, newdata, ...) {
  
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  pred <- predict(object = object$object,
                  newdata = newdata)
  pred

}


