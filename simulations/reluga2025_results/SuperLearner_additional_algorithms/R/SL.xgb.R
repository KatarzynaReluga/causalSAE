#' @title Wrapper for xgb
#' @description Wrapper for generalised random forests with tuning via xgboost().
#'
#' Note that for outcomes bounded by [0, 1] the binomial family can be used in
#' addition to gaussian.
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
#' @importFrom xgboost xgboost xgb.cv
#' 
SL.xgb <- function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
  
  if (!is.matrix(X)) {
    X = data.matrix(X)
  }
  

  
  fit <- xgboost(data = X, label = Y,
                 nrounds = 100,
                 verbose = FALSE)
  
  

  if (!is.matrix(newX)) {
    newX = data.matrix(newX)
  }
  
  
  pred <- unname(unlist(predict(fit, newdata = newX)))
  
  # For binomial family restrict predicted probability to [0, 1].
  #if (family$family == "binomial") {
  #  pred = pmin(pmax(pred, 0), 1)
  #}
  
  fit <- list(object = fit, family = family)
  class(fit) <- "SL.xgb"
  
  out <- list(pred = pred, fit = fit)
  
  return(out)
}

#' @title Prediction for SL.xgb
#' @description Prediction for SL.xgb
#'
#' @param object SL.xgb object
#' @param newdata Dataframe to generate predictions
#' @param ... Unused additional arguments
#'
#'
#' @export
predict.SL.xgb <- function(object, newdata, ...) {
  
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  
  pred <- unname(unlist(predict(object = object$object, newdata = newdata)))
  
  pred
}

