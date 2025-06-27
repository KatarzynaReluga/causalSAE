#' @title Wrapper for xgb with tunning
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
SL.xgbt <- function(Y, X, newX, family, obsWeights, model = TRUE, ...) {
  
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
#  fit <- xgboost(data = X, label = Y,
#                 nrounds = 50,
#                 verbose = FALSE)
  
  
  xgboost_cv <- xgb.cv(data = X, label = Y,
                       nfold = 5,
                       nrounds = 100,
                       verbose = FALSE)
  best_iter_xgb = which.min(xgboost_cv$evaluation_log$test_rmse_mean)
  
  fit <- xgboost(data = X, label = Y,
                         nrounds = best_iter_xgb,
                         verbose = FALSE)
  
 
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  
  
  pred <- unname(unlist(predict(fit, newdata = newX)))
  
  # For binomial family restrict predicted probability to [0, 1].
  #if (family$family == "binomial") {
  #  pred = pmin(pmax(pred, 0), 1)
  #}
  
  fit <- list(object = fit, family = family)
  class(fit) <- "SL.xgbt"
  
  out <- list(pred = pred, fit = fit)
  
  return(out)
}

#' @title Prediction for SL.xgbt
#' @description Prediction for SL.xgbt
#'
#' @param object SL.xgbt object
#' @param newdata Dataframe to generate predictions
#' @param ... Unused additional arguments
#'
#'
#' @export
predict.SL.xgbt <- function(object, newdata, ...) {
  
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  
  pred <- unname(unlist(predict(object = object$object, newdata = newdata)))
  
  pred
}

