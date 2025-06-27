#' @title Wrapper for grf
#' @description Wrapper for generalised random forests with tuning via grf().
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
#' @importFrom grf regression_forest
#' 
SL.grf_t <- function(Y, X, newX, family, obsWeights, model = FALSE, ...) {
  
  # X must be a dataframe, not a matrix.
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }

  clusters <- X$group
    
  X <- X[ ,!names(X) %in% c("group")]

  
  fit <- regression_forest(X, Y, clusters = clusters,
                           tune.parameters = "all")
  

   if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  
  newX <- newX[ ,!names(newX) %in% c("group")]
  
  pred <- unname(unlist(predict(fit, newdata = newX)))
  
  # For binomial family restrict predicted probability to [0, 1].
  if (family$family == "binomial") {
    pred = pmin(pmax(pred, 0), 1)
  }
  
  fit <- list(object = fit, family = family)
  class(fit) <- "SL.grf_t"
  
  out <- list(pred = pred, fit = fit)
  
  return(out)
}

#' @title Prediction for SL.grf
#' @description Prediction for SL.grf
#'
#' @param object SL.grf object
#' @param newdata Dataframe to generate predictions
#' @param ... Unused additional arguments
#'
#'
#' @export
predict.SL.grf_t <- function(object, newdata, ...) {
  
  newdata <- newdata[ ,!names(newdata) %in% c("group")]
  
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  
  pred <- unname(unlist(predict(object = object$object, newdata = newdata)))
  
  pred
}

