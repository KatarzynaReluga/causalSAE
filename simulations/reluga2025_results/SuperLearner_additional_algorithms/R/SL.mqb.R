#' @title Wrapper for mq binom
#' @description Wrapper for MQ regression via mquantreg() for binary outcomes.
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
SL.mqb <- function(Y, X, newX, family, obsWeights, model = FALSE, ...) {
  
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
                  method = "binom")
  
  # newX must be a dataframe, not a matrix.
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  
  pred <- predict(fit,
                  newdata = newX,
                  regression_type = "binary")
  
  fit <- list(object = fit, family = family)
  class(fit) <- "SL.mqb"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @title Prediction for SL.mqb
#' @description Prediction for SL.mqb
#'
#' @param object SL.mqb object
#' @param newdata Dataframe to generate predictions
#' @param ... Unused additional arguments
#'
# ' @seealso \code{\link{SL.glm}} \code{\link[stats]{glm}}
# '   \code{\link[stats]{predict.glm}}  \code{\link{SL.speedglm}}
# '
#' @export
#' 
#' 
predict.SL.mqb <- function(object, newdata, ...) {
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  #  pred <- predict(object = object$object, newdata = newdata, type = "response")
  
  pred <- predict(object = object$object,
                  newdata = newdata, 
                  regression_type = "binary")
  
  pred
  
  
}


