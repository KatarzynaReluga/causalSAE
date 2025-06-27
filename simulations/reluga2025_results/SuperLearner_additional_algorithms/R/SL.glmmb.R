#' @title Wrapper for glmmb
#' @description Wrapper for glmm via lm4().
#'
#' @param Y Outcome variable
#' @param X Training dataframe
#' @param newX Test dataframe
#' @param family Gaussian or binomial
#' @param obsWeights Observation-level weights
#' @param model Whether to save model.matrix of data in fit object. Set to FALSE
#'   to save memory.
#' @param ... Any remaining arguments, not used.
#'
#'
#' @export
#' 
#' @importFrom lme4 glmer lmer
#' 
SL.glmmb <- function(Y, X, newX, family, obsWeights, model = FALSE, ...) {
  
  # X must be a dataframe, not a matrix.
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  
  varsfe <- names(X[ ,!names(X) %in% c("group")])
  varsre <- "(1|group)"
  vars <- c(varsfe, varsre)
  predictors <- paste(vars, collapse = " + ")
  model_formula <- paste("Y", " ~ ", predictors)
  
  df <- data.frame(Y, X)
  
  
  fit <- glmer(formula = model_formula, 
                 data = df, 
                 family = binomial)
 
  
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  
  pred0 <- unname(predict(fit, newdata = newX, 
                          allow.new.levels = TRUE))
  pred <- pmin(pmax(pred0, 0), 1)

  fit <- list(object = fit, family = family)
  class(fit) <- "SL.glmmb"
  
  out <- list(pred = pred, fit = fit)
  
  return(out)
}

#' @title Prediction for SL.glmmb
#' @description Prediction for SL.glmmb
#'
#' @param object SL.glmmb object
#' @param newdata Dataframe to generate predictions
#' @param ... Unused additional arguments
#'
#'
#' @export
#' 
#' @importFrom lme4 glmer
#' 
predict.SL.glmmb <- function(object, newdata, ...) {
  
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  
  pred <-  unname(predict(object = object$object, newdata = newdata, 
                          allow.new.levels = TRUE))
  
  pred
}
