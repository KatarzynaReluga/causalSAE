#' @title Wrapper for glmm
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
#'
#' @export
#' 
#' @importFrom lme4 glmer lmer
#' 
SL.glmm <- function(Y, X, newX, family, obsWeights, model = FALSE, ...) {
  
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
 
  fit <- lmer(formula = model_formula, 
              data = df)
  
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  
  pred <- unname(predict(fit, newdata = newX, 
                         allow.new.levels = TRUE))
  

  fit <- list(object = fit, family = family)
  class(fit) <- "SL.glmm"
  
  out <- list(pred = pred, fit = fit)
  
  return(out)
}

#' @title Prediction for SL.glmm
#' @description Prediction for SL.glmm
#'
#' @param object SL.glmm object
#' @param newdata Dataframe to generate predictions
#' @param ... Unused additional arguments
#'
#'
#' @export
#' 
#' @importFrom lme4 glmer
#' 
predict.SL.glmm <- function(object, newdata, ...) {
  
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  
  pred <-  unname(predict(object = object$object, newdata = newdata, 
                  allow.new.levels = TRUE))

  pred
}
