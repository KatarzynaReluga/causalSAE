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
SL.mqchb <- function(Y, X, newX, family, obsWeights, model = FALSE, ...) {
  
  # X must be a dataframe, not a matrix.
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  
  group <- X[ ,names(X) %in% c("group")]
  X <- X[ ,!names(X) %in% c("group")]
  group_names <- as.data.frame(table(group))$group
  
  Q <-c(0.25, 0.30, 0.4, seq(from = 0.45, to = 0.55, by = 0.005), 0.60, 0.65, 0.75)
  
  # Control -----------------------------------------
  M_p0 <- QRLM(y = Y,
                x = cbind(1, as.matrix(X)), 
                q = Q,
                k = 1.6, 
                maxit = 100) 
  
  tmp.scores <- QSCORE(Y, 
                       M_p0$fitted, 
                       Q)
  scores = (cbind(group, 
                  tmp.scores))

#  tmp.scores <- QSCORE(data_fullp$A, 
#                       M_p0$fitted, 
#                       Q)
#  scores = (cbind(data_fullp$group, 
#                  tmp.scores))
  
  
  MQ_p = rep(0, length(group_names))
  for (i in 1:length(group_names)){
    MQ_p[i] = (mean(scores[, 2][scores[, 1]==group_names[i]]))
  }
  
  
  fit <- glm.mq.binom(y = Y,
                      x = cbind(1, as.matrix(X)),
                      q = MQ_p,
                      k = 1.6,
                      maxit = 100) 
  

#  fit <- list(object = fit, family = family)
  class(fit) <- "SL.mqchb"
  
  # newX must be a dataframe, not a matrix.
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }

  
  pred <- predict(fit,
                  newdata = newX)
  
  out <- list(pred = pred, fit = fit)
  return(out)
 
}
#' 
#' @title Prediction for SL.mq
#' @description Prediction for SL.mqch
#'
#' @param object SL.mq object
#' @param newdata Dataframe to generate predictions
#' @param ... Unused additional arguments
#'
#' @export
#' 
#' 
predict.SL.mqchb <- function(object, newdata, ...) {
  
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  
  group <- newdata[ ,names(newdata) %in% c("group")]
  newdata <- as.matrix(newdata[ ,!names(newdata) %in% c("group")])
  group_names <- as.data.frame(table(group))$group
  

  pred <- NULL
  
#  for (i in 1:length(group_names))
#  { 
#    pred <- c(pred, (cbind(1, as.matrix(newdata[group == group_names[i], ])) %*%  object$coef[, i]) ) 
#    
#  }
#  p_hat_Mch <- NULL
  
  # Check this against the old example
  for (i in 1:length(group_names)){
    pred <- c(pred, exp(cbind(1, as.matrix(newdata[group == group_names[i], ])) %*% object$coefficients[, i])/
                     (1 + exp(cbind(1, as.matrix(newdata[group == group_names[i], ])) %*% object$coefficients[, i])))  
  }

  
  pred

}


