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
SL.mqch <- function(Y, X, newX, family, obsWeights, model = FALSE, ...) {
  
  # X must be a dataframe, not a matrix.
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  
  group <- X[ ,names(X) %in% c("group")]
  X <- X[ ,!names(X) %in% c("group")]
  
  Q = sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))  
  
  # Control -----------------------------------------
  M_y00 <- QRLM(y = Y,
                x = cbind(1, X), 
                q = Q,
                k = 1.345, 
                maxit = 100) 
  
  q0 <- matrix(c(gridfitinter(y  = Y,
                              expectile = M_y00$fitted.values,
                              Q  = M_y00$q.values)),
               nrow = length(Y), 
               ncol = 1)
  
  qmat0 <- matrix(c(q0, group), nrow = length(Y),
                  ncol = 2)
  
  Qi0 <- tapply(qmat0[,1], qmat0[,2], mean)
  
  fit <- QRLM(y = Y,
               x = cbind(1, X),
               q = Qi0, 
               k = 1.345, 
               maxit = 100)
  
#  fit <- list(object = fit, family = family)
  class(fit) <- "SL.mqch"
  
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
predict.SL.mqch <- function(object, newdata, ...) {
  
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  
  group <- newdata[ ,names(newdata) %in% c("group")]
  newdata <- newdata[ ,!names(newdata) %in% c("group")]
  group_names <- as.data.frame(table(group))$group
  

  pred <- NULL
  
  for (i in 1:length(group_names))
  { 
    pred <- c(pred, (cbind(1, as.matrix(newdata[group == group_names[i], ])) %*%  object$coef[, i]) ) 
    
  }

#  pred <- predict(object = object$object,
#                  newdata = newdata)
  pred

}


