#' Heterogenous treatment effect
#' 
#' Calculates heterogrenous treatment effects for subpopulations
#'
#' @param formula_y Formula for fitting outcome regression
#' @param formula_p_score Formula for fitting outcome regression
#' @param type_tau Type of computation, \code{regular}: constant weights, 
#' \code{adaptive}: adaptive weights
#' @param data_sample Data from the sample
#' @param data_out_of_sample Data from the population 
#' @param method_y Method to predict/fit the outcome
#' @param method_p_score Method to predict/fit propnsity score
#' @param tune_RF Tune random forest parameters
#' @param p_score_estimate Vector of estimated propensity scores
#' @param ... Additional parameters
#' 
#' @importFrom lme4 lmer glmer
#' @importFrom grf regression_forest 
#' @importFrom dplyr arrange select
#'
## #' @examples 
#' 
#' 

hte <- function(formula_y, 
                formula_p_score, 
                type_tau, 
                data_sample, 
                data_out_of_sample,
                method_y = c("EBLUP", "MQ", "RF"),
                method_p_score =  c("EBLUP", "MQ", "RF"), 
                p_score_estimate = NULL, 
                tune_RF = F, ...){
  
  method_y <- match.arg(method_y)
  method_p_score <- match.arg(method_p_score)
  
  y_hat <- outcome_regression(formula_y, 
                              data_sample, 
                              data_out_of_sample, 
                              method_y, 
                              family = "gaussian", 
                              method = "continuous", 
                              tune_RF)
 
  if (is.null(p_score_estimate)) {
    
    data_sample_p_score <- select(data_sample, -y)

    if (sum(names(data_out_of_sample) == "y") == 1) {
      data_out_of_sample_p_score <- select(data_out_of_sample, -y)
    } else {
      data_out_of_sample_p_score <- data_out_of_sample
    }
    
    data_p_score <- rbind(data_sample_p_score, data_out_of_sample_p_score)
    ps_hat <-  p_score(formula_p_score, 
                       data_sample = data_p_score, 
                       method_p_score = method_p_score)
  } 

  
  y_full <- c(data_sample$y, unlist(y_hat))
  group_full <- c(data_sample$group, data_out_of_sample$group)
  A_full <- c(data_sample$A, data_out_of_sample$A)
 
  pop_est <- data.frame(y = y_full, A = A_full, 
                        group = group_full, p_score = ps_hat)
 
 
#  pop_full_a <- arrange(pop_est, group)
  pop_full <- list(pop_est)
  
#  tau_hat <- calculate_tau(list(pop_full_a))
  tau_hat <- calculate_tau(list(pop_full), type_tau)
  
  output <- list(tau_hat = tau_hat, 
                  pop_full = pop_full) 
  output
}


#' Outcome regression
#' 
#' Fit outcome regression model
#' 
#'
#' @inheritParams hte
#' @param family Specification for the model link function in glmer()
#' @param method Specification for the model link function in mquantreg()
#' 
#' @importFrom lme4 lmer glmer
#' @importFrom stats terms
#' @importFrom grf regression_forest 
#' @importFrom mquantreg mquantreg
#' 
#' @export
#' 



outcome_regression <- function(formula_y, 
                               data_sample, 
                               data_out_of_sample, 
                               method_y, 
                               family = "gaussian", 
                               method = "continuous", 
                               tune_RF) {
  
  
  if (method_y ==  "EBLUP") {
    if (family == "gaussian") {
      
      outcome_fit <- lmer(formula_y, data = data_sample)
      y_hat <- predict(outcome_fit, newdata = data_out_of_sample, allow.new.levels = TRUE)
      
    } else {
      outcome_fit <- glmer(formula_y, data = data_sample, family = family) 
      y_hat <- as.vector(predict(outcome_fit, newdata = data_out_of_sample , allow.new.levels = TRUE))
    }
  } else {
    
    # Get the (predictor) variables
    vars <- attr(terms(formula_y), which = "term.labels")
    
    # Get the response
    response <- as.character(attr(terms(formula_y), which = "variables")[[2]]) # The response variable name
    
    # Get the predictors without group variable
    predictor <- paste(vars[-length(vars)], collapse = " + ")
    
    if (method_y ==  "MQ") {
      # Build a formula
      formula_yMQ <- paste(response, " ~ ", predictor)
      outcome_fit = mquantreg(formula = formula_yMQ, data = data_sample, q  = 0.5, method = method)
      y_hat <- as.vector(predict(outcome_fit, newdata = data_out_of_sample, 
                       regression_type = "continuous"))
      
    }  else {
      
      X = data_sample[, vars[ - length(vars)],  drop = F]
      X_newdata = data_out_of_sample[, vars[ - length(vars)],  drop = F]
      Y = data_sample[, response]
      clusters  = as.numeric(data_sample$group)
      
      if (tune_RF) {
        outcome_fit <- regression_forest(X, Y, 
                                         clusters = clusters, tune.parameters = "all")
      } else { 
        outcome_fit <- regression_forest(X, Y, clusters = clusters)
      }

      y_hat <- predict(outcome_fit, newdata = X_newdata)
      
    }
  }
  y_hat
  
 }

#'
#' Propensity score
#' 
#' Fit propensity score model
#' 
#'
#' @inheritParams hte
#' 
#' @importFrom lme4 lmer glmer
#' @importFrom stats terms binomial
#' @importFrom grf regression_forest 
#' @importFrom mquantreg mquantreg
#' 
#' @export
#' 


p_score <- function(formula_p_score, 
                    data_sample, 
                    data_out_of_sample = NULL,
                    method_p_score) {
  
  
  if (method_p_score == "EBLUP") {

      ps_fit <- glmer(formula_p_score, data = data_sample, family = binomial(link = "logit")) 
      ps_hat_sample <- as.vector(predict(ps_fit, newdata = data_sample, 
                                      type = "response", allow.new.levels = TRUE))
#      ps_hat_out <- as.vector(predict(ps_fit, newdata = data_out_of_sample, 
#                                  type = "response", allow.new.levels = TRUE))
#      ps_hat <- c(ps_hat_sample, ps_hat_out)
      
      if (is.null(data_out_of_sample)) {
        ps_hat <- c(ps_hat_sample)
      } else {
        ps_hat_out <- as.vector(predict(ps_fit, newdata = data_out_of_sample, 
                               type = "response", allow.new.levels = TRUE))
        ps_hat <- c(ps_hat_sample, ps_hat_out)
      }
      
  } else {
    
    # Get the (predictor) variables
    vars <- attr(terms(formula_p_score), which = "term.labels")
    
    # Get the response
    response <- as.character(attr(terms(formula_p_score), which = "variables")[[2]]) # The response variable name
    
    # Get the predictors without group variable
    predictor <- paste(vars[-length(vars)], collapse = " + ")
    
    if (method_p_score ==  "MQ") {
      # Build a formula
      formula_p_scoreMQ <- paste(response, " ~ ", predictor)
      ps_fit <- mquantreg(formula = formula_p_scoreMQ, data = data_sample, 
                         q  = 0.5, method = "binom")
      ps_hat_sample <- unlist(predict(ps_fit, newdata = data_sample, 
                            regression_type = "binary"))
#      ps_hat_out <- unlist(predict(ps_fit, newdata = data_out_of_sample, 
#                       regression_type = "binary"))
#      ps_hat <- c(ps_hat_sample, ps_hat_out)
      
      if (is.null(data_out_of_sample)) {
        ps_hat <- c(ps_hat_sample)
      } else {
        ps_hat_out <- unlist(predict(ps_fit, newdata = data_out_of_sample, 
                             regression_type = "binary"))
        ps_hat <- c(ps_hat_sample, ps_hat_out)
      }
      
    }  else {
      
      X = data_sample[, vars[ - length(vars)], drop = F]
      Y = data_sample[, response]
      clusters  = as.numeric(data_sample$group)
      ps_fit <- regression_forest(X, Y, clusters = clusters)
      
      ps_hat_sample <- unlist(predict(ps_fit, newdata = X), use.names = FALSE)
#      ps_hat_out <- predict(ps_fit, newdata = X_newdata)
#      ps_hat <- c(ps_hat_sample, ps_hat_out)
      
      if (is.null(data_out_of_sample)) {
        ps_hat <- c(ps_hat_sample)
      } else {
        X_newdata = data_out_of_sample[, vars[ - length(vars)], drop = F]
        ps_hat_out <- unlist(predict(ps_fit, newdata = X_newdata), use.names = FALSE)
        ps_hat <- c(ps_hat_sample, ps_hat_out)
      }
      
    }
  }
  ps_hat
}