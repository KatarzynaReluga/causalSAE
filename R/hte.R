#' Heterogenous treatment effect
#'
#' Calculates heterogrenous treatment effects for subpopulations
#'
#' @inheritParams impute_y
#' @param type_hte Type of estimation method, choose between: \code{OR}, \code{IPW},
#' \code{NIPW}, \code{AIPW}.
#' @param p_score_estimate Vector of estimated propensity scores. Default: \code{p_score_estimate = NULL}.
#' @param params_p_score List with parameters to fit propensity score:
#'  \itemize{
#'  \item model_formula - model formula,
#'  \item method - estimation method, choose between: \code{EBLUP},
#'  \code{MQ}, \code{RF}, \code{XGB},
#'  \item tune_RF - tune random forest parameters,
#'  \item xgboost_params - list of parameters to obtain predictions
#'  using gradient boosting, default: \code{CV_XGB = TRUE},
#'  \code{nfolds = 5}, \code{nrounds = 50}.
#'  }
#' @param params_impute_y List with parameters to fit imputation model:
#'  \itemize{
#'  \item model_formula - model formula,
#'  \item method - estimation method, choose between: \code{EBLUP},
#'  \code{MQ}, \code{RF}, \code{XGB},
#'  \item tune_RF - tune random forest parameters,
#'  \item xgboost_params - list of parameters to obtain predictions
#'  using gradient boosting, default: \code{CV_XGB = TRUE},
#'  \code{nfolds = 5}, \code{nrounds = 50}.
#'  \item type_model - type of outcome.
#'  }
#' @param params_OR List with parameters to fit outcome regression model:
#'  \itemize{
#'  \item model_formula - model formula,
#'  \item method - estimation method, choose between: \code{EBLUP},
#'  \code{MQ}, \code{RF}, \code{XGB},
#'  \item tune_RF - tune random forest parameters,
#'  \item xgboost_params - list of parameters to obtain predictions
#'  using gradient boosting, default: \code{CV_XGB = TRUE},
#'  \code{nfolds = 5}, \code{nrounds = 50}.
#'  \item type_model - type of outcome.
#' }
#' @param ... Additional parameters
#'
#' @importFrom lme4 lmer glmer
#' @importFrom grf regression_forest
#' @importFrom dplyr arrange select
#'
#' @examples
#'
#' #' @examples
#'
#' m = 50
#' ni = rep(10, m)
#' Ni = rep(200, m)
#' N = sum(Ni)
#' n = sum(ni)
#'
#' X <- generate_X(
#'  n = N,
#'  p = 1,
#'  covariance_norm = NULL,
#'  cov_type = "unif",
#'  seed = 1
#' )
#'
#' X_outcome <- generate_X(
#'  n = N,
#'  p = 1,
#'  covariance_norm = NULL,
#'  cov_type = "lognorm",
#'  seed = 1
#' )
#'
#' populations <- generate_pop(X, X_outcome,
#' coeffs = get_default_coeffs(),
#' errors_outcome = get_default_errors_outcome(),
#' rand_eff_outcome = get_default_rand_eff_outcome(),
#' rand_eff_p_score = get_default_rand_eff_p_score(),
#' regression_type = "continuous",
#' Ni_size  = 200,
#' m = 50,
#' no_sim = 1,
#' seed = 10)
#'
#' samples <- generate_sample(populations, ni_size = 10,
#'                            sample_part = "sampled",
#'                            get_index = TRUE)
#'
#' data_sample <- data.frame(samples[[1]]$samp_data)
#' index_sample <- samples[[1]]$index_s
#' data_out_of_sample <- populations[-index_sample, ]
#'
#' model_formula_OR = y ~ X1 + Xo1 + (1|group)
#'
#'

hte <- function(type_hte = c("OR", "IPW", "NIPW", "AIPW"),
                data_sample,
                data_out_of_sample,
                p_score_estimate = NULL,
                params_p_score = list(model_formula = A ~ X1 + (1|group),
                                      method = "EBLUP",
                                      tune_RF = FALSE,
                                      xgboost_params = list(CV_XGB = TRUE,
                                                            nfolds = 5,
                                                            nrounds = 50)),
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                 method = "EBLUP",
                                 tune_RF = FALSE,
                                 xgboost_params = list(CV_XGB = TRUE,
                                                       nfolds = 5,
                                                       nrounds = 50),
                                 type_model = "gaussian"),
                params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                       method = "EBLUP",
                                       tune_RF = FALSE,
                                       xgboost_params = list(CV_XGB = TRUE,
                                                             nfolds = 5,
                                                             nrounds = 50),
                                       type_model = "gaussian"), ...) {

  type_hte <- match.arg(type_hte)
  obj_hte <- list(data_sample = data_sample,
                  data_out_of_sample = data_out_of_sample)
  class(obj_hte) <- type_hte



}


#' Estimate hte
#'
#' Internal generic function to estimate hte
#' @inheritParams hte
#' @param obj_hte Object to estimate hte
#'
#' @return List with following parameters:
#' \item{y_hat_out_of_sample}{imputed out of sample outcomes}
#' \item{y_hat_sample}{predicted sample data using fitted model}
#' \item{outcome_fit}{fitted model to impute the data}
#'
#' @importFrom lme4 lmer glmer
#' @importFrom stats terms
#' @importFrom grf regression_forest
#' @importFrom xgboost xgboost xgb.cv
#' @importFrom mquantreg mquantreg
#'


estimate_hte <- function(...)
  UseMethod("estimate_hte")


#'
#' @describeIn imputation_y Function to estimate imputation_y using EBLUP
#' @export
#'


estimate_hte.OR <- function(obj_hte,
                            params_OR,
                            ...) {
  # Full data --------------------------
  data_full <- rbind(obj_hte$data_sample, obj_hte$data_out_of_sample)
#  data_sample  = obj_hte$data_sample
#  data_out_of_sample = obj_hte$data_out_of_sample

  # Controls -------------------------------------
  data_sample0 = data_sample[data_sample$A == 0, ]

  OR0 <- impute_y(model_formula = params_OR$model_formula,
                  data_sample = data_sample0,
                  data_out_of_sample = data_out_of_sample,
                  method = params_OR$method,
                  type_model = params_OR$type_model,
                  tune_RF = params_OR$tune_RF,
                  xgboost_params = params_OR$xgboost_params)


  mu0_y <- predict(object = OR0$outcome_fit, newdata = data_full,
                   allow.new.levels = TRUE)

  # Treated --------------------------------------
  data_sample1 = data_sample[data_sample$A == 1, ]

  OR1 <- impute_y(model_formula = params_OR$model_formula,
                  data_sample = data_sample1,
                  data_out_of_sample = data_out_of_sample,
                  method = params_OR$method,
                  type_model = params_OR$type_model,
                  tune_RF = params_OR$tune_RF,
                  xgboost_params = params_OR$xgboost_params)

  mu1_y <- predict(OR1$outcome_fit, newdata = data_full,
                   allow.new.levels = TRUE)

  data_OR <- data.frame(mu1_y = mu1_y,
                        mu0_y = mu0_y,
                        group = c(data_sample$group, data_out_of_sample$group))

  tau_1 = aggregate(data_OR$mu1_y, list(data_OR$group), FUN = mean)$x
  tau_0 = aggregate(data_OR$mu0_y, list(data_OR$group), FUN = mean)$x
  tau = tau_1 - tau_0

  output <- list(tau_1 = tau_1,
                 tau_0 = tau_0,
                 tau = tau,
                 group = )
  return(output)

  }

#'
#' @describeIn imputation_y Function to estimate imputation_y using EBLUP
#' @export
#'


estimate_hte.IPW <- function(obj_hte,
                             params_p_score,
                             params_impute_y,
                             ...) {


  # Obtain fitted propensity score -------------------------------
  fitted_p_score <- fit_p_score(obj_hte = obj_hte,
                                params_p_score = params_p_score)


  # Obtain imputation model -------------------------------
  data_sample  = obj_hte$data_sample
  data_out_of_sample = obj_hte$data_out_of_sample

  imputed_y <- impute_y(model_formula = params_impute_y$model_formula,
                  data_sample = data_sample,
                  data_out_of_sample = data_out_of_sample,
                  method = params_impute_y$method,
                  type_model = params_impute_y$type_model,
                  tune_RF = params_impute_y$tune_RF,
                  xgboost_params = params_impute_y$xgboost_params)

}

#' Fit propensity score
#'
#' Internal function to obtain subpopulation heterogenous treatment effects
#'
#' @inheritParams hte
#'
#' @importFrom dplyr select
#'
#' @return Estimates of propensity score.
#'


fit_p_score <- function(obj_hte, params_p_score) {

  data_out_of_sample <- obj_hte$data_out_of_sample
  names_out_of_sample <- names(data_out_of_sample)

  data_sample <- obj_hte$data_sample


  if ("y" %in% names_out_of_sample) {
    #    y <- data_out_of_sample$y
    data_out_of_sample_p_score <- select(data_out_of_sample, -y)
    data_sample_p_score <- select(data_sample, -y)
    data_p_score <- rbind(data_sample_p_score, data_out_of_sample_p_score)
    obj_p_score <- list(data_p_score = data_p_score)
  } else {
    data_sample_p_score <- select(data_sample, -y)
    data_p_score <- rbind(data_sample_p_score, data_out_of_sample_p_score)
    obj_p_score <- list(data_p_score = data_p_score)
  }

  class(obj_p_score) <- params_p_score$method

  ps_hat <- p_score(obj_p_score = obj_p_score,
                    model_formula = params_p_score$model_formula,
                    xgboost_params = params_p_score$xgboost_params,
                    tune_RF = params_p_score$tune_RF)
  return(ps_hat)
}



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
