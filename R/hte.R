#' Heterogenous treatment effect
#'
#' Calculates heterogrenous treatment effects for subpopulations
#'
#' @inheritParams impute_y
#' @param type_hte Type of estimation method, choose between: \code{OR}, \code{IPW},
#' \code{NIPW}, \code{AIPW}.
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
#' @param boot_var Compute bootstrap variance? Default: \code{boot_var = FALSE}
#' @param ... Additional parameters
#'
#' @importFrom lme4 lmer glmer
#' @importFrom grf regression_forest
#' @importFrom dplyr arrange select
#'
#' @export
#'
#' @examples
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
#' hte_OR <- hte(type_hte = "OR",
#'               data_sample,
#'               data_out_of_sample,
#'               params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
#'                                method = "EBLUP",
#'                                tune_RF = FALSE,
#'                                xgboost_params = list(CV_XGB = TRUE,
#'                                                      nfolds = 5,
#'                                                      nrounds = 50),
#'                                                      type_model = "gaussian"),
#'                                      boot_var = FALSE)
#'
#' hte_IPW <- hte(type_hte = "NIPW",
#'               data_sample,
#'               data_out_of_sample,
#'               params_p_score = list(model_formula = A ~ X1 + (1|group),
#'                           method = "RF",
#'                           tune_RF = FALSE,
#'                           xgboost_params = list(CV_XGB = TRUE,
#'                                                 nfolds = 5,
#'                                                 nrounds = 50)),
#'               params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
#'                                      method = "EBLUP",
#'                                      tune_RF = FALSE,
#'                                      xgboost_params = list(CV_XGB = TRUE,
#'                                                            nfolds = 5,
#'                                                            nrounds = 50),
#'                                      type_model = "gaussian"),
#'                                      boot_var = FALSE)
#'
#'
#' hte_AIPW <- hte(type_hte = "NIPW",
#'               data_sample,
#'               data_out_of_sample,
#'               params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
#'                                method = "EBLUP",
#'                                tune_RF = FALSE,
#'                                xgboost_params = list(CV_XGB = TRUE,
#'                                                      nfolds = 5,
#'                                                      nrounds = 50),
#'                                                      type_model = "gaussian"),
#'               params_p_score = list(model_formula = A ~ X1 + (1|group),
#'                           method = "EBLUP",
#'                           tune_RF = FALSE,
#'                           xgboost_params = list(CV_XGB = TRUE,
#'                                                 nfolds = 5,
#'                                                 nrounds = 50)),
#'               params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
#'                                      method = "RF",
#'                                      tune_RF = FALSE,
#'                                      xgboost_params = list(CV_XGB = TRUE,
#'                                                            nfolds = 5,
#'                                                            nrounds = 50),
#'                                      type_model = "gaussian"),
#'               boot_var = FALSE)
#'
#'
hte <- function(type_hte = c("OR", "IPW", "NIPW", "AIPW"),
                data_sample,
                data_out_of_sample,
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
                                       type_model = "gaussian"),
                boot_var = FALSE,
                ...) {

  type_hte <- match.arg(type_hte)
  obj_hte <- list(data_sample = data_sample,
                  data_out_of_sample = data_out_of_sample)
  class(obj_hte) <- type_hte

  # Estimate heterogenous treatment effects ---------------
  estimate_hte <- estimate_hte(obj_hte = obj_hte,
                               params_p_score = params_p_score,
                               params_impute_y = params_impute_y,
                               params_OR = params_OR)

  if (boot_var) {

  }

  return(estimate_hte)
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
#' @describeIn estimate_hte Estimate heterogeneous treatment effects using OR
#' @export
#'


estimate_hte.OR <- function(obj_hte,
                            params_OR,
                            ...) {


  data_OR <- fit_OR(obj_hte, params_OR)

  tau_treat = aggregate(data_OR$mu1_y, list(data_OR$group), FUN = mean)$x
  tau_untreat = aggregate(data_OR$mu0_y, list(data_OR$group), FUN = mean)$x
  tau = tau_treat - tau_untreat

  output <- data.frame(tau_treat = tau_treat,
                 tau_untreat = tau_untreat,
                 tau = tau,
                 group_name = unique(data_OR$group))
  return(output)

  }

#'
#' @describeIn estimate_hte Estimate heterogeneous treatment effects using IPW
#' @export
#'


estimate_hte.IPW <- function(obj_hte,
                             params_p_score,
                             params_impute_y,
                             ...) {

  # Obtain IPW data -------------------------------------------
a = Sys.time()
  IPW_data <- obtain_IPW_data(obj_hte,
                              params_p_score,
                              params_impute_y)
b = Sys.time()
  tau_data <- as.data.frame(calculate_tau(IPW_data, type_tau = "HT"))
  return(tau_data)
}

#'
#' @describeIn estimate_hte Estimate heterogeneous treatment effects using NIPW
#' @export
#'


estimate_hte.NIPW <- function(obj_hte,
                             params_p_score,
                             params_impute_y,
                             ...) {


  # Obtain IPW data  -------------------------------
  IPW_data <- obtain_IPW_data(obj_hte,
                              params_p_score,
                              params_impute_y)

  tau_data <- as.data.frame(calculate_tau(IPW_data, type_tau = "H"))
  return(tau_data)
}


#'
#' @describeIn estimate_hte Estimate heterogeneous treatment effects using AIPW
#' @export
#'


estimate_hte.AIPW <- function(obj_hte,
                              params_p_score,
                              params_impute_y,
                              params_OR,
                              ...) {


  if (params_OR$method == params_impute_y$method) {
    stop("In AIPW, the method to impute outcomes cannot be the same as the method for fitting outcome regression model. Choose another estimation method.")
  }

  # Obtain IPW data ----------------------------------------------
  IPW_data <- obtain_IPW_data(obj_hte,
                              params_p_score,
                              params_impute_y)
  # Obtain outcome regression data -----------------------------------------
  data_OR <- fit_OR(obj_hte, params_OR)

  AIPW_data <- IPW_data
  AIPW_data$mu0_y <- data_OR$mu0_y
  AIPW_data$mu1_y <- data_OR$mu1_y


  tau_data <- as.data.frame(calculate_tau(AIPW_data, type_tau = "AIPW"))
  return(tau_data)

}



#' Fit propensity score
#'
#' Internal function to obtain subpopulation heterogenous treatment effects
#'
#' @inheritParams hte
#' @param obj_hte Object to estimate hte
#'
#' @importFrom dplyr select
#'
#' @return Estimates of propensity score.
#'


fit_p_score <- function(obj_hte, params_p_score) {

  a = Sys.time()
  data_out_of_sample <- obj_hte$data_out_of_sample
  names_out_of_sample <- names(data_out_of_sample)

  data_sample <- obj_hte$data_sample
  names_data_sample <- names(data_sample)

  if (length(names_out_of_sample) == length(names_data_sample)) {
    data_p_score <- rbind(data_sample_p_score, data_out_of_sample_p_score)
  } else {
    y <- data_sample$y
    data_sample_p_score <- select(data_sample, -y)
    data_p_score <- rbind(data_sample_p_score, data_out_of_sample_p_score)
  }

#  if ("y" %in% names_out_of_sample) {
#    #    y <- data_out_of_sample$y
#    y <- data_out_of_sample$y
#    data_out_of_sample_p_score <- select(data_out_of_sample, -y)

#    y <- data_sample$y
#    data_sample_p_score <- select(data_sample, -y)

#    data_p_score <- rbind(data_sample_p_score, data_out_of_sample_p_score)
#    obj_p_score <- list(data_p_score = data_p_score)
#  } else {
#    y <- data_sample$y
#    data_sample_p_score <- select(data_sample, -y)
#    data_p_score <- rbind(data_sample_p_score, data_out_of_sample_p_score)
#    obj_p_score <- list(data_p_score = data_p_score)
#  }
  b = Sys.time()
  class(obj_p_score) <- params_p_score$method

  ps_hat <- p_score(obj_p_score = obj_p_score,
                    model_formula = params_p_score$model_formula,
                    xgboost_params = params_p_score$xgboost_params,
                    tune_RF = params_p_score$tune_RF)

  return(ps_hat)
}

#' Obtain IPW data
#'
#' Internal function to obtain data to fit IPW and NIPW estimators
#'
#' @inheritParams hte
#' @param obj_hte Object to estimate hte
#'
#' @importFrom dplyr select
#'
#' @return Estimates of propensity score.
#'

obtain_IPW_data <- function(obj_hte,
                         params_p_score,
                         params_impute_y) {
  # Obtain fitted propensity score -------------------------------
  a = Sys.time()
  fitted_p_score <- fit_p_score(obj_hte = obj_hte,
                                params_p_score = params_p_score)

  b = Sys.time()
  b -a
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

  y_full_impute <- imputed_y$y_full_imputed

  IPW_data <- data.frame(y = y_full_impute,
                                 A = c(data_sample$A, data_out_of_sample$A),
                                 group = c(data_sample$group, data_out_of_sample$group),
                                 p_score  = fitted_p_score)
  return(IPW_data)

}
