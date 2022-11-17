#' Impute y
#'
#' Impute missing values of y in the remaining part of population.
#'
#' @param data_sample Data from the sample.
#' @param data_out_of_sample Data from the remaining part of population.
#' @param model_formula Model formula.
#' @param method Estimation method, choose between: \code{EBLUP}, \code{MQ}, \code{RF}, \code{XGB}.
#' @param type_model Type of outcome.
#' @param tune_RF Tune parameters in random forest? Default = FALSE.
#' @param xgboost_params List with parameters to run xgboost:
#' \itemize{
#'  \item CV_XGB - logical variable, use cross-validation for gradient boosting? Default: TRUE.
#'  \item nfolds - number of folds in cross-validation, default: nfolds = 5.
#'  \item nrounds - the max number of iterations, default: nrounds = 50.
#' }
#' @param ... Additional parameters.
#'
#' @export
#'
#' @return List with following parameters:
#' \item{y_hat_out_of_sample}{imputed out of sample outcomes}
#' \item{y_hat_sample}{predicted sample data using fitted model}
#' \item{outcome_fit}{fitted model to impute the data}
#' \item{y_full_imputed}{original outcomes and imputed out of sample outcomes in one vector}
#' \item{group_sample}{group labels in the sample}
#' \item{group_out_of_sample}{group lables out of sample}
#'
#'
#' @examples
#'
#' m = 50
#' ni = rep(5, m)
#' Ni = rep(100, m)
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
#' Ni_size  = 100,
#' m = 50,
#' no_sim = 1,
#' seed = 10)
#'
#' samples <- generate_sample(populations, ni_size = 5,
#'                            sample_part = "sampled",
#'                            get_index = TRUE)
#'
#' data_sample <- data.frame(samples[[1]]$samp_data)
#' index_sample <- samples[[1]]$index_s
#' data_out_of_sample <- populations[-index_sample, ]
#'
#' model_formula = y ~ X1 + Xo1 + A + (1 + A||group)
#'
#' impute_EBLUP <- impute_y(model_formula,
#'                    data_sample,
#'                    data_out_of_sample,
#'                    method = "EBLUP",
#'                    type_model = "gaussian")
#'
#' impute_MQ <- impute_y(model_formula,
#'                    data_sample,
#'                    data_out_of_sample,
#'                    method = "MQ",
#'                    type_model = "continuous")
#'
#' impute_RF <- impute_y(model_formula,
#'                    data_sample,
#'                    data_out_of_sample,
#'                    method = "RF",
#'                    tune_RF = TRUE,
#'                    xgboost_params = list(CV_XGB = TRUE,
#'                                          nfolds = 5,
#'                                          nrounds = 50))
#'
#' impute_XGB <- impute_y(model_formula,
#'                    data_sample,
#'                    data_out_of_sample,
#'                    method = "XGB",
#'                    xgboost_params = list(CV_XGB = TRUE,
#'                                          nfolds = 5,
#'                                          nrounds = 50))
#'
#'


impute_y <- function(model_formula,
                     data_sample,
                     data_out_of_sample,
                     method = c("EBLUP", "MQ", "RF", "XGB"),
                     type_model = "gaussian",
                     tune_RF = FALSE,
                     xgboost_params = list(CV_XGB = TRUE,
                                           nfolds = 5,
                                           nrounds = 50),
                     ...) {

  # Check and format data
  method <- match.arg(method)
  obj_imputation_y <- list(data_sample = data_sample,
                           data_out_of_sample = data_out_of_sample)
  class(obj_imputation_y) <- method

  y_hat <- imputation_y(obj_imputation_y,
                        model_formula = model_formula,
                        type_model = type_model,
                        tune_RF = tune_RF,
                        xgboost_params = xgboost_params)


  y_full_imputed <- c(data_sample$y, unlist(y_hat$y_hat_out_of_sample))
  y_hat$y_full_imputed <- y_full_imputed
  y_hat$group_sample <- data_sample$group
  y_hat$group_out_of_sample <- data_out_of_sample$group

  output <- y_hat

  return(output)
}


#'
#' Internal generic function to impute values of y
#'
#' @inheritParams impute_y
#' @param obj_imputation_y Object to impute data
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


imputation_y <- function(...)
  UseMethod("imputation_y")


#'
#' @describeIn imputation_y Function to estimate imputation_y using EBLUP
#' @export
#'


imputation_y.EBLUP <- function(obj_imputation_y,
                               model_formula,
                               type_model,
                               ...) {

  # Add additional parameter for sample

  data_sample = obj_imputation_y$data_sample
  data_out_of_sample = obj_imputation_y$data_out_of_sample

  if (type_model == "gaussian") {

    outcome_fit <- lmer(model_formula, data = data_sample)
  } else {
    outcome_fit <- glmer(model_formula, data = data_sample, family = type_model)
  }

  y_hat_out_of_sample <- predict(outcome_fit,
                                 newdata = data_out_of_sample,
                                 allow.new.levels = TRUE)

  y_hat_sample <- predict(outcome_fit,
                          newdata = data_sample,
                          allow.new.levels = TRUE)

  output <- list(y_hat_out_of_sample = y_hat_out_of_sample,
                 y_hat_sample = y_hat_sample,
                 outcome_fit = outcome_fit)


  return(output)
}

#'
#' @describeIn imputation_y Function to estimate imputation_y using MQ
#' @export
#'


imputation_y.MQ <- function(obj_imputation_y,
                            model_formula,
                            type_model,
                            ...) {

  data_sample = obj_imputation_y$data_sample
  data_out_of_sample = obj_imputation_y$data_out_of_sample

  # Get the (predictor) variables
  vars <- attr(terms(model_formula), which = "term.labels")

  # Get the response
  response <- as.character(attr(terms(model_formula), which = "variables")[[2]]) # The response variable name

  # Get the predictors without group variable
  predictor <- paste(vars[-length(vars)], collapse = " + ")

  # Build a formula
  model_formulaMQ <- paste(response, " ~ ", predictor)

  outcome_fit = mquantreg(formula = model_formulaMQ,
                          data = data_sample,
                          q  = 0.5,
                          method = type_model)

  y_hat_out_of_sample <- unname(unlist(predict(outcome_fit,
                             newdata = data_out_of_sample,
                             regression_type = type_model)))

  y_hat_sample <- unname(unlist(predict(outcome_fit,
                            newdata = data_sample,
                            regression_type = type_model)))

  output <- list(y_hat_out_of_sample = y_hat_out_of_sample,
                 y_hat_sample = y_hat_sample,
                 outcome_fit = outcome_fit)

  return(output)

}


#'
#' @describeIn imputation_y Function to estimate imputation_y using RF
#' @export
#'



imputation_y.RF <- function(obj_imputation_y,
                            model_formula,
                            type_model,
                            tune_RF,
                            ...) {

  data_sample = obj_imputation_y$data_sample
  data_out_of_sample = obj_imputation_y$data_out_of_sample

  # Get the (predictor) variables

  formatted_data <- format_data(model_formula = model_formula,
                                data_sample = data_sample,
                                data_out_of_sample = data_out_of_sample)
  X = formatted_data$X
  X_newdata = formatted_data$X_newdata
  Y = formatted_data$Y
  clusters = formatted_data$clusters_sample

  if (tune_RF) {
    outcome_fit <- regression_forest(X, Y,
                                     clusters = clusters,
                                     tune.parameters = "all")
  } else {
    outcome_fit <- regression_forest(X, Y, clusters = clusters)
  }

  y_hat_out_of_sample <- unname(unlist(predict(outcome_fit, newdata = X_newdata)))

  y_hat_sample <- unname(unlist(predict(outcome_fit, newdata = X)))

  output <- list(y_hat_out_of_sample = y_hat_out_of_sample,
                 y_hat_sample = y_hat_sample,
                 outcome_fit = outcome_fit)

  return(output)

}

#'
#' @describeIn imputation_y Function to estimate imputation_y using RF
#' @export
#'


imputation_y.XGB <- function(obj_imputation_y,
                             model_formula,
                             type_model,
                             xgboost_params,
                             ...) {

  data_sample = obj_imputation_y$data_sample
  data_out_of_sample = obj_imputation_y$data_out_of_sample

  # Format data

  formatted_data <- format_data(model_formula = model_formula,
                                 data_sample = data_sample,
                                 data_out_of_sample = data_out_of_sample)
  X = as.matrix(formatted_data$X)
  X_newdata = as.matrix(formatted_data$X_newdata)
  Y = unlist(formatted_data$Y)

  CV_XGB = xgboost_params$CV_XGB
  nfolds = xgboost_params$nfolds
  nrounds = xgboost_params$nrounds

  #clusters  = as.numeric(data_sample$group)
  if (CV_XGB) {

    xgboost_cv <- xgb.cv(data = X, label = Y,
                         nfold = nfolds,
                         nrounds = nrounds,
                         verbose = FALSE)
    best_iter_xgb = which.min(xgboost_cv$evaluation_log$test_rmse_mean)

    outcome_fit <- xgboost(data = X, label = Y,
                           nrounds = best_iter_xgb,
                           verbose = FALSE)

    y_hat_out_of_sample <- predict(outcome_fit,
                                   newdata = X_newdata,
                                   iteration_range =  best_iter_xgb)

    y_hat_sample <- predict(outcome_fit,
                            newdata = X,
                            iteration_range =  best_iter_xgb)

  } else {
    outcome_fit <- xgboost(data = X, label = Y,
                           nrounds = nrounds,
                           verbose = FALSE)

    y_hat_out_of_sample <- predict(outcome_fit,
                                   newdata = X_newdata)

    y_hat_sample <- predict(outcome_fit,
                            newdata = X)
  }

  output <- list(y_hat_out_of_sample = y_hat_out_of_sample,
                 y_hat_sample = y_hat_sample,
                 outcome_fit = outcome_fit)


  return(output)
}
