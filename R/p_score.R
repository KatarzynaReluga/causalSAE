#' Propensity score
#'
#' Fit propensity score model
#'
#' @inheritParams impute_y
#' @param obj_p_score Object defining propensity score fit.
#' @param tune_RF Tune parameters for fitting random forest? Default: \code{tune_RF = FALSE}.
#' @param xgboost_params List of parameters to obtain predictions using gradient boosting:
#'  \itemize{
#'  \item CV_XGB - logical variable, use cross-validation for gradient boosting? Default: \code{CV_XGB = TRUE}.
#'  \item nfolds - number of folds in cross-validation. Default: \code{nfolds = 5}.
#'  \item nrounds - the max number of iterations. Default: \code{nrounds = 50}.
#'  }
#'
#' @importFrom lme4 lmer glmer
#' @importFrom stats terms binomial
#' @importFrom grf regression_forest
#' @importFrom mquantreg mquantreg
#'
#' @export
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
#'                             coeffs = get_default_coeffs(),
#'                             errors_outcome = get_default_errors_outcome(),
#'                             rand_eff_outcome = get_default_rand_eff_outcome(),
#'                             rand_eff_p_score = get_default_rand_eff_p_score(),
#'                             regression_type = "continuous",
#'                             Ni_size  = 100,
#'                             m = 50,
#'                             no_sim = 1,
#'                             seed = 10)
#'
#'
#' # EBLUP
#'
#' obj_p_score_EBLUP <- list(data_p_score = populations)
#' class(obj_p_score_EBLUP) <- "EBLUP"
#'
#' ps_hat_EBLUP <-  p_score(obj_p_score = obj_p_score_EBLUP,
#'                          model_formula = A ~ X1 + (1|group))
#'
#' # MQ
#'
#' obj_p_score_MQ <- list(data_p_score = populations)
#' class(obj_p_score_MQ) <- "MQ"
#'
#' ps_hat_MQ <-  p_score(obj_p_score = obj_p_score_MQ,
#'                       model_formula = A ~ X1 + (1|group))
#'
#' # RF
#'
#' obj_p_score_RF <- list(data_p_score = populations)
#' class(obj_p_score_RF) <- "RF"
#'
#' ps_hat_RF <-  p_score(obj_p_score = obj_p_score_RF,
#'                      model_formula = A ~ X1 + (1|group),
#'                      tune_RF = FALSE)
#'
#' # XGB
#'
#' obj_p_score_XGB <- list(data_p_score = populations)
#' class(obj_p_score_XGB) <- "XGB"
#'
#' ps_hat_XGB <-  p_score(obj_p_score = obj_p_score_XGB,
#'                        model_formula = A ~ X1 + (1|group),
#'                        xgboost_params = list(CV_XGB = FALSE,
#'                                              nfolds = 5,
#'                                              nrounds = 50))
#'
#'
#'

p_score <- function(...)
  UseMethod("p_score")

#'
#' @describeIn p_score Function to estimate p_score using EBLUP
#' @export
#'


p_score.EBLUP <- function(obj_p_score,
                          model_formula,
                          ...) {

  data_p_score <- obj_p_score$data_p_score

  ps_fit <- glmer(model_formula, data = data.frame(data_p_score),
                  family = binomial(link = "logit"))

  ps_hat <- as.vector(predict(ps_fit,
                              newdata = data_p_score,
                              type = "response",
                              allow.new.levels = TRUE))

  return(ps_hat)

}


#'
#' @describeIn p_score Function to estimate p_score using MQ
#' @export
#'

p_score.MQ <- function(obj_p_score,
                       model_formula,
                       ...) {

  data_p_score <- obj_p_score$data_p_score

  # Get the (predictor) variables
  vars <- attr(terms(model_formula), which = "term.labels")

  # Get the response
  response <- as.character(attr(terms(model_formula), which = "variables")[[2]]) # The response variable name

  # Get the predictors without group variable
  predictor <- paste(vars[-length(vars)], collapse = " + ")

  # Build a formula

  ps_fit <- mquantreg(formula = model_formulaMQ,
                      data = data_p_score,
                      q  = 0.5, method = "binom")

  ps_hat <- unname(unlist(predict(ps_fit, newdata = data_p_score,
                                  regression_type = "binary")))

  return(ps_hat)

}


#'
#' @describeIn p_score Function to estimate p_score using RF
#' @export
#'

p_score.RF <- function(obj_p_score,
                       model_formula,
                       tune_RF = FALSE, ...) {

  data_p_score <- obj_p_score$data_p_score

  # Formatted data
  formatted_data <- format_data(model_formula = model_formula,
                                data_sample = obj_p_score$data_p_score,
                                data_out_of_sample = NULL)

  X = formatted_data$X
  Y = formatted_data$Y

  clusters = formatted_data$cluster_sample

  if (tune_RF) {
    ps_fit <- regression_forest(X, Y,
                                clusters = clusters,
                                tune.parameters = "all")
  } else {
    ps_fit <- regression_forest(X, Y, clusters = clusters)
  }

  ps_hat <- unname(unlist(predict(ps_fit, newdata = X)))

  return(ps_hat)

}

#'
#' @describeIn p_score Function to estimate p_score using RF
#' @export
#'

p_score.XGB <- function(obj_p_score,
                        model_formula,
                        xgboost_params,
                        ...) {

  data_p_score <- obj_p_score$data_p_score

  # Formatted data
  formatted_data <- format_data(model_formula = model_formula,
                                data_sample = obj_p_score$data_p_score,
                                data_out_of_sample = NULL)

  X = formatted_data$X
  Y = formatted_data$Y

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

    ps_fit <- xgboost(data = X, label = Y,
                      nrounds = best_iter_xgb,
                      verbose = FALSE)

    ps_hat <- predict(ps_fit, newdata = X,
                      iteration_range =  best_iter_xgb)

  } else {
    ps_fit <- xgboost(data = X, label = Y,
                      nrounds = nrounds,
                      verbose = FALSE)

    ps_hat <- predict(ps_fit, newdata = X)
  }

  return(ps_hat)

}
