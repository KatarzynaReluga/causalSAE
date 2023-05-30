#' Geeneric function to fit values of y
#'
#' @inheritParams impute_y
#' @param obj_fit_y Object to impute data
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
#' @export
#'
#'
#'
#'


fit_y <- function(...)
  UseMethod("fit_y")


#'
#' @describeIn fit_y Function to estimate fit_y using EBLUP
#' @export
#'


fit_y.EBLUP <- function(obj_fit_y,
                        #                               model_formula,
                        type_model,
                        params_bootstrap,
                        ...) {

  data_sample = obj_fit_y$data_sample
  model_formula  = obj_fit_y$model_formula

  if (type_model == "gaussian") {
    outcome_fit <- lmer(formula = model_formula, data = data_sample)
  } else {
    outcome_fit <- glmer(formula = model_formula, data = data_sample, family = type_model)
  }

#  y_hat_sample <- unname(predict(outcome_fit,
#                                 newdata = data_sample))

  y_hat_sample <- unname(predict(outcome_fit))



  output <- list(y_hat_sample = y_hat_sample,
                 outcome_fit = outcome_fit)



  return(output)
}

#'
#' @describeIn fit_y Function to estimate fit_y using MQ
#' @export
#'


fit_y.MQ <- function(obj_fit_y,
                     type_model,
                     ...) {

  data_sample = obj_fit_y$data_sample
  model_formulaMQ = obj_fit_y$model_formula

  outcome_fit = mquantreg(formula = model_formulaMQ,
                          data = data_sample,
                          q  = 0.5,
                          method = type_model)

  y_hat_sample <- predict(outcome_fit,
                          newdata = data_sample,
                          regression_type = type_model)

  output <- list(y_hat_sample = y_hat_sample,
                 outcome_fit = outcome_fit)

  return(output)

}


#'
#' @describeIn fit_y Function to estimate fit_y using RF
#' @export
#'



fit_y.RF <- function(obj_fit_y,
                     type_model,
                     tune_RF,
                     ...) {
  X = obj_fit_y$X
  Y = obj_fit_y$Y
  clusters = obj_fit_y$clusters_sample
   if (tune_RF) {
       test  = "try-error"
    while(test == "try-error") {
      outcome_fit <-  try(regression_forest(X, Y,
                                            clusters = clusters,
                                            tune.parameters = "all"), silent = TRUE)
      test = class(outcome_fit)[1]
    }
  } else {
    outcome_fit <- regression_forest(X, Y, clusters = clusters)
  }
  y_hat_sample <- c(outcome_fit$predictions)


  output <- list(y_hat_sample = y_hat_sample,
                 outcome_fit = outcome_fit)

  return(output)

}

#'
#' @describeIn fit_y Function to estimate fit_y using RF
#' @export
#'


fit_y.XGB <- function(obj_fit_y,
                      #                            model_formula,
                      type_model,
                      xgboost_params,
                      ...) {

  X = as.matrix(obj_fit_y$X)
  Y = unlist(obj_fit_y$Y)

  CV_XGB = xgboost_params$CV_XGB
  nfolds = xgboost_params$nfolds
  nrounds = xgboost_params$nrounds

  if (CV_XGB) {

    xgboost_cv <- xgb.cv(data = X, label = Y,
                         nfold = nfolds,
                         nrounds = nrounds,
                         verbose = FALSE)
    best_iter_xgb = which.min(xgboost_cv$evaluation_log$test_rmse_mean)

    outcome_fit <- xgboost(data = X, label = Y,
                           nrounds = best_iter_xgb,
                           verbose = FALSE)
    y_hat_sample <- predict(outcome_fit,
                            newdata = X,
                            iteration_range =  best_iter_xgb)

  } else {
    outcome_fit <- xgboost(data = X, label = Y,
                           nrounds = nrounds,
                           verbose = FALSE)

    y_hat_sample <- predict(outcome_fit,
                            newdata = X)
  }

  output <- list(y_hat_sample = y_hat_sample,
                 outcome_fit = outcome_fit)


  return(output)
}

# #'
# #' @describeIn fit_y Function to estimate fit_y using RF
# #' @export
# #'



#fit_y.SL <- function(obj_fit_y,
#                     type_model,
#                     tune_RF,
#                     ...) {
#  X = obj_fit_y$X
#  Y = obj_fit_y$Y

#  outcome_fit = SuperLearner(Y = obj_fit_y$Y, X = obj_fit_y$X, family = gaussian(),
#                    SL.library = c("SL.mq","SL.grf", "SL.xgboost", "SL.glmm"))

#  y_hat_sample <-
#  clusters = obj_fit_y$clusters_sample
#  if (tune_RF) {
#    test  = "try-error"
#    while(test == "try-error") {
#      outcome_fit <-  try(regression_forest(X, Y,
#                                            clusters = clusters,
#                                            tune.parameters = "all"), silent = TRUE)
#      test = class(outcome_fit)[1]
#    }
#  } else {
#    outcome_fit <- regression_forest(X, Y, clusters = clusters)
#  }
#  y_hat_sample <- c(outcome_fit$predictions)


#  output <- list(y_hat_sample = y_hat_sample,
#                 outcome_fit = outcome_fit)

#  return(output)

#}

#' Generic function to format data
#'
#' @inheritParams impute_y
#' @param obj_fit_y Object to fit outcome y
#'
#' @export
#'


mutate_obj_fit <- function(...)
  UseMethod("mutate_obj_fit")


#'
#' @describeIn mutate_obj_fit Function to estimate fit_y using EBLUP
#' @export
#'

mutate_obj_fit.EBLUP <- function(obj_fit_y, model_formula, ...) {

  output <- obj_fit_y
  output$model_formula <- model_formula
  class(output) <- class(obj_fit_y)

  return(output)

}

#'
#' @describeIn mutate_obj_fit Function to estimate fit_y using MQ
#' @export
#'

mutate_obj_fit.MQ <- function(obj_fit_y, model_formula, ...) {


  # Get the (predictor) variables
  vars <- attr(terms(model_formula), which = "term.labels")

  # Get the response
  response <- as.character(attr(terms(model_formula), which = "variables")[[2]])

  # Get the predictors without group variable
  predictor <- paste(vars[-length(vars)], collapse = " + ")

  # Build a formula
  model_formulaMQ <- paste(response, " ~ ", predictor)

  output <- obj_fit_y
  output$model_formula <- model_formulaMQ

  class(output) <- class(obj_fit_y)

  return(output)

}


#'
#' @describeIn mutate_obj_fit Function to estimate fit_y using RF
#' @export
#'

mutate_obj_fit.RF <- function(obj_fit_y, model_formula, ...) {


  # Get the (predictor) variables

  formatted_data <- format_data(model_formula = model_formula,
                                data_sample = obj_fit_y$data_sample,
                                data_out_of_sample = obj_fit_y$data_out_of_sample)

  output <- formatted_data
  #  output$model_formula <- model_formula

  class(output) <- class(obj_fit_y)

  return(output)

}


#'
#' @describeIn mutate_obj_fit Function to estimate fit_y using XGB
#' @export
#'

mutate_obj_fit.XGB <- function(obj_fit_y, model_formula, ...) {


  # Get the (predictor) variables

  formatted_data <- format_data(model_formula = model_formula,
                                data_sample = obj_fit_y$data_sample,
                                data_out_of_sample = obj_fit_y$data_out_of_sample)

  output <- formatted_data
  #  output$model_formula <- model_formula

  class(output) <- class(obj_fit_y)

  return(output)

}

# '
# ' @describeIn mutate_obj_fit Function to estimate fit_y using SuperLearner
# ' @export
# '

#mutate_obj_fit.SL <- function(obj_fit_y, model_formula, ...) {


  # Get the (predictor) variables

#  formatted_data <- format_data(model_formula = model_formula,
#                                data_sample = obj_fit_y$data_sample,
#                                data_out_of_sample = obj_fit_y$data_out_of_sample)

#  output <- formatted_data
  #  output$model_formula <- model_formula

#  class(output) <- class(obj_fit_y)

#  return(output)

#}
