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
#' @param params_bootstrap List with parameters to obtain bootstrap variance:
#' \itemize{
#'  \item boot_var = FALSE,
#'  \item n_boot = 500,
#'  \item boot_seed = 10,
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
#'                                          nfolds = 5))
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
                     params_bootstrap  = list(boot_var = FALSE,
                                              n_boot = 100,
                                              boot_seed = 10,
                                              type_boot = "br1",
#                                              method_scale = "mad",
#                                              method_center = "median",
                                              rand_clust = TRUE),
                     ...) {

  # Check and format data
  method <- match.arg(method)
  obj_fit_y <- list(data_sample = data_sample,
                    data_out_of_sample = data_out_of_sample)

  class(obj_fit_y) <- method
  mutated_obj <- mutate_obj_fit(obj_fit_y, model_formula)

  y_hat <- fit_y(mutated_obj,
                 type_model = type_model,
                 tune_RF = tune_RF,
                 xgboost_params = xgboost_params)

  output <- y_hat

  if (method == "EBLUP") {
    y_hat_out_of_sample <- predict(y_hat$outcome_fit,
                                   newdata = data_out_of_sample)
  } else if (method == "MQ") {
    y_hat_out_of_sample <- predict(y_hat$outcome_fit,
                                   newdata = data_out_of_sample,
                                   regression_type = type_model)
  } else if (method == "RF") {
    y_hat_out_of_sample <- unname(unlist(predict(y_hat$outcome_fit,
                                   newdata = mutated_obj$X_newdata,
                                   clusters = mutated_obj$clusters_out_of_sample)))
  } else {
    y_hat_out_of_sample <- unname(unlist(predict(y_hat$outcome_fit,
                                   newdata = mutated_obj$X_newdata)))
  }

  output$y_full_imputed <- c(data_sample$y, y_hat_out_of_sample)
  output$group_sample <- data_sample$group
  output$group_out_of_sample <- data_out_of_sample$group

  if (params_bootstrap$boot_var) {


    y_boot <- residual_bootstrap(y = data_sample$y,
                       y_hat = output$y_hat,
                       data_sample = data_sample,
                       n_boot  = params_bootstrap$n_boot,
                       type_boot = params_bootstrap$type_boot,
 #                      method_scale  = params_bootstrap$method_scale,
 #                      method_center  = params_bootstrap$method_center,
                       boot_seed = params_bootstrap$boot_seed,
                       rand_clust = params_bootstrap$rand_clust)



    y_full_b <- list()
    obj_fit_y_b <- obj_fit_y

    for (i in 1 : params_bootstrap$n_boot) {
      #print(i)
      # Modulus operation
      if(i %% 10 == 0) {
        # Print on the screen some message
        cat(paste0("Bootstrap iteration: ", i, "\n"))
      }
#      data_sample_b <- data_sample

      obj_fit_y_b$data_sample$y <- y_boot[[i]]
      mutated_obj_b <- mutate_obj_fit(obj_fit_y_b, model_formula)

      y_hatb <- fit_y(mutated_obj_b,
                      type_model = type_model,
                      tune_RF = tune_RF,
                      xgboost_params = xgboost_params)


      if (method == "EBLUP") {
        y_hat_out_of_sample_b <- unname(unlist(predict(y_hatb$outcome_fit,
                                       newdata = data_out_of_sample,
                                       allow.new.levels = TRUE)))
      } else if (method == "MQ") {
        y_hat_out_of_sample_b <- unname(unlist(predict(y_hatb$outcome_fit,
                                       newdata = data_out_of_sample,
                                       regression_type = type_model)))
      } else if (method == "RF") {
        y_hat_out_of_sample_b <- unname(unlist(predict(y_hatb$outcome_fit,
                                                     newdata = mutated_obj$X_newdata,
                                                     clusters = mutated_obj$clusters_out_of_sample)))
      } else {
        y_hat_out_of_sample_b <- unname(unlist(predict(y_hatb$outcome_fit,
                                                     newdata = mutated_obj$X_newdata)))
      }

      y_full_boot <- c(y_boot[[i]], y_hat_out_of_sample_b)

      y_full_b[[i]] <- y_full_boot
    }
    output$y_full_b <- y_full_b
  }

  return(output)
}
