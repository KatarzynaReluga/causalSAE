#' Fit outcome regression models
#'
#' Function to obtain outcome regession models
#'
#' @inheritParams hte
#' @param obj_hte Object to estimate hte
#'
#' @return Data frame with outcome regression models.
#'
#' @export
#'
#'

fit_OR <- function(obj_hte,
                   params_OR) {

  # Full data --------------------------
 if (params_OR$method == "RF"|params_OR$method == "XGB") {
  formatted_data <- format_data(model_formula = params_OR$model_formula,
                                  data_sample = rbind(obj_hte$data_sample, obj_hte$data_out_of_sample))
  data_full <- formatted_data$X
 } else {
   data_full <- rbind(obj_hte$data_sample, obj_hte$data_out_of_sample)
 }
  data_sample  = obj_hte$data_sample
  data_out_of_sample = obj_hte$data_out_of_sample

  # Controls -------------------------------------
  data_sample0 = data_sample[data_sample$A == 0, ]

  OR0 <- impute_y(model_formula = params_OR$model_formula,
                  data_sample = data_sample0,
                  data_out_of_sample = data_out_of_sample,
                  method = params_OR$method,
                  type_model = params_OR$type_model,
                  tune_RF = params_OR$tune_RF,
                  xgboost_params = params_OR$xgboost_params)


  mu0_y <-  unname(unlist(predict(object = OR0$outcome_fit, newdata = data_full,
                   allow.new.levels = TRUE)))

  # Treated --------------------------------------
  data_sample1 = data_sample[data_sample$A == 1, ]

  OR1 <- impute_y(model_formula = params_OR$model_formula,
                  data_sample = data_sample1,
                  data_out_of_sample = data_out_of_sample,
                  method = params_OR$method,
                  type_model = params_OR$type_model,
                  tune_RF = params_OR$tune_RF,
                  xgboost_params = params_OR$xgboost_params)

  mu1_y <- unname(unlist(predict(OR1$outcome_fit, newdata = data_full,
                   allow.new.levels = TRUE)))
  data_OR <- data.frame(mu1_y = mu1_y,
                        mu0_y = mu0_y,
                        group = c(data_sample$group, data_out_of_sample$group))
#  fit0 <- OR0$outcome_fit
#  fit1 <- OR1$outcome_fit

#  output <-list(data_OR = data_OR,
#                fit0 = fit0,
#                fit1 = fit1)

  return(data_OR)
}
