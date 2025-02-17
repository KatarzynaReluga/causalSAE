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
                   params_OR,
                   params_bootstrap) {


  if (params_OR$method == "RF"|params_OR$method == "XGB") {
   formatted_data <- format_data(model_formula = params_OR$model_formula,
                                 data_sample = obj_hte$data_sample,
                                 data_out_of_sample  = obj_hte$data_out_of_sample)
    data_full <- rbind(formatted_data$X, formatted_data$X_newdata)
  } else {
    names_sample <- names(obj_hte$data_sample)
    names_out_of_sample <- names(obj_hte$data_out_of_sample)
    data_full <- rbind(obj_hte$data_sample[, !names_sample %in% "y"],
                       obj_hte$data_out_of_sample[, !names_sample %in% "y"])



  }

  # Check and format data
  data_sample <- obj_hte$data_sample

  # Controls -------------------------------------
  data_sample0 = data_sample[data_sample$A == 0, ]

  obj_fit_y0 <- list(data_sample = data_sample0)
  class(obj_fit_y0) <- params_OR$method
  mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                                  model_formula = params_OR$model_formula)


  OR0 <- fit_y(mutated_obj0,
                 type_model = params_OR$type_model,
                 tune_RF = params_OR$tune_RF,
                 clust_RF = params_OR$clust_RF,
                 xgboost_params = params_OR$xgboost_params,
                 params_bootstrap = params_bootstrap)

  mu0_y <-  unname(unlist(predict(object = OR0$outcome_fit,
                                  newdata = data_full,
                                  allow.new.levels = TRUE)))


  # Treated --------------------------------------
  data_sample1 = data_sample[data_sample$A == 1, ]

  obj_fit_y1 <- list(data_sample = data_sample1)
  class(obj_fit_y1) <- params_OR$method
  mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                                 model_formula = params_OR$model_formula)


  OR1 <- fit_y(mutated_obj1,
               type_model = params_OR$type_model,
               tune_RF = params_OR$tune_RF,
               clust_RF = params_OR$clust_RF,
               xgboost_params = params_OR$xgboost_params,
               params_bootstrap = params_bootstrap)

  mu1_y <-  unname(unlist(predict(object = OR1$outcome_fit,
                                  newdata = data_full,
                                  allow.new.levels = TRUE)))


  # Here we should sample y only once becasue we have only one sample
  data_OR <- data.frame(mu1_y = mu1_y,
                        mu0_y = mu0_y,
                        group = c(data_sample$group, obj_hte$data_out_of_sample$group))

  output <- list(data_OR = data_OR)

a = Sys.time()
  if (params_bootstrap$boot_var) {


    y_boot0 <- residual_bootstrap(y = data_sample0$y,
                                  y_hat = OR0$y_hat_sample,
                                  data_sample = data_sample0,
                                  n_boot  = params_bootstrap$n_boot,
                                  type_boot = params_bootstrap$type_boot,
                                  boot_seed = params_bootstrap$boot_seed,
                                  rand_clust = params_bootstrap$rand_clust)


    y_boot1 <- residual_bootstrap(y = data_sample1$y,
                                  y_hat = OR1$y_hat_sample,
                                  data_sample = data_sample1,
                                  n_boot  = params_bootstrap$n_boot,
                                  type_boot = params_bootstrap$type_boot,
                                  boot_seed = params_bootstrap$boot_seed,
                                  rand_clust = params_bootstrap$rand_clust)



    data_OR_boot <- list()
    obj_fit_y_b0 <- obj_fit_y0
    obj_fit_y_b1 <- obj_fit_y1

    for (i in 1 : params_bootstrap$n_boot) {
      #print(i)
      # Modulus operation
      if(i %% 10 == 0) {
        # Print on the screen some message
        cat(paste0("Bootstrap iteration: ", i, "\n"))
      }

      # Treated cases
      obj_fit_y_b1$data_sample$y <- y_boot1[[i]]
      mutated_obj_b1 <- mutate_obj_fit(obj_fit_y_b1, model_formula  = params_OR$model_formula)

      OR1b <- fit_y(mutated_obj_b1,
                      type_model = params_OR$type_model,
                      tune_RF = params_OR$tune_RF,
                      clust_RF = params_OR$clust_RF,
                      xgboost_params = params_OR$xgboost_params)
      mu1b_y <-  unname(unlist(predict(object = OR1b$outcome_fit,
                                      newdata = data_full,
                                      allow.new.levels = TRUE)))


      # Control cases
      obj_fit_y_b0$data_sample$y <- y_boot0[[i]]
      mutated_obj_b0 <- mutate_obj_fit(obj_fit_y_b0, model_formula  = params_OR$model_formula)

      OR0b <- fit_y(mutated_obj_b0,
                       type_model = params_OR$type_model,
                       tune_RF = params_OR$tune_RF,
                       clust_RF = params_OR$clust_RF,
                       xgboost_params = params_OR$xgboost_params)
      mu0b_y <-  unname(unlist(predict(object = OR0b$outcome_fit,
                                      newdata = data_full,
                                      allow.new.levels = TRUE)))

      # Save data
      data_OR_boot[[i]] <-  data.frame(mu1_y = mu1b_y,
                                       mu0_y = mu0b_y,
                                       group = c(data_sample$group,  obj_hte$data_out_of_sample$group))
    }
    output$data_OR_boot <- data_OR_boot
  }
b = Sys.time()


  return(output)
}
