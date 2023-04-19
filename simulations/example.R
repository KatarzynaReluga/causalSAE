m = 50
ni = rep(5, m)
Ni = rep(100, m)
N = sum(Ni)
n = sum(ni)

X <- generate_X(
  n = N,
  p = 1,
  covariance_norm = NULL,
  cov_type = "unif",
  seed = 1
 )

 X_outcome <- generate_X(
  n = N,
  p = 1,
  covariance_norm = NULL,
  cov_type = "lognorm",
  seed = 1
 )

 populations <- generate_pop(X, X_outcome,
 coeffs = get_default_coeffs(),
 errors_outcome = get_default_errors_outcome(),
 rand_eff_outcome = get_default_rand_eff_outcome(),
 rand_eff_p_score = get_default_rand_eff_p_score(),
 regression_type = "continuous",
 Ni_size  = 100,
 m = 50,
 no_sim = 1,
 seed = 10)

# samples <- generate_sample(populations, ni_size = 5,
#                            sample_part = "sampled",
#                            get_index = TRUE)


 subpopulation <- sample_subpopulations(populations,
 frac_nc = 0.1, frac_nt = 0.1, seed = 1)


data_sample <- populations[subpopulation, ]

#data_sample <- data.frame(samples[[1]]$samp_data)
#index_sample <- samples[[1]]$index_s
data_out_of_sample <- populations[-subpopulation, ]

model_formula = y ~ X1 + Xo1 + A + (1 + A||group)

hte_OR <- hte(type_hte = "OR",
             data_sample,
             data_out_of_sample,
             params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                              method = "EBLUP",
                              type_model = "gaussian"),
             params_bootstrap  = list(
               boot_var = FALSE,
               n_boot = 250,
               boot_seed = 10,
               type_boot = "br1",
               method_scale = "sd",
               method_center = "mean"
             ))

hte_ORb <- hte(type_hte = "OR",
              data_sample,
              data_out_of_sample,
              params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                               method = "EBLUP",
                               type_model = "gaussian"),
              params_bootstrap  = list(
                boot_var = TRUE,
                n_boot = 250,
                boot_seed = 10,
                type_boot = "br1",
                method_scale = "sd",
                method_center = "mean"
              ))


plot(hte_ORb$tau, type = "l", lwd = 2)
lines(sqrt(hte_ORb$var_tau) + hte_ORb$tau, col = 1, lty = 2)
lines(-sqrt(hte_ORb$var_tau) + hte_ORb$tau, col = 1, lty = 2)

lines(hte_NIPWb$tau, type = "l", col = 2, lwd = 2)
lines(sqrt(hte_NIPWb$var_tau) + hte_NIPWb$tau, col = 2, lty = 2)
lines(-sqrt(hte_NIPWb$var_tau) + hte_NIPWb$tau, col = 2, lty = 2)

lines(hte_AIPWb$tau, type = "l", col = 3, lwd = 2)
lines(sqrt(hte_AIPWb$var_tau) + hte_AIPWb$tau, col = 3, lty = 2)
lines(-sqrt(hte_AIPWb$var_tau) + hte_AIPWb$tau, col = 3, lty = 2)





hte_NIPWb <- hte(type_hte = "NIPW",
                 data_sample,
                 data_out_of_sample,
                 params_p_score = list(
                   model_formula = A ~ X1 + (1 | group),
                   method = "EBLUP",
                   tune_RF = FALSE,
                   xgboost_params = list(
                     CV_XGB = TRUE,
                     nfolds = 5,
                     nrounds = 50
                   )
                 ),
                 params_impute_y = list(
                   model_formula = y ~ X1 + Xo1 + A + (1 + A || group),
                   method = "EBLUP",
                   tune_RF = FALSE,
                   xgboost_params = list(
                     CV_XGB = TRUE,
                     nfolds = 5,
                     nrounds = 50
                   ),
                   type_model = "gaussian"
                 ),
                 params_bootstrap  = list(
                   boot_var = TRUE,
                   n_boot = 250,
                   boot_seed = 10,
                   type_boot = "br1",
                   method_scale = "sd",
                   method_center = "mean"
                 ))



hte_AIPW <- hte(type_hte = "AIPW",
              data_sample,
              data_out_of_sample,
              params_p_score = list(
                model_formula = A ~ X1 + (1 | group),
                method = "EBLUP",
                tune_RF = FALSE,
                xgboost_params = list(
                  CV_XGB = TRUE,
                  nfolds = 5,
                  nrounds = 50
                )
              ),
              params_OR = list(
                model_formula = y ~ X1 + Xo1 + (1 | group),
                method = "EBLUP",
                tune_RF = FALSE,
                xgboost_params = list(
                  CV_XGB = TRUE,
                  nfolds = 5,
                  nrounds = 50
                ),
                type_model = "gaussian"
              ),
              params_impute_y = list(
                model_formula = y ~ X1 + Xo1 + A + (1 + A || group),
                method = "MQ",
                tune_RF = FALSE,
                xgboost_params = list(
                  CV_XGB = TRUE,
                  nfolds = 5,
                  nrounds = 50
                ),
                type_model = "continuous"
              ),
              params_bootstrap  = list(
                boot_var = FALSE,
                n_boot = 250,
                boot_seed = 10,
                type_boot = "br1",
                method_scale = "sd",
                method_center = "mean"
              ))

hte_AIPWb <- hte(type_hte = "AIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(
                  model_formula = A ~ X1 + (1 | group),
                  method = "EBLUP",
                  tune_RF = FALSE,
                  xgboost_params = list(
                    CV_XGB = TRUE,
                    nfolds = 5,
                    nrounds = 50
                  )
                ),
                params_OR = list(
                  model_formula = y ~ X1 + Xo1 + (1 | group),
                  method = "EBLUP",
                  tune_RF = FALSE,
                  xgboost_params = list(
                    CV_XGB = TRUE,
                    nfolds = 5,
                    nrounds = 50
                  ),
                  type_model = "gaussian"
                ),
                params_impute_y = list(
                  model_formula = y ~ X1 + Xo1 + A + (1 + A || group),
                  method = "MQ",
                  tune_RF = FALSE,
                  xgboost_params = list(
                    CV_XGB = TRUE,
                    nfolds = 5,
                    nrounds = 50
                  ),
                  type_model = "continuous"
                ),
                params_bootstrap  = list(
                  boot_var = TRUE,
                  n_boot = 250,
                  boot_seed = 10,
                  type_boot = "br1",
                  method_scale = "sd",
                  method_center = "mean"
                ))




hte_NIPW <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(
                  model_formula = A ~ X1 + (1 | group),
                  method = "EBLUP",
                  tune_RF = FALSE,
                  xgboost_params = list(
                    CV_XGB = TRUE,
                    nfolds = 5,
                    nrounds = 50
                  )
                ),
                params_impute_y = list(
                  model_formula = y ~ X1 + Xo1 + A + (1 + A || group),
                  method = "EBLUP",
                  tune_RF = FALSE,
                  xgboost_params = list(
                    CV_XGB = TRUE,
                    nfolds = 5,
                    nrounds = 50
                  ),
                  type_model = "gaussian"
                ),
                params_bootstrap  = list(
                  boot_var = FALSE,
                  n_boot = 250,
                  boot_seed = 10,
                  type_boot = "br1",
                  method_scale = "sd",
                  method_center = "mean"
                ))
