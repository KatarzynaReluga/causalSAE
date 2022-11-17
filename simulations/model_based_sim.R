# Model based simulations

m = 50
ni = rep(5, m)
Ni = rep(50, m)
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
                            Ni_size = 50,
                            m = 50,
                            no_sim = 1,
                            seed = 10)

NoSim <- 100

# OR ----------------------------
# HTE
EBLUP_OR <- matrix(NA, NoSim, m)
MQ_OR <- matrix(NA, NoSim, m)
RF_OR <- matrix(NA, NoSim, m)
XGB_OR <- matrix(NA, NoSim, m)
# var
EBLUP_OR_var <- matrix(NA, NoSim, m)
MQ_OR_var <- matrix(NA, NoSim, m)
RF_OR_var <- matrix(NA, NoSim, m)
XGB_OR_var <- matrix(NA, NoSim, m)


# NIPW ---------------------------------------------
EBLUP_NIPW <- matrix(NA, NoSim, m)
MQ_NIPW <- matrix(NA, NoSim, m)
RF_NIPW <- matrix(NA, NoSim, m)
XGB_NIPW <- matrix(NA, NoSim, m)

EBLUP_RF_NIPW <- matrix(NA, NoSim, m)
EBLUP_XGB_NIPW <- matrix(NA, NoSim, m)
MQ_RF_NIPW <- matrix(NA, NoSim, m)
MQ_XGB_NIPW <- matrix(NA, NoSim, m)
# var
EBLUP_NIPW_var <- matrix(NA, NoSim, m)
EBLUP_NIPW_var <- matrix(NA, NoSim, m)
EBLUP_NIPW_var <- matrix(NA, NoSim, m)
EBLUP_NIPW_var <- matrix(NA, NoSim, m)

EBLUP_RF_NIPW_var <- matrix(NA, NoSim, m)
EBLUP_XGB_NIPW_var <- matrix(NA, NoSim, m)
MQ_RF_NIPW_var <- matrix(NA, NoSim, m)
MQ_XGB_NIPW_var <- matrix(NA, NoSim, m)

# AIPW --------------------------------------------------
EBLUP_RF_AIPW <- matrix(NA, NoSim, m)
MQ_RF_AIPW <- matrix(NA, NoSim, m)
RF_EBLUP_AIPW <- matrix(NA, NoSim, m)
XGB_EBLUP_AIPW <- matrix(NA, NoSim, m)

# Direct estimators -----------------------------------------
Dir_EBLUP_tau <- matrix(NA, NoSim, m)
Dir_MQ_tau <- matrix(NA, NoSim, m)
Dir_RF_tau <- matrix(NA, NoSim, m)
Dir_XGB_tau <- matrix(NA, NoSim, m)




# Simple checks of the code ------------------------------------------------------------------
for (i in 1:NoSim) {

  print(i)
  set.seed(i * 2022)

  samples <- generate_sample(populations, ni_size = 5,
                             sample_part = "sampled",
                             get_index = TRUE)

  data_sample <- data.frame(samples[[1]]$samp_data)
  index_sample <- samples[[1]]$index_s
  data_out_of_sample <- populations[-index_sample, ]


  # OR -------------------------------------------------------------------
  # EBLUP OR
  EBLUP_ORf <- hte(type_hte = "OR",
                data_sample,
                data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  tune_RF = FALSE,
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50),
                                                        type_model = "gaussian"),
                 params_bootstrap  = list(boot_var = TRUE,
                                          n_boot = 200,
                                          boot_seed = 10))
  EBLUP_OR[i, ] <- EBLUP_ORf$tau
  EBLUP_OR_var[i, ] <- EBLUP_ORf$var_tau

  # MQ OR
  MQ_ORf <- hte(type_hte = "OR",
                   data_sample,
                   data_out_of_sample,
                   params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                    method = "MQ",
                                    tune_RF = FALSE,
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50),
                                    type_model = "continuous"),
                   params_bootstrap  = list(boot_var = TRUE,
                                            n_boot = 200,
                                            boot_seed = 10))
  MQ_OR[i, ] <- MQ_ORf$tau
  MQ_OR_var[i, ] <- MQ_ORf$var_tau

  # RF OR

  RF_ORf <- hte(type_hte = "OR",
                data_sample,
                data_out_of_sample,
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                 method = "RF",
                                 tune_RF = FALSE,
                                 xgboost_params = list(CV_XGB = TRUE,
                                                       nfolds = 5,
                                                       nrounds = 50),
                                 type_model = "continuous"),
                params_bootstrap  = list(boot_var = TRUE,
                                         n_boot = 200,
                                         boot_seed = 10))
  RF_OR[i, ] <- RF_ORf$tau
  RF_OR_var[i, ] <- RF_ORf$var_tau

  # EBLUP XGB

  XGB_ORf <- hte(type_hte = "OR",
                data_sample,
                data_out_of_sample,
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                 method = "XGB",
                                 tune_RF = FALSE,
                                 xgboost_params = list(CV_XGB = FALSE,
                                                       nfolds = 5,
                                                       nrounds = 50),
                                 type_model = "continuous"),
                params_bootstrap  = list(boot_var = TRUE,
                                         n_boot = 200,
                                         boot_seed = 10))
  XGB_OR[i, ] <- XGB_ORf$tau
  XGB_OR_var[i, ] <- XGB_ORf$var_tau

  # NIPW -----------------------------------------------------------------------
  # EBLUP OR
  a = Sys.time()
  EBLUP_NIPWf <- hte(type_hte = "NIPW",
                   data_sample,
                   data_out_of_sample,
                   params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                         method =  "EBLUP",
                                         tune_RF = FALSE,
                                         xgboost_params = list(CV_XGB = FALSE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                          method = "EBLUP",
                                          tune_RF = FALSE,
                                          xgboost_params = list(CV_XGB = FALSE,
                                                                nfolds = 5,
                                                                nrounds = 50),
                                          type_model = "gaussian"),
                   params_bootstrap  = list(boot_var = TRUE,
                                            n_boot = 500,
                                            boot_seed = 10))
  b = Sys.time()
  EBLUP_NIPW[i, ] <- EBLUP_NIPWf$tau
  EBLUP_NIPW_var[i, ] <- EBLUP_NIPWf$var_tau

  # MQ NIPWf
  a = Sys.time()
  MQ_NIPWf <- hte(type_hte = "NIPW",
                  data_sample,
                  data_out_of_sample,
                  params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                        method =  "MQ",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)),
                  params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                          method = "MQ",
                                          tune_RF = FALSE,
                                          xgboost_params = list(CV_XGB = FALSE,
                                                                nfolds = 5,
                                                               nrounds = 50),
                                              type_model = "continuous"),
                  params_bootstrap  = list(boot_var = TRUE,
                                                n_boot = 500,
                                                boot_seed = 10))
  b = Sys.time()
  b-a
  MQ_NIPW[i, ] <- MQ_NIPWf$tau
  MQ_NIPW_var[i, ] <- MQ_NIPWf$var_tau

  # RF OR

  a = Sys.time()
  RF_NIPWf <- hte(type_hte = "NIPW",
                  data_sample,
                  data_out_of_sample,
                  params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                        method =  "RF",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)),
                  params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                         method = "RF",
                                         tune_RF = FALSE,
                                         xgboost_params = list(CV_XGB = FALSE,
                                                               nfolds = 5,
                                                               nrounds = 50),
                                         type_model = "continuous"),
                  params_bootstrap  = list(boot_var = TRUE,
                                           n_boot = 200,
                                           boot_seed = 10))
  b = Sys.time()
  RF_NIPW[i, ] <- RF_NIPWf$tau
  RF_NIPW_var[i, ] <- RF_NIPWf$var_tau

  # EBLUP XGB
  a = Sys.time()
  XGB_NIPWf <- hte(type_hte = "NIPW",
                  data_sample,
                  data_out_of_sample,
                  params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                        method =  "XGB",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)),
                  params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                         method = "XGB",
                                         tune_RF = FALSE,
                                         xgboost_params = list(CV_XGB = FALSE,
                                                               nfolds = 5,
                                                               nrounds = 50),
                                         type_model = "continuous"),
                  params_bootstrap  = list(boot_var = TRUE,
                                           n_boot = 200,
                                           boot_seed = 10))
  b = Sys.time()
  XGB_NIPW[i, ] <- XGB_NIPWf$tau
  XGB_NIPW_var[i, ] <- XGN_NIPWf$var_tau


  # EBLUP RF
  a = Sys.time()
  EBLUP_RF_NIPWf <- hte(type_hte = "NIPW",
                     data_sample,
                     data_out_of_sample,
                     params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                           method =  "EBLUP",
                                           tune_RF = FALSE,
                                           xgboost_params = list(CV_XGB = FALSE,
                                                                 nfolds = 5,
                                                                 nrounds = 50)),
                     params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                            method = "RF",
                                            tune_RF = FALSE,
                                            xgboost_params = list(CV_XGB = FALSE,
                                                                  nfolds = 5,
                                                                  nrounds = 50),
                                            type_model = "gaussian"),
                     params_bootstrap  = list(boot_var = TRUE,
                                              n_boot = 500,
                                              boot_seed = 10))
  b = Sys.time()
  EBLUP_RF_NIPW[i, ] <- EBLUP_RF_NIPWf$tau
  EBLUP_RF_NIPW_var[i, ] <- EBLUP_NIPWf$var_tau

  # EBLUP XGB
  a = Sys.time()
  EBLUP_XGB_NIPWf <- hte(type_hte = "NIPW",
                        data_sample,
                        data_out_of_sample,
                        params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                              method =  "EBLUP",
                                              tune_RF = FALSE,
                                              xgboost_params = list(CV_XGB = FALSE,
                                                                    nfolds = 5,
                                                                    nrounds = 50)),
                        params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                               method = "XGB",
                                               tune_RF = FALSE,
                                               xgboost_params = list(CV_XGB = FALSE,
                                                                     nfolds = 5,
                                                                     nrounds = 50),
                                               type_model = "gaussian"),
                        params_bootstrap  = list(boot_var = TRUE,
                                                 n_boot = 500,
                                                 boot_seed = 10))
  b = Sys.time()
  EBLUP_XGB_NIPW[i, ] <- EBLUP_XGB_NIPWf$tau
  EBLUP_XGB_NIPW_var[i, ] <- EBLUP_XGB_NIPWf$var_tau

  # MQ
  a = Sys.time()
  MQ_RF_NIPWf <- hte(type_hte = "NIPW",
                  data_sample,
                  data_out_of_sample,
                  params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                        method =  "MQ",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)),
                  params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                         method = "RF",
                                         tune_RF = FALSE,
                                         xgboost_params = list(CV_XGB = FALSE,
                                                               nfolds = 5,
                                                               nrounds = 50),
                                         type_model = "continuous"),
                  params_bootstrap  = list(boot_var = TRUE,
                                           n_boot = 500,
                                           boot_seed = 10))
  b = Sys.time()
  b-a
  MQ_RF_NIPW[i, ] <- MQ_RF_NIPWf$tau
  MQ_NIPW_RF_var[i, ] <- MQ_RF_NIPWf$var_tau

  MQ_RF_NIPWf <- hte(type_hte = "NIPW",
                     data_sample,
                     data_out_of_sample,
                     params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                           method =  "MQ",
                                           tune_RF = FALSE,
                                           xgboost_params = list(CV_XGB = FALSE,
                                                                 nfolds = 5,
                                                                 nrounds = 50)),
                     params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                            method = "XGB",
                                            tune_RF = FALSE,
                                            xgboost_params = list(CV_XGB = FALSE,
                                                                  nfolds = 5,
                                                                  nrounds = 50),
                                            type_model = "continuous"),
                     params_bootstrap  = list(boot_var = TRUE,
                                              n_boot = 500,
                                              boot_seed = 10))
  b = Sys.time()
  b-a
  MQ_XGB_NIPW[i, ] <- MQ_XGB_NIPWf$tau
  MQ_NIPW_XGB_var[i, ] <- MQ_XGB_NIPWf$var_tau

  # AIPW -----------------------------------------------------
  a = Sys.time()
  EBLUP_RF_AIPWf <- hte(type_hte = "AIPW",
                  data_sample,
                  data_out_of_sample,
                  params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                   method = "EBLUP",
                                   tune_RF = FALSE,
                                   xgboost_params = list(CV_XGB = TRUE,
                                                         nfolds = 5,
                                                         nrounds = 50),
                                                         type_model = "gaussian"),
                  params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                        method =  "EBLUP",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)),
                  params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                         method = "RF",
                                         tune_RF = FALSE,
                                         xgboost_params = list(CV_XGB = FALSE,
                                                               nfolds = 5,
                                                               nrounds = 50),
                                         type_model = "continuous"),
                  params_bootstrap  = list(boot_var = TRUE,
                                           n_boot = 200,
                                           boot_seed = 10))
  b = Sys.time()
  EBLUP_RF_AIPW[i, ] <- EBLUP_RF_AIPWf$tau
  EBLUP_RF_AIPW_var[i, ] <- EBLUP_RF_AIPWf$var_tau

  # MQ EBLUP --------------------------------------------------------------------------------------
  a = Sys.time()
  MQ_RF_AIPWf <- hte(type_hte = "AIPW",
                        data_sample,
                        data_out_of_sample,
                        params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                         method = "MQ",
                                         tune_RF = FALSE,
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50),
                                         type_model = "continuous"),
                        params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                              method =  "MQ",
                                              tune_RF = FALSE,
                                              xgboost_params = list(CV_XGB = FALSE,
                                                                    nfolds = 5,
                                                                    nrounds = 50)),
                        params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                               method = "RF",
                                               tune_RF = FALSE,
                                               xgboost_params = list(CV_XGB = FALSE,
                                                                     nfolds = 5,
                                                                     nrounds = 50),
                                               type_model = "continuous"),
                        params_bootstrap  = list(boot_var = TRUE,
                                                 n_boot = 200,
                                                 boot_seed = 10))
  b = Sys.time()
  MQ_RF_AIPW[i, ] <- EBLUP_RF_AIPWf$tau
  MQ_RF_AIPW_var[i, ] <- EBLUP_RF_AIPWf$var_tau
  # RF EBLUP --------------------------------------------------------------------------------------
  a = Sys.time()
  RF_EBLUP_AIPWf <- hte(type_hte = "AIPW",
                     data_sample,
                     data_out_of_sample,
                     params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                      method = "RF",
                                      tune_RF = FALSE,
                                      xgboost_params = list(CV_XGB = TRUE,
                                                            nfolds = 5,
                                                            nrounds = 50),
                                      type_model = "continuous"),
                     params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                           method =  "RF",
                                           tune_RF = FALSE,
                                           xgboost_params = list(CV_XGB = FALSE,
                                                                 nfolds = 5,
                                                                 nrounds = 50)),
                     params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                            method = "RF",
                                            tune_RF = FALSE,
                                            xgboost_params = list(CV_XGB = FALSE,
                                                                  nfolds = 5,
                                                                  nrounds = 50),
                                            type_model = "continuous"),
                     params_bootstrap  = list(boot_var = TRUE,
                                              n_boot = 200,
                                              boot_seed = 10))
  b = Sys.time()
  EBLUP_RF_AIPW[i, ] <- EBLUP_RF_AIPWf$tau
  EBLUP_RF_AIPW_var[i, ] <- EBLUP_RF_AIPWf$var_tau

  # Direct estimators
  # Direct EBLUP

  Dir_EBLUP_tau[i, ] <- (calculate_tau(list(data_sample),
                                       type_tau = "H"))[[1]]$tau

  # Direct MQ
  samp_MQ <- data_pop_MQ[samp_index, ]
  Dir_MQ_tau[i, ] <- (calculate_tau(list(samp_MQ),
                                    type_tau = "H"))[[1]]$tau

  # Direct RF
  samp_RF <- data_pop_RF[samp_index, ]
  Dir_RF_tau[i, ] <- (calculate_tau(list(samp_RF),
                                    type_tau = "H"))[[1]]$tau

  # Direct XGB
  samp_XGB <- data_pop_XGB[samp_index, ]
  Dir_XGB_tau[i, ] <- (calculate_tau(list(samp_XGB),
                                     type_tau = "H"))[[1]]$tau


}

 #' hte_NIPW <- hte(type_hte = "NIPW",
#'               data_sample,
#'               data_out_of_sample,
#'               params_p_score = list(model_formula = A ~ X1 + (1|group),
#'                           method = "EBLUP",
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
#'               params_bootstrap  = list(boot_var = FALSE,
#'                                        n_boot = 500,
#'                                        boot_seed = 10))
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
