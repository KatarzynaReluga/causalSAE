# Model based simulations
setwd("./causalSAE")
devtools::load_all()

#aa <- Sys.time()

# Set seed
set.seed(100)

m = 50
ni = rep(10, m)
Ni = rep(100, m)
N = sum(Ni)
n = sum(ni)

n_boot = 500

# Generate covariates
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

# Generate populations
populations <- generate_pop(X, X_outcome,
                            coeffs = get_default_coeffs(),
                            errors_outcome = get_default_errors_outcome(),
                            rand_eff_outcome = get_default_rand_eff_outcome(),
                            rand_eff_p_score = get_default_rand_eff_p_score(),
                            regression_type = "continuous",
                            Ni_size = 100,
                            m = 50,
                            no_sim = 1,
                            seed = 10)
tau_true <- calculate_tau(list(populations), type_tau = "H")

#a = 6
a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# Simple checks of the code ------------------------------------------------------------------
#for (i in 1:NoSim) {
#a  = Sys.time()
#  print(i)
set.seed(a * 2022)

subpopulation <- sample_subpopulations(populations, frac_nc = 0.05, frac_nt = 0.05)
data_sample <- data.frame(populations[subpopulation, ])
data_out_of_sample <- populations[-subpopulation, ]

# AIPW -------------------------------------------------------------------
# EER
EEM_AIPWf <- hte(type_hte = "NIPW",
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
                                        method = "MQ",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50),
                                        type_model = "continuous"))

EEM_AIPW <- EEM_AIPWf$tau
# EER
EER_AIPWf <- hte(type_hte = "NIPW",
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
                                       type_model = "gaussian"))

EER_AIPW <- EER_AIPWf$tau

# EEX ----------------------
EEX_AIPWf <- hte(type_hte = "NIPW",
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
                                        method = "XGB",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50),
                                        type_model = "gaussian"))

EEX_AIPW <- EEX_AIPWf$tau

# EMM ----------------------------------------------------------------------
EMM_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

EMM_AIPW <- EMM_AIPWf$tau

# EMR --------------------------------------------
EMR_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "gaussian"))

EMR_AIPW <- EMR_AIPWf$tau

# EMX -----------------------------------------
EMX_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "gaussian"))

EMX_AIPW <- EMX_AIPWf$tau

# ERM ------------------------------------------------------
ERM_AIPWf <- hte(type_hte = "NIPW",
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
                                       method =  "RF",
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
                                        type_model = "continuous"))

ERM_AIPW <- ERM_AIPWf$tau

# ERR ------------------------------------------------------------------
ERR_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "gaussian"))

ERR_AIPW <- ERR_AIPWf$tau

# ERX -------------------------------------------------------------------
ERX_AIPWf <- hte(type_hte = "NIPW",
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
                                       method =  "RF",
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
                                        type_model = "gaussian"))

ERX_AIPW <- ERX_AIPWf$tau

# EXM --------------------------------------------------
EXM_AIPWf <- hte(type_hte = "NIPW",
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
                                        method = "MQ",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50),
                                        type_model = "continuous"))

EXM_AIPW <- EXM_AIPWf$tau

# EXR --------------------------------------------------
EXR_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

EXR_AIPW <- EXR_AIPWf$tau

# EXX ----------------------------------------------------------------
EXX_AIPWf <- hte(type_hte = "NIPW",
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
                                        method = "XGB",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50),
                                        type_model = "gaussian"))

EXX_AIPW <- EXX_AIPWf$tau

# MM -----------------------------------------------------------
#######################################################################
# MEE ------------------------------------------------------
MEE_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

MEE_AIPW <- MEE_AIPWf$tau

# MER ------------------------------------------------------
MER_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

MER_AIPW <- MER_AIPWf$tau

# MER ------------------------------------------------------
MEX_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

MEX_AIPW <- MEX_AIPWf$tau

# MME ------------------------------------------------------
MME_AIPWf <- hte(type_hte = "NIPW",
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
                                        method = "EBLUP",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50),
                                        type_model = "continuous"))

MME_AIPW <- MME_AIPWf$tau

# MER ------------------------------------------------------
MMR_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

MMR_AIPW <- MMR_AIPWf$tau

# MER ------------------------------------------------------
MMX_AIPWf <- hte(type_hte = "NIPW",
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
                                        method = "XGB",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50),
                                        type_model = "continuous"))

MMX_AIPW <- MMX_AIPWf$tau

# MME ------------------------------------------------------
MME_AIPWf <- hte(type_hte = "NIPW",
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
                                        method = "EBLUP",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50),
                                        type_model = "continuous"))

MME_AIPW <- MME_AIPWf$tau

# MER ------------------------------------------------------
MMR_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

MMR_AIPW <- MMR_AIPWf$tau

# MER ------------------------------------------------------
MMX_AIPWf <- hte(type_hte = "NIPW",
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
                                        method = "XGB",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50),
                                        type_model = "continuous"))

MMX_AIPW <- MMX_AIPWf$tau

# MRE ------------------------------------------------------
MRE_AIPWf <- hte(type_hte = "NIPW",
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
                                       method =  "RF",
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
                                        type_model = "continuous"))

MRE_AIPW <- MRE_AIPWf$tau

# MRR ------------------------------------------------------
MRR_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

MRR_AIPW <- MRR_AIPWf$tau

# MRX ------------------------------------------------------
MRX_AIPWf <- hte(type_hte = "NIPW",
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
                                       method =  "RF",
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
                                        type_model = "continuous"))

MRX_AIPW <- MRX_AIPWf$tau

# MXE ------------------------------------------------------
MXE_AIPWf <- hte(type_hte = "NIPW",
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
                                       method =  "XGB",
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
                                        type_model = "continuous"))

MXE_AIPW <- MXE_AIPWf$tau

# MRR ------------------------------------------------------
MXR_AIPWf <- hte(type_hte = "NIPW",
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
                                       method =  "XGB",
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
                                        type_model = "continuous"))

MXR_AIPW <- MXR_AIPWf$tau

# MRX ------------------------------------------------------
MXX_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

MXX_AIPW <- MXX_AIPWf$tau

# MM -----------------------------------------------------------
#######################################################################
# REE ------------------------------------------------------
REE_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

REE_AIPW <- REE_AIPWf$tau

# REM ------------------------------------------------------
REM_AIPWf <- hte(type_hte = "NIPW",
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
                                       method =  "EBLUP",
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
                                        type_model = "continuous"))

REM_AIPW <- REM_AIPWf$tau
# MER ------------------------------------------------------
REX_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

REX_AIPW <- REX_AIPWf$tau


# R -----------------------------------------------------------
#######################################################################
# REE ------------------------------------------------------
REE_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

REE_AIPW <- REE_AIPWf$tau

# RMM ------------------------------------------------------
REM_AIPWf <- hte(type_hte = "NIPW",
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
                                       method =  "EBLUP",
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
                                        type_model = "continuous"))

REM_AIPW <- REM_AIPWf$tau
# REX ------------------------------------------------------
REX_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

REX_AIPW <- REX_AIPWf$tau

# REE ------------------------------------------------------
RME_AIPWf <- hte(type_hte = "NIPW",
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
                                       method =  "MQ",
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
                                        type_model = "continuous"))

RME_AIPW <- RME_AIPWf$tau

# RMM ------------------------------------------------------
RMM_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

RMM_AIPW <- RMM_AIPWf$tau
# REX ------------------------------------------------------
REX_AIPWf <- hte(type_hte = "NIPW",
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
                                        type_model = "continuous"))

REX_AIPW <- REX_AIPWf$tau

# Direct estimator
Dir_tau <- (calculate_tau(list(data_sample), type_tau = "H"))[[1]]$tau

###########################################################################
# Store results in the list - standard for baobab.                         #
############################################################################
Results = list(tau_true = tau_true[[1]],

               EBLUP_OR = EBLUP_OR,
               MQ_OR = MQ_OR,
               RF_OR = RF_OR,
               XGB_OR = XGB_OR,

               EE_NIPW = EE_NIPW,
               MM_NIPW = MM_NIPW,
               RR_NIPW = RR_NIPW,
               XX_NIPW = XX_NIPW,

               EER_AIPW = EER_AIPW,
               EEX_AIPW = EEX_AIPW,
               MMR_AIPW = MMR_AIPW,
               MMX_AIPW = MMX_AIPW,

               Dir_tau = Dir_tau)

outputName = paste("sim_nb", a, ".RData",sep="")
outputPath = file.path("/home/reluga/Comp", outputName)
#outputPath = file.path( "C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE",outputName)
save("Results", file = outputPath)

#bb <- Sys.time()
