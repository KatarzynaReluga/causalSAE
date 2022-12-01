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
# EEM
EEM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

EEM_AIPW <- EEM_AIPWf$tau
# EER --------------------------------------------------
EER_AIPWf <- hte(type_hte = "AIPW",
                data_sample,
                data_out_of_sample,
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                 method = "EBLUP",
                                 tune_RF = FALSE,
                                 type_model = "gaussian"),
                params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                      method =  "EBLUP"),
                params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                       method = "RF",
                                       tune_RF = FALSE))

EER_AIPW <- EER_AIPWf$tau

# EEX ----------------------
EEX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

EEX_AIPW <- EEX_AIPWf$tau

# EMM ----------------------------------------------------------------------
EMM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

EMM_AIPW <- EMM_AIPWf$tau

# EMR --------------------------------------------
EMR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

EMR_AIPW <- EMR_AIPWf$tau

# EMX -----------------------------------------
EMX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

EMX_AIPW <- EMX_AIPWf$tau

# ERM ------------------------------------------------------
ERM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        tune_RF = FALSE,
                                        type_model = "continuous"))

ERM_AIPW <- ERM_AIPWf$tau

# ERR ------------------------------------------------------------------
ERR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

ERR_AIPW <- ERR_AIPWf$tau

# ERX -------------------------------------------------------------------
ERX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

ERX_AIPW <- ERX_AIPWf$tau

# EXM --------------------------------------------------
EXM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = A ~ X1 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

EXM_AIPW <- EXM_AIPWf$tau

# EXR --------------------------------------------------
EXR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

EXR_AIPW <- EXR_AIPWf$tau

# EXX ----------------------------------------------------------------
EXX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

EXX_AIPW <- EXX_AIPWf$tau

# MM -----------------------------------------------------------
#######################################################################
# MEE ------------------------------------------------------
MEE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "EBLUP",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MEE_AIPW <- MEE_AIPWf$tau

# MER ------------------------------------------------------
MER_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

MER_AIPW <- MER_AIPWf$tau

# MER ------------------------------------------------------
MEX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

MEX_AIPW <- MEX_AIPWf$tau

# MME ------------------------------------------------------
MME_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MME_AIPW <- MME_AIPWf$tau

# MER ------------------------------------------------------
MMR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

MMR_AIPW <- MMR_AIPWf$tau

# MMX ------------------------------------------------------
MMX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

MMX_AIPW <- MMX_AIPWf$tau


# MRE ------------------------------------------------------
MRE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = A ~ X1 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MRE_AIPW <- MRE_AIPWf$tau

# MRR ------------------------------------------------------
MRR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method = "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

MRR_AIPW <- MRR_AIPWf$tau

# MRX ------------------------------------------------------
MRX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

MRX_AIPW <- MRX_AIPWf$tau

# MXE ------------------------------------------------------
MXE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MXE_AIPW <- MXE_AIPWf$tau

# MRR ------------------------------------------------------
MXR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

MXR_AIPW <- MXR_AIPWf$tau

# MRX ------------------------------------------------------
MXX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method = "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

MXX_AIPW <- MXX_AIPWf$tau

# RR -----------------------------------------------------------
#######################################################################
# REE ------------------------------------------------------
REE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

REE_AIPW <- REE_AIPWf$tau

# REM ------------------------------------------------------
REM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = A ~ X1 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

REM_AIPW <- REM_AIPWf$tau
# REX ------------------------------------------------------
REX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

REX_AIPW <- REX_AIPWf$tau


# RME ------------------------------------------------------
RME_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        tune_RF = FALSE,
                                        type_model = "gaussian"))

RME_AIPW <- RME_AIPWf$tau

# RMM ------------------------------------------------------
RMM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

RMM_AIPW <- RMM_AIPWf$tau
# RMX ------------------------------------------------------
RMX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

RMX_AIPW <- RMX_AIPWf$tau

# RRE ------------------------------------------------------
RRE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

RRE_AIPW <- RRE_AIPWf$tau

# RRM ------------------------------------------------------
RRM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

RRM_AIPW <- RRM_AIPWf$tau
# RRX ------------------------------------------------------
RRX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

RRX_AIPW <- RRX_AIPWf$tau

# RXE ------------------------------------------------------
RXE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

RXE_AIPW <- RXE_AIPWf$tau

# RXM ------------------------------------------------------
RXM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

RXM_AIPW <- RXM_AIPWf$tau

# RXX ------------------------------------------------------
RXX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

RXX_AIPW <- RXX_AIPWf$tau

# X -----------------------------------------------------------
#######################################################################
# XEE ------------------------------------------------------
XEE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XEE_AIPW <- XEE_AIPWf$tau

# XEM ------------------------------------------------------
XEM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "EBLUP",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

XEM_AIPW <- XEM_AIPWf$tau

# XER ------------------------------------------------------
XER_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

XER_AIPW <- XER_AIPWf$tau

# XME ------------------------------------------------------
XME_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XME_AIPW <- XME_AIPWf$tau

# XMM ------------------------------------------------------
XMM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

XMM_AIPW <- XMM_AIPWf$tau


# XER ------------------------------------------------------
XMR_AIPWf <- hte(type_hte = "AIPW",
                data_sample,
                data_out_of_sample,
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                 method = "XGB",
                                 xgboost_params = list(CV_XGB = FALSE,
                                                       nfolds = 5,
                                                       nrounds = 50)),
                params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                      method =  "MQ"),
                params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                       method = "RF",
                                       tune_RF = FALSE))

XMR_AIPW <- XMR_AIPWf$tau

# XRE ------------------------------------------------------
XRE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XRE_AIPW <- XRE_AIPWf$tau

# XRM ------------------------------------------------------
XRM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

XRM_AIPW <- XRM_AIPWf$tau


# XRR ------------------------------------------------------
XRR_AIPWf <- hte(type_hte = "AIPW",
                data_sample,
                data_out_of_sample,
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                 method = "XGB",
                                 xgboost_params = list(CV_XGB = TRUE,
                                                       nfolds = 5,
                                                       nrounds = 50)),
                params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                      method =  "RF",
                                      tune_RF = FALSE),
                params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                       method = "RF",
                                       tune_RF = FALSE))

XRR_AIPW <- XRR_AIPWf$tau

# XRE ------------------------------------------------------
XXE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XXE_AIPW <- XXE_AIPWf$tau

# XRM ------------------------------------------------------
XXM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

XXM_AIPW <- XXM_AIPWf$tau


# XER ------------------------------------------------------
XXR_AIPWf <- hte(type_hte = "AIPW",
                data_sample,
                data_out_of_sample,
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                 method = "XGB",
                                 xgboost_params = list(CV_XGB = TRUE,
                                                       nfolds = 5,
                                                       nrounds = 50)),
                params_p_score = list(model_formula =   A ~ X1 + (1|group),
                                      method =  "XGB",
                                      xgboost_params = list(CV_XGB = FALSE,
                                                            nfolds = 5,
                                                            nrounds = 50)),
                params_impute_y = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
                                       method = "RF",
                                       tune_RF = FALSE))

XXR_AIPW <- XXR_AIPWf$tau


# Direct estimator
Dir_tau <- (calculate_tau(list(data_sample), type_tau = "H"))[[1]]$tau

###########################################################################
# Store results in the list - standard for baobab.                         #
############################################################################
Results = list(tau_true = tau_true[[1]],

               # EBLUP -------
               EEM_AIPW = EEM_AIPW,
               EER_AIPW = EER_AIPW,
               EEX_AIPW = EEX_AIPW,

               EMM_AIPW = EMM_AIPW,
               EMR_AIPW = EMR_AIPW,
               EMX_AIPW = EMX_AIPW,

               ERM_AIPW = ERM_AIPW,
               ERR_AIPW = ERR_AIPW,
               ERX_AIPW = ERX_AIPW,

               EXM_AIPW = EXM_AIPW,
               EXR_AIPW = EXR_AIPW,
               EXX_AIPW = EXX_AIPW,

               # MQ ----------------
               MEE_AIPW = MEE_AIPW,
               MER_AIPW = MER_AIPW,
               MEX_AIPW = MEX_AIPW,

               MME_AIPW = MME_AIPW,
               MMR_AIPW = MMR_AIPW,
               MMX_AIPW = MMX_AIPW,

               MRE_AIPW = MRE_AIPW,
               MRR_AIPW = MRR_AIPW,
               MRX_AIPW = MRX_AIPW,

               MXE_AIPW = MXE_AIPW,
               MXR_AIPW = MXR_AIPW,
               MXX_AIPW = MXX_AIPW,

               # R --------------
               REE_AIPW = REE_AIPW,
               REM_AIPW = REM_AIPW,
               REX_AIPW = REX_AIPW,

               RME_AIPW = RME_AIPW,
               RMM_AIPW = RMM_AIPW,
               RMX_AIPW = RMX_AIPW,

               RRE_AIPW = RRE_AIPW,
               RRM_AIPW = RRM_AIPW,
               RRX_AIPW = RRX_AIPW,

               RXE_AIPW = RXE_AIPW,
               RXM_AIPW = RXM_AIPW,
               RXX_AIPW = RXX_AIPW,

               # X --------------
               XEE_AIPW = XEE_AIPW,
               XEM_AIPW = XEM_AIPW,
               XER_AIPW = XER_AIPW,

               XME_AIPW = XME_AIPW,
               XMM_AIPW = XMM_AIPW,
               XMR_AIPW = XMR_AIPW,

               XRE_AIPW = XRE_AIPW,
               XRM_AIPW = XRM_AIPW,
               XRR_AIPW = XRR_AIPW,

               XXE_AIPW = XXE_AIPW,
               XXM_AIPW = XXM_AIPW,
               XXR_AIPW = XXR_AIPW,

               Dir_tau = Dir_tau)

outputName = paste("sim_nb_AIPW", a, ".RData",sep="")
outputPath = file.path("/home/reluga/Comp", outputName)
#outputPath = file.path( "C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE",outputName)
save("Results", file = outputPath)

#bb <- Sys.time()
