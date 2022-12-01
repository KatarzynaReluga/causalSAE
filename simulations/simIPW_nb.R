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


a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# Simple checks of the code ------------------------------------------------------------------
#for (i in 1:NoSim) {
#a  = Sys.time()
#  print(i)
set.seed(a * 2022)

subpopulation <- sample_subpopulations(populations, frac_nc = 0.05, frac_nt = 0.05)
data_sample <- data.frame(populations[subpopulation, ])
data_out_of_sample <- populations[-subpopulation, ]


# NIPW -------------------------------------------------------------------
# EBLUP NIPW
#a = Sys.time()
EE_NIPWf <- hte(type_hte = "NIPW",
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
                                          type_model = "gaussian"))

EE_NIPW <- EE_NIPWf$tau

# MM NIPW
MM_NIPWf <- hte(type_hte = "NIPW",
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
                                       type_model = "continuous"))
MM_NIPW <- MM_NIPWf$tau

# RR NIPW
RR_NIPWf <- hte(type_hte = "NIPW",
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
                                       type_model = "gaussian"))
RR_NIPW <- RR_NIPWf$tau

# OR -------------------------------------------------------------------
# EBLUP NIPW
#a = Sys.time()
XX_NIPWf <- hte(type_hte = "NIPW",
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
                                       type_model = "gaussian"))

XX_NIPW <- XX_NIPWf$tau

# Direct estimator
Dir_tau <- (calculate_tau(list(data_sample), type_tau = "H"))[[1]]$tau


###########################################################################
# Store results in the list - standard for baobab.                         #
############################################################################
Results = list(tau_true = tau_true[[1]],

               EE_NIPW = EE_NIPW,
               MM_NIPW = MM_NIPW,
               RR_NIPW = RR_NIPW,
               XX_NIPW = XX_NIPW,

               Dir_tau = Dir_tau)

outputName = paste("sim_IPW_nb", a, ".RData",sep="")
outputPath = file.path("/home/reluga/Comp", outputName)
#outputPath = file.path( "C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE",outputName)
save("Results", file = outputPath)

#bb <- Sys.time()
