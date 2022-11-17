#library(causalSAE)

set.seed(123456)

m = 50
ni = rep(5, m)
Ni = rep(1000, m)
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

populations <- generate_pop(
  X,
  X_outcome,
  coeffs = get_default_coeffs(),
  errors_outcome = get_default_errors_outcome(),
  rand_eff_outcome = get_default_rand_eff_outcome(),
  rand_eff_p_score = get_default_rand_eff_p_score(),
  regression_type = "continuous",
  Ni_size  = 1000,
  m = 50,
  no_sim = 1,
  seed = 10
)


# EBLUP

obj_p_score_EBLUP <- list(data_p_score = populations)
class(obj_p_score_EBLUP) <- "EBLUP"

ps_hat_EBLUP <-  p_score(obj_p_score = obj_p_score_EBLUP,
                         model_formula = A ~ X1 + (1 | group))
# MQ

obj_p_score_MQ <- list(data_p_score = populations)
class(obj_p_score_MQ) <- "MQ"

ps_hat_MQ <-  p_score(obj_p_score = obj_p_score_MQ,
                      model_formula = A ~ X1 + (1 | group))

# RF

obj_p_score_RF <- list(data_p_score = populations)
class(obj_p_score_RF) <- "RF"

ps_hat_RF <-  p_score(
  obj_p_score = obj_p_score_RF,
  model_formula = A ~ X1 + (1 | group),
  tune_RF = FALSE
)

# XGB

obj_p_score_XGB <- list(data_p_score = populations)
class(obj_p_score_XGB) <- "XGB"

ps_hat_XGB <-  p_score(
  obj_p_score = obj_p_score_XGB,
  model_formula = A ~ X1 + (1 | group),
  xgboost_params = list(
    CV_XGB = FALSE,
    nfolds = 5,
    nrounds = 50
  )
)
