
m = 50
ni = rep(10, m)
Ni = rep(200, m)
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
 Ni_size  = 200,
 m = 50,
 no_sim = 1,
 seed = 10)

 samples <- generate_sample(populations, ni_size = 10,
                            sample_part = "sampled",
                            get_index = TRUE)

 data_sample <- data.frame(samples[[1]]$samp_data)
 index_sample <- samples[[1]]$index_s
 data_out_of_sample <- populations[-index_sample, ]

 model_formula_OR = y ~ X1 + Xo1 + (1|group)


 # EBLUP

 obj_p_score_EBLUP <- list(data_p_score = populations)
 class(obj_p_score_EBLUP) <- "EBLUP"

 a  = Sys.time()
 ps_hat_EBLUP <-  p_score(obj_p_score = obj_p_score_EBLUP,
                          model_formula = A ~ X1 + (1 | group))
 b = Sys.time()

 # MQ

 obj_p_score_MQ <- list(data_p_score = populations)
 class(obj_p_score_MQ) <- "MQ"
 a  = Sys.time()
 ps_hat_MQ <-  p_score(obj_p_score = obj_p_score_MQ,
                       model_formula = A ~ X1 + (1 | group))
 b = Sys.time()
 # RF

 obj_p_score_RF <- list(data_p_score = populations)
 class(obj_p_score_RF) <- "RF"

 a = Sys.time()
 ps_hat_RF <-  p_score(
   obj_p_score = obj_p_score_RF,
   model_formula = A ~ X1 + (1 | group),
   tune_RF = TRUE
 )
 b = Sys.time()
 # XGB

 obj_p_score_XGB <- list(data_p_score = populations)
 class(obj_p_score_XGB) <- "XGB"
 a = Sys.time()
 ps_hat_XGB <-  p_score(
   obj_p_score = obj_p_score_XGB,
   model_formula = A ~ X1 + (1 | group),
   xgboost_params = list(
     CV_XGB = TRUE,
     nfolds = 5,
     nrounds = 50
   )
 )
 b = Sys.time()
