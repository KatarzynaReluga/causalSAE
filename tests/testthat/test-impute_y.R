#library(causalSAE)

set.seed(123456)

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


samples <- generate_sample(
  populations,
  ni_size = 5,
  sample_part = "sampled",
  get_index = TRUE
)

data_sample <- data.frame(samples[[1]]$samp_data)
index_sample <- samples[[1]]$index_s
data_out_of_sample <- populations[-index_sample,]

model_formula = y ~ X1 + Xo1 + A + (1 + A || group)

# Impute EBLUP ----------------------------------------------------

impute_EBLUP <- impute_y(
  model_formula,
  data_sample,
  data_out_of_sample,
  method = "EBLUP",
  type_model = "gaussian"
)

MSE_EBLUP <- mean((populations$y - impute_EBLUP$y_full_imputed)^2)
RMSE_EBLUP <- mean((populations$y - impute_EBLUP$y_full_imputed)^2/populations$y)

test_that("Output is correct", {
  expect_output(str(impute_EBLUP), "List")
  expect_equal(length(impute_EBLUP), 6)
})

# Impute MQ ----------------------------------------------------

impute_MQ <- impute_y(
  model_formula,
  data_sample,
  data_out_of_sample,
  method = "MQ",
  type_model = "continuous"
)

MSE_MQ <- mean((populations$y - impute_MQ$y_full_imputed)^2)
RMSE_MQ <- mean((populations$y - impute_MQ$y_full_imputed)^2/populations$y)

test_that("Output is correct", {
  expect_output(str(impute_MQ), "List")
  expect_equal(length(impute_MQ), 6)
})

# Impute RF --------------------------------------------------------

impute_RF <- impute_y(
  model_formula,
  data_sample,
  data_out_of_sample,
  method = "RF",
  tune_RF = TRUE,
  xgboost_params = list(
    CV_XGB = TRUE,
    nfolds = 5,
    nrounds = 50
  )
)

MSE_RF <- mean((populations$y - impute_RF$y_full_imputed)^2)
RMSE_RF <- mean((populations$y - impute_RF$y_full_imputed)^2/populations$y)

test_that("Output is correct", {
  expect_output(str(impute_RF), "List")
  expect_equal(length(impute_RF), 6)
})

# Impute XGB ---------------------------------------------------------

impute_XGB <- impute_y(
  model_formula,
  data_sample,
  data_out_of_sample,
  method = "XGB",
  xgboost_params = list(
    CV_XGB = TRUE,
    nfolds = 5,
    nrounds = 50
  )
)

MSE_XGB <- mean((populations$y - impute_XGB$y_full_imputed)^2)
RMSE_XGB <- mean((populations$y - impute_XGB$y_full_imputed)^2/populations$y)

test_that("Output is correct", {
  expect_output(str(impute_XGB), "List")
  expect_equal(length(impute_XGB), 6)
})
