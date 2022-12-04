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
 no_sim = 100,
 seed = 1)
