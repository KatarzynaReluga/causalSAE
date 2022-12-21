#
# LMM with sinexp --------------------------------------
# No random effects
#

m = 50
ni = rep(10, m)
Ni = rep(200, m)
N = sum(Ni)
n = sum(ni)

# Exp works but it needs to have small coef

X_unif <- generate_X(
  n = N,
  p = 5,
  covariance_norm = NULL,
  cov_type = "unif",
  seed = 3
)

var_norm <- matrix(c(1, 0.3, 0.3, 0.3, 0.3,
                     0.3, 1, 0.3, 0.3, 0.3,
                     0.3, 0.3, 1, 0.3, 0.3,
                     0.3, 0.3, 0.3, 1, 0.3,
                     0.3, 0.3, 0.3, 0.3, 1), nrow = 5, ncol = 5, byrow = T)

X_norm <- generate_X(
  n = N,
  p = 5,
  covariance_norm = var_norm,
  cov_type = "norm",
  seed = 3
)

X = cbind(X_unif, X_norm)


coeffs = list(intercept_outcome = 100,
              intercept_p_score = - 0.3,

              coef_outcome = rep(c(1, 2), 5),
              coef_p_score = rep(0.5, 10),

              mean_A = 10,
              var_A = 1)

#
errors_outcome <- list(var_e = 1,
                       mean_e = 0,

                       frac_out = 0.1,
                       var_e_out = 10,

                       mean_e_out = 0,
                       disturbance_outcome = 5)

re_outcome <- list(var_re = 3,
                   mean_re = 0,

                   frac_out = 0,
                   var_re_out = 0,

                   mean_re_out = 0)

re_p_score <- list(var_re = 0.25,
                   mean_re = 0,

                   frac_out = 0,
                   var_re_out = 0,

                   mean_re_out = 0)

populations <- generate_pop(X = X,
                            coeffs = coeffs,
                            errors_outcome = errors_outcome,
                            rand_eff_outcome = re_outcome,
                            rand_eff_p_score = re_p_score,
                            regression_type = "continuous",
                            Ni_size  = 200,
                            m = 50,
                            no_sim = 1,
                            fct_coef = "sin_exp",
                            seed = 1)
#head(populations)
# Retrieve y and groups -------------------------------------------
y1 <- populations$y1
y0 <- populations$y0
group <- populations$group
# Compute true tau -------------------------------------------
tau_treat <- aggregate(y1, list(group), FUN = mean)$x
tau_untreat <- aggregate(y0, list(group), FUN = mean)$x
tau = tau_treat - tau_untreat

a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# Simple checks of the code ------------------------------------------------------------------
#for (i in 1:NoSim) {
#a  = Sys.time()
#  print(i)
set.seed(a * 2022)

subpopulation <- sample_subpopulations(populations,
                                       frac_nc = 0.1, frac_nt = 0.1,
                                       seed = set.seed(a * 2022))
data_sample <- data.frame(populations[subpopulation, ])
data_out_of_sample <- populations[-subpopulation, ]
