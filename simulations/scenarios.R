# LMM
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

                       frac_out = 0,
                       var_e_out = 0,

                       mean_e_out = 0,
                       disturbance_outcome = 5)

populations <- generate_pop(X = X,
                            coeffs = coeffs,
                            errors_outcome = errors_outcome,
                            rand_eff_outcome = get_default_rand_eff_outcome(),
                            rand_eff_p_score = get_default_rand_eff_p_score(),
                            regression_type = "continuous",
                            Ni_size  = 200,
                            m = 50,
                            no_sim = 1,
                            seed = 1)
y1 <- populations$y1
y0 <- populations$y0
tau_treat <- aggregate(y1, list(group), FUN = mean)$x
tau_untreat <- aggregate(y0, list(group), FUN = mean)$x
tau = tau_treat - tau_untreat

hajek <- calculate_tau(populations, type_tau = "H")
ht <- calculate_tau(populations, type_tau = "HT")

#plot(1:50, hajek[[1]]$tau, type = "l")
#lines(tau, col = 2)

# Sin

# Exp

# Non linear scenario


x = sort(X[, 1])
nu_x <- sin(10*x)
f_x <- cos(x) + sin(10*x2)

 # coef

X_outcome <- generate_X(
 n = N,
 p = 5,
 covariance_norm = NULL,
 cov_type = "lognorm",
 seed = 1
)

x1 = X[, 1]
x2 = X[, 2]
x3 = X[, 3]
x4 = X[, 4]
x5 = X[, 5]
# Maybe scenario with normal distribution (?)
f_x <- (10*sin(pi*x1*x2) + 20*(x3 - 0.5)^2 + 10*x4 + 5*x5 - 14.3)/4.9
f_x <- (0.1*exp(4*x1)+4/(1 + exp(-20*(x2-0.5))) + 3*x3 +
          2*x4 + x5 - 6.3)/2.5



 nu_x <- nu(x1, x2, x3, x4, x5)
 f_x <- 5*(x1 - 0.5)^2 + 5*(x2 - 0.5)^2 +
   6*(x3 - 0.5)^2 + 6.7*(x4 - 0.5)^2 + 7.5*(x5 - 0.5)^2
 pi_score <- plogis(0.5 * x1 + 0.5 * x2 + 0.5 * x3 + 0.5 * x4 + 0.5 * x5 - 0.5)
 a <- rbinom(n, 1, pi_score)
 y <- nu_x + a * f_x + rnorm(n)


 nu_x <- nu(x1, x2, x3, x4, x5)
 f_x <- (0.1*exp(4*x1)+4/(1 + exp(-20*(x2-0.5))) + 3*x3 +
           2*x4 + x5 - 6.3)/2.5
 pi_score <- plogis(0.5 * x1 + 0.5 * x2 + 0.5 * x3 + 0.5 * x4 + 0.5 * x5 - 0.5)
 a <- rbinom(n, 1, pi_score)
 y <- nu_x + a * f_x + rnorm(n)


 y = coeffs$intercept_outcome + Xreg_outcome + coef_A * A + re_repeat + e

 #Zmień X_reg function, dodaj funkcje
