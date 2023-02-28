# rm(list=ls())
# .rs.restartR()
#
#

setwd("./causalSAE")
devtools::load_all()

m = 100
ni = rep(10, m)
Ni = rep(200, m)
N = sum(Ni)
n = sum(ni)

# Exp works but it needs to have small coef

var_norm <- matrix(c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5,
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5,
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5,
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1), nrow = 10, ncol = 10, byrow = T)

## Retrieve group indicator ---------------------------------------------------
group = rep(1:m, times = Ni)

X <- generate_X(
  n = N,
  p = 10,
  covariance_norm = var_norm,
  cov_type = "norm",
  seed = 3
)

untreat1 <- list()
untreat2 <- list()

treat1 <- list()
treat2 <- list()

for (i in 1:m) {
  untreat1[[i]] <- matrix(rep(rt(10, 3), 200), nrow = 200, ncol = 10, byrow = T)
  untreat2[[i]] <- 0.1 *  matrix(rep(rt(10, 3), 200), nrow = 200, ncol = 10, byrow = T)

  treat1[[i]] <- untreat1[[i]] + 100 * matrix(rep(rt(10, 3), 200), nrow = 200, ncol = 10, byrow = T)
  treat2[[i]] <- untreat2[[i]] + 0.1 * matrix(rep(rt(10, 3), 200), nrow = 200, ncol = 10, byrow = T)
}

coef_01 <- do.call("rbind", untreat1)
coef_02 <- do.call("rbind", untreat2)

coef_11 <- do.call("rbind", treat1)
coef_12 <- do.call("rbind", treat2)

set.seed(39)
y0ne <- rowMeans(10 + coef_01 * X + exp(coef_02 * X))
y0 <- y0ne + rnorm(200 * m, 0, mean(y0ne))

y1ne <- rowMeans(20 + coef_11 * X + exp(coef_12 * X))
y1 <- y1ne + rnorm(200 * m, 0, mean(y1ne))


tau_treat <- aggregate(y1, list(group), FUN = mean)$x
tau_untreat <- aggregate(y0, list(group), FUN = mean)$x
tau_true = tau_treat - tau_untreat
#sort(tau_true)
#plot(1:100, tau_true)

# Propensity score
intercept_p_score = - 0.1
coef_p_score = rep(0.5, 10)
Xreg_p_score = rowMeans(coef_p_score * X)

rand_eff_p_score <- list(var_re = 0.25,
                         mean_re = 0,

                         frac_out = 0,
                         var_re_out = 0,

                         mean_re_out = 0)

re_treat <- generate_re(
  n = m,
  type_re = "random_effects",
  var_re = rand_eff_p_score$var_re,
  mean_re = rand_eff_p_score$mean_re,
  frac_out = rand_eff_p_score$frac_out,
  var_re_out = rand_eff_p_score$var_re_out,
  mean_re_out = rand_eff_p_score$mean_re_out,
  seed = 1
)

re_treat_repeat <- rep(re_treat, times = Ni)
exp_p_score = exp(intercept_p_score + Xreg_p_score + re_treat_repeat)

p_score = exp_p_score * (1 + exp_p_score)^(-1)
A <- rbinom(N, 1, p_score)
#mean(A)
A_group <- aggregate(A, list(group), FUN = mean)$x
#A_group
names(X) <-  paste0("X", 1:ncol(X))
y = A * y1 + (1 - A) * y0

populations <- data.frame(X, A, group, p_score, y)


a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# Simple checks of the code ------------------------------------------------------------------
#for (i in 1:NoSim) {
#b  = Sys.time()
#  print(i)
set.seed(a * 2022)

subpopulation <- sample_subpopulations(populations,
                                       frac_nc = 0.1, frac_nt = 0.1,
                                       seed = set.seed(a * 2022))
data_sample <- data.frame(populations[subpopulation, ])
data_out_of_sample <- populations[-subpopulation, ]
#######################################################################################################
# X -----------------------------------------------------------
#######################################################################
# XEE ------------------------------------------------------
XEE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XEE_AIPW <- XEE_AIPWf$tau

# XEM ------------------------------------------------------
XEM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

XEM_AIPW <- XEM_AIPWf$tau

# XER ------------------------------------------------------
XER_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

XER_AIPW <- XER_AIPWf$tau

# XME ------------------------------------------------------
XME_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XME_AIPW <- XME_AIPWf$tau

# XMM ------------------------------------------------------
XMM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

XMM_AIPW <- XMM_AIPWf$tau


# XER ------------------------------------------------------
XMR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

XMR_AIPW <- XMR_AIPWf$tau

###########################################################################
# Store results in the list - standard for baobab.                         #
############################################################################
Results = list(tau_true = tau_true,
               tau_treat = tau_treat,
               tau_untreat = tau_untreat,

               XEE_AIPW = XEE_AIPW,
               XEM_AIPW = XEM_AIPW,
               XER_AIPW = XER_AIPW,

               XME_AIPW = XME_AIPW,
               XMM_AIPW = XMM_AIPW,
               XMR_AIPW = XMR_AIPW)

outputName = paste("Mplus100_2", a, ".RData", sep = "")
outputPath = file.path("/home/reluga/Comp", outputName)
#outputPath = file.path("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE",outputName)
save("Results", file = outputPath)

#c = Sys.time()
