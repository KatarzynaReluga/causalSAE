# rm(list=ls())
# .rs.restartR()
#
#

setwd("./causalSAE")
devtools::load_all()

#aaa <- Sys.time()

m = 25
Nii = 1000
Ni = rep(Nii, m)
N = sum(Ni)

#ni = rep(10, m)
#n = sum(ni)

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
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1),
                   nrow = 10, ncol = 10, byrow = T)

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
  untreat1[[i]] <- matrix(rep(rt(10, 3), Nii), nrow = Nii, ncol = 10, byrow = T)
  untreat2[[i]] <- 0.1 *  matrix(rep(rt(10, 3), Nii), nrow = Nii, ncol = 10, byrow = T)

  treat1[[i]] <- untreat1[[i]] + 100 * matrix(rep(rt(10, 3), Nii), nrow = Nii, ncol = 10, byrow = T)
  treat2[[i]] <- untreat2[[i]] + 0.1 * matrix(rep(rt(10, 3), Nii), nrow = Nii, ncol = 10, byrow = T)
}

coef_01 <- do.call("rbind", untreat1)
coef_02 <- do.call("rbind", untreat2)

coef_11 <- do.call("rbind", treat1)
coef_12 <- do.call("rbind", treat2)

set.seed(39)
y0ne <- rowMeans(10 + coef_01 * X + exp(coef_02 * X))
y0 <- y0ne + rnorm(Nii * m, 0, mean(y0ne))

y1ne <- rowMeans(20 + coef_11 * X + exp(coef_12 * X))
y1 <- y1ne + rnorm(Nii * m, 0, mean(y1ne))


tau_treat <- aggregate(y1, list(group), FUN = mean)$x
tau_untreat <- aggregate(y0, list(group), FUN = mean)$x
tau_true = tau_treat - tau_untreat
#sort(tau_true)
#tau_true
#plot(1:25, tau_true)

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

A_group <- aggregate(A, list(group), FUN = mean)$x

names(X) <-  paste0("X", 1:ncol(X))
y = A * y1 + (1 - A) * y0

populations <- data.frame(X, A, group, p_score, y)

Nc = as.numeric(table(populations$group[populations$A == 0]))
Nt = as.numeric(table(populations$group[populations$A == 1]))

frac_nt = 0.02
nt <- round(frac_nt * Nt)
frac_nc0 <- c(runif(3 * m/5, 0.06, 0.5), runif(2 * m/5, 0.51, 1))
nc <- ceiling(frac_nc0 * nt)
frac_nc <- nc/Nc

#bbb <- Sys.time()

a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# Simple checks of the code ------------------------------------------------------------------
#for (i in 1:NoSim) {
#ccc  = Sys.time()

set.seed(a * 2022)
subpopulation <- sample_subpopulations(populations,
                                       frac_nc = frac_nc, frac_nt = frac_nt,
                                       seed = set.seed(a * 2022))
data_sample <- data.frame(populations[subpopulation, ])
data_out_of_sample <- populations[-subpopulation, ]

ni = as.numeric(table(data_sample$group))
nc = as.numeric(table(data_sample$group[data_sample$A == 0]))
nt = as.numeric(table(data_sample$group[data_sample$A == 1]))
frac_nN <- dim(data_sample)[1]/dim(populations)[1]

#nc
#nt
#frac_nN
#######################################################################################################
######
# OR #
######
#######################################################################################################
# EBLUP OR ----------------------------------------------------------------------------------------------
E_ORf <- hte(type_hte = "OR",
             data_sample,
             data_out_of_sample,
             params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                              method = "EBLUP",
                              type_model = "gaussian"))
E_OR <- E_ORf$tau

# MQ OR --------------------------------------------------------------------------------------------------------
M_ORf <- hte(type_hte = "OR",
             data_sample,
             data_out_of_sample,
             params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                              method = "MQ",
                              type_model = "continuous"))
M_OR <- M_ORf$tau

# RF OR ------------------------------------------------------------------------------------------------------------

R_ORf <- hte(type_hte = "OR",
             data_sample,
             data_out_of_sample,
             params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                              method = "RF",
                              tune_RF = FALSE))
R_OR <- R_ORf$tau

# EBLUP XGB -----------------------------------------------------------------------------------------------------------

X_ORf <- hte(type_hte = "OR",
             data_sample,
             data_out_of_sample,
             params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                              method = "XGB",
                              xgboost_params = list(CV_XGB = FALSE,
                                                    nfolds = 5,
                                                    nrounds = 50)))
X_OR <- X_ORf$tau
##############################################################################################################
########
# NIPW #
########
###############################################################################################################
# E
# EE
EE_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula =   A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "EBLUP"),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "EBLUP",
                                       type_model = "gaussian"))

EE_NIPW <- EE_NIPWf$tau

#EM ------------------------------------------------------------------------------------------------
EM_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula =   A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "EBLUP"),
                params_impute_y = list(model_formula = y ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "MQ",
                                       type_model = "continuous"))

EM_NIPW <- EM_NIPWf$tau

# ER ------------------------------------------------------------------------------------------------
ER_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula =  A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "EBLUP"),
                params_impute_y = list(model_formula =  y ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "RF",
                                       tune_RF = FALSE))

ER_NIPW <- ER_NIPWf$tau

#EX --------------------------------------------------------------------------------------------
EX_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "EBLUP"),
                params_impute_y = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)))

EX_NIPW <- EX_NIPWf$tau
##############################################################################################################
# M
# MM -------------------------------------------------------------------------------------------
MM_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "MQ"),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "MQ",
                                       type_model = "continuous"))
MM_NIPW <- MM_NIPWf$tau


# ME --------------------------------------------------------
ME_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "MQ"),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "EBLUP",
                                       type_model = "gaussian"))
ME_NIPW <- ME_NIPWf$tau


# MR --------------------------------------------------------
MR_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "MQ"),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "RF",
                                       tune_RF = FALSE))
MR_NIPW <- MR_NIPWf$tau

# MX --------------------------------------------------------
MX_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "MQ"),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)))
MX_NIPW <- MX_NIPWf$tau

######################################################################################################
# R #
# RR ------------------------------------------------------------------
RR_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "RF",
                                      tune_RF = FALSE),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "RF",
                                       tune_RF = FALSE))
RR_NIPW <- RR_NIPWf$tau
# RE ------------------------------------------------------------------
RE_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "RF",
                                      tune_RF = FALSE),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "EBLUP",
                                       type_model = "gaussian"))
RE_NIPW <- RE_NIPWf$tau

# RM ------------------------------------------------------------------
RM_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "RF",
                                      tune_RF = FALSE),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "MQ",
                                       type_model = "continuous"))
RM_NIPW <- RM_NIPWf$tau
# RX ------------------------------------------------------------------
RX_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "RF",
                                      tune_RF = FALSE),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)))
RX_NIPW <- RX_NIPWf$tau
######################################################################################################
# X #
# XX -------------------------------------------------------------------
XX_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "XGB",
                                      xgboost_params = list(CV_XGB = FALSE,
                                                            nfolds = 5,
                                                            nrounds = 50)),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)))

XX_NIPW <- XX_NIPWf$tau

# XE -------------------------------------------------------------------
XE_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "XGB",
                                      xgboost_params = list(CV_XGB = FALSE,
                                                            nfolds = 5,
                                                            nrounds = 50)),
                params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "EBLUP",
                                       type_model = "gaussian"))

XE_NIPW <- XE_NIPWf$tau
# XM -------------------------------------------------------------------
XM_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "XGB",
                                      xgboost_params = list(CV_XGB = FALSE,
                                                            nfolds = 5,
                                                            nrounds = 50)),
                params_impute_y = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "MQ",
                                       type_model = "continuous"))

XM_NIPW <- XM_NIPWf$tau

# XR -------------------------------------------------------------------
XR_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                      method =  "XGB",
                                      xgboost_params = list(CV_XGB = FALSE,
                                                            nfolds = 5,
                                                            nrounds = 50)),
                params_impute_y = list(model_formula = y ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                       method = "RF",
                                       tune_RF = FALSE))

XR_NIPW <- XR_NIPWf$tau

##############################################################################################################
########
# AIPW #
########
###############################################################################################################
# AIPW -------------------------------------------------------------------
# EEM
EEM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =   A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

EEM_AIPW <- EEM_AIPWf$tau
# EER --------------------------------------------------
EER_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

EER_AIPW <- EER_AIPWf$tau

# EEX ----------------------
EEX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

EEX_AIPW <- EEX_AIPWf$tau

# EMM ----------------------------------------------------------------------
EMM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

EMM_AIPW <- EMM_AIPWf$tau

# EMR --------------------------------------------
EMR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

EMR_AIPW <- EMR_AIPWf$tau

# EMX -----------------------------------------
EMX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

EMX_AIPW <- EMX_AIPWf$tau

# ERM ------------------------------------------------------
ERM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

ERM_AIPW <- ERM_AIPWf$tau

# ERR ------------------------------------------------------------------
ERR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

ERR_AIPW <- ERR_AIPWf$tau

# ERX -------------------------------------------------------------------
ERX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

ERX_AIPW <- ERX_AIPWf$tau

# EXM --------------------------------------------------
EXM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

EXM_AIPW <- EXM_AIPWf$tau

# EXR --------------------------------------------------
EXR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

EXR_AIPW <- EXR_AIPWf$tau

# EXX ----------------------------------------------------------------
EXX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
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
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MEE_AIPW <- MEE_AIPWf$tau

# MER ------------------------------------------------------
MER_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

MER_AIPW <- MER_AIPWf$tau

# MER ------------------------------------------------------
MEX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

MEX_AIPW <- MEX_AIPWf$tau

# MME ------------------------------------------------------
MME_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MME_AIPW <- MME_AIPWf$tau

# MER ------------------------------------------------------
MMR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

MMR_AIPW <- MMR_AIPWf$tau

# MMX ------------------------------------------------------
MMX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

MMX_AIPW <- MMX_AIPWf$tau


# MRE ------------------------------------------------------
MRE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MRE_AIPW <- MRE_AIPWf$tau

# MRR ------------------------------------------------------
MRR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

MRR_AIPW <- MRR_AIPWf$tau

# MRX ------------------------------------------------------
MRX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

MRX_AIPW <- MRX_AIPWf$tau

# MXE ------------------------------------------------------
MXE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MXE_AIPW <- MXE_AIPWf$tau

# MRR ------------------------------------------------------
MXR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula =   A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

MXR_AIPW <- MXR_AIPWf$tau

# MRX ------------------------------------------------------
MXX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
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
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

REE_AIPW <- REE_AIPWf$tau

# REM ------------------------------------------------------
REM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

REM_AIPW <- REM_AIPWf$tau
# REX ------------------------------------------------------
REX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

REX_AIPW <- REX_AIPWf$tau


# RME ------------------------------------------------------
RME_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

RME_AIPW <- RME_AIPWf$tau

# RMM ------------------------------------------------------
RMM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

RMM_AIPW <- RMM_AIPWf$tau
# RMX ------------------------------------------------------
RMX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

RMX_AIPW <- RMX_AIPWf$tau

# RRE ------------------------------------------------------
RRE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

RRE_AIPW <- RRE_AIPWf$tau

# RRM ------------------------------------------------------
RRM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

RRM_AIPW <- RRM_AIPWf$tau
# RRX ------------------------------------------------------
RRX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

RRX_AIPW <- RRX_AIPWf$tau

# RXE ------------------------------------------------------
RXE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

RXE_AIPW <- RXE_AIPWf$tau

# RXM ------------------------------------------------------
RXM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

RXM_AIPW <- RXM_AIPWf$tau

# RXX ------------------------------------------------------
RXX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula =  A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
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
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
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
                                  xgboost_params = list(CV_XGB = FALSE,
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
                                  xgboost_params = list(CV_XGB = FALSE,
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
                                  xgboost_params = list(CV_XGB = FALSE,
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
                                  xgboost_params = list(CV_XGB = FALSE,
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
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

XMR_AIPW <- XMR_AIPWf$tau

# XRE ------------------------------------------------------
XRE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XRE_AIPW <- XRE_AIPWf$tau

# XRM ------------------------------------------------------
XRM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

XRM_AIPW <- XRM_AIPWf$tau


# XRR ------------------------------------------------------
XRR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

XRR_AIPW <- XRR_AIPWf$tau

# XXE ------------------------------------------------------
XXE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XXE_AIPW <- XXE_AIPWf$tau

# XRM ------------------------------------------------------
XXM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "MQ",
                                        type_model = "continuous"))

XXM_AIPW <- XXM_AIPWf$tau


# XER ------------------------------------------------------
XXR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + A + (1 + A||group),
                                        method = "RF",
                                        tune_RF = FALSE))

XXR_AIPW <- XXR_AIPWf$tau
####################
# Direct estimator #
####################
Dir_tau <- (calculate_tau(list(data_sample), type_tau = "H"))[[1]]$tau

###########################################################################
# Store results in the list - standard for baobab.                         #
############################################################################
Results = list(ni = ni,
               nc = nc,
               nt = nt,
               frac_nN = frac_nN,

               n = sum(ni),
               ncsum = sum(nc),
               ntsum = sum(nt),

               N = sum(Ni),
               Ncsum = sum(Nc),
               Ntsum = sum(Nt),

               #Tau
               tau_true = tau_true,
               tau_treat = tau_treat,
               tau_untreat = tau_untreat,

               # OR
               E_OR = E_OR,
               M_OR = M_OR,
               R_OR = R_OR,
               X_OR = X_OR,

               # NIPW
               EE_NIPW = EE_NIPW,
               EM_NIPW = EM_NIPW,
               ER_NIPW = ER_NIPW,
               EX_NIPW = EX_NIPW,

               MM_NIPW = MM_NIPW,
               ME_NIPW = ME_NIPW,
               MR_NIPW = MR_NIPW,
               MX_NIPW = MX_NIPW,

               RR_NIPW = RR_NIPW,
               RE_NIPW = RE_NIPW,
               RM_NIPW = RM_NIPW,
               RX_NIPW = RX_NIPW,

               XX_NIPW = XX_NIPW,
               XE_NIPW = XE_NIPW,
               XM_NIPW = XM_NIPW,
               XR_NIPW = XR_NIPW,

               # AIPW
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

#ddd <- Sys.time()

outputName = paste("DB25", a, ".RData", sep = "")
outputPath = file.path("/home/reluga/Comp", outputName)
#outputPath = file.path("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE",outputName)
save("Results", file = outputPath)

#c = Sys.time()
