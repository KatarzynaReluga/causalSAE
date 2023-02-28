# rm(list=ls())
# .rs.restartR()
#

library(dplyr)
library(sampling)
library(xtable)
library(ggplot2)
library(grid)
library(gtable)
library(skewt)

setwd("./causalSAE")
devtools::load_all()

# Read and pre_process data
data_survey <- read.csv2("./simulations/data_filtered_survey.csv")
data_survey_inc <- data_survey[, -c(17)]
data_survey_inc <- data_survey_inc %>% filter(eq_house_disp_inc > 0)
y = log(data_survey_inc$eq_house_disp_inc)
data_survey_inc$y <- y

# Standardize continuous variables
data_survey_inc$age <- scale(data_survey_inc$age,
                             center = T, scale = T)

data_survey_inc$house_size <- scale(data_survey_inc$house_size,
                                    center = T, scale = T)

data_pop <- data_survey_inc %>%
  rename(A  = type_contract,
         group  = province) %>%
  select(-eq_house_disp_inc)
plot(density(data_pop$y))

####################################################################################
# Estimate propensity score
####################################################################################
#formula_p_score = A ~ X1 + (1|group)
formula_p_score = A ~ sex + nationality + age + house_size + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + (1|group)

# EBLUP

obj_p_score_EBLUP <- list(data_p_score = data_pop)
class(obj_p_score_EBLUP) <- "EBLUP"

ps_hat_EBLUP <-  p_score(obj_p_score = obj_p_score_EBLUP,
                         model_formula = formula_p_score)

# MQ

obj_p_score_MQ <- list(data_p_score = data_pop)
class(obj_p_score_MQ) <- "MQ"

ps_hat_MQ <-  p_score(obj_p_score = obj_p_score_MQ,
                      model_formula = formula_p_score)

# RF

obj_p_score_RF <- list(data_p_score = data_pop)
class(obj_p_score_RF) <- "RF"

ps_hat_RF <-  p_score(obj_p_score = obj_p_score_RF,
                      model_formula = formula_p_score,
                      tune_RF = FALSE)

# XGB

obj_p_score_XGB <- list(data_p_score = data_pop)
class(obj_p_score_XGB) <- "XGB"

ps_hat_XGB <-  p_score(obj_p_score = obj_p_score_XGB,
                       model_formula = formula_p_score,
                       xgboost_params = list(CV_XGB = FALSE,
                                             nfolds = 5,
                                             nrounds = 50))

# RF tune

#obj_p_score_RF <- list(data_p_score = data_pop)
#class(obj_p_score_RF) <- "RF"

ps_hat_RFt <-  p_score(obj_p_score = obj_p_score_RF,
                      model_formula = formula_p_score,
                      tune_RF = TRUE)

# XGB tune

#obj_p_score_XGB <- list(data_p_score = data_pop)
#class(obj_p_score_XGB) <- "XGB"

ps_hat_XGBt <- p_score(obj_p_score = obj_p_score_XGB,
                      model_formula = formula_p_score,
                      xgboost_params = list(CV_XGB = TRUE,
                                             nfolds = 5,
                                             nrounds = 50))

# Some checks
# EBLUP data
data_pop_EB <- data_pop
data_pop_EB$p_score <- ps_hat_EBLUP
ps_hat_EBLUP_sort <- sort(ps_hat_EBLUP, index.return = T)

#MQ data
data_pop_MQ <- data_pop
data_pop_MQ$p_score <- ps_hat_MQ
ps_hat_MQ_sort <- sort(ps_hat_MQ, index.return = T)

#RF data
data_pop_RF <- data_pop
data_pop_RF$p_score <- ps_hat_RF
ps_hat_RF_sort <- sort(unlist(ps_hat_RF), index.return = T)

#XGBoost data
data_pop_XGB <- data_pop
data_pop_XGB$p_score <- ps_hat_XGB
ps_hat_XGB_sort <- sort(unlist(ps_hat_XGB), index.return = T)

#RF tune data
data_pop_RFt <- data_pop
data_pop_RFt$p_score <- ps_hat_RFt
ps_hat_RFt_sort <- sort(unlist(ps_hat_RFt), index.return = T)

#XGBoost tune data
data_pop_XGBt <- data_pop
data_pop_XGBt$p_score <- ps_hat_XGBt
ps_hat_XGBt_sort <- sort(unlist(ps_hat_XGBt), index.return = T)


# Compute tau true
tau_EB <- calculate_tau(list(data_pop_EB), type_tau = "H")
tau_MQ <- calculate_tau(list(data_pop_MQ), type_tau = "H")
tau_RF <- calculate_tau(list(data_pop_RF), type_tau = "H")
tau_XGB <- calculate_tau(list(data_pop_XGB), type_tau = "H")
tau_RFt <- calculate_tau(list(data_pop_RFt), type_tau = "H")
tau_XGBt <- calculate_tau(list(data_pop_XGBt), type_tau = "H")


# Load provinces
source("./simulations/Province.R")

plot(1:41, tau_EB[[1]]$tau, col = 1, main = "Heterogeneous treatment effect",
     type = "l", ylab  = expression(tau), xlab = "Province",
     lwd = 2)
lines(tau_MQ[[1]]$tau, col = 3, lwd = 2)
lines(tau_RF[[1]]$tau, col = 2, lwd = 2)
lines(tau_XGB[[1]]$tau, col = 4, lwd = 2)
lines(tau_RFt[[1]]$tau, col = 5, lwd = 2)
lines(tau_XGBt[[1]]$tau, col = 6, lwd = 2)
legend(28, 1.5, legend=c("EBLUP", "RF"), col = 1:2, lwd = 2, cex=0.8)

df <- data.frame(effect = c(tau_EB[[1]]$tau, tau_RF[[1]]$tau),
                 method = c(rep("EBLUP", 41), rep("RF", 41)),
                 province = c(Province, Province))

################################################
## Formulas
###################################################################
# Drop: single, edu0
#
#formula_y = y ~ sex + nationality + age + house_size  + married +
#  separated + widowed + divorced  + edu1 + edu2 + edu3 + A + (1 + A||group)

formula_y_OR = y ~ sex + nationality + age + house_size  + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + (1|group)

formula_y_impute = y ~ sex + nationality + age + house_size  + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + (1 + A||group)

formula_p_score = A ~ sex + nationality + age + house_size + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + (1|group)

####################################################################################
## Design based simulations
####################################################################################

Ni = as.numeric(table(data_pop$group))
Nc = as.numeric(table(data_pop$group[data_pop$A == 0]))
Nt = as.numeric(table(data_pop$group[data_pop$A == 1]))

N = sum(Ni)
m = length(unique(data_pop$group))

a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(a * 2022)
# Here remember to change sth not to repeat that in simulations

data_pop <- data_pop_EB
subpopulation <- sample_subpopulations(data_pop,
                                       frac_nc = 0.1, frac_nt = 0.1,
                                       seed = set.seed(a * 2022))
data_sample <- data.frame(data_pop[subpopulation, ])
data_out_of_sample <- data_pop[-subpopulation, ]

#######################################################################################################
######
# OR #
######
#######################################################################################################
# EBLUP OR ----------------------------------------------------------------------------------------------
E_ORf <- hte(type_hte = "OR",
             data_sample,
             data_out_of_sample,
             params_OR = list(model_formula = formula_y_OR,
                              method = "EBLUP",
                              type_model = "gaussian"))
E_OR <- E_ORf$tau

# MQ OR --------------------------------------------------------------------------------------------------------
M_ORf <- hte(type_hte = "OR",
             data_sample,
             data_out_of_sample,
             params_OR = list(model_formula = formula_y_OR,
                              method = "MQ",
                              type_model = "continuous"))
M_OR <- M_ORf$tau

# RF OR ------------------------------------------------------------------------------------------------------------

R_ORf <- hte(type_hte = "OR",
             data_sample,
             data_out_of_sample,
             params_OR = list(model_formula = formula_y_OR,
                              method = "RF",
                              tune_RF = FALSE))
R_OR <- R_ORf$tau

# EBLUP XGB -----------------------------------------------------------------------------------------------------------

X_ORf <- hte(type_hte = "OR",
             data_sample,
             data_out_of_sample,
             params_OR = list(model_formula = formula_y_OR,
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
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "EBLUP"),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "EBLUP",
                                       type_model = "gaussian"),
                estimated_p_score = c(data_sample$p_score, data_out_of_sample$p_score))

EE_NIPW <- EE_NIPWf$tau

#EM ------------------------------------------------------------------------------------------------
EM_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "EBLUP"),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "MQ",
                                       type_model = "continuous"))

EM_NIPW <- EM_NIPWf$tau

# ER ------------------------------------------------------------------------------------------------
ER_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "EBLUP"),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "RF",
                                       tune_RF = FALSE))

ER_NIPW <- ER_NIPWf$tau

#EX --------------------------------------------------------------------------------------------
EX_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "EBLUP"),
                params_impute_y = list(model_formula = formula_y_impute,
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
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "MQ"),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "MQ",
                                       type_model = "continuous"))
MM_NIPW <- MM_NIPWf$tau


# ME --------------------------------------------------------
ME_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "MQ"),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "EBLUP",
                                       type_model = "gaussian"))
ME_NIPW <- ME_NIPWf$tau


# MR --------------------------------------------------------
MR_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "MQ"),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "RF",
                                       tune_RF = FALSE))
MR_NIPW <- MR_NIPWf$tau

# MX --------------------------------------------------------
MX_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "MQ"),
                params_impute_y = list(model_formula = formula_y_impute,
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
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "RF",
                                      tune_RF = FALSE),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "RF",
                                       tune_RF = FALSE))
RR_NIPW <- RR_NIPWf$tau
# RE ------------------------------------------------------------------
RE_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "RF",
                                      tune_RF = FALSE),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "EBLUP",
                                       type_model = "gaussian"))
RE_NIPW <- RE_NIPWf$tau

# RM ------------------------------------------------------------------
RM_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "RF",
                                      tune_RF = FALSE),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "MQ",
                                       type_model = "continuous"))
RM_NIPW <- RM_NIPWf$tau
# RX ------------------------------------------------------------------
RX_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "RF",
                                      tune_RF = FALSE),
                params_impute_y = list(model_formula = formula_y_impute,
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
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "XGB",
                                      xgboost_params = list(CV_XGB = FALSE,
                                                            nfolds = 5,
                                                            nrounds = 50)),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)))

XX_NIPW <- XX_NIPWf$tau

# XE -------------------------------------------------------------------
XE_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "XGB",
                                      xgboost_params = list(CV_XGB = FALSE,
                                                            nfolds = 5,
                                                            nrounds = 50)),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "EBLUP",
                                       type_model = "gaussian"))

XE_NIPW <- XE_NIPWf$tau
# XM -------------------------------------------------------------------
XM_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "XGB",
                                      xgboost_params = list(CV_XGB = FALSE,
                                                            nfolds = 5,
                                                            nrounds = 50)),
                params_impute_y = list(model_formula = formula_y_impute,
                                       method = "MQ",
                                       type_model = "continuous"))

XM_NIPW <- XM_NIPWf$tau

# XR -------------------------------------------------------------------
XR_NIPWf <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_p_score = list(model_formula = formula_p_score,
                                      method =  "XGB",
                                      xgboost_params = list(CV_XGB = FALSE,
                                                            nfolds = 5,
                                                            nrounds = 50)),
                params_impute_y = list(model_formula = formula_y_impute,
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
                 params_OR = list(model_formula =  formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

EEM_AIPW <- EEM_AIPWf$tau
# EER --------------------------------------------------
EER_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

EER_AIPW <- EER_AIPWf$tau

# EEX ----------------------
EEX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

EEX_AIPW <- EEX_AIPWf$tau

# EMM ----------------------------------------------------------------------
EMM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

EMM_AIPW <- EMM_AIPWf$tau

# EMR --------------------------------------------
EMR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

EMR_AIPW <- EMR_AIPWf$tau

# EMX -----------------------------------------
EMX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

EMX_AIPW <- EMX_AIPWf$tau

# ERM ------------------------------------------------------
ERM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

ERM_AIPW <- ERM_AIPWf$tau

# ERR ------------------------------------------------------------------
ERR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

ERR_AIPW <- ERR_AIPWf$tau

# ERX -------------------------------------------------------------------
ERX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

ERX_AIPW <- ERX_AIPWf$tau

# EXM --------------------------------------------------
EXM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

EXM_AIPW <- EXM_AIPWf$tau

# EXR --------------------------------------------------
EXR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

EXR_AIPW <- EXR_AIPWf$tau

# EXX ----------------------------------------------------------------
EXX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "EBLUP",
                                  type_model = "gaussian"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
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
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MEE_AIPW <- MEE_AIPWf$tau

# MER ------------------------------------------------------
MER_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

MER_AIPW <- MER_AIPWf$tau

# MER ------------------------------------------------------
MEX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

MEX_AIPW <- MEX_AIPWf$tau

# MME ------------------------------------------------------
MME_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MME_AIPW <- MME_AIPWf$tau

# MER ------------------------------------------------------
MMR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

MMR_AIPW <- MMR_AIPWf$tau

# MMX ------------------------------------------------------
MMX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

MMX_AIPW <- MMX_AIPWf$tau


# MRE ------------------------------------------------------
MRE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MRE_AIPW <- MRE_AIPWf$tau

# MRR ------------------------------------------------------
MRR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

MRR_AIPW <- MRR_AIPWf$tau

# MRX ------------------------------------------------------
MRX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

MRX_AIPW <- MRX_AIPWf$tau

# MXE ------------------------------------------------------
MXE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

MXE_AIPW <- MXE_AIPWf$tau

# MRR ------------------------------------------------------
MXR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

MXR_AIPW <- MXR_AIPWf$tau

# MRX ------------------------------------------------------
MXX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "MQ",
                                  type_model = "continuous"),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
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
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

REE_AIPW <- REE_AIPWf$tau

# REM ------------------------------------------------------
REM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

REM_AIPW <- REM_AIPWf$tau
# REX ------------------------------------------------------
REX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

REX_AIPW <- REX_AIPWf$tau


# RME ------------------------------------------------------
RME_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

RME_AIPW <- RME_AIPWf$tau

# RMM ------------------------------------------------------
RMM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

RMM_AIPW <- RMM_AIPWf$tau
# RMX ------------------------------------------------------
RMX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

RMX_AIPW <- RMX_AIPWf$tau

# RRE ------------------------------------------------------
RRE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

RRE_AIPW <- RRE_AIPWf$tau

# RRM ------------------------------------------------------
RRM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

RRM_AIPW <- RRM_AIPWf$tau
# RRX ------------------------------------------------------
RRX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "XGB",
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50)))

RRX_AIPW <- RRX_AIPWf$tau

# RXE ------------------------------------------------------
RXE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

RXE_AIPW <- RXE_AIPWf$tau

# RXM ------------------------------------------------------
RXM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

RXM_AIPW <- RXM_AIPWf$tau

# RXX ------------------------------------------------------
RXX_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = FALSE),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
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
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XEE_AIPW <- XEE_AIPWf$tau

# XEM ------------------------------------------------------
XEM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

XEM_AIPW <- XEM_AIPWf$tau

# XER ------------------------------------------------------
XER_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "EBLUP"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

XER_AIPW <- XER_AIPWf$tau

# XME ------------------------------------------------------
XME_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XME_AIPW <- XME_AIPWf$tau

# XMM ------------------------------------------------------
XMM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

XMM_AIPW <- XMM_AIPWf$tau


# XER ------------------------------------------------------
XMR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "MQ"),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

XMR_AIPW <- XMR_AIPWf$tau

# XRE ------------------------------------------------------
XRE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XRE_AIPW <- XRE_AIPWf$tau

# XRM ------------------------------------------------------
XRM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

XRM_AIPW <- XRM_AIPWf$tau


# XRR ------------------------------------------------------
XRR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "RF",
                                       tune_RF = FALSE),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "RF",
                                        tune_RF = FALSE))

XRR_AIPW <- XRR_AIPWf$tau

# XXE ------------------------------------------------------
XXE_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        type_model = "gaussian"))

XXE_AIPW <- XXE_AIPWf$tau

# XRM ------------------------------------------------------
XXM_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "MQ",
                                        type_model = "continuous"))

XXM_AIPW <- XXM_AIPWf$tau


# XER ------------------------------------------------------
XXR_AIPWf <- hte(type_hte = "AIPW",
                 data_sample,
                 data_out_of_sample,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "XGB",
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 50)),
                 params_p_score = list(model_formula = formula_p_score,
                                       method =  "XGB",
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
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
Results = list(tau_true = tau_true,
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

outputName = paste("LMMerrorssim", a, ".RData",sep="")
outputPath = file.path("/home/reluga/Comp", outputName)
#outputPath = file.path("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE",outputName)
save("Results", file = outputPath)
