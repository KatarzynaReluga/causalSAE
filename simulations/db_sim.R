library(dplyr)
#library(sampling)
#library(xtable)
#library(ggplot2)
#library(grid)
#library(gtable)
#library(skewt)

# Set path
# setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2021_2022/Causality/Comparative_HTE/EU-SILC_Censimento_SAE_Causal_Inference/Censimento/Filtered")

# Read and pre_process data
data_survey <- read.csv2("data_filtered_survey.csv")
#data_survey <- read.csv2("./data_filtered_survey.csv")
#setwd("./causalSAE")

data_survey_inc <- data_survey[, -c(17)]
data_survey_inc <- data_survey_inc %>% filter(eq_house_disp_inc > 0)
y = log(data_survey_inc$eq_house_disp_inc)
data_survey_inc$y <- y

# Standardize count variables
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
# Estimate the propensity score
####################################################################################
# (parametric estimation!!!!)
# Drop: single, edu0

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

######################################################################
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

# Very different depending on the method
# To the presentation
plot(1:length(ps_hat_EBLUP), ps_hat_EBLUP_sort$x, type = "l")
lines(ps_hat_MQ[ps_hat_EBLUP_sort$ix], col = 2)
lines(unlist(ps_hat_RF)[ps_hat_EBLUP_sort$ix], col = 3)
lines(ps_hat_XGB[ps_hat_EBLUP_sort$ix], col = 4)

# Compute tau true
tau_EB <- calculate_tau(list(data_pop_EB), type_tau = "H")
tau_MQ <- calculate_tau(list(data_pop_MQ), type_tau = "H")
tau_RF <- calculate_tau(list(data_pop_RF), type_tau = "H")
tau_XGB <- calculate_tau(list(data_pop_XGB), type_tau = "H")

source("Province.R")

plot(1:41, tau_EB[[1]]$tau, col = 1, main = "Heterogeneous treatment effect",
     type = "l", ylab  = expression(tau), xlab = "Province",
     lwd = 2)
lines(tau_RF[[1]]$tau, col = 2, lwd = 2)
lines(tau_MQ[[1]]$tau, col = 3, lwd = 2)
lines(tau_XGB[[1]]$tau, col = 4, lwd = 2)
legend(28, 1.5, legend = c("EBLUP", "RF"), col = 1:2, lwd = 2, cex=0.8)
df <- data.frame(effect = c(tau_EB[[1]]$tau, tau_RF[[1]]$tau),
                 method = c(rep("EBLUP", 41), rep("RF", 41)),
                 province = c(Province, Province))

####################################################################################
## Design based simulations
####################################################################################

Ni = as.numeric(table(data_pop$group))
Nc = as.numeric(table(data_pop$group[data_pop$A == 0]))
Nt = as.numeric(table(data_pop$group[data_pop$A == 1]))
#Nt/Ni
N = sum(Ni)
m = length(unique(data_pop$group))

NoSim <- 250

EBLUP_OR <- matrix(NA, NoSim, m)
MQ_OR <- matrix(NA, NoSim, m)
RF_OR <- matrix(NA, NoSim, m)
XGB_OR <- matrix(NA, NoSim, m)

EBLUP_NIPW <- matrix(NA, NoSim, m)
MQ_NIPW <- matrix(NA, NoSim, m)
RF_NIPW <- matrix(NA, NoSim, m)
XGB_NIPW <- matrix(NA, NoSim, m)

EBLUP_AIPW <- matrix(NA, NoSim, m)
MQ_AIPW <- matrix(NA, NoSim, m)

Dir_EBLUP_tau <- matrix(NA, NoSim, m)
Dir_MQ_tau <- matrix(NA, NoSim, m)
Dir_RF_tau <- matrix(NA, NoSim, m)
Dir_XGB_tau <- matrix(NA, NoSim, m)

nic <- round(0.1 * Nc)
nic[which(nic == 0)] <- 1
nit <- round(0.1 * Nt)
nit[which(nit == 0)] <- 1
ni = nic + nit
# round(ni/Ni, digits = 3)
#[1] 0.101 0.098 0.143 0.099 0.104 0.100 0.100 0.114 0.222 0.102 0.098 0.101 0.108 0.099 0.096 0.089 0.098
#[18] 0.109 0.101 0.099 0.117 0.113 0.101 0.097 0.102 0.125 0.100 0.100 0.097 0.118 0.100 0.122 0.116 0.122
#[35] 0.116 0.099 0.106 0.120 0.115 0.106 0.125
n <- sum(ni)

# Drop: single, edu0
#
formula_y = y ~ sex + nationality + age + house_size  + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + A + (1 + A||group)

formula_y_OR = y ~ sex + nationality + age + house_size  + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + (1|group)

formula_y_impute = y ~ sex + nationality + age + house_size  + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + (1 + A||group)

formula_p_score = A ~ sex + nationality + age + house_size + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + (1|group)

# 6 8 12 14 19
# 20 22 46 47 48
# 49 53 55 56 60
# 65 71 75 78 84
# 93 97 109 112 129
# 134 136 137 138 139
# 140 146 147 150 165
# 170 181 187 190 194
# 200 204 207 210 213
# 214 217 225 228 232
# 239 242 245 246 249
seqSim <- seq(1:250)
seqSim <- seqSim[-index_NA]
for (i in seqSim) {

  set.seed(i * 2022)

  samp <- NULL
  samp_index <- sample_subpopulations(data_pop, frac_nc = 0.1, frac_nt = 0.1)
  samp <- data_pop[samp_index, ]
  non_samp <- data_pop[ - samp_index, ]

  # OR type estimators -----------------------------------------------
  # EBLUP
  hte_OR_EBLUP <- hte(type_hte = "OR",
                      data_sample = samp,
                      data_out_of_sample = non_samp,
                      params_p_estimate = NULL,
                      params_OR = list(model_formula = formula_y_OR,
                                       method = "EBLUP",
                                       tune_RF = FALSE,
                                       xgboost_params = list(CV_XGB = TRUE,
                                                             nfolds = 5,
                                                             nrounds = 50),
                                        type_model = "gaussian"))
  EBLUP_OR[i, ] <- hte_OR_EBLUP$tau

  # MQ
  hte_OR_MQ <- hte(type_hte = "OR",
                      data_sample = samp,
                      data_out_of_sample = non_samp,
                      params_p_estimate = NULL,
                      params_OR = list(model_formula = formula_y_OR,
                                       method = "MQ",
                                       tune_RF = FALSE,
                                       xgboost_params = list(CV_XGB = TRUE,
                                                             nfolds = 5,
                                                             nrounds = 50),
                                       type_model = "continuous"))
  MQ_OR[i, ] <- hte_OR_MQ$tau

  # RF
  hte_OR_RF <- hte(type_hte = "OR",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_p_estimate = NULL,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = FALSE,
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50),
                                    type_model = "continuous"))
  RF_OR[i, ] <- hte_OR_RF$tau

  # XGB
  hte_OR_XGB <- hte(type_hte = "OR",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_p_estimate = NULL,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    tune_RF = FALSE,
                                    xgboost_params = list(CV_XGB = FALSE,
                                                          nfolds = 5,
                                                          nrounds = 50),
                                    type_model = "continuous"))

  XGB_OR[i, ] <- hte_OR_XGB$tau

  # IPW type estimators -----------------------------------------------
  # EBLUP

  hte_NIPW_EBLUP <- hte(type_hte = "NIPW",
                 data_sample = samp,
                 data_out_of_sample = non_samp,
                 params_p_estimate = NULL,
                 params_p_score = list(model_formula =  formula_p_score,
                                       method =  "EBLUP",
                                       tune_RF = FALSE,
                                       xgboost_params = list(CV_XGB = FALSE,
                                                             nfolds = 5,
                                                             nrounds = 50)),
                 params_impute_y = list(model_formula = formula_y_impute,
                                        method = "EBLUP",
                                        tune_RF = FALSE,
                                        xgboost_params = list(CV_XGB = FALSE,
                                                              nfolds = 5,
                                                              nrounds = 50),
                                        type_model = "gaussian"))

   EBLUP_NIPW[i, ] <- hte_NIPW_EBLUP$tau

   # MQ NIPW
   hte_NIPW_MQ <- hte(type_hte = "NIPW",
                         data_sample = samp,
                         data_out_of_sample = non_samp,
                         params_p_estimate = NULL,
                         params_p_score = list(model_formula =  formula_p_score,
                                               method =  "MQ",
                                               tune_RF = FALSE,
                                               xgboost_params = list(CV_XGB = FALSE,
                                                                     nfolds = 5,
                                                                     nrounds = 50)),
                         params_impute_y = list(model_formula = formula_y_impute,
                                                method = "MQ",
                                                tune_RF = FALSE,
                                                xgboost_params = list(CV_XGB = FALSE,
                                                                      nfolds = 5,
                                                                      nrounds = 50),
                                                type_model = "continuous"))

   MQ_NIPW[i, ] <- hte_NIPW_MQ$tau

   # RF NIPW
   hte_NIPW_RF <- hte(type_hte = "NIPW",
                      data_sample = samp,
                      data_out_of_sample = non_samp,
                      params_p_estimate = NULL,
                      params_p_score = list(model_formula =  formula_p_score,
                                            method =  "RF",
                                            tune_RF = FALSE,
                                            xgboost_params = list(CV_XGB = FALSE,
                                                                  nfolds = 5,
                                                                  nrounds = 50)),
                      params_impute_y = list(model_formula = formula_y_impute,
                                             method = "RF",
                                             tune_RF = FALSE,
                                             xgboost_params = list(CV_XGB = FALSE,
                                                                   nfolds = 5,
                                                                   nrounds = 50),
                                             type_model = "continuous"))

   RF_NIPW[i, ] <- hte_NIPW_RF$tau

   # XGB NIPW
   hte_NIPW_XGB <- hte(type_hte = "NIPW",
                      data_sample = samp,
                      data_out_of_sample = non_samp,
                      params_p_estimate = NULL,
                      params_p_score = list(model_formula =  formula_p_score,
                                            method =  "XGB",
                                            tune_RF = FALSE,
                                            xgboost_params = list(CV_XGB = FALSE,
                                                                  nfolds = 5,
                                                                  nrounds = 50)),
                      params_impute_y = list(model_formula = formula_y_impute,
                                             method = "XGB",
                                             tune_RF = FALSE,
                                             xgboost_params = list(CV_XGB = FALSE,
                                                                   nfolds = 5,
                                                                   nrounds = 50),
                                             type_model = "continuous"))

   XGB_NIPW[i, ] <- hte_NIPW_XGB$tau

   # AIPW type estimators -----------------------------------------------
   # EBLUP

   hte_AIPW_EBLUP <- hte(type_hte = "AIPW",
                         data_sample = samp,
                         data_out_of_sample = non_samp,
                         params_p_estimate = NULL,
                         params_p_score = list(model_formula =  formula_p_score,
                                               method =  "EBLUP",
                                               tune_RF = FALSE,
                                               xgboost_params = list(CV_XGB = FALSE,
                                                                     nfolds = 5,
                                                                     nrounds = 50)),
                         params_impute_y = list(model_formula = formula_y_impute,
                                                method = "RF",
                                                tune_RF = TRUE,
                                                xgboost_params = list(CV_XGB = FALSE,
                                                                      nfolds = 5,
                                                                      nrounds = 50),
                                                type_model = "continuous"),
                        params_OR = list(model_formula = formula_y_OR,
                                             method = "EBLUP",
                                             tune_RF = FALSE,
                                             xgboost_params = list(CV_XGB = FALSE,
                                                                   nfolds = 5,
                                                                   nrounds = 50),
                                             type_model = "gaussian"))

   EBLUP_AIPW[i, ] <- hte_AIPW_EBLUP$tau


   hte_AIPW_MQ <- hte(type_hte = "AIPW",
                         data_sample = samp,
                         data_out_of_sample = non_samp,
                         params_p_estimate = NULL,
                         params_p_score = list(model_formula =  formula_p_score,
                                               method =  "EBLUP",
                                               tune_RF = FALSE,
                                               xgboost_params = list(CV_XGB = FALSE,
                                                                     nfolds = 5,
                                                                     nrounds = 50)),
                         params_impute_y = list(model_formula = formula_y_impute,
                                                method = "RF",
                                                tune_RF = FALSE,
                                                xgboost_params = list(CV_XGB = FALSE,
                                                                      nfolds = 5,
                                                                      nrounds = 50),
                                                type_model = "continuous"),
                         params_OR = list(model_formula = formula_y_OR,
                                          method = "EBLUP",
                                          tune_RF = FALSE,
                                          xgboost_params = list(CV_XGB = FALSE,
                                                                nfolds = 5,
                                                                nrounds = 50),
                                          type_model = "gaussian"))

   MQ_AIPW[i, ] <- hte_AIPW_MQ$tau

  # Direct EBLUP

  samp_EBLUP <- data_pop_EB[samp_index, ]
  Dir_EBLUP_tau[i, ] <- (calculate_tau(list(samp_EBLUP),
                              type_tau = "Hajek"))[[1]]$tau

  # Direct MQ
  samp_MQ <- data_pop_MQ[samp_index, ]
  Dir_MQ_tau[i, ] <- (calculate_tau(list(samp_MQ),
                                       type_tau = "Hajek"))[[1]]$tau

  # Direct RF
  samp_RF <- data_pop_RF[samp_index, ]
  Dir_RF_tau[i, ] <- (calculate_tau(list(samp_RF),
                              type_tau = "Hajek"))[[1]]$tau

  # Direct XGB
  samp_XGB <- data_pop_XGB[samp_index, ]
  Dir_XGB_tau[i, ] <- (calculate_tau(list(samp_XGB),
                                    type_tau = "Hajek"))[[1]]$tau

}
#########################################################################################
## Save the simulations
#########################################################################################
index_NA <- c(6, 8, 12, 14, 19, 20, 22, 46, 47, 48,
              49, 53, 55, 56, 60, 65, 71, 75, 78, 84,
              93, 97, 109, 112, 129, 134, 136, 137, 138, 139,
              140, 146, 147, 150, 165, 170, 181, 187, 190, 194,
              200, 204, 207, 210, 213, 214, 217, 225, 228, 232,
              239, 242, 245, 246, 249)

EBLUP_ORf = EBLUP_OR[-index_NA, ]
MQ_ORf = MQ_OR[-index_NA, ]
RF_ORf = RF_OR[-index_NA, ]
XGB_ORf = XGB_OR[-index_NA, ]

EBLUP_NIPWf = EBLUP_NIPW[-index_NA, ]
MQ_NIPWf = MQ_NIPW[-index_NA, ]
RF_NIPWf = RF_NIPW[-index_NA, ]
XGB_NIPWf = XGB_NIPW[-index_NA, ]

EBLUP_AIPWf = EBLUP_AIPW[-index_NA, ]
MQ_AIPWf = MQ_AIPW[-index_NA, ]

Dir_EBLUP_tauf = Dir_EBLUP_tau[-index_NA, ]
Dir_MQ_tauf = Dir_MQ_tau[-index_NA, ]
Dir_RF_tauf = Dir_RF_tau[-index_NA, ]
Dir_XGB_tauf = Dir_XGB_tau[-index_NA, ]


# Save as RData
#save(data_pop_EB, data_pop_MQ,
#     data_pop_RF, data_pop_XGB,

#     EBLUP_ORf, MQ_ORf,
#     RF_ORf, XGB_ORf,

#     EBLUP_NIPWf, MQ_NIPWf,
#     RF_NIPWf, XGB_NIPWf,

#     EBLUP_AIPWf, MQ_AIPWf,
#     Dir_EBLUP_tauf, Dir_MQ_tauf,
#     Dir_RF_tauf, Dir_XGB_tauf,
#     file = "Design_based_simulations32.RData")
# Save as csf
#write.csv2(data_pop_EB, "data_pop_EB_sim32.csv", row.names = F)
#write.csv2(data_pop_MQ, "data_pop_MQ_sim32.csv", row.names = F)
#write.csv2(data_pop_RF, "data_pop_RF_sim32.csv", row.names = F)
#write.csv2(data_pop_XGB, "data_pop_XGB_sim32.csv", row.names = F)

load("Design_based_sim1cor2.RData")

#write.csv2(EBLUP_ORf, "EBLUP_tau_sim32.csv", row.names = F)
#write.csv2(MQ_ORf, "MQ_tau_sim32.csv", row.names = F)
#write.csv2(RF_ORf, "RF_tau_sim32.csv", row.names = F)
#write.csv2(XGB_ORf, "XGB_tau_sim32.csv", row.names = F)

#write.csv2(EBLUP_NIPWf, "EBLUP_NIPW_sim32.csv", row.names = F)
#write.csv2(MQ_NIPWf, "MQ_NIPW_sim32.csv", row.names = F)
#write.csv2(RF_NIPWf, "RF_NIPW_sim32.csv", row.names = F)
#write.csv2(XGB_NIPWf, "XGB_NIPW_sim32.csv", row.names = F)

#write.csv2(EBLUP_AIPWf, "EBLUP_AIPW_sim32.csv", row.names = F)
#write.csv2(MQ_AIPWf, "MQ_AIPW_sim32.csv", row.names = F)

#write.csv2(Dir_EBLUP_tauf, "Dir_EBLUP_sim32.csv", row.names = F)
#write.csv2(Dir_MQ_tauf, "Dir_MQ_sim32.csv", row.names = F)
#write.csv2(Dir_RF_tauf, "Dir_RF_sim32.csv", row.names = F)
#write.csv2(Dir_XGB_tauf, "Dir_XGB_sim32.csv", row.names = F)

plot(1:41, tau_EB[[1]]$tau, col = 1, lwd = 2.5,
     type = "l", ylab = expression(tau), ylim = c(-1, 1.5))
lines(tau_RF[[1]]$tau, col = "grey", lwd = 2.5)

lines(colMeans(EBLUP_ORf), col = 2, lwd = 2)
lines(colMeans(MQ_ORf), col = 3, lwd = 2)
lines(colMeans(RF_ORf), col = 4, lwd = 2)
lines(colMeans(XGB_ORf), col = 5, lwd = 2)


plot(1:41, tau_EB[[1]]$tau, col = 1, lwd = 2.5,
     type = "l", ylab = expression(tau), ylim = c(-1, 1.5))
lines(tau_RF[[1]]$tau, col = "grey", lwd = 2.5)

lines(colMeans(EBLUP_NIPWf), col = 4, lwd = 2)
lines(colMeans(MQ_NIPWf), col = 3, lwd = 2)
lines(colMeans(RF_NIPWf), col = 4, lwd = 2)
lines(colMeans(XGB_NIPWf), col = 5, lwd = 2)


plot(1:41, tau_EB[[1]]$tau, col = 1, lwd = 2.5,
     type = "l", ylab = expression(tau), ylim = c(-1, 1.5))
lines(tau_RF[[1]]$tau, col = "grey", lwd = 2.5)

lines(colMeans(EBLUP_AIPWf), col = 2, lwd = 2)
lines(colMeans(MQ_AIPWf), col = 3, lwd = 2)

plot(1:41, tau_EB[[1]]$tau, col = 1, lwd = 2.5,
     type = "l", ylab = expression(tau), ylim = c(-1, 1.5))
lines(tau_RF[[1]]$tau, col = "grey", lwd = 2.5)
lines(colMeans(Dir_EBLUP_tauf), col = 6, lwd = 2)
lines(colMeans(Dir_MQ_tauf), col = 7, lwd = 2)
lines(colMeans(Dir_RF_tauf), col = 3, lwd = 2)
lines(colMeans(Dir_XGB_tauf), col = 4, lwd = 2)

legend(33, 1.44,
       legend=c("Tau_EBLUP", "Tau_RF",
                "EBLUP", "MQ", "RF", "RFt", "Dir_EBLUP", "Dir_RF"),
       col = c(1, "grey", 1:7), lwd = 2, cex=0.8)

############################################################################
# RBias
#############################################################################
RBias_percent <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RB.EBLUP <- (apply(EBLUP_tau, 2, mean) - True) / abs(True) * 100
  value <- (apply(estimator, 2, mean) - true_value) / abs(true_value) * 100
  value

}

RBias <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RB.EBLUP <- (apply(EBLUP_tau, 2, mean) - True) / abs(True) * 100
  value <- (apply(estimator, 2, mean) - true_value) / abs(true_value)
  value

}

#EBLUP_ORf
#MQ_ORf
#RF_ORf
#XGB_ORf

#EBLUP_NIPWf
#MQ_NIPWf
#RF_NIPWf
#XGB_NIPWf

#EBLUP_AIPWf
#MQ_AIPWf

#Dir_EBLUP_tauf
#Dir_MQ_tauf
#Dir_RF_tauf
#Dir_XGB_tauf


list_estimators <- list()

list_estimators[[1]] <- EBLUP_ORf
list_estimators[[2]] <- MQ_ORf
list_estimators[[3]] <- RF_ORf
list_estimators[[4]] <- XGB_ORf

list_estimators[[5]] <- EBLUP_NIPWf
list_estimators[[6]] <- MQ_NIPWf
list_estimators[[7]] <- RF_NIPWf
list_estimators[[8]] <- XGB_NIPWf

list_estimators[[9]] <- EBLUP_AIPWf
list_estimators[[10]] <- MQ_AIPWf

list_estimators[[11]] <- Dir_EBLUP_tauf
list_estimators[[12]] <- Dir_MQ_tauf
list_estimators[[13]] <- Dir_RF_tauf
list_estimators[[14]] <- Dir_XGB_tauf

true_EB = tau_EB[[1]]$tau
true_MQ = tau_MQ[[1]]$tau
true_RF = tau_RF[[1]]$tau
true_XGB = tau_XGB[[1]]$tau

RBper_estEB <- lapply(list_estimators, RBias_percent, true_value = true_EB)
RBper_estMQ <- lapply(list_estimators, RBias_percent, true_value = true_MQ)
RBper_estRF <- lapply(list_estimators, RBias_percent, true_value = true_RF)
RBper_estXGB <- lapply(list_estimators, RBias_percent, true_value = true_XGB)

plot(1:41, RBper_estEB[[1]], col = 1,  type = "l",
     lwd = 2, xlab = "Provinces", ylab = expression(tau))
lines(RBper_estEB[[2]], col = 2, lwd = 2)
lines(RBper_estEB[[3]], col = 3, lwd = 2)
lines(RBper_estEB[[4]], col = 4, lwd = 2)

lines(RBper_estEB[[5]], col = 5, lwd = 2)
lines(RBper_estEB[[6]], col = 6, lwd = 2)
lines(RBper_estEB[[7]], col = 5, lwd = 2)
lines(RBper_estEB[[8]], col = 6, lwd = 2)

summary(RBper_estEB[[1]])
summary(RBper_estEB[[2]])
summary(RBper_estEB[[3]])
summary(RBper_estEB[[4]])

summary(RBper_estEB[[5]])
summary(RBper_estEB[[6]])

RB_estEB <- lapply(list_estimators, RBias, true_value = true_EB)
RB_estRF <- lapply(list_estimators, RBias, true_value = true_RF)

plot(1:41, RB_estEB[[1]], col = 1, type = "l",
     lwd = 2, xlab = "Provinces", ylab = expression(tau))
lines(RB_estEB[[2]], col = 2, lwd = 2)
lines(RB_estEB[[3]], col = 3, lwd = 2)
lines(RB_estEB[[4]], col = 4, lwd = 2)
lines(RB_estEB[[5]], col = 5, lwd = 2)
lines(RB_estEB[[6]], col = 6, lwd = 2)

summary(RB_estEB[[1]])
summary(RB_estEB[[2]])
summary(RB_estEB[[3]])
summary(RB_estEB[[4]])

summary(RB_estEB[[5]])
summary(RB_estEB[[6]])
summary(RB_estEB[[7]])
summary(RB_estEB[[8]])

#####################################################################
# RRMSE
####################################################################

RBias_percent <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RB.EBLUP <- (apply(EBLUP_tau, 2, mean) - True) / abs(True) * 100
  value <- (apply(estimator, 2, mean) - true_value) / abs(true_value) * 100
  value

}

RRMSE_percent <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RRMSE.Direct <-
  #  (sqrt(apply(
  #    apply(Direct, 1, function(x) {
  #      (x - True) ^ 2
  #    }), 1, mean, na.rm = TRUE
  #  )) / abs(True)) * 100
  value <- (sqrt(apply(apply(estimator, 1, function(x) {
        (x - true_value) ^ 2
      }), 1, mean, na.rm = TRUE
    )) / abs(true_value)) * 100
  value
}

RRMSE <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RRMSE.Direct <-
  #  (sqrt(apply(
  #    apply(Direct, 1, function(x) {
  #      (x - True) ^ 2
  #    }), 1, mean, na.rm = TRUE
  #  )) / abs(True)) * 100
  value <- (sqrt(apply(apply(estimator, 1, function(x) {
    (x - true_value) ^ 2
  }), 1, mean, na.rm = TRUE
  )) / abs(true_value))
  value
}

RRMSEper_estEB <- lapply(list_estimators, RRMSE_percent, true_value = true_EB)
RRMSEper_estMQ <- lapply(list_estimators, RRMSE_percent, true_value = true_MQ)
RRMSEper_estRF <- lapply(list_estimators, RRMSE_percent, true_value = true_RF)
RRMSEper_estXGB <- lapply(list_estimators, RRMSE_percent, true_value = true_XGB)

#RRMSE_estEB <- lapply(list_estimators, RRMSE, true_value = true_EB)
#RRMSE_estRF <- lapply(list_estimators, RRMSE, true_value = true_RF)

# RBias
# EBLUP
# Dir
boxplot(
  cbind(Dir1 = RBper_estEB[[11]],
        Dir2 = RBper_estEB[[12]],
        Dir3 = RBper_estEB[[13]],
        Dir4 =  RBper_estEB[[14]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# OR
boxplot(
  cbind(OR1 = RBper_estEB[[1]],
        OR2 = RBper_estEB[[2]],
        OR3 = RBper_estEB[[3]],
        OR4 =  RBper_estEB[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# NIPW
boxplot(
  cbind(NIPW1 = RBper_estEB[[5]],
        NIPW2 = RBper_estEB[[6]],
        NIPW3 = RBper_estEB[[7]],
        NIPW4 =  RBper_estEB[[8]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

boxplot(
  cbind(AIPW1 = RBper_estEB[[9]],
        AIPW2 = RBper_estEB[[10]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# MQ
# Dir
boxplot(
  cbind(Dir1 = RBper_estEB[[11]],
        Dir2 = RBper_estEB[[12]],
        Dir3 = RBper_estEB[[13]],
        Dir4 =  RBper_estEB[[14]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# OR
boxplot(
  cbind(OR1 = RBper_estEB[[1]],
        OR2 = RBper_estEB[[2]],
        OR3 = RBper_estEB[[3]],
        OR4 =  RBper_estEB[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# NIPW
boxplot(
  cbind(NIPW1 = RBper_estEB[[5]],
        NIPW2 = RBper_estEB[[6]],
        NIPW3 = RBper_estEB[[7]],
        NIPW4 =  RBper_estEB[[8]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

boxplot(
  cbind(AIPW1 = RBper_estEB[[9]],
        AIPW2 = RBper_estEB[[10]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# RF
boxplot(
  cbind(Dir1 = RBper_estRF[[11]],
        Dir2 = RBper_estRF[[12]],
        Dir3 = RBper_estRF[[13]],
        Dir4 =  RBper_estRF[[14]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# RF
# OR
boxplot(
  cbind(OR1 = RBper_estRF[[1]],
        OR2 = RBper_estRF[[2]],
        OR3 = RBper_estRF[[3]],
        OR4 =  RBper_estRF[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# NIPW
boxplot(
  cbind(NIPW1 = RBper_estRF[[5]],
        NIPW2 = RBper_estRF[[6]],
        NIPW3 = RBper_estRF[[7]],
        NIPW4 =  RBper_estRF[[8]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

boxplot(
  cbind(AIPW1 = RBper_estRF[[9]],
        AIPW2 = RBper_estRF[[10]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)


# XGB
# OR
boxplot(
  cbind(OR1 = RBper_estXGB[[1]],
        OR2 = RBper_estXGB[[2]],
        OR3 = RBper_estXGB[[3]],
        OR4 =  RBper_estXGB[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# NIPW
boxplot(
  cbind(NIPW1 = RBper_estXGB[[5]],
        NIPW2 = RBper_estXGB[[6]],
        NIPW3 = RBper_estXGB[[7]],
        NIPW4 =  RBper_estXGB[[8]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

boxplot(
  cbind(AIPW1 = RBper_estXGB[[9]],
        AIPW2 = RBper_estXGB[[10]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)


# RMSE
pdf(file = "rb_rrmse_per_EB_cor2.pdf", width = 20)
par(mfrow = c(1, 2))
boxplot(
  cbind(Dir1 = RBper_estEB[[5]],
        Dir2 = RBper_estEB[[6]],
        EBLUP = RBper_estEB[[1]],
        MQ =  RBper_estEB[[2]],
        RF = RBper_estEB[[3]],
        RFt = RBper_estEB[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
boxplot(
  cbind(Dir1 = RRMSEper_estEB[[5]],
        Dir2 = RRMSEper_estEB[[6]],
        EBLUP = RRMSEper_estEB[[1]],
        MQ =  RRMSEper_estEB[[2]],
        RF = RRMSEper_estEB[[3]],
        RFt = RRMSEper_estEB[[4]]),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
dev.off()

pdf(file = "rb_rrmse_EB_cor2.pdf", width = 20)
par(mfrow = c(1, 2))
boxplot(
  cbind(Dir1 = RB_estEB[[5]],
        Dir2 = RB_estEB[[6]],
        EBLUP = RB_estEB[[1]],
        MQ =  RB_estEB[[2]],
        RF = RB_estEB[[3]],
        RFt = RB_estEB[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
boxplot(
  cbind(Dir1 = RRMSE_estEB[[5]],
        Dir2 = RRMSE_estEB[[6]],
        EBLUP = RRMSE_estEB[[1]],
        MQ =  RRMSE_estEB[[2]],
        RF = RRMSE_estEB[[3]],
        RFt = RRMSE_estEB[[4]]),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
dev.off()

pdf(file = "rb_rrmse_per_RF_cor2.pdf", width = 20)
par(mfrow = c(1, 2))
boxplot(
  cbind(Dir1 = RBper_estRF[[5]],
        Dir2 = RBper_estRF[[6]],
        EBLUP = RBper_estRF[[1]],
        MQ =  RBper_estRF[[2]],
        RF = RBper_estRF[[3]],
        RFt = RBper_estRF[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
boxplot(
  cbind(Dir1 = RRMSEper_estRF[[5]],
        Dir2 = RRMSEper_estRF[[6]],
        RFLUP = RRMSEper_estRF[[1]],
        MQ =  RRMSEper_estRF[[2]],
        RF = RRMSEper_estRF[[3]],
        RFt = RRMSEper_estRF[[4]]),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
dev.off()

pdf(file = "rb_rrmse_RF_cor2.pdf", width = 20)
par(mfrow = c(1, 2))
boxplot(
  cbind(Dir1 = RB_estRF[[5]],
        Dir2 = RB_estRF[[6]],
        EBLUP = RB_estRF[[1]],
        MQ =  RB_estRF[[2]],
        RF = RB_estRF[[3]],
        RFt = RB_estRF[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
boxplot(
  cbind(Dir1 = RRMSE_estRF[[5]],
        Dir2 = RRMSE_estRF[[6]],
        EBLUP = RRMSE_estRF[[1]],
        MQ =  RRMSE_estRF[[2]],
        RF = RRMSE_estRF[[3]],
        RFt = RRMSE_estRF[[4]]),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
dev.off()
######################################################################
#RMSE
# EBLUP
# Dir
boxplot(
  cbind(Dir1 = RRMSEper_estEB[[11]],
        Dir2 = RRMSEper_estEB[[12]],
        Dir3 = RRMSEper_estEB[[13]],
        Dir4 =  RRMSEper_estEB[[14]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# OR
boxplot(
  cbind(OR1 = RRMSEper_estEB[[1]],
        OR2 = RRMSEper_estEB[[2]],
        OR3 = RRMSEper_estEB[[3]],
        OR4 =  RRMSEper_estEB[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# NIPW
boxplot(
  cbind(NIPW1 = RRMSEper_estEB[[5]],
        NIPW2 = RRMSEper_estEB[[6]],
        NIPW3 = RRMSEper_estEB[[7]],
        NIPW4 =  RRMSEper_estEB[[8]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

boxplot(
  cbind(AIPW1 = RRMSEper_estEB[[9]],
        AIPW2 = RRMSEper_estEB[[10]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# MQ
# Dir
boxplot(
  cbind(Dir1 = RRMSEper_estMQ[[11]],
        Dir2 = RRMSEper_estMQ[[12]],
        Dir3 = RRMSEper_estMQ[[13]],
        Dir4 =  RRMSEper_estMQ[[14]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# OR
boxplot(
  cbind(OR1 = RRMSEper_estMQ[[1]],
        OR2 = RRMSEper_estMQ[[2]],
        OR3 = RRMSEper_estMQ[[3]],
        OR4 =  RRMSEper_estMQ[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# NIPW
boxplot(
  cbind(NIPW1 = RRMSEper_estMQ[[5]],
        NIPW2 = RRMSEper_estMQ[[6]],
        NIPW3 = RRMSEper_estMQ[[7]],
        NIPW4 =  RRMSEper_estMQ[[8]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

boxplot(
  cbind(AIPW1 = RRMSEper_estMQ[[9]],
        AIPW2 = RRMSEper_estMQ[[10]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# RF
boxplot(
  cbind(Dir1 = RRMSEper_estRF[[11]],
        Dir2 = RRMSEper_estRF[[12]],
        Dir3 = RRMSEper_estRF[[13]],
        Dir4 =  RRMSEper_estRF[[14]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# RF
# OR
boxplot(
  cbind(OR1 = RRMSEper_estRF[[1]],
        OR2 = RRMSEper_estRF[[2]],
        OR3 = RRMSEper_estRF[[3]],
        OR4 =  RRMSEper_estRF[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# NIPW
boxplot(
  cbind(NIPW1 = RRMSEper_estRF[[5]],
        NIPW2 = RRMSEper_estRF[[6]],
        NIPW3 = RRMSEper_estRF[[7]],
        NIPW4 =  RRMSEper_estRF[[8]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

boxplot(
  cbind(AIPW1 = RRMSEper_estRF[[9]],
        AIPW2 = RRMSEper_estRF[[10]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)


# XGB
# DIR
boxplot(
  cbind(Dir1 = RRMSEper_estXGB[[11]],
        Dir2 = RRMSEper_estXGB[[12]],
        Dir3 = RRMSEper_estXGB[[13]],
        Dir4 =  RRMSEper_estXGB[[14]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
# OR
boxplot(
  cbind(OR1 = RRMSEper_estXGB[[1]],
        OR2 = RRMSEper_estXGB[[2]],
        OR3 = RRMSEper_estXGB[[3]],
        OR4 =  RRMSEper_estXGB[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# NIPW
boxplot(
  cbind(NIPW1 = RRMSEper_estXGB[[5]],
        NIPW2 = RRMSEper_estXGB[[6]],
        NIPW3 = RRMSEper_estXGB[[7]],
        NIPW4 =  RRMSEper_estXGB[[8]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

boxplot(
  cbind(AIPW1 = RBper_estXGB[[9]],
        AIPW2 = RBper_estXGB[[10]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)


# RMSE
pdf(file = "rb_rrmse_per_EB_cor2.pdf", width = 20)
par(mfrow = c(1, 2))
boxplot(
  cbind(Dir1 = RBper_estEB[[5]],
        Dir2 = RBper_estEB[[6]],
        EBLUP = RBper_estEB[[1]],
        MQ =  RBper_estEB[[2]],
        RF = RBper_estEB[[3]],
        RFt = RBper_estEB[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
boxplot(
  cbind(Dir1 = RRMSEper_estEB[[5]],
        Dir2 = RRMSEper_estEB[[6]],
        EBLUP = RRMSEper_estEB[[1]],
        MQ =  RRMSEper_estEB[[2]],
        RF = RRMSEper_estEB[[3]],
        RFt = RRMSEper_estEB[[4]]),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
dev.off()

pdf(file = "rb_rrmse_EB_cor2.pdf", width = 20)
par(mfrow = c(1, 2))
boxplot(
  cbind(Dir1 = RB_estEB[[5]],
        Dir2 = RB_estEB[[6]],
        EBLUP = RB_estEB[[1]],
        MQ =  RB_estEB[[2]],
        RF = RB_estEB[[3]],
        RFt = RB_estEB[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
boxplot(
  cbind(Dir1 = RRMSE_estEB[[5]],
        Dir2 = RRMSE_estEB[[6]],
        EBLUP = RRMSE_estEB[[1]],
        MQ =  RRMSE_estEB[[2]],
        RF = RRMSE_estEB[[3]],
        RFt = RRMSE_estEB[[4]]),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
dev.off()

pdf(file = "rb_rrmse_per_RF_cor2.pdf", width = 20)
par(mfrow = c(1, 2))
boxplot(
  cbind(Dir1 = RBper_estRF[[5]],
        Dir2 = RBper_estRF[[6]],
        EBLUP = RBper_estRF[[1]],
        MQ =  RBper_estRF[[2]],
        RF = RBper_estRF[[3]],
        RFt = RBper_estRF[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
boxplot(
  cbind(Dir1 = RRMSEper_estRF[[5]],
        Dir2 = RRMSEper_estRF[[6]],
        RFLUP = RRMSEper_estRF[[1]],
        MQ =  RRMSEper_estRF[[2]],
        RF = RRMSEper_estRF[[3]],
        RFt = RRMSEper_estRF[[4]]),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
dev.off()

pdf(file = "rb_rrmse_RF_cor2.pdf", width = 20)
par(mfrow = c(1, 2))
boxplot(
  cbind(Dir1 = RB_estRF[[5]],
        Dir2 = RB_estRF[[6]],
        EBLUP = RB_estRF[[1]],
        MQ =  RB_estRF[[2]],
        RF = RB_estRF[[3]],
        RFt = RB_estRF[[4]]),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
boxplot(
  cbind(Dir1 = RRMSE_estRF[[5]],
        Dir2 = RRMSE_estRF[[6]],
        EBLUP = RRMSE_estRF[[1]],
        MQ =  RRMSE_estRF[[2]],
        RF = RRMSE_estRF[[3]],
        RFt = RRMSE_estRF[[4]]),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
dev.off()

## table of estimates and their s.e.
results <- NULL
results <-
  data.frame(
    "Province" = Province,

    "tEB" = true_EB,
    "tRF" = true_RF,

    "Dir1" = apply(Dir_EBLUP_tau, 2, mean, na.rm = TRUE),
    "seDir1" = sqrt(apply(
      apply(Dir_EBLUP_tau, 1, function(x) {
        (x - true_EB) ^ 2
      }), 1, mean, na.rm = TRUE
    )),

    "Dir2" = apply(Dir_RF_tau, 2, mean, na.rm = TRUE),
    "seDir2" = sqrt(apply(
      apply(Dir_RF_tau, 1, function(x) {
        (x - true_EB) ^ 2
      }), 1, mean, na.rm = TRUE
    )),

    "EB" = apply(EBLUP_tau, 2, mean, na.rm = TRUE),
    "seEB" = sqrt(apply(
      apply(EBLUP_tau, 1, function(x) {
        (x - true_EB) ^ 2
      }), 1, mean, na.rm = TRUE
    )),

    "MQ" = apply(MQ_tau, 2, mean, na.rm = TRUE),
    "seMQ" = sqrt(apply(
      apply(MQ_tau, 1, function(x) {
        (x - true_EB) ^ 2
      }), 1, mean, na.rm = TRUE
    )),

    "RF" = apply(RF_tau, 2, mean, na.rm = TRUE),
    "seRF" = sqrt(apply(
      apply(RF_tau, 1, function(x) {
        (x - true_EB) ^ 2
      }), 1, mean, na.rm = TRUE
    )),

    "RFt" = apply(RFt_tau, 2, mean, na.rm = TRUE),
    "seRF" = sqrt(apply(
      apply(RFt_tau, 1, function(x) {
        (x - true_EB) ^ 2
      }), 1, mean, na.rm = TRUE
    ))
  )
est.results <- xtable(results, digits = 2)
print(est.results, include.rownames = FALSE)
##############################################################

## table of 95% confidence intervals of the estimates
results_CI <- NULL
CI_Dir1 <-
  apply(Dir_EBLUP_tau, 2, function(x)
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
CI_Dir2 <-
  apply(Dir_RF_tau, 2, function(x)
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
CI_EBLUP <-
  apply(EBLUP_tau, 2, function(x)
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
CI_MQ <-
  apply(MQ_tau, 2, function(x)
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
CI_RF <-
  apply(RF_tau, 2, function(x)
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
CI_RFt <-
  apply(RFt_tau, 2, function(x)
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
results_CI <-
  data.frame(
    "Province" = Province,
    "tEB" = true_EB,
    "tRF" = true_RF,
    "LL_D1" = CI_Dir1[1, ],
    "UL_D1" = CI_Dir1[2, ],

    "LL_D2" = CI_Dir2[1, ],
    "UL_D2" = CI_Dir2[2, ],

    "LL_E" = CI_EBLUP[1, ],
    "UL_E" = CI_EBLUP[2, ],

    "LL_MQ" = CI_MQ[1, ],
    "UL_MQ" = CI_MQ[2, ],

    "LL_RF" = CI_RF[1, ],
    "UL_RF" = CI_RF[2, ],

    "LL_RFt" = CI_RFt[1, ],
    "UL_RFt" = CI_RFt[2, ]

  )
est.results_CI <- xtable(results_CI, digits = 2)
print(est.results_CI, include.rownames = FALSE)



pdf(file = "CI_compare2_cor2.pdf",
    width = 20,
    height = 10)
ggplot(data = results_CI) + theme_bw()+
  geom_hline(yintercept = 0,
             color = "red",
             lty = "dashed",
             lwd = 1.1) +
  scale_color_manual(name = " ", values = c("Dir1" = "chartreuse4",
                                            "Dir2" = "chartreuse3",
                                            "EBLUP" = "darkred",
                                            "MQ" = "blue",
                                            "RF" = "green",
                                            "RFt" = "magenta",
                                            "tau1" = "grey40",
                                            "tau2" = "black")) +
  geom_errorbar(mapping = aes(x = Province, ymin = LL_D1, ymax = UL_D1, col = "Dir1"),
                lwd = 1.1, width = 0.75) +
  geom_errorbar(mapping = aes(x = Province, ymin = LL_D2, ymax = UL_D2, col = "Dir2"),
                lwd = 1.1, width = 0.75) +
  geom_errorbar(mapping = aes(x = Province, ymin = LL_E, ymax = UL_E,  color = "EBLUP"),
                lwd  = 1.1, width = 0.75) +
  geom_errorbar(mapping = aes(x = Province, ymin = LL_MQ, ymax = UL_MQ, color = "MQ"),
                lwd  = 1.1, width = 0.75) +
  geom_errorbar(mapping = aes(x = Province, ymin = LL_RF, ymax = UL_RF,  color = "RF"),
               lwd  = 1.1, width = 0.75) +
  geom_errorbar(mapping = aes(x = Province, ymin = LL_RFt, ymax = UL_RFt, color = "RFt"),
                lwd  = 1.1, width = 0.75) +
  geom_point(mapping = aes(x = Province, y = true_EB, col = "tau1"), size = 3) +
  geom_point(mapping = aes(x = Province, y = true_RF, col = "tau2"), size = 3) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.y = element_text(size = rel(1.5), angle = 00)) +
  ylab("Treatment effects") + theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 16, face = "bold")) +
#  theme(legend.title = element_blank())
  theme(axis.title.y = element_text(size = 24)) +
  theme(axis.title.x = element_text(size = 24)) +
  theme(axis.title.x = element_text(vjust = 1))  +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.title = element_blank(), legend.position = c(0.935, 0.85))
dev.off()

pdf(file = "App_ATE_log_cor2.pdf", width = 20)
res <- NULL
res <-
  data.frame(
    "True_EBLUP" = true_EB,
    "True_RF" = true_RF,
    "Dir1" = apply(Dir_EBLUP_tau, 2, mean),
    "Dir2" = apply(Dir_RF_tau, 2, mean),
    "EBLUP" = apply(EBLUP_tau, 2, mean),
    "MQ" = apply(MQ_tau, 2, mean),
    "RF" = apply(RF_tau, 2, mean),
    "RFt" = apply(RFt_tau, 2, mean))
boxplot(res,
        cex.axis = 1.5,
        font.lab = 2,
        lwd = 1.5)
dev.off()


# Capturing the heterogeniety of the effects
True_EBLUP = true_EB
True_RF = true_RF
Dir1 = apply(Dir_EBLUP_tau, 2, mean)
Dir2 = apply(Dir_RF_tau, 2, mean)
EBLUP = apply(EBLUP_tau, 2, mean)
MQ = apply(MQ_tau, 2, mean)
RF = apply(RF_tau, 2, mean)
RFt = apply(RFt_tau, 2, mean)

min.x <- min(c(True_EBLUP, True_RF,
               Dir1, Dir2,
               EBLUP, MQ,
               RF, RFt))

max.x <- max(c(True_EBLUP, True_RF,
               Dir1, Dir2,
               EBLUP, MQ,
               RF, RFt))

max.y <- 3.7
par(mfrow = c(1, 1))
pdf(file = "hetero-est-log2.pdf", width = 15)
plot(
  density(True_EBLUP),
  xlim = c(-1.1, max.x),
  ylim = c(0, max.y),
  lty = 2,
  col = "black",
  xlab = "Group-specific effect",
  main = "Heteregenous treatment effects - different estimators",
  lwd = 2,
  cex.lab=1.5
)
lines(density(True_RF, na.rm = T), col = "darkgrey", lty = 2, lwd = 2)

lines(density(Dir1, na.rm = T), col = 2, lwd = 2)
lines(density(Dir2, na.rm = T), col = 3, lwd = 2)

lines(density(EBLUP, na.rm = T), col = 4, lwd = 2)
lines(density(MQ, na.rm = T), co = 5, lwd = 2)

lines(density(RF, na.rm = T), col = 6, lwd = 2)
lines(density(RFt, na.rm = T), co = 7, lwd = 2)

legend(
  -1.2,
  max.y,
  legend = c("True_EBLUP", "True_RF", "Dir1", "Dir2", "EBLUP", "MQ", "RF", "RFt"),
  col = c("black", "darkgrey", 2:7),
  lty = c(2, 2, rep(1, 5)),
  lwd = 2,
  bty = "n"
)
dev.off()












