library(dplyr)

data_survey <- read.csv2("data_filtered_survey.csv")
setwd("./causalSAE")
devtools::load_all()

# Read and pre_process data

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
#plot(density(data_pop$y))
####################################################################################
# Estimate the propensity score
####################################################################################
# Drop: single, edu0

#formula_p_score = A ~ X1 + (1|group)
formula_p_score = A ~ sex + nationality + age + house_size + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + (1|group)

# EBLUP

obj_p_score_EBLUP <- list(data_p_score = data_pop)
class(obj_p_score_EBLUP) <- "EBLUP"

ps_hat_EBLUP <-  p_score(obj_p_score = obj_p_score_EBLUP,
                         model_formula = formula_p_score)

data_pop$p_score <- ps_hat_EBLUP

tauL <- calculate_tau(list(data_pop), type_tau = "H")

tau_true = tauL[[1]]$tau
tau_treat = tauL[[1]]$tau_treat
tau_untreat = tauL[[1]]$tau_untreat
####################################################################################
## Design based simulations
####################################################################################

Ni = as.numeric(table(data_pop$group))
Nc = as.numeric(table(data_pop$group[data_pop$A == 0]))
Nt = as.numeric(table(data_pop$group[data_pop$A == 1]))

N = sum(Ni)
m = length(unique(data_pop$group))

nic <- round(0.1 * Nc)
nic[which(nic == 0)] <- 1
nit <- round(0.1 * Nt)
nit[which(nit == 0)] <- 1
ni = nic + nit
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

a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# Simple checks of the code ------------------------------------------------------------------
#for (i in 1:NoSim) {
#b  = Sys.time()
#  print(i)
set.seed(a * 2022)

#  aa <- Sys.time()

  samp <- NULL
  samp_index <- sample_subpopulations(data_pop,
                                      frac_nc = 0.1, frac_nt = 0.1,
                                      seed = a * 2022)
  samp <- data_pop[samp_index, ]
  non_samp <- data_pop[ - samp_index, ]


  #######################################################################################################
  # OR #
  #######################################################################################################
  # EBLUP OR ----------------------------------------------------------------------------------------------
  E_ORf <- hte(type_hte = "OR",
               data_sample = samp,
               data_out_of_sample = non_samp,
               params_OR = list(model_formula = formula_y_OR,
                                method = "EBLUP",
                                type_model = "gaussian"))
  E_OR <- E_ORf$tau

  # MQ OR --------------------------------------------------------------------------------------------------------
  M_ORf <- hte(type_hte = "OR",
               data_sample = samp,
               data_out_of_sample = non_samp,
               params_OR = list(model_formula = formula_y_OR,
                                method = "MQ",
                                type_model = "continuous"))
  M_OR <- M_ORf$tau

  # RF OR ------------------------------------------------------------------------------------------------------------

#  R_ORf <- hte(type_hte = "OR",
#               data_sample = samp,
#               data_out_of_sample = non_samp,
#               params_OR = list(model_formula = formula_y_OR,
#                                method = "RF",
#                                tune_RF = TRUE))
#  R_OR <- R_ORf$tau


  test  = "try-error"
  while(test == "try-error") {
    R_ORf <- try(hte(type_hte = "OR",
                 data_sample = samp,
                 data_out_of_sample = non_samp,
                 params_OR = list(model_formula = formula_y_OR,
                                  method = "RF",
                                  tune_RF = TRUE)), silent=TRUE)
    test = class(R_ORf)
  }

  R_OR <- R_ORf$tau
  # EBLUP XGB -----------------------------------------------------------------------------------------------------------

  X_ORf <- hte(type_hte = "OR",
               data_sample = samp,
               data_out_of_sample = non_samp,
               params_OR = list(model_formula = formula_y_OR,
                                method = "XGB",
                                xgboost_params = list(CV_XGB = TRUE,
                                                      nfolds = 5,
                                                      nrounds = 50)))
  X_OR <- X_ORf$tau

  ##############################################################################################################
  # NIPW #
  ##############################################################################################################
  # E
  # EE
  EE_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "EBLUP"),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "EBLUP",
                                         type_model = "gaussian"))

  EE_NIPW <- EE_NIPWf$tau

  #EM ------------------------------------------------------------------------------------------------
  EM_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "EBLUP"),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "MQ",
                                         type_model = "continuous"))

  EM_NIPW <- EM_NIPWf$tau

  # ER ------------------------------------------------------------------------------------------------
  ER_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "EBLUP"),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "RF",
                                         tune_RF = TRUE))
  ER_NIPW <- ER_NIPWf$tau

  #EX --------------------------------------------------------------------------------------------
  EX_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "EBLUP"),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)))

  EX_NIPW <- EX_NIPWf$tau
  ##############################################################################################################
  # M
  # MM -------------------------------------------------------------------------------------------
  MM_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "MQ"),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "MQ",
                                         type_model = "continuous"))
  MM_NIPW <- MM_NIPWf$tau

  # ME --------------------------------------------------------
  ME_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "MQ"),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "EBLUP",
                                         type_model = "gaussian"))
  ME_NIPW <- ME_NIPWf$tau


  # MR --------------------------------------------------------
  MR_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "MQ"),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "RF",
                                         tune_RF = TRUE))

  MR_NIPW <- MR_NIPWf$tau

  # MX --------------------------------------------------------
  MX_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "MQ"),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)))
  MX_NIPW <- MX_NIPWf$tau

  ######################################################################################################
  # R #
  # RR ------------------------------------------------------------------
 # test  = "try-error"
 # while(test == "try-error") {
  RR_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "RF",
                                        tune_RF = TRUE),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "RF",
                                         tune_RF = TRUE))
#  }
  RR_NIPW <- RR_NIPWf$tau
  # RE ------------------------------------------------------------------
  RE_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "RF",
                                        tune_RF = TRUE),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "EBLUP",
                                         type_model = "gaussian"))

  RE_NIPW <- RE_NIPWf$tau

  # RM ------------------------------------------------------------------
  RM_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "RF",
                                        tune_RF = TRUE),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "MQ",
                                         type_model = "continuous"))
  RM_NIPW <- RM_NIPWf$tau
  # RX ------------------------------------------------------------------
#  test  = "try-error"
#  while(test == "try-error") {
  RX_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "RF",
                                        tune_RF = TRUE),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)))
#  }
  RX_NIPW <- RX_NIPWf$tau
  ######################################################################################################
  # X #
  # XX -------------------------------------------------------------------
  XX_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "XGB",
                                        xgboost_params = list(CV_XGB = TRUE,
                                                              nfolds = 5,
                                                              nrounds = 50)),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)))

  XX_NIPW <- XX_NIPWf$tau

  # XE -------------------------------------------------------------------
  XE_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "XGB",
                                        xgboost_params = list(CV_XGB = TRUE,
                                                              nfolds = 5,
                                                              nrounds = 50)),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "EBLUP",
                                         type_model = "gaussian"))

  XE_NIPW <- XE_NIPWf$tau
  # XM -------------------------------------------------------------------
  XM_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "XGB",
                                        xgboost_params = list(CV_XGB = TRUE,
                                                              nfolds = 5,
                                                              nrounds = 50)),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "MQ",
                                         type_model = "continuous"))

  XM_NIPW <- XM_NIPWf$tau

  # XR -------------------------------------------------------------------
#  test  = "try-error"
#  while(test == "try-error") {
  XR_NIPWf <- hte(type_hte = "NIPW",
                  data_sample = samp,
                  data_out_of_sample = non_samp,
                  params_p_score = list(model_formula = formula_p_score,
                                        method =  "XGB",
                                        xgboost_params = list(CV_XGB = TRUE,
                                                              nfolds = 5,
                                                              nrounds = 50)),
                  params_impute_y = list(model_formula = formula_y_impute,
                                         method = "RF",
                                         tune_RF = TRUE))
#  }
  XR_NIPW <- XR_NIPWf$tau


  ##############################################################################################################
  # AIPW #
  ###############################################################################################################
  # AIPW -------------------------------------------------------------------
  # EEM
  EEM_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "EBLUP"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "MQ",
                                          type_model = "continuous"))

  EEM_AIPW <- EEM_AIPWf$tau
  # EER --------------------------------------------------
#  test  = "try-error"
#  while(test == "try-error") {
  EER_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "EBLUP"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))

  EER_AIPW <- EER_AIPWf$tau
#}
  # EEX ----------------------
  EEX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "EBLUP"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))

  EEX_AIPW <- EEX_AIPWf$tau

  # EMM ----------------------------------------------------------------------
  EMM_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
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
 # test  = "try-error"
#  while(test == "try-error") {
  EMR_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "MQ"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
#  }

  EMR_AIPW <- EMR_AIPWf$tau

  # EMX -----------------------------------------
  EMX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "MQ"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))

  EMX_AIPW <- EMX_AIPWf$tau

  # ERM ------------------------------------------------------
#  test  = "try-error"
#  while(test == "try-error") {
  ERM_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "MQ",
                                          type_model = "continuous"))
#}
  ERM_AIPW <- ERM_AIPWf$tau

  # ERR ------------------------------------------------------------------
#  test  = "try-error"
#  while(test == "try-error") {
  ERR_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
#}
  ERR_AIPW <- ERR_AIPWf$tau

  # ERX -------------------------------------------------------------------
#  test  = "try-error"
#  while(test == "try-error") {
  ERX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))
#}
  ERX_AIPW <- ERX_AIPWf$tau

  # EXM --------------------------------------------------
  EXM_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "MQ",
                                          type_model = "continuous"))

  EXM_AIPW <- EXM_AIPWf$tau

  # EXR --------------------------------------------------
#  test  = "try-error"
#  while(test == "try-error") {
  EXR_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
#}
  EXR_AIPW <- EXR_AIPWf$tau

  # EXX ----------------------------------------------------------------
  EXX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "EBLUP",
                                    type_model = "gaussian"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))

  EXX_AIPW <- EXX_AIPWf$tau

  # MM -----------------------------------------------------------
  #######################################################################
  # MEE ------------------------------------------------------
  MEE_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
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
 # test  = "try-error"
#  while(test == "try-error") {
  MER_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "MQ",
                                    type_model = "continuous"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "EBLUP"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
#}
  MER_AIPW <- MER_AIPWf$tau

  # MER ------------------------------------------------------
  MEX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "MQ",
                                    type_model = "continuous"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "EBLUP"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))

  MEX_AIPW <- MEX_AIPWf$tau

  # MME ------------------------------------------------------
  MME_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
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
#  test  = "try-error"
#  while(test == "try-error") {
  MMR_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "MQ",
                                    type_model = "continuous"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "MQ"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
#}
  MMR_AIPW <- MMR_AIPWf$tau

  # MMX ------------------------------------------------------
  MMX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "MQ",
                                    type_model = "continuous"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "MQ"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))

  MMX_AIPW <- MMX_AIPWf$tau

  # MRE ------------------------------------------------------
#  test  = "try-error"
#  while(test == "try-error") {
  MRE_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "MQ",
                                    type_model = "continuous"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "EBLUP",
                                          type_model = "gaussian"))
#}
  MRE_AIPW <- MRE_AIPWf$tau

  # MRR ------------------------------------------------------
#  test  = "try-error"
#  while(test == "try-error") {
  MRR_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "MQ",
                                    type_model = "continuous"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
#}
  MRR_AIPW <- MRR_AIPWf$tau

  # MRX ------------------------------------------------------
#  test  = "try-error"
#  while(test == "try-error") {
  MRX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "MQ",
                                    type_model = "continuous"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))
#}
  MRX_AIPW <- MRX_AIPWf$tau

  # MXE ------------------------------------------------------
  MXE_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "MQ",
                                    type_model = "continuous"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "EBLUP",
                                          type_model = "gaussian"))

  MXE_AIPW <- MXE_AIPWf$tau

  # MRR ------------------------------------------------------
#  test  = "try-error"
##  while(test == "try-error") {
  MXR_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "MQ",
                                    type_model = "continuous"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
#}
  MXR_AIPW <- MXR_AIPWf$tau

  # MRX ------------------------------------------------------
  MXX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR ,
                                    method = "MQ",
                                    type_model = "continuous"),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))

  MXX_AIPW <- MXX_AIPWf$tau

  # RR -----------------------------------------------------------
  #######################################################################
  # REE ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  REE_AIPWf <- try(hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "EBLUP"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "EBLUP",
                                          type_model = "gaussian")), silent = TRUE)
  test = class(REE_AIPWf)
}
  REE_AIPW <- REE_AIPWf$tau

  # REM ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  REM_AIPWf <- try(hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "EBLUP"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "MQ",
                                          type_model = "continuous")), silent = TRUE)
  test = class(REM_AIPWf)
}
  REM_AIPW <- REM_AIPWf$tau
  # REX ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  REX_AIPWf <- try(hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "EBLUP"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50))), silent = TRUE)
  test = class(REX_AIPWf)
}
  REX_AIPW <- REE_AIPWf$tau

  # RME ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  RME_AIPWf <- try(hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "MQ"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "EBLUP",
                                          type_model = "gaussian")), silent = TRUE)
  test = class(RME_AIPWf)
  }
  RME_AIPW <- RME_AIPWf$tau

  # RMM ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  RMM_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "MQ"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "MQ",
                                          type_model = "continuous"))
  }
  RMM_AIPW <- RMM_AIPWf$tau
  # RMX ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  RMX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "MQ"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))
  }
  RMX_AIPW <- RMX_AIPWf$tau

  # RRE ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  RRE_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "EBLUP",
                                          type_model = "gaussian"))
  }
  RRE_AIPW <- RRE_AIPWf$tau

  # RRM ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  RRM_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "MQ",
                                          type_model = "continuous"))
  }
  RRM_AIPW <- RRM_AIPWf$tau
  # RRX ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  RRX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))
  }
  RRX_AIPW <- RRX_AIPWf$tau

  # RXE ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  RXE_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "EBLUP",
                                          type_model = "gaussian"))
  }
  RXE_AIPW <- RXE_AIPWf$tau
  # RXM ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  RXM_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "MQ",
                                          type_model = "continuous"))
  }
  RXM_AIPW <- RXM_AIPWf$tau

  # RXX ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  RXX_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "RF",
                                    tune_RF = TRUE),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "XGB",
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50)))

  }
  RXX_AIPW <- RXX_AIPWf$tau

  # X -----------------------------------------------------------
  #######################################################################
  # XEE ------------------------------------------------------
  XEE_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
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
                   data_sample = samp,
                   data_out_of_sample = non_samp,
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
  test  = "try-error"
  while(test == "try-error") {
  XER_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50)),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "EBLUP"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
  }
  XER_AIPW <- XER_AIPWf$tau

  # XME ------------------------------------------------------
  XME_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
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
                   data_sample = samp,
                   data_out_of_sample = non_samp,
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
  test  = "try-error"
  while(test == "try-error") {
  XMR_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50)),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "MQ"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
  }
  XMR_AIPW <- XMR_AIPWf$tau

  # XRE ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  XRE_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50)),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "EBLUP",
                                          type_model = "gaussian"))
  }
  XRE_AIPW <- XRE_AIPWf$tau

  # XRM ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  XRM_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50)),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "MQ",
                                          type_model = "continuous"))
  }
  XRM_AIPW <- XRM_AIPWf$tau


  # XRR ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  XRR_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50)),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "RF",
                                         tune_RF = TRUE),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
  }
  XRR_AIPW <- XRR_AIPWf$tau

  # XXE ------------------------------------------------------
  XXE_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50)),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "EBLUP",
                                          type_model = "gaussian"))

  XXE_AIPW <- XXE_AIPWf$tau

  # XRM ------------------------------------------------------
  XXM_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50)),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "MQ",
                                          type_model = "continuous"))

  XXM_AIPW <- XXM_AIPWf$tau


  # XER ------------------------------------------------------
  test  = "try-error"
  while(test == "try-error") {
  XXR_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50)),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "XGB",
                                         xgboost_params = list(CV_XGB = TRUE,
                                                               nfolds = 5,
                                                               nrounds = 50)),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = TRUE))
  }
  XXR_AIPW <- XXR_AIPWf$tau
  ####################
  # Direct estimator #
  ####################
  Dir_tau <- (calculate_tau(list(samp), type_tau = "H"))[[1]]$tau

#  bb <- Sys.time()
#}
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

  outputName = paste("db_psEt", a, ".RData", sep = "")
  outputPath = file.path("/home/reluga/Comp", outputName)
  #outputPath = file.path("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE",outputName)
  save("Results", file = outputPath)