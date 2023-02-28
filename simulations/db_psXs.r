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
plot(density(data_pop$y))
####################################################################################
# Estimate the propensity score
####################################################################################
# Drop: single, edu0

formula_p_score = A ~ sex + nationality + age + house_size + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + (1|group)


# XGB

obj_p_score_XGB <- list(data_p_score = data_pop)
class(obj_p_score_XGB) <- "XGB"

ps_hat_XGB <-  p_score(obj_p_score = obj_p_score_XGB,
                       model_formula = formula_p_score,
                       xgboost_params = list(CV_XGB = FALSE,
                                             nfolds = 5,
                                             nrounds = 50))

data_pop$p_score <- ps_hat_XGB

tauL <- calculate_tau(list(data_pop), type_tau = "H")

tau_true <- tauL[[1]]$tau
tau_treat <- tauL[[1]]$tau_treat
tau_untreat <- tauL[[1]]$tau_untreat
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
  samp_index <- sample_subpopulations(data_pop, frac_nc = 0.1, frac_nt = 0.1)
  samp <- data_pop[samp_index, ]
  non_samp <- data_pop[ - samp_index, ]


  #######################################################################
  # X -----------------------------------------------------------
  #######################################################################
  # XEE ------------------------------------------------------
  XEE_AIPWf <- hte(type_hte = "AIPW",
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = FALSE,
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
                                    xgboost_params = list(CV_XGB = FALSE,
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
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = FALSE,
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
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = FALSE,
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
                                    xgboost_params = list(CV_XGB = FALSE,
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
                   data_sample = samp,
                   data_out_of_sample = non_samp,
                   params_OR = list(model_formula = formula_y_OR,
                                    method = "XGB",
                                    xgboost_params = list(CV_XGB = FALSE,
                                                          nfolds = 5,
                                                          nrounds = 50)),
                   params_p_score = list(model_formula = formula_p_score,
                                         method =  "MQ"),
                   params_impute_y = list(model_formula = formula_y_impute,
                                          method = "RF",
                                          tune_RF = FALSE))

  XMR_AIPW <- XMR_AIPWf$tau

  #  bb <- Sys.time()
  #}
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

  outputName = paste("db_psX", a, ".RData", sep = "")
  outputPath = file.path("/home/reluga/Comp", outputName)
  #outputPath = file.path("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE",outputName)
  save("Results", file = outputPath)
