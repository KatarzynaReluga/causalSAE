library("MASS")

library(mlr3)
library(mlr3learners)
library(DoubleML)

setwd("/user/work/ze23696/Comp/causalSAE")
devtools::load_all()

setwd("/user/work/ze23696/Comp/BinaryMQ")
devtools::load_all()

source("/user/work/ze23696/Comp/QRLM.r")

populations <- read.csv('/user/work/ze23696/Comp/population_sae4.csv')
taus <- read.csv('/user/work/ze23696/Comp/tau_sae4.csv')
tau_true <- taus$tau
tau_treat <- taus$tau_treat
tau_untreat <- taus$tau_untreat



m = 41

Nii = 1000
Ni = rep(Nii, m)
N = sum(Ni)
group_names <- as.data.frame(table(populations$group))$Var1
Nc = as.numeric(table(populations$group[populations$A == 0]))
Nt = as.numeric(table(populations$group[populations$A == 1]))

frac_nt = 0.02
nt <- round(frac_nt * Nt)
frac_nc0 <- c(runif(25, 0.01, 0.5), runif(10, 0.51, 1), runif(6, 0.9, 1))
nc <- ceiling(frac_nc0 * nt)
frac_nc <- nc/Nc

SN = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(SN * 2022)
subpopulation <- sample_subpopulations(populations,
                                       frac_nc = frac_nc, frac_nt = frac_nt,
                                       seed = set.seed(SN * 2022))
data_sample <- data.frame(populations[subpopulation, c(1:15)])
data_out_of_sample <- populations[-subpopulation, c(1:15)]

data_sample_ps <-  data.frame(populations[subpopulation, -c(1:15)])
data_out_of_sample_ps <-  data.frame(populations[-subpopulation, -c(1:15)])

ni = as.numeric(table(data_sample$group))
nc = as.numeric(table(data_sample$group[data_sample$A == 0]))
nt = as.numeric(table(data_sample$group[data_sample$A == 1]))
frac_nN <- dim(data_sample)[1]/dim(populations)[1]


########################################################################################
########################################################################################
# Super Learner ------------------------------------------------------------------

#y_train <- data_sample$y
#x_train <- data_sample[, c(1:11, 13)]

names_sample <- names(data_sample)
#names_out_of_sample <- names(data_out_of_sample)
data_full_prelim <- rbind(data_sample[, !names_sample %in% "y"],
                          data_out_of_sample[, !names_sample %in% "y"])
data_out_of_sampley <- data_out_of_sample[, !names_sample %in% "y"]

data_full <- data_full_prelim[, c(1:11, 13)]

# Imputation --------------------------------------------------------------------

################################################################
# Linear regression -------------------------------------------
################################################################

impute_E <- impute_y(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A + (1|group),
                     data_sample,
                     data_out_of_sample,
                     method = "EBLUP",
                     type_model = "gaussian")

y_full_imputeE <- impute_E$y_full_imputed

data_full_i <- data_full_prelim
data_full_i$y <- y_full_imputeE

########################################################################################
########################################################################################

LL <- NULL
LR <- NULL
LX <- NULL

RL <- NULL
RR <- NULL
RX <- NULL

XL <- NULL
XR <- NULL
XX <- NULL

LL_up <- NULL
LR_up <- NULL
LX_up <- NULL

RL_up <- NULL
RR_up <- NULL
RX_up <- NULL

XL_up <- NULL
XR_up <- NULL
XX_up <- NULL


LL_do <- NULL
LR_do <- NULL
LX_do <- NULL

RL_do <- NULL
RR_do <- NULL
RX_do <- NULL

XL_do <- NULL
XR_do <- NULL
XX_do <- NULL

LL_se <- NULL
LR_se <- NULL
LX_se <- NULL

RL_se <- NULL
RR_se <- NULL
RX_se <- NULL

XL_se <- NULL
XR_se <- NULL
XX_se <- NULL




start_time = Sys.time()
for (i in 1:m) {
  print(i)
  data_full_g <- data_full_i[data_full_i$group == group_names[i], ]
  
  obj_dml = DoubleMLData$new(data_full_g,
                             y_col = "y",
                             d_cols = "A",
                             x_cols = c("X1", "X2", "X3", "X4", "X5",  
                                        "X6", "X7", "X8", "X9", "X10"))
  
  set.seed(SN * m)
  
  learner_lL = lrn("regr.lm")
  learner_lR = lrn("regr.ranger", num.trees = 500,
                   min.node.size = 2, max.depth = 5)
  learner_lX = lrn("regr.xgboost", nrounds = 100)
  
  
  learner_mL = lrn("classif.log_reg") # For prop score  
  learner_mR = lrn("regr.ranger", num.trees = 500,
                   min.node.size = 2, max.depth = 5) # For prop score
  learner_mX = lrn("regr.xgboost", nrounds = 100) # For prop score
  
  # L
  # LL
  doubleml_LL = DoubleMLPLR$new(obj_dml,
                                ml_l = learner_lL,
                                ml_m = learner_mL,
                                score = "partialling out",
                                dml_procedure = "dml1",
                                n_rep = 5)
  
  fit_LL <- doubleml_LL$fit()
  
  LL <- c(LL, fit_LL$coef)
  LL_up <- c(LL_up, fit_LL$confint()[1])
  LL_do <- c(LL_do, fit_LL$confint()[2])
  LL_se <- unname(c(LL_se, fit_LL$se))
  
  # LR
  doubleml_LR = DoubleMLPLR$new(obj_dml,
                                ml_l = learner_lL,
                                ml_m = learner_mR,
                                score = "partialling out",
                                dml_procedure = "dml1",
                                n_rep = 5)
  fit_LR <- doubleml_LR$fit()
  
  LR <- c(LR, fit_LR$coef)
  LR_up <- c(LR_up, fit_LR$confint()[1])
  LR_do <- c(LR_do, fit_LR$confint()[2])
  LR_se <-  unname(c(LR_se, fit_LR$se))
  
  # LX
  doubleml_LX = DoubleMLPLR$new(obj_dml,
                                ml_l = learner_lL,
                                ml_m = learner_mX,
                                score = "partialling out",
                                dml_procedure = "dml1",
                                n_rep = 5)
  fit_LX <- doubleml_LX$fit()
  
  LX <- c(LX, fit_LX$coef)
  LX_up <- c(LX_up, fit_LX$confint()[1])
  LX_do <- c(LX_do, fit_LX$confint()[2])
  LX_se <-  unname(c(LX_se, fit_LX$se))
  
  # R
  # RL
  doubleml_RL = DoubleMLPLR$new(obj_dml,
                                ml_l = learner_lR,
                                ml_m = learner_mL,
                                score = "partialling out",
                                dml_procedure = "dml1",
                                n_rep = 5)
  
  fit_RL <- doubleml_RL$fit()
  
  RL <- c(RL, fit_RL$coef)
  RL_up <- c(RL_up, fit_RL$confint()[1])
  RL_do <- c(RL_do, fit_RL$confint()[2])
  RL_se <- c(RL_se, fit_RL$se)
  
  # RR
  doubleml_RR = DoubleMLPLR$new(obj_dml,
                                ml_l = learner_lR,
                                ml_m = learner_mR,
                                score = "partialling out",
                                dml_procedure = "dml1",
                                n_rep = 5)
  fit_RR <- doubleml_RR$fit()
  
  RR <- c(RR, fit_RR$coef)
  RR_up <- c(RR_up, fit_RR$confint()[1])
  RR_do <- c(RR_do, fit_RR$confint()[2])
  RR_se <- c(RR_se, fit_RR$se)
  
  # RX
  doubleml_RX = DoubleMLPLR$new(obj_dml,
                                ml_l = learner_lR,
                                ml_m = learner_mX,
                                score = "partialling out",
                                dml_procedure = "dml1",
                                n_rep = 5)
  fit_RX <- doubleml_RX$fit()
  
  RX <- c(RX, fit_RX$coef)
  RX_up <- c(RX_up, fit_RX$confint()[1])
  RX_do <- c(RX_do, fit_RX$confint()[2])
  RX_se <- c(RX_se, fit_RX$se)
  
  
  # X
  # XL
  doubleml_XL = DoubleMLPLR$new(obj_dml,
                                ml_l = learner_lX,
                                ml_m = learner_mL,
                                score = "partialling out",
                                dml_procedure = "dml1",
                                n_rep = 5)
  
  fit_XL <- doubleml_XL$fit()
  
  XL <- c(XL, fit_XL$coef)
  XL_up <- c(XL_up, fit_XL$confint()[1])
  XL_do <- c(XL_do, fit_XL$confint()[2])
  XL_se <- c(XL_se, fit_XL$se)
  
  # XR
  doubleml_XR = DoubleMLPLR$new(obj_dml,
                                ml_l = learner_lX,
                                ml_m = learner_mR,
                                score = "partialling out",
                                dml_procedure = "dml1",
                                n_rep = 5)
  fit_XR <- doubleml_XR$fit()
  
  XR <- c(XR, fit_XR$coef)
  XR_up <- c(XR_up, fit_XR$confint()[1])
  XR_do <- c(XR_do, fit_XR$confint()[2])
  XR_se <- c(XR_se, fit_XR$se)
  
  # XX
  doubleml_XX = DoubleMLPLR$new(obj_dml,
                                ml_l = learner_lX,
                                ml_m = learner_mX,
                                score = "partialling out",
                                dml_procedure = "dml1",
                                n_rep = 5)
  fit_XX <- doubleml_XX$fit()
  
  XX <- c(XX, fit_XX$coef)
  XX_up <- c(XX_up, fit_XX$confint()[1])
  XX_do <- c(XX_do, fit_XX$confint()[2])
  XX_se <- c(XX_se, fit_XX$se)
  
  
  
}
#Estimation
stop_time = Sys.time() #Time difference of 2.317492 hours

Results = list(ni = ni,
               nc = nc,
               nt = nt,
               n = sum(ni),
               
               Ni = Ni,
               Nc = Nc,
               Nt = Nt,
               N = sum(Ni),
               frac_nN = frac_nN,
               
               
               #Tau
               tau_true = tau_true,
               tau_treat = tau_treat,
               tau_untreat = tau_untreat,
               
               # LL
               LLE_DML = LL,
               LL_up = LL_up,
               LL_do = LL_do,
               LL_se = LL_se,
               
               # LR
               LRE_DML = LR, 
               LR_up = LR_up, 
               LR_do = LR_do, 
               LR_se = LR_se, 
               
               # LX
               LXE_DML  = LX, 
               LX_up = LX_up, 
               LX_do = LX_do,
               LX_se = LX_se,
               
               # R
               # RL
               RLE_DML = RL, 
               RL_up = RL_up, 
               RL_do  = RL_do, 
               RL_se = RL_se, 
               
               # RR
               RRE_DML = RR, 
               RR_up = RR_up, 
               RR_do = RR_do, 
               RR_se = RR_se, 
               
               # RX
               RXE_DML = RX, 
               RX_up = RX_up, 
               RX_do = RX_do, 
               RX_se = RX_se, 
               
               
               # X
               # XL
               XLE_DML = XL, 
               XL_up = XL_up, 
               XL_do = XL_do, 
               XL_se = XL_se, 
               
               # XR
               XRE_DML = XR, 
               XR_up = XR_up, 
               XR_do = XR_do, 
               XR_se = XR_se, 
               
               # XX
               XXE_DML = XX, 
               XX_up = XX_up, 
               XX_do = XX_do, 
               XX_se = XX_se)
stop_time = Sys.time() #Time difference of 8.192187 mins
outputName = paste("sae4_E_DML_", SN, ".RData", sep = "")
outputPath = file.path("/user/work/ze23696/Comp", outputName)
save("Results", file = outputPath)