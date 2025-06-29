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

impute_M <- impute_y(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A + (1|group),
                     data_sample,
                     data_out_of_sample,
                     method = "MQ",
                     type_model = "continuous")
y_full_imputeM <- impute_M$y_full_imputed

data_full_i <- data_full_prelim
data_full_i$y <- y_full_imputeM

########################################################################################
########################################################################################


LLM <- NULL
LRM <- NULL
LXM <- NULL

RLM <- NULL
RRM <- NULL
RXM <- NULL

XLM <- NULL
XRM <- NULL
XXM <- NULL

LLM_up <- NULL
LRM_up <- NULL
LXM_up <- NULL

RLM_up <- NULL
RRM_up <- NULL
RXM_up <- NULL

XLM_up <- NULL
XRM_up <- NULL
XXM_up <- NULL


LLM_do <- NULL
LRM_do <- NULL
LXM_do <- NULL

RLM_do <- NULL
RRM_do <- NULL
RXM_do <- NULL

XLM_do <- NULL
XRM_do <- NULL
XXM_do <- NULL

LLM_se <- NULL
LRM_se <- NULL
LXM_se <- NULL

RLM_se <- NULL
RRM_se <- NULL
RXM_se <- NULL

XLM_se <- NULL
XRM_se <- NULL
XXM_se <- NULL




start_time = Sys.time()
for (i in 1:m) {
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
  doubleml_LLM = DoubleMLPLR$new(obj_dml,
                                 ml_l = learner_lL,
                                 ml_m = learner_mL,
                                 score = "partialling out",
                                 dml_procedure = "dml1",
                                 n_rep = 5)
  
  fit_LLM <- doubleml_LLM$fit()
  
  LLM <- c(LLM, fit_LLM$coef)
  LLM_up <- c(LLM_up, fit_LLM$confint()[1])
  LLM_do <- c(LLM_do, fit_LLM$confint()[2])
  LLM_se <- unname(c(LLM_se, fit_LLM$se))
  
  # LRM
  doubleml_LRM = DoubleMLPLR$new(obj_dml,
                                 ml_l = learner_lL,
                                 ml_m = learner_mR,
                                 score = "partialling out",
                                 dml_procedure = "dml1",
                                 n_rep = 5)
  fit_LRM <- doubleml_LRM$fit()
  
  LRM <- c(LRM, fit_LRM$coef)
  LRM_up <- c(LRM_up, fit_LRM$confint()[1])
  LRM_do <- c(LRM_do, fit_LRM$confint()[2])
  LRM_se <-  unname(c(LRM_se, fit_LRM$se))
  
  # LXM
  doubleml_LXM = DoubleMLPLR$new(obj_dml,
                                 ml_l = learner_lL,
                                 ml_m = learner_mX,
                                 score = "partialling out",
                                 dml_procedure = "dml1",
                                 n_rep = 5)
  fit_LXM <- doubleml_LXM$fit()
  
  LXM <- c(LXM, fit_LXM$coef)
  LXM_up <- c(LXM_up, fit_LXM$confint()[1])
  LXM_do <- c(LXM_do, fit_LXM$confint()[2])
  LXM_se <-  unname(c(LXM_se, fit_LXM$se))
  
  # R
  # RLM
  doubleml_RLM = DoubleMLPLR$new(obj_dml,
                                 ml_l = learner_lR,
                                 ml_m = learner_mL,
                                 score = "partialling out",
                                 dml_procedure = "dml1",
                                 n_rep = 5)
  
  fit_RLM <- doubleml_RLM$fit()
  
  RLM <- c(RLM, fit_RLM$coef)
  RLM_up <- c(RLM_up, fit_RLM$confint()[1])
  RLM_do <- c(RLM_do, fit_RLM$confint()[2])
  RLM_se <- c(RLM_se, fit_RLM$se)
  
  # RRM
  doubleml_RRM = DoubleMLPLR$new(obj_dml,
                                 ml_l = learner_lR,
                                 ml_m = learner_mR,
                                 score = "partialling out",
                                 dml_procedure = "dml1",
                                 n_rep = 5)
  fit_RRM <- doubleml_RRM$fit()
  
  RRM <- c(RRM, fit_RRM$coef)
  RRM_up <- c(RRM_up, fit_RRM$confint()[1])
  RRM_do <- c(RRM_do, fit_RRM$confint()[2])
  RRM_se <- c(RRM_se, fit_RRM$se)
  
  # RXM
  doubleml_RXM = DoubleMLPLR$new(obj_dml,
                                 ml_l = learner_lR,
                                 ml_m = learner_mX,
                                 score = "partialling out",
                                 dml_procedure = "dml1",
                                 n_rep = 5)
  fit_RXM <- doubleml_RXM$fit()
  
  RXM <- c(RXM, fit_RXM$coef)
  RXM_up <- c(RXM_up, fit_RXM$confint()[1])
  RXM_do <- c(RXM_do, fit_RXM$confint()[2])
  RXM_se <- c(RXM_se, fit_RXM$se)
  
  
  # X
  # XLM
  doubleml_XLM = DoubleMLPLR$new(obj_dml,
                                 ml_l = learner_lX,
                                 ml_m = learner_mL,
                                 score = "partialling out",
                                 dml_procedure = "dml1",
                                 n_rep = 5)
  
  fit_XLM <- doubleml_XLM$fit()
  
  XLM <- c(XLM, fit_XLM$coef)
  XLM_up <- c(XLM_up, fit_XLM$confint()[1])
  XLM_do <- c(XLM_do, fit_XLM$confint()[2])
  XLM_se <- c(XLM_se, fit_XLM$se)
  
  # XRM
  doubleml_XRM = DoubleMLPLR$new(obj_dml,
                                 ml_l = learner_lX,
                                 ml_m = learner_mR,
                                 score = "partialling out",
                                 dml_procedure = "dml1",
                                 n_rep = 5)
  fit_XRM <- doubleml_XRM$fit()
  
  XRM <- c(XRM, fit_XRM$coef)
  XRM_up <- c(XRM_up, fit_XRM$confint()[1])
  XRM_do <- c(XRM_do, fit_XRM$confint()[2])
  XRM_se <- c(XRM_se, fit_XRM$se)
  
  # XXM
  doubleml_XXM = DoubleMLPLR$new(obj_dml,
                                 ml_l = learner_lX,
                                 ml_m = learner_mX,
                                 score = "partialling out",
                                 dml_procedure = "dml1",
                                 n_rep = 5)
  fit_XXM <- doubleml_XXM$fit()
  
  XXM <- c(XXM, fit_XXM$coef)
  XXM_up <- c(XXM_up, fit_XXM$confint()[1])
  XXM_do <- c(XXM_do, fit_XXM$confint()[2])
  XXM_se <- c(XXM_se, fit_XXM$se)
  
  
  
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
               
               # LLM
               LLM_DML = LLM,
               LLM_up = LLM_up,
               LLM_do = LLM_do,
               LLM_se = LLM_se,
               
               # LRM
               LRM_DML = LRM, 
               LRM_up = LRM_up, 
               LRM_do = LRM_do, 
               LRM_se = LRM_se, 
               
               # LXM
               LXM_DML  = LXM, 
               LXM_up = LXM_up, 
               LXM_do = LXM_do,
               LXM_se = LXM_se,
               
               # R
               # RLM
               RLM_DML = RLM, 
               RLM_up = RLM_up, 
               RLM_do  = RLM_do, 
               RLM_se = RLM_se, 
               
               # RRM
               RRM_DML = RRM, 
               RRM_up = RRM_up, 
               RRM_do = RRM_do, 
               RRM_se = RRM_se, 
               
               # RXM
               RXM_DML = RXM, 
               RXM_up = RXM_up, 
               RXM_do = RXM_do, 
               RXM_se = RXM_se, 
               
               
               # X
               # XLM
               XLM_DML = XLM, 
               XLM_up = XLM_up, 
               XLM_do = XLM_do, 
               XLM_se = XLM_se, 
               
               # XRM
               XRM_DML = XRM, 
               XRM_up = XRM_up, 
               XRM_do = XRM_do, 
               XRM_se = XRM_se, 
               
               # XXM
               XXM_DML = XXM, 
               XXM_up = XXM_up, 
               XXM_do = XXM_do, 
               XXM_se = XXM_se)
C = Sys.time()
stop_time = Sys.time() #Time difference of 8.192187 mins
outputName = paste("sae4_M_DML_", SN, ".RData", sep = "")
outputPath = file.path("/user/work/ze23696/Comp", outputName)
save("Results", file = outputPath)