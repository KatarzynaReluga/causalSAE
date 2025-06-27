library("MASS")

setwd("/user/work/ze23696/Comp/causalSAE")
devtools::load_all()

setwd("/user/work/ze23696/Comp/SuperLearner")
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
start_time = Sys.time()

y_train <- data_sample$y
x_train <- data_sample[, c(1:11, 13)]

names_sample <- names(data_sample)
#names_out_of_sample <- names(data_out_of_sample)
data_full_prelim <- rbind(data_sample[, !names_sample %in% "y"],
                          data_out_of_sample[, !names_sample %in% "y"])
data_out_of_sampley <- data_out_of_sample[, !names_sample %in% "y"]

data_full <- data_full_prelim[, c(1:11, 13)]

# Imputation ------------------------------------------------------------------------------

################################################################
# 2 models
################################################################
# # Check and format data
model_formula_OR <- y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + (1|group)

# Controls -------------------------------------
data_sample0 = data_sample[data_sample$A == 0, ]
obj_fit_y0 <- list(data_sample = data_sample0)

class(obj_fit_y0) <- "EBLUP"
mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                               model_formula = model_formula_OR)

OR0_E <- fit_y(mutated_obj0,
               type_model =  "gaussian")

mu0_y_E <-  unname(unlist(predict(object = OR0_E$outcome_fit,
                                  newdata = data_full_prelim,
                                  allow.new.levels = TRUE)))


# Treated --------------------------------------
data_sample1 = data_sample[data_sample$A == 1, ]
obj_fit_y1 <- list(data_sample = data_sample1)


obj_fit_y1 <- list(data_sample = data_sample1)
class(obj_fit_y1) <- "EBLUP"
mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                               model_formula = model_formula_OR)


OR1_E <- fit_y(mutated_obj1,
               type_model = "gaussian")

mu1_y_E <-  unname(unlist(predict(object = OR1_E$outcome_fit,
                                  newdata = data_full_prelim,
                                  allow.new.levels = TRUE)))

y_full_imputeE2 <- numeric()
A0_full <- which(data_full_prelim$A == 0)
A1_full <- which(data_full_prelim$A == 1)
y_full_imputeE2[A0_full] <- mu0_y_E[A0_full]
y_full_imputeE2[A1_full] <- mu1_y_E[A1_full]

#######################################################################
# csda
####################################################################

impute_Ecsda <-
  impute_y(
    model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A + (1+A|group),
    data_sample,
    data_out_of_sample,
    method = "EBLUP",
    type_model = "gaussian"
  )

smcsda = summary(impute_Ecsda$outcome_fit)

if (is.null(smcsda$optinfo$conv$lme4$messages)) {
  error_csda <- 0
} else {
  error_csda <- 1
  stop("Convergence issue")
}

y_full_imputecsda <- impute_Ecsda$y_full_imputed


df_imputation <- data.frame(E2 = y_full_imputeE2, 
                            csda = y_full_imputecsda)



# Controls -------------------------------------
data_sample0 = data_sample[data_sample$A == 0, ]

y_train0 <- data_sample0$y
x_train0 <- data_sample0[, c(1:11, 13)]

SL0 = SuperLearner(Y = y_train0,
                   X = x_train0,
                   id = data_sample0$group,
                   newX = data_full,
                   family = gaussian(),
                   SL.library = c("SL.lm",
                                  "SL.glmm",
                                  "SL.mq",  "SL.mqch", 
                                  "SL.grf", "SL.grf_nc", "SL.grf_t", "SL.grf_nct",
                                  "SL.xgboost", "SL.xgboost_t"),
                   cvControl = list(V = 5))
mu0_y_S <- SL0$SL.predict[, 1]
#y0_pred_samp <- predict(SL0, newdata = x_train0)$pred

# Treat -----------------------------------------
data_sample1 = data_sample[data_sample$A == 1, ]

y_train1 <- data_sample1$y
x_train1 <- data_sample1[, c(1:11, 13)]
A = Sys.time()
SL1 = SuperLearner(Y = y_train1,
                   X = x_train1,
                   id = data_sample1$group,
                   newX = data_full,
                   family = gaussian(),
                   SL.library =c("SL.lm",
                                 "SL.glmm",
                                 "SL.mq",  "SL.mqch", 
                                 "SL.grf", "SL.grf_nc", "SL.grf_t", "SL.grf_nct",
                                 "SL.xgboost", "SL.xgboost_t"),
                   cvControl = list(V = 5))

mu1_y_S <- SL1$SL.predict[,1]
B = Sys.time() #1 minute
### Other methods
# Check and format data
model_formula_OR <- y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +X11 + (1|group)
# Controls -------------------------------------
data_sample0 = data_sample[data_sample$A == 0, ]
obj_fit_y0 <- list(data_sample = data_sample0)


# Treated --------------------------------------
data_sample1 = data_sample[data_sample$A == 1, ]
obj_fit_y1 <- list(data_sample = data_sample1)

C = Sys.time()
################################################################
# Linear regression -------------------------------------------
################################################################

OR1_L <- lm(y ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11, data = data_sample1)
mu1_y_L <- predict(OR1_L, newdata = data_full_prelim)

OR0_L <- lm(y ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11, data = data_sample0)
mu0_y_L <- predict(OR0_L, newdata = data_full_prelim)

################################################################
# EBLUP ------------------------------------------------------------------------------
###############################################################

class(obj_fit_y0) <- "EBLUP"
mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                               model_formula = model_formula_OR)

OR0_E <- fit_y(mutated_obj0,
               type_model =  "gaussian")

mu0_y_E <-  unname(unlist(predict(object = OR0_E$outcome_fit,
                                  newdata = data_full_prelim,
                                  allow.new.levels = TRUE)))


# Treated --------------------------------------
#data_sample1 = data_sample[data_sample$A == 1, ]

obj_fit_y1 <- list(data_sample = data_sample1)
class(obj_fit_y1) <- "EBLUP"
mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                               model_formula = model_formula_OR)


OR1_E <- fit_y(mutated_obj1,
             type_model = "gaussian")

mu1_y_E <-  unname(unlist(predict(object = OR1_E$outcome_fit,
                                newdata = data_full_prelim,
                                allow.new.levels = TRUE)))
###########################################################
# MQ -----------------------------------------------------------------------------
#############################################################
start_timeM <- Sys.time()
class(obj_fit_y0) <- "MQ"
mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                               model_formula = model_formula_OR)


OR0_M <- fit_y(mutated_obj0,
               type_model =  "continuous")

mu0_y_M <-  unname(unlist(predict(object = OR0_M$outcome_fit,
                                  newdata = data_full_prelim,
                                  allow.new.levels = TRUE)))


# Treated --------------------------------------
class(obj_fit_y1) <- "MQ"
mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                               model_formula = model_formula_OR)


OR1_M <- fit_y(mutated_obj1,
               type_model = "continuous")

mu1_y_M <-  unname(unlist(predict(object = OR1_M$outcome_fit,
                                  newdata = data_full_prelim,
                                  allow.new.levels = TRUE)))
stop_timeM <- Sys.time()
stop_timeM - start_timeM
# Hierarchical M Quantile  --------------------------------------------------------------------------------------------------------

Q = sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))  

# Control -----------------------------------------
M_y00 <- QRLM(y = data_sample0$y,
              x = cbind(1, data_sample0[, -c(12:15)]), 
              q = Q,
              k = 1.345, 
              maxit = 100) 

q0 <- matrix(c(gridfitinter(y  = data_sample0$y,
                            expectile = M_y00$fitted.values,
                            Q  = M_y00$q.values)),
             nrow = length(data_sample0$y), 
             ncol = 1)

qmat0 <- matrix(c(q0, data_sample0$group), 
                nrow = length(data_sample0$y),
                ncol = 2)

Qi0 <- tapply(qmat0[,1], qmat0[,2], mean)

M_y0 <- QRLM(y = data_sample0$y,
             x = cbind(1, data_sample0[, -c(12:15)]),
             q = Qi0, 
             k = 1.345, 
             maxit = 100)

mu0_y_Mch <- NULL

for (i in 1:m) { 
  mu0_y_Mch <- c(mu0_y_Mch, (cbind(1, as.matrix(data_full_prelim[data_full_prelim$group == group_names[i], -c(12:14)])) %*%  M_y0$coef[, i]) ) 
  
}

# Treated --------------------------------------------------
M_y10 <- QRLM(y = data_sample1$y,
              x = cbind(1, data_sample1[,  -c(12:15)]), 
              q = Q,
              k = 1.345, 
              maxit = 100) 

q1 <- matrix(c(gridfitinter(y  = data_sample1$y,
                            expectile = M_y10$fitted.values,
                            Q  = M_y10$q.values)),
             nrow = length(data_sample1$y), 
             ncol = 1)

qmat1 <- matrix(c(q1, data_sample1$group), nrow = length(data_sample1$y),
                ncol = 2)

Qi1 <- tapply(qmat1[,1], qmat1[,2], mean)

M_y1 <- QRLM(y = data_sample1$y,
             x = cbind(1, data_sample1[, -c(12:15)]),
             q = Qi1, 
             k = 1.345, 
             maxit = 100)

mu1_y_Mch <- NULL

for (i in 1:m) { 
  mu1_y_Mch <- c(mu1_y_Mch, (cbind(1, as.matrix(data_full_prelim[data_full_prelim$group == group_names[i], -c(12:15)])) %*%  M_y1$coef[, i]) ) 
  
}


# RF -----------------------------------------------------------------------------
formatted_data <- format_data(model_formula = model_formula_OR,
                              data_sample = data_sample,
                              data_out_of_sample  = data_out_of_sample)
data_full_RX <- rbind(formatted_data$X, formatted_data$X_newdata)
###################################
# No clustering, no tuning
#############################################
startR <- Sys.time()
class(obj_fit_y0) <- "RF"
mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                               model_formula = model_formula_OR)


OR0_R <- fit_y(mutated_obj0,
               type_model = "continuous", 
               tune_RF = FALSE, 
               clust_RF = FALSE)

mu0_y_R <-  unname(unlist(predict(object = OR0_R$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))


# Treated --------------------------------------
class(obj_fit_y1) <- "RF"
mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                               model_formula = model_formula_OR)


OR1_R <- fit_y(mutated_obj1,
               type_model = "continuous", 
               tune_RF = FALSE,
               clust_RF = FALSE)

mu1_y_R <-  unname(unlist(predict(object = OR1_R$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))
stoopR <- Sys.time()
#############################################
# RF no clustering, tunning --------------------------------------------------------------------------
#######################################
OR0_Rt <- fit_y(mutated_obj0,
               type_model = "continuous", 
               tune_RF = TRUE, 
               clust_RF = FALSE)



mu0_y_Rt <-  unname(unlist(predict(object = OR0_Rt$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))


# Treated --------------------------------------
OR1_Rt <- fit_y(mutated_obj1,
               type_model = "continuous", 
               tune_RF = TRUE, 
               clust_RF = FALSE)

mu1_y_Rt <-  unname(unlist(predict(object = OR1_Rt$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))


###################################
# Clustering, no tuning
#############################################
startRc <- Sys.time()
OR0_Rc <- fit_y(mutated_obj0,
               type_model = "continuous", 
               tune_RF = FALSE, 
               clust_RF = TRUE)



mu0_y_Rc <-  unname(unlist(predict(object = OR0_Rc$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))


# Treated --------------------------------------
OR1_Rc <- fit_y(mutated_obj1,
               type_model = "continuous", 
               tune_RF = FALSE,
               clust_RF = TRUE)

mu1_y_Rc <-  unname(unlist(predict(object = OR1_Rc$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))
stopRc <- Sys.time()
stopRc - startRc
#############################################
# RF clustering, tunning --------------------------------------------------------------------------
#######################################
startRct <- Sys.time()
OR0_Rct <- fit_y(mutated_obj0,
                type_model = "continuous", 
                tune_RF = TRUE, 
                clust_RF = TRUE)



mu0_y_Rct <-  unname(unlist(predict(object = OR0_Rct$outcome_fit,
                                   newdata = data_full_RX,
                                   allow.new.levels = TRUE)))


# Treated --------------------------------------
OR1_Rct <- fit_y(mutated_obj1,
                type_model = "continuous", 
                tune_RF = TRUE, 
                clust_RF = TRUE)

mu1_y_Rct <-  unname(unlist(predict(object = OR1_Rct$outcome_fit,
                                   newdata = data_full_RX,
                                   allow.new.levels = TRUE)))

stopRct <- Sys.time()
######################################################
# XGboost
#####################################################
# X --------------------------------------------------------------------------------
class(obj_fit_y0) <- "XGB"
mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                               model_formula = model_formula_OR)


OR0_X <- fit_y(mutated_obj0,
               type_model = "continuous", 
               xgboost_params = list(CV_XGB = FALSE,
                                    nfolds = 5,
                                    nrounds = 100))



mu0_y_X <-  unname(unlist(predict(object = OR0_X$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))


# Treated --------------------------------------
class(obj_fit_y1) <- "XGB"
mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                               model_formula = model_formula_OR)


OR1_X <- fit_y(mutated_obj1,
               type_model = "continuous", 
               xgboost_params = list(CV_XGB = FALSE,
                                     nfolds = 5,
                                     nrounds = 100))

mu1_y_X <-  unname(unlist(predict(object = OR1_X$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))

# X t --------------------------------------------------------------------------------
OR0_Xt <- fit_y(mutated_obj0,
               type_model = "continuous", 
               xgboost_params = list(CV_XGB = TRUE,
                                     nfolds = 5,
                                     nrounds = 100))



mu0_y_Xt <-  unname(unlist(predict(object = OR0_Xt$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))


# Treated --------------------------------------
OR1_Xt <- fit_y(mutated_obj1,
               type_model = "continuous", 
               xgboost_params = list(CV_XGB = TRUE,
                                     nfolds = 5,
                                     nrounds = 100))

mu1_y_Xt <-  unname(unlist(predict(object = OR1_Xt$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))

df_mu1 <- data.frame(L = mu1_y_L, 
                     E = mu1_y_E, 
                    
                     M = mu1_y_M, 
                     Mch = mu1_y_Mch,
                    
                     R = mu1_y_R,
                     Rt = mu1_y_Rt,
                     Rc = mu1_y_Rc,
                     Rct = mu1_y_Rct,
                    
                     X = mu1_y_X,
                     Xt = mu1_y_Xt, 
                     
                     S = mu1_y_S)

df_mu0 <- data.frame(L = mu0_y_L, 
                     E = mu0_y_E, 
                     
                     M = mu0_y_M, 
                     Mch = mu0_y_Mch,
                     
                     R = mu0_y_R,
                     Rt = mu0_y_Rt,
                     Rc = mu0_y_Rc,
                     Rct = mu0_y_Rct,
                     
                     X = mu0_y_X,
                     Xt = mu0_y_Xt, 
                     
                     S = mu0_y_S)

df_ps <- data.frame(L = c(data_sample_ps$ps_hat_L, data_out_of_sample_ps$ps_hat_L), 
                    E = c(data_sample_ps$ps_hat_E, data_out_of_sample_ps$ps_hat_E), 
                            
                    M = c(data_sample_ps$ps_hat_M, data_out_of_sample_ps$ps_hat_M), 
                    Mch = c(data_sample_ps$ps_hat_Mch, data_out_of_sample_ps$ps_hat_Mch),
                          
                    R = c(data_sample_ps$ps_hat_R, data_out_of_sample_ps$ps_hat_R),
                    Rt = c(data_sample_ps$ps_hat_Rt, data_out_of_sample_ps$ps_hat_Rt),
                    Rc = c(data_sample_ps$ps_hat_Rc, data_out_of_sample_ps$ps_hat_Rc),
                    Rct = c(data_sample_ps$ps_hat_Rct, data_out_of_sample_ps$ps_hat_Rct),
                            
                    X = c(data_sample_ps$ps_hat_X, data_out_of_sample_ps$ps_hat_X),
                    Xt = c(data_sample_ps$ps_hat_Xt, data_out_of_sample_ps$ps_hat_Xt),
                    
                    S = c(data_sample_ps$ps_hat_S, data_out_of_sample_ps$ps_hat_S))

E = Sys.time() #Time difference of 4.381926 mins

# Compute AIPW -----------------------------------------------
names_method <- c("L", "E", "M", "Mch", "R", "Rt", "Rc", "Rct", "X", "Xt", "S")
names_method_imp <- c("E2", "csda") 

AIPW_names <- expand.grid(names_method, names_method, names_method_imp)

#names_methods <- numeric(dim(AIPW_names)[1])


C = Sys.time()
for (k in 1:dim(AIPW_names)[1]) {
#for (k in 1:2) {
  
#  group <- X[ ,names(X) %in% c("group")]
  
  mu1_hat <- df_mu1[, names(df_mu1) %in% unlist(unname(AIPW_names[k, 1]))]
  mu0_hat <- df_mu0[, names(df_mu0) %in% unlist(unname(AIPW_names[k, 1]))]
  
  ps_hat <- df_ps[, names(df_ps) %in% unlist(unname(AIPW_names[k, 2]))]
  
  y_hat <- df_imputation[, names(df_imputation) %in% unlist(unname(AIPW_names[k, 3]))]

  
  AIPW_data <- data.frame(
    y = y_hat,
    A = c(data_sample$A, data_out_of_sample$A),
    group = c(data_sample$group, data_out_of_sample$group),
    p_score  = ps_hat, 
    mu0_y = mu0_hat, 
    mu1_y = mu1_hat
  )
  
  AIPWf <-
    as.data.frame(calculate_tau(AIPW_data, type_tau = "AIPW"))

  assign(paste0(unlist(unname(AIPW_names[k, ])), collapse = ""), AIPWf$tau)
#  names_methods[k] <- paste0(paste0(unlist(unname(AIPW_names[k, ])), collapse = ""), " = ", 
#                             paste0(unlist(unname(AIPW_names[k, ])), collapse = ""), ",")
}

D = Sys.time() #Time difference of 16.61508 mins

# Compute IPW -----------------------------------------------
#names_method <- c("L", "E", "M", "Mch", "R", "Rt", "Rc", "Rct", "X", "Xt", "S")

IPW_names <- expand.grid(names_method, names_method_imp)

#names_methods <- numeric(dim(IPW_names)[1])


C = Sys.time()
for (k in 1:dim(IPW_names)[1]) {
  
  #  group <- X[ ,names(X) %in% c("group")]
  
  ps_hat <- df_ps[, names(df_ps) %in% unlist(unname(IPW_names[k, 1]))]
  
  y_hat <- df_imputation[, names(df_imputation) %in% unlist(unname(IPW_names[k, 2]))]
  
  
  IPW_data <- data.frame(
    y = y_hat,
    A = c(data_sample$A, data_out_of_sample$A),
    group = c(data_sample$group, data_out_of_sample$group),
    p_score  = ps_hat
  )
  
  IPWf <-
    as.data.frame(calculate_tau(IPW_data, type_tau = "H"))
  
  assign(paste0(unlist(unname(IPW_names[k, ])), collapse = ""), IPWf$tau)
  #  names_methods[k] <- paste0(paste0(unlist(unname(IPW_names[k, ])), collapse = ""), " = ", 
  #                             paste0(unlist(unname(IPW_names[k, ])), collapse = ""), ",")
  
  
}

###########################################################################
# Store results in the list - standard for baobab.                         #
############################################################################
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
               
               # 1 - 100
               LLE2 = LLE2,
               ELE2 = ELE2,
               MLE2 = MLE2,        
               MchLE2 = MchLE2,
               RLE2 = RLE2,
               RtLE2 = RtLE2,      
               RcLE2 = RcLE2,
               RctLE2 = RctLE2,     
               XLE2 = XLE2,        
               XtLE2 = XtLE2,
               SLE2 = SLE2,
               LEE2 = LEE2,
               EEE2 = EEE2,
               MEE2 = MEE2,
               MchEE2 = MchEE2,    
               REE2 = REE2,
               RtEE2 = RtEE2,     
               RcEE2 = RcEE2,      
               RctEE2 = RctEE2,   
               XEE2 = XEE2,       
               XtEE2 = XtEE2,      
               SEE2 = SEE2,        
               LME2 = LME2,      
               EME2 = EME2,        
               MME2 = MME2,       
               MchME2 = MchME2,  
               RME2 = RME2,        
               RtME2 = RtME2,      
               RcME2 = RcME2,     
               RctME2 = RctME2,    
               XME2 = XME2,       
               XtME2 = XtME2,    
               SME2 = SME2,        
               LMchE2 = LMchE2,    
               EMchE2 = EMchE2,  
               MMchE2 = MMchE2,    
               MchMchE2 = MchMchE2, 
               RMchE2 = RMchE2,   
               RtMchE2 = RtMchE2,  
               RcMchE2 = RcMchE2,  
               RctMchE2 = RctMchE2,
               XMchE2 = XMchE2,    
               XtMchE2 = XtMchE2, 
               SMchE2 = SMchE2, 
               LRE2 = LRE2,        
               ERE2 = ERE2,
               MRE2 = MRE2, 
               MchRE2 = MchRE2,    
               RRE2 = RRE2,       
               RtRE2 = RtRE2, 
               RcRE2 = RcRE2,      
               RctRE2 = RctRE2,
               XRE2 = XRE2, 
               XtRE2 = XtRE2, 
               
               SRE2 = SRE2,      
               LRtE2 = LRtE2,
               ERtE2 = ERtE2,      
               MRtE2 = MRtE2,
               MchRtE2 = MchRtE2,
               RRtE2 = RRtE2,      
               RtRtE2 = RtRtE2,  
               RcRtE2 = RcRtE2,  
               RctRtE2 = RctRtE2,  
               XRtE2 = XRtE2,   
               XtRtE2 = XtRtE2, 
               SRtE2 = SRtE2,      
               LRcE2 = LRcE2,     
               ERcE2 = ERcE2,      
               MRcE2 = MRcE2,      
               MchRcE2 = MchRcE2,
               RRcE2 = RRcE2,      
               RtRcE2 = RtRcE2,    
               RcRcE2 = RcRcE2, 
               RctRcE2 = RctRcE2,  
               XRcE2 = XRcE2,      
               XtRcE2 = XtRcE2,   
               SRcE2 = SRcE2,    
               LRctE2 = LRctE2,    
               ERctE2 = ERctE2, 
               MRctE2 = MRctE2,  
               MchRctE2 = MchRctE2,
               RRctE2 = RRctE2,   
               RtRctE2 = RtRctE2,  
               RcRctE2 = RcRctE2,  
               RctRctE2 = RctRctE2,
               XRctE2 = XRctE2,    
               XtRctE2 = XtRctE2,  
               SRctE2 = SRctE2, 
               LXE2 = LXE2,      
               EXE2 = EXE2,        
               MXE2 = MXE2,     
               MchXE2 = MchXE2,  
               RXE2 = RXE2,        
               RtXE2 = RtXE2,   
               RcXE2 = RcXE2,  
               RctXE2 = RctXE2,    
               XXE2 = XXE2,     
               XtXE2 = XtXE2,    
               SXE2 = SXE2,        
               LXtE2 = LXtE2,
               
               #101 -200
               EXtE2 = EXtE2, 
               MXtE2 = MXtE2, 
               MchXtE2 = MchXtE2,  
               RXtE2 = RXtE2, 
               RtXtE2 = RtXtE2,
               RcXtE2 = RcXtE2,    
               RctXtE2 = RctXtE2,
               XXtE2 = XXtE2, 
               XtXtE2 = XtXtE2,    
               SXtE2 = SXtE2,   
               LSE2 = LSE2, 
               ESE2 = ESE2,        
               MSE2 = MSE2, 
               MchSE2 = MchSE2,
               RSE2 = RSE2,        
               RtSE2 = RtSE2,
               RcSE2 = RcSE2,
               RctSE2 = RctSE2,    
               XSE2 = XSE2,     
               XtSE2 = XtSE2, 
               SSE2 = SSE2,   
               
               # csda
               LLcsda = LLcsda, 
               ELcsda = ELcsda,      
               MLcsda = MLcsda,       
               MchLcsda = MchLcsda,
               RLcsda = RLcsda,  
               RtLcsda = RtLcsda,     
               RcLcsda = RcLcsda,  
               RctLcsda = RctLcsda,  
               XLcsda = XLcsda,       
               XtLcsda = XtLcsda,
               SLcsda = SLcsda,    
               LEcsda = LEcsda,       
               EEcsda = EEcsda,  
               MEcsda = MEcsda,  
               MchEcsda = MchEcsda,   
               REcsda = REcsda,  
               RtEcsda = RtEcsda,  
               RcEcsda = RcEcsda,     
               RctEcsda = RctEcsda,
               XEcsda = XEcsda,
               XtEcsda = XtEcsda,     
               SEcsda = SEcsda, 
               LMcsda = LMcsda,
               EMcsda = EMcsda,       
               MMcsda = MMcsda,
               MchMcsda = MchMcsda,    
               RMcsda = RMcsda,       
               RtMcsda = RtMcsda,
               RcMcsda = RcMcsda,    
               RctMcsda = RctMcsda,   
               XMcsda = XMcsda,    
               XtMcsda = XtMcsda,
               SMcsda = SMcsda,       
               LMchcsda = LMchcsda,
               EMchcsda = EMchcsda,
               MMchcsda = MMchcsda,   
               MchMchcsda = MchMchcsda,
               RMchcsda = RMchcsda,
               RtMchcsda = RtMchcsda, 
               RcMchcsda = RcMchcsda,
               RctMchcsda = RctMchcsda,
               XMchcsda = XMchcsda,   
               XtMchcsda = XtMchcsda,
               SMchcsda = SMchcsda,
               LRcsda = LRcsda,       
               ERcsda = ERcsda,
               MRcsda = MRcsda,
               MchRcsda = MchRcsda,   
               RRcsda = RRcsda,   
               RtRcsda = RtRcsda,  
               RcRcsda = RcRcsda,     
               RctRcsda = RctRcsda,
               XRcsda = XRcsda,  
               XtRcsda = XtRcsda,     
               SRcsda = SRcsda, 
               LRtcsda = LRtcsda,  
               ERtcsda = ERtcsda,     
               MRtcsda = MRtcsda,  
               MchRtcsda = MchRtcsda,
               RRtcsda = RRtcsda,     
               RtRtcsda = RtRtcsda,
               RcRtcsda = RcRtcsda,
               RctRtcsda = RctRtcsda, 
               XRtcsda = XRtcsda,    
               XtRtcsda = XtRtcsda,  
               SRtcsda = SRtcsda,     
               LRccsda = LRccsda,      
               ERccsda = ERccsda,      
               MRccsda = MRccsda,     
               MchRccsda = MchRccsda,
               RRccsda = RRccsda,    
               RtRccsda = RtRccsda,   
               RcRccsda = RcRccsda,
               RctRccsda = RctRccsda,
               XRccsda = XRccsda,     
               XtRccsda = XtRccsda,
               SRccsda = SRccsda,    
               LRctcsda = LRctcsda,   
               ERctcsda = ERctcsda, 
               
               #201 - 300
               
               MRctcsda = MRctcsda,  
               MchRctcsda = MchRctcsda,
               RRctcsda = RRctcsda,   
               RtRctcsda = RtRctcsda,
               RcRctcsda = RcRctcsda,
               RctRctcsda = RctRctcsda,
               XRctcsda = XRctcsda,
               XtRctcsda = XtRctcsda,
               SRctcsda = SRctcsda,   
               LXcsda = LXcsda,    
               EXcsda = EXcsda,  
               MXcsda = MXcsda,       
               MchXcsda = MchXcsda,
               RXcsda = RXcsda,    
               RtXcsda = RtXcsda,     
               RcXcsda = RcXcsda,
               RctXcsda = RctXcsda,  
               XXcsda = XXcsda,       
               XtXcsda = XtXcsda,
               SXcsda = SXcsda,
               LXtcsda = LXtcsda,     
               EXtcsda = EXtcsda,  
               MXtcsda = MXtcsda,
               MchXtcsda = MchXtcsda, 
               RXtcsda = RXtcsda,  
               RtXtcsda = RtXtcsda,
               RcXtcsda = RcXtcsda,   
               RctXtcsda = RctXtcsda,
               XXtcsda = XXtcsda,    
               XtXtcsda = XtXtcsda,   
               SXtcsda = SXtcsda,    
               LScsda = LScsda,    
               EScsda = EScsda,       
               MScsda = MScsda,  
               MchScsda = MchScsda,  
               RScsda = RScsda,       
               RtScsda = RtScsda,
               RcScsda = RcScsda,  
               RctScsda = RctScsda,   
               XScsda = XScsda,  
               XtScsda = XtScsda,    
               SScsda = SScsda, 
               
               
 
 # IPW
 # 1 - 100
 LE2 = LE2,
 EE2 = EE2,       
 ME2 = ME2,
 MchE2 = MchE2,   
 RE2 = RE2,    
 RtE2 = RtE2,     
 RcE2 = RcE2,   
 RctE2 = RctE2,   
 XE2 = XE2,  
 XtE2 = XtE2,     
 SE2 = SE2, 
 
 
 Lcsda = Lcsda,       
 Ecsda = Ecsda,   
 Mcsda = Mcsda,       
 Mchcsda = Mchcsda,
 Rcsda = Rcsda,       
 Rtcsda = Rtcsda,     
 Rccsda = Rccsda,     
 Rctcsda = Rctcsda,    
 Xcsda = Xcsda,       
 Xtcsda = Xtcsda,   
 Scsda = Scsda)

stop_time = Sys.time() #Time difference of 8.192187 mins
outputName = paste("sae4_csda_e2_", SN, ".RData", sep = "")
outputPath = file.path("/user/work/ze23696/Comp", outputName)
save("Results", file = outputPath)