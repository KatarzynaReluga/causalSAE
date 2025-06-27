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
B = Sys.time() 
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

#############################################
# RF clustering, tunning --------------------------------------------------------------------------
#######################################
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
# Imputation ------------------------------------------------------------------------------


# SL ---------------------------------------------------------------------
y_train <- data_sample$y
x_train <- data_sample[, c(1:13)]


SL = SuperLearner(Y = y_train,
                  X = x_train,
                  id = data_sample$group,
                  newX = data_out_of_sampley[, c(1:13)],
                  family = gaussian(),
                  SL.library =c("SL.lm",
                                "SL.glmm",
                                "SL.mq",  "SL.mqch", 
                                "SL.grf", "SL.grf_nc", "SL.grf_t", "SL.grf_nct",
                                "SL.xgboost", "SL.xgboost_t"),
                  cvControl = list(V = 5))
y_full_imputeS0 <- SL$SL.predict[, 1]
y_full_imputeS <- c(data_sample$y, y_full_imputeS0)

################################################################
# Linear regression -------------------------------------------
################################################################

impute_L <- lm(y ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A, data = data_sample)

y_full_imputeL0 <- predict(impute_L, newdata = data_out_of_sampley)
y_full_imputeL <- c(data_sample$y, y_full_imputeL0)


# EBLUP -----------------------------------------------------------------

impute_E <- impute_y(model_formula =  y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A + (1|group),
                         data_sample,
                         data_out_of_sample,
                         method = "EBLUP",
                         type_model = "gaussian")

y_full_imputeE <- impute_E$y_full_imputed


# MQ ------------------------------------------------------------------------

impute_M <- impute_y(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A + (1|group),
                     data_sample,
                     data_out_of_sample,
                     method = "MQ",
                     type_model = "continuous")
y_full_imputeM <- impute_M$y_full_imputed

#  Mch --------------------------------------------------------------------------
Q = sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))  

# Control -----------------------------------------
M_y0 <- QRLM(y = data_sample$y,
              x = cbind(1, data_sample[, -c(13:15)]), 
              q = Q,
              k = 1.345, 
              maxit = 100) 

q <- matrix(c(gridfitinter(y  = data_sample$y,
                            expectile = M_y0$fitted.values,
                            Q  = M_y0$q.values)),
             nrow = length(data_sample$y), 
             ncol = 1)

qmat <- matrix(c(q, data_sample$group), 
                nrow = length(data_sample$y),
                ncol = 2)

Qi <- tapply(qmat[,1], qmat[,2], mean)

M_y <- QRLM(y = data_sample$y,
             x = cbind(1, data_sample[, -c(13:15)]),
             q = Qi, 
             k = 1.345, 
             maxit = 100)

y_full_imputeMch0 <- NULL

for (i in 1:m) { 
  y_full_imputeMch0 <- c(y_full_imputeMch0, (cbind(1, as.matrix(data_out_of_sampley[data_out_of_sampley$group == group_names[i], -c(13:14)])) %*%  M_y$coef[, i]) ) 
  
}

y_full_imputeMch <- c(data_sample$y, y_full_imputeMch0)


# RF ----------------------------------------------------------------------------
# No tuning, no cluster
impute_R <- impute_y(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ A + (1|group),
                      data_sample,
                      data_out_of_sample,
                      method = "RF",
                      tune_RF = FALSE, 
                      clust_RF  = FALSE)

y_full_imputeR <- impute_R$y_full_imputed

# No tuning, cluster
impute_Rc <- impute_y(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ A + (1|group),
                     data_sample,
                     data_out_of_sample,
                     method = "RF",
                     tune_RF = FALSE, 
                     clust_RF  = TRUE)

y_full_imputeRc <- impute_Rc$y_full_imputed

# Tuning, no cluster
impute_Rt <- impute_y(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ A + (1|group),
                      data_sample,
                      data_out_of_sample,
                      method = "RF",
                      tune_RF = TRUE, 
                      clust_RF  = FALSE)

y_full_imputeRt <- impute_Rt$y_full_imputed

# Tuning, cluster
impute_Rct <- impute_y(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ A + (1|group),
                      data_sample,
                      data_out_of_sample,
                      method = "RF",
                      tune_RF = TRUE, 
                      clust_RF  = TRUE)

y_full_imputeRct <- impute_Rct$y_full_imputed


# XGB ----------------------------------------------------------------------------

impute_X <- impute_y(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ A + (1|group),
                     data_sample,
                     data_out_of_sample,
                     method = "XGB",
                     xgboost_params = list(CV_XGB = FALSE,
                                           nfolds = 5,
                                           nrounds = 100))

y_full_imputeX <- impute_X$y_full_imputed


# XGB t----------------------------------------------------------------------------

impute_Xt <- impute_y(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ A + (1|group),
                      data_sample,
                      data_out_of_sample,
                      method = "XGB",
                      xgboost_params = list(CV_XGB = TRUE,
                                            nfolds = 5,
                                            nrounds = 100))

y_full_imputeXt <- impute_Xt$y_full_imputed


df_imputation <- data.frame(L = y_full_imputeL, 
                            E = y_full_imputeE, 
                            
                            M = y_full_imputeM, 
                            Mch = y_full_imputeMch,
                            
                            R = y_full_imputeR,
                            Rt = y_full_imputeRt,
                            Rc = y_full_imputeRc,
                            Rct = y_full_imputeRct,
                            
                            X = y_full_imputeX,
                            Xt = y_full_imputeXt, 
                            
                            S = y_full_imputeS)

E = Sys.time() #Time difference of 4.381926 mins

# Compute AIPW -----------------------------------------------
names_method <- c("L", "E", "M", "Mch", "R", "Rt", "Rc", "Rct", "X", "Xt", "S")
  
AIPW_names <- expand.grid(names_method, names_method, names_method)

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

IPW_names <- expand.grid(names_method, names_method)

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

# OR S
df_S <- data.frame(mu1_y = mu1_y_S,
                   mu0_y = mu0_y_S,
                   group = c(data_sample$group, 
                             data_out_of_sample$group))

tau_treatS = aggregate(df_S$mu1_y, list(df_S$group), FUN = mean)$x
tau_untreatS = aggregate(df_S$mu0_y, list(df_S$group), FUN = mean)$x
S = tau_treatS - tau_untreatS


# OR L 
df_L <- data.frame(mu1_y = mu1_y_L,
                   mu0_y = mu0_y_L,
                   group = c(data_sample$group, 
                             data_out_of_sample$group))

tau_treatL = aggregate(df_L$mu1_y, list(df_L$group), FUN = mean)$x
tau_untreatL = aggregate(df_L$mu0_y, list(df_L$group), FUN = mean)$x
L = tau_treatL - tau_untreatL

# OR E 
df_E <- data.frame(mu1_y = mu1_y_E,
                   mu0_y = mu0_y_E,
                   group = c(data_sample$group, 
                             data_out_of_sample$group))

tau_treatE = aggregate(df_E$mu1_y, list(df_E$group), FUN = mean)$x
tau_untreatE = aggregate(df_E$mu0_y, list(df_E$group), FUN = mean)$x
E = tau_treatE - tau_untreatE


# OR M 
df_M <- data.frame(mu1_y = mu1_y_M,
                   mu0_y = mu0_y_M,
                   group = c(data_sample$group, 
                             data_out_of_sample$group))

tau_treatM = aggregate(df_M$mu1_y, list(df_M$group), FUN = mean)$x
tau_untreatM = aggregate(df_M$mu0_y, list(df_M$group), FUN = mean)$x
M = tau_treatM - tau_untreatM

# OR Mch 
df_Mch <- data.frame(mu1_y = mu1_y_Mch,
                     mu0_y = mu0_y_Mch,
                     group = c(data_sample$group, 
                               data_out_of_sample$group))

tau_treatMch = aggregate(df_Mch$mu1_y, list(df_Mch$group), FUN = mean)$x
tau_untreatMch = aggregate(df_Mch$mu0_y, list(df_Mch$group), FUN = mean)$x
Mch = tau_treatMch - tau_untreatMch


# OR R
df_R <- data.frame(mu1_y = mu1_y_R,
                   mu0_y = mu0_y_R,
                   group = c(data_sample$group, 
                             data_out_of_sample$group))

tau_treatR = aggregate(df_R$mu1_y, list(df_R$group), FUN = mean)$x
tau_untreatR = aggregate(df_R$mu0_y, list(df_R$group), FUN = mean)$x
R = tau_treatR - tau_untreatR

# OR Rc
df_Rc <- data.frame(mu1_y = mu1_y_Rc,
                    mu0_y = mu0_y_Rc,
                    group = c(data_sample$group, 
                              data_out_of_sample$group))

tau_treatRc = aggregate(df_Rc$mu1_y, list(df_Rc$group), FUN = mean)$x
tau_untreatRc = aggregate(df_Rc$mu0_y, list(df_Rc$group), FUN = mean)$x
Rc = tau_treatRc - tau_untreatRc

# OR Rt
df_Rt <- data.frame(mu1_y = mu1_y_Rt,
                    mu0_y = mu0_y_Rt,
                    group = c(data_sample$group, 
                              data_out_of_sample$group))

tau_treatRt = aggregate(df_Rt$mu1_y, list(df_Rt$group), FUN = mean)$x
tau_untreatRt = aggregate(df_Rt$mu0_y, list(df_Rt$group), FUN = mean)$x
Rt = tau_treatRt - tau_untreatRt

# OR Rct
df_Rct <- data.frame(mu1_y = mu1_y_Rct,
                     mu0_y = mu0_y_Rct,
                     group = c(data_sample$group, 
                               data_out_of_sample$group))

tau_treatRct = aggregate(df_Rct$mu1_y, list(df_Rct$group), FUN = mean)$x
tau_untreatRct = aggregate(df_Rct$mu0_y, list(df_Rct$group), FUN = mean)$x
Rct = tau_treatRct - tau_untreatRct

# OR X
df_X <- data.frame(mu1_y = mu1_y_X,
                   mu0_y = mu0_y_X,
                   group = c(data_sample$group, 
                             data_out_of_sample$group))

tau_treatX = aggregate(df_X$mu1_y, list(df_X$group), FUN = mean)$x
tau_untreatX = aggregate(df_X$mu0_y, list(df_X$group), FUN = mean)$x
X = tau_treatX - tau_untreatX

# OR Xt
df_Xt <- data.frame(mu1_y = mu1_y_Xt,
                    mu0_y = mu0_y_Xt,
                    group = c(data_sample$group, 
                              data_out_of_sample$group))

tau_treatXt = aggregate(df_Xt$mu1_y, list(df_Xt$group), FUN = mean)$x
tau_untreatXt = aggregate(df_Xt$mu0_y, list(df_Xt$group), FUN = mean)$x
Xt = tau_treatXt - tau_untreatXt

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
               LLL = LLL,
               ELL = ELL,
               MLL = MLL,        
               MchLL = MchLL,
               RLL = RLL,
               RtLL = RtLL,      
               RcLL = RcLL,
               RctLL = RctLL,     
               XLL = XLL,        
               XtLL = XtLL,
               SLL = SLL,
               LEL = LEL,
               EEL = EEL,
               MEL = MEL,
               MchEL = MchEL,    
               REL = REL,
               RtEL = RtEL,     
               RcEL = RcEL,      
               RctEL = RctEL,   
               XEL = XEL,       
               XtEL = XtEL,      
               SEL = SEL,        
               LML = LML,      
               EML = EML,        
               MML = MML,       
               MchML = MchML,  
               RML = RML,        
               RtML = RtML,      
               RcML = RcML,     
               RctML = RctML,    
               XML = XML,       
               XtML = XtML,    
               SML = SML,        
               LMchL = LMchL,    
               EMchL = EMchL,  
               MMchL = MMchL,    
               MchMchL = MchMchL, 
               RMchL = RMchL,   
               RtMchL = RtMchL,  
               RcMchL = RcMchL,  
               RctMchL = RctMchL,
               XMchL = XMchL,    
               XtMchL = XtMchL, 
               SMchL = SMchL, 
               LRL = LRL,        
               ERL = ERL,
               MRL = MRL, 
               MchRL = MchRL,    
               RRL = RRL,       
               RtRL = RtRL, 
               RcRL = RcRL,      
               RctRL = RctRL,
               XRL = XRL, 
               XtRL = XtRL, 
               
               SRL = SRL,      
               LRtL = LRtL,
               ERtL = ERtL,      
               MRtL = MRtL,
               MchRtL = MchRtL,
               RRtL = RRtL,      
               RtRtL = RtRtL,  
               RcRtL = RcRtL,  
               RctRtL = RctRtL,  
               XRtL = XRtL,   
               XtRtL = XtRtL, 
               SRtL = SRtL,      
               LRcL = LRcL,     
               ERcL = ERcL,      
               MRcL = MRcL,      
               MchRcL = MchRcL,
               RRcL = RRcL,      
               RtRcL = RtRcL,    
               RcRcL = RcRcL, 
               RctRcL = RctRcL,  
               XRcL = XRcL,      
               XtRcL = XtRcL,   
               SRcL = SRcL,    
               LRctL = LRctL,    
               ERctL = ERctL, 
               MRctL = MRctL,  
               MchRctL = MchRctL,
               RRctL = RRctL,   
               RtRctL = RtRctL,  
               RcRctL = RcRctL,  
               RctRctL = RctRctL,
               XRctL = XRctL,    
               XtRctL = XtRctL,  
               SRctL = SRctL, 
               LXL = LXL,      
               EXL = EXL,        
               MXL = MXL,     
               MchXL = MchXL,  
               RXL = RXL,        
               RtXL = RtXL,   
               RcXL = RcXL,  
               RctXL = RctXL,    
               XXL = XXL,     
               XtXL = XtXL,    
               SXL = SXL,        
               LXtL = LXtL,
               
               #101 -200
               EXtL = EXtL, 
               MXtL = MXtL, 
               MchXtL = MchXtL,  
               RXtL = RXtL, 
               RtXtL = RtXtL,
               RcXtL = RcXtL,    
               RctXtL = RctXtL,
               XXtL = XXtL, 
               XtXtL = XtXtL,    
               SXtL = SXtL,   
               LSL = LSL, 
               ESL = ESL,        
               MSL = MSL, 
               MchSL = MchSL,
               RSL = RSL,        
               RtSL = RtSL,
               RcSL = RcSL,
               RctSL = RctSL,    
               XSL = XSL,     
               XtSL = XtSL, 
               SSL = SSL,        
               LLE = LLE,  
               ELE = ELE,       
               MLE = MLE,        
               MchLE = MchLE,
               RLE = RLE,   
               RtLE = RtLE,      
               RcLE = RcLE,   
               RctLE = RctLE,   
               XLE = XLE,        
               XtLE = XtLE, 
               SLE = SLE,     
               LEE = LEE,        
               EEE = EEE,   
               MEE = MEE,   
               MchEE = MchEE,    
               REE = REE,   
               RtEE = RtEE,   
               RcEE = RcEE,      
               RctEE = RctEE, 
               XEE = XEE, 
               XtEE = XtEE,      
               SEE = SEE,  
               LME = LME, 
               EME = EME,        
               MME = MME, 
               MchME = MchME,     
               RME = RME,        
               RtME = RtME, 
               RcME = RcME,     
               RctME = RctME,    
               XME = XME,     
               XtME = XtME, 
               SME = SME,        
               LMchE = LMchE, 
               EMchE = EMchE,
               MMchE = MMchE,    
               MchMchE = MchMchE,
               RMchE = RMchE, 
               RtMchE = RtMchE,  
               RcMchE = RcMchE, 
               RctMchE = RctMchE,
               XMchE = XMchE,    
               XtMchE = XtMchE, 
               SMchE = SMchE, 
               LRE = LRE,        
               ERE = ERE, 
               MRE = MRE, 
               MchRE = MchRE,    
               RRE = RRE,    
               RtRE = RtRE,   
               RcRE = RcRE,      
               RctRE = RctRE, 
               XRE = XRE,   
               XtRE = XtRE,      
               SRE = SRE,  
               LRtE = LRtE,   
               ERtE = ERtE,      
               MRtE = MRtE,   
               MchRtE = MchRtE, 
               RRtE = RRtE,      
               RtRtE = RtRtE, 
               RcRtE = RcRtE, 
               RctRtE = RctRtE,  
               XRtE = XRtE,     
               XtRtE = XtRtE,   
               SRtE = SRtE,      
               LRcE = LRcE,       
               ERcE = ERcE,       
               MRcE = MRcE,      
               MchRcE = MchRcE, 
               RRcE = RRcE,     
               RtRcE = RtRcE,    
               RcRcE = RcRcE, 
               RctRcE = RctRcE,
               XRcE = XRcE,      
               XtRcE = XtRcE, 
               SRcE = SRcE,     
               LRctE = LRctE,    
               ERctE = ERctE,  
               
               #201 - 300
               
               MRctE = MRctE,   
               MchRctE = MchRctE,
               RRctE = RRctE,    
               RtRctE = RtRctE,
               RcRctE = RcRctE,
               RctRctE = RctRctE,
               XRctE = XRctE,
               XtRctE = XtRctE,
               SRctE = SRctE,    
               LXE = LXE,     
               EXE = EXE,   
               MXE = MXE,        
               MchXE = MchXE, 
               RXE = RXE,     
               RtXE = RtXE,      
               RcXE = RcXE,
               RctXE = RctXE,   
               XXE = XXE,        
               XtXE = XtXE, 
               SXE = SXE, 
               LXtE = LXtE,      
               EXtE = EXtE,   
               MXtE = MXtE, 
               MchXtE = MchXtE,  
               RXtE = RXtE,   
               RtXtE = RtXtE, 
               RcXtE = RcXtE,    
               RctXtE = RctXtE,
               XXtE = XXtE,     
               XtXtE = XtXtE,    
               SXtE = SXtE,     
               LSE = LSE,     
               ESE = ESE,        
               MSE = MSE,   
               MchSE = MchSE,   
               RSE = RSE,        
               RtSE = RtSE, 
               RcSE = RcSE,   
               RctSE = RctSE,    
               XSE = XSE,   
               XtSE = XtSE,     
               SSE = SSE,        
               LLM = LLM,     
               ELM = ELM,     
               MLM = MLM,        
               MchLM = MchLM,   
               RLM = RLM,     
               RtLM = RtLM,      
               RcLM = RcLM,
               RctLM = RctLM,
               XLM = XLM,        
               XtLM = XtLM, 
               SLM = SLM, 
               LEM = LEM,        
               EEM = EEM,   
               MEM = MEM, 
               MchEM = MchEM,    
               REM = REM,      
               RtEM = RtEM,   
               RcEM = RcEM,      
               RctEM = RctEM, 
               XEM = XEM,       
               XtEM = XtEM,      
               SEM = SEM, 
               LMM = LMM,     
               EMM = EMM,        
               MMM = MMM,   
               MchMM = MchMM,     
               RMM = RMM,        
               RtMM = RtMM,     
               RcMM = RcMM, 
               RctMM = RctMM,    
               XMM = XMM,   
               XtMM = XtMM,     
               SMM = SMM,        
               LMchM = LMchM,     
               EMchM = EMchM, 
               MMchM = MMchM,    
               MchMchM = MchMchM, 
               RMchM = RMchM,
               RtMchM = RtMchM,  
               RcMchM = RcMchM,
               RctMchM = RctMchM,
               XMchM = XMchM,    
               XtMchM = XtMchM,  
               SMchM = SMchM, 
               LRM = LRM,        
               ERM = ERM, 
               MRM = MRM,     
               MchRM = MchRM,    
               RRM = RRM, 
               RtRM = RtRM,    
               RcRM = RcRM,      
               RctRM = RctRM,   
               XRM = XRM,       
               XtRM = XtRM,      
               SRM = SRM,      
               LRtM = LRtM,    
               ERtM = ERtM,      
               MRtM = MRtM,  
               
               #301 - 400
               MchRtM = MchRtM, 
               RRtM = RRtM,          
               RtRtM = RtRtM,    
               RcRtM = RcRtM,        
               RctRtM = RctRtM,  
               XRtM = XRtM,          
               XtRtM = XtRtM,    
               SRtM = SRtM,          
               LRcM = LRcM,     
               ERcM = ERcM,          
               MRcM = MRcM,     
               MchRcM = MchRcM,      
               RRcM = RRcM,       
               RtRcM = RtRcM,        
               RcRcM = RcRcM,     
               RctRcM = RctRcM,      
               XRcM = XRcM,       
               XtRcM = XtRcM,        
               SRcM = SRcM,        
               LRctM = LRctM,        
               ERctM = ERctM,      
               MRctM = MRctM,        
               MchRctM = MchRctM,   
               RRctM = RRctM,        
               RtRctM = RtRctM,    
               RcRctM = RcRctM,      
               RctRctM = RctRctM,  
               XRctM = XRctM,        
               XtRctM = XtRctM,    
               SRctM = SRctM,        
               LXM = LXM,        
               EXM = EXM,            
               MXM = MXM,         
               MchXM = MchXM,        
               RXM = RXM,         
               RtXM = RtXM,          
               RcXM = RcXM,        
               RctXM = RctXM,        
               XXM = XXM,         
               XtXM = XtXM,          
               SXM = SXM,         
               LXtM = LXtM,          
               EXtM = EXtM,        
               MXtM = MXtM,          
               MchXtM = MchXtM,    
               RXtM = RXtM,          
               RtXtM = RtXtM,     
               RcXtM = RcXtM,        
               RctXtM = RctXtM,    
               XXtM = XXtM,          
               XtXtM = XtXtM,      
               SXtM = SXtM,          
               LSM = LSM,         
               ESM = ESM,            
               MSM = MSM,          
               MchSM = MchSM,        
               RSM = RSM,          
               RtSM = RtSM,          
               RcSM = RcSM,        
               RctSM = RctSM,        
               XSM = XSM,          
               XtSM = XtSM,          
               SSM = SSM,           
               LLMch = LLMch,        
               ELMch = ELMch,      
               MLMch = MLMch,        
               MchLMch = MchLMch,    
               RLMch = RLMch,        
               RtLMch = RtLMch,     
               RcLMch = RcLMch,      
               RctLMch = RctLMch,   
               XLMch = XLMch,        
               XtLMch = XtLMch,    
               SLMch = SLMch,        
               LEMch = LEMch,      
               EEMch = EEMch,        
               MEMch = MEMch,      
               MchEMch = MchEMch,    
               REMch = REMch,       
               RtEMch = RtEMch,      
               RcEMch = RcEMch,     
               RctEMch = RctEMch,    
               XEMch = XEMch,       
               XtEMch = XtEMch,      
               SEMch = SEMch,       
               LMMch = LMMch,        
               EMMch = EMMch,       
               MMMch = MMMch,        
               MchMMch = MchMMch,  
               RMMch = RMMch,        
               RtMMch = RtMMch,    
               RcMMch = RcMMch,      
               RctMMch = RctMMch,  
               XMMch = XMMch,        
               XtMMch = XtMMch,   
               SMMch = SMMch,        
               LMchMch = LMchMch, 
               EMchMch = EMchMch,    
               MMchMch = MMchMch,  
               MchMchMch = MchMchMch,
               
               #401 - 500
               RMchMch = RMchMch,    
               RtMchMch = RtMchMch,  
               RcMchMch = RcMchMch,  
               RctMchMch = RctMchMch,
               XMchMch = XMchMch,  
               XtMchMch = XtMchMch,  
               SMchMch = SMchMch,  
               LRMch = LRMch,        
               ERMch = ERMch,      
               MRMch = MRMch,        
               MchRMch = MchRMch,  
               RRMch = RRMch,        
               RtRMch = RtRMch,    
               RcRMch = RcRMch,      
               RctRMch = RctRMch,  
               XRMch = XRMch,        
               XtRMch = XtRMch,   
               SRMch = SRMch,        
               LRtMch = LRtMch,   
               ERtMch = ERtMch,      
               MRtMch = MRtMch,    
               MchRtMch = MchRtMch,  
               RRtMch = RRtMch,     
               RtRtMch = RtRtMch,    
               RcRtMch = RcRtMch,  
               RctRtMch = RctRtMch,  
               XRtMch = XRtMch,     
               XtRtMch = XtRtMch,    
               SRtMch = SRtMch,   
               LRcMch = LRcMch,      
               ERcMch = ERcMch,      
               MRcMch = MRcMch,      
               MchRcMch = MchRcMch, 
               RRcMch = RRcMch,      
               RtRcMch = RtRcMch,   
               RcRcMch = RcRcMch,    
               RctRcMch = RctRcMch, 
               XRcMch = XRcMch,      
               XtRcMch = XtRcMch,   
               SRcMch = SRcMch,      
               LRctMch = LRctMch,  
               ERctMch = ERctMch,    
               MRctMch = MRctMch,   
               MchRctMch = MchRctMch,
               RRctMch = RRctMch,    
               RtRctMch = RtRctMch,  
               RcRctMch = RcRctMch, 
               RctRctMch = RctRctMch,
               XRctMch = XRctMch,   
               XtRctMch = XtRctMch,  
               SRctMch = SRctMch,  
               LXMch = LXMch,        
               EXMch = EXMch,      
               MXMch = MXMch,        
               MchXMch = MchXMch,  
               RXMch = RXMch,        
               RtXMch = RtXMch,    
               RcXMch = RcXMch,      
               RctXMch = RctXMch,  
               XXMch = XXMch,        
               XtXMch = XtXMch,     
               SXMch = SXMch,        
               LXtMch = LXtMch,    
               EXtMch = EXtMch,      
               MXtMch = MXtMch,    
               MchXtMch = MchXtMch,  
               RXtMch = RXtMch,    
               RtXtMch = RtXtMch,    
               RcXtMch = RcXtMch,  
               RctXtMch = RctXtMch,  
               XXtMch = XXtMch,    
               XtXtMch = XtXtMch,    
               SXtMch = SXtMch,    
               LSMch = LSMch,        
               ESMch = ESMch,      
               MSMch = MSMch,        
               MchSMch = MchSMch,  
               RSMch = RSMch,        
               RtSMch = RtSMch,     
               RcSMch = RcSMch,      
               RctSMch = RctSMch,   
               XSMch = XSMch,       
               XtSMch = XtSMch,    
               SSMch = SSMch,        
               LLR = LLR,         
               ELR = ELR,            
               MLR = MLR,          
               MchLR = MchLR,        
               RLR = RLR,          
               RtLR = RtLR,          
               RcLR = RcLR,       
               RctLR = RctLR,        
               XLR = XLR,        
               XtLR = XtLR,          
               SLR = SLR,       
               LER = LER,            
               EER = EER,       
               MER = MER,            
               MchER = MchER,   
               RER = RER,
               
               #501 - 600
               RtER = RtER,  
               RcER = RcER,  
               RctER = RctER,    
               XER = XER,      
               XtER = XtER,   
               SER = SER,        
               LMR = LMR,     
               EMR = EMR,      
               MMR = MMR,        
               MchMR = MchMR,   
               RMR = RMR,      
               RtMR = RtMR,      
               RcMR = RcMR,    
               RctMR = RctMR, 
               XMR = XMR,        
               XtMR = XtMR,    
               SMR = SMR,     
               LMchR = LMchR,    
               EMchR = EMchR,  
               MMchR = MMchR,  
               MchMchR = MchMchR,
               RMchR = RMchR,   
               RtMchR = RtMchR, 
               RcMchR = RcMchR,  
               RctMchR = RctMchR,
               XMchR = XMchR,   
               XtMchR = XtMchR,  
               SMchR = SMchR,   
               LRR = LRR,      
               ERR = ERR,        
               MRR = MRR,       
               MchRR = MchRR,  
               RRR = RRR,        
               RtRR = RtRR,    
               RcRR = RcRR,    
               RctRR = RctRR,    
               XRR = XRR,      
               XtRR = XtRR,    
               SRR = SRR,      
               LRtR = LRtR,     
               ERtR = ERtR,     
               MRtR = MRtR,      
               MchRtR = MchRtR, 
               RRtR = RRtR,    
               RtRtR = RtRtR,    
               RcRtR = RcRtR,   
               RctRtR = RctRtR, 
               XRtR = XRtR,      
               XtRtR = XtRtR,   
               SRtR = SRtR,     
               LRcR = LRcR,      
               ERcR = ERcR,     
               MRcR = MRcR,    
               MchRcR = MchRcR,  
               RRcR = RRcR,    
               RtRcR = RtRcR,   
               RcRcR = RcRcR,    
               RctRcR = RctRcR, 
               XRcR = XRcR,     
               XtRcR = XtRcR,    
               SRcR = SRcR,     
               LRctR = LRctR,   
               ERctR = ERctR,    
               MRctR = MRctR,   
               MchRctR = MchRctR,
               RRctR = RRctR,    
                RtRctR = RtRctR, 
                RcRctR = RcRctR,  
                RctRctR = RctRctR,
                XRctR = XRctR, 
                XtRctR = XtRctR, 
                SRctR = SRctR,    
                LXR = LXR,       
                EXR = EXR,       
                MXR = MXR,        
                MchXR = MchXR,   
                RXR = RXR,       
                RtXR = RtXR,      
                RcXR = RcXR,     
                RctXR = RctXR,    
                XXR = XXR,        
                XtXR = XtXR,     
                SXR = SXR,        
                LXtR = LXtR,      
                EXtR = EXtR,   
                MXtR = MXtR,   
                MchXtR = MchXtR,  
                RXtR = RXtR,    
                RtXtR = RtXtR,   
                RcXtR = RcXtR,    
                RctXtR = RctXtR, 
                XXtR = XXtR,     
                XtXtR = XtXtR,    
                SXtR = SXtR,   
                LSR = LSR,      
                ESR = ESR,        
                MSR = MSR,      
                MchSR = MchSR,  
                RSR = RSR,        
                RtSR = RtSR, 
               
               #601 - 700
               RcSR = RcSR,     
               RctSR = RctSR,   
               XSR = XSR,          
               XtSR = XtSR,     
               SSR = SSR,        
               LLRt = LLRt,        
               ELRt = ELRt,      
               MLRt = MLRt,      
               MchLRt = MchLRt,    
               RLRt = RLRt,      
               RtLRt = RtLRt,    
               RcLRt = RcLRt,      
               RctLRt = RctLRt,   
               XLRt = XLRt,     
               XtLRt = XtLRt,      
               SLRt = SLRt,      
               LERt = LERt,      
               EERt = EERt,        
               MERt = MERt,     
               MchERt = MchERt,  
               RERt = RERt,        
               RtERt = RtERt,     
               RcERt = RcERt,     
               RctERt = RctERt,    
               XERt = XERt,     
                XtERt = XtERt,   
               SERt = SERt,        
                LMRt = LMRt,     
                EMRt = EMRt,      
                MMRt = MMRt,        
               MchMRt = MchMRt,  
               RMRt = RMRt,       
               RtMRt = RtMRt,      
               RcMRt = RcMRt,     
               RctMRt = RctMRt,   
               XMRt = XMRt,        
               XtMRt = XtMRt,     
               SMRt = SMRt,       
               LMchRt = LMchRt,    
               EMchRt = EMchRt,   
               MMchRt = MMchRt,   
               MchMchRt = MchMchRt,
               RMchRt = RMchRt,    
               RtMchRt = RtMchRt,  
               RcMchRt = RcMchRt,  
               RctMchRt = RctMchRt,
               XMchRt = XMchRt,   
               XtMchRt = XtMchRt,  
               SMchRt = SMchRt,    
               LRRt = LRRt,       
               ERRt = ERRt,        
               MRRt = MRRt,      
               MchRRt = MchRRt,   
               RRRt = RRRt,        
               RtRRt = RtRRt,    
               RcRRt = RcRRt,    
               RctRRt = RctRRt,    
               XRRt = XRRt,      
               XtRRt = XtRRt,    
               SRRt = SRRt,        
               LRtRt = LRtRt,     
               ERtRt = ERtRt,      
               MRtRt = MRtRt,      
               MchRtRt = MchRtRt, 
               RRtRt = RRtRt,     
               RtRtRt = RtRtRt,    
               RcRtRt = RcRtRt,   
               RctRtRt = RctRtRt,  
               XRtRt = XRtRt,      
               XtRtRt = XtRtRt,   
               SRtRt = SRtRt,     
               LRcRt = LRcRt,      
               ERcRt = ERcRt,    
               MRcRt = MRcRt,     
               MchRcRt = MchRcRt,  
               RRcRt = RRcRt,     
               RtRcRt = RtRcRt,   
               RcRcRt = RcRcRt,    
               RctRcRt = RctRcRt, 
               XRcRt = XRcRt,     
               XtRcRt = XtRcRt,    
               SRcRt = SRcRt,     
               LRctRt = LRctRt,  
               ERctRt = ERctRt,    
               MRctRt = MRctRt,  
               MchRctRt = MchRctRt,
               RRctRt = RRctRt,    
               RtRctRt = RtRctRt,  
               RcRctRt = RcRctRt,  
               RctRctRt = RctRctRt,
               XRctRt = XRctRt,    
               XtRctRt = XtRctRt, 
               SRctRt = SRctRt,    
               LXRt = LXRt,       
               EXRt = EXRt,       
               MXRt = MXRt,        
               MchXRt = MchXRt,   
               RXRt = RXRt,      
               RtXRt = RtXRt,      
               RcXRt = RcXRt,    
               
               #701-800
               RctXRt = RctXRt,   
               XXRt = XXRt,      
               XtXRt = XtXRt,      
               SXRt = SXRt,       
               LXtRt = LXtRt,    
               EXtRt = EXtRt,      
               MXtRt = MXtRt,    
               MchXtRt = MchXtRt, 
               RXtRt = RXtRt,      
               RtXtRt = RtXtRt,   
               RcXtRt = RcXtRt,   
               RctXtRt = RctXtRt,  
               XXtRt = XXtRt,     
               XtXtRt = XtXtRt,   
               SXtRt = SXtRt,      
               LSRt = LSRt,       
               ESRt = ESRt,       
               MSRt = MSRt,        
               MchSRt = MchSRt,   
               RSRt = RSRt,       
               RtSRt = RtSRt,      
                RcSRt = RcSRt,    
                RctSRt = RctSRt,   
                XSRt = XSRt,        
                XtSRt = XtSRt,    
                SSRt = SSRt,      
                LLRc = LLRc,        
                ELRc = ELRc,      
                MLRc = MLRc,    
                MchLRc = MchLRc,    
                RLRc = RLRc,     
                RtLRc = RtLRc,     
                RcLRc = RcLRc,      
                RctLRc = RctLRc,    
                XLRc = XLRc,       
                XtLRc = XtLRc,      
                SLRc = SLRc,     
                LERc = LERc,      
                EERc = EERc,        
                MERc = MERc,      
                MchERc = MchERc,    
                RERc = RERc,        
                RtERc = RtERc,    
                RcERc = RcERc,     
                RctERc = RctERc,    
                XERc = XERc,       
                XtERc = XtERc,     
                SERc = SERc,        
                LMRc = LMRc,      
                EMRc = EMRc,       
                MMRc = MMRc,        
                MchMRc = MchMRc,   
                RMRc = RMRc,       
                RtMRc = RtMRc,      
                RcMRc = RcMRc,     
                RctMRc = RctMRc,   
                XMRc = XMRc,        
                XtMRc = XtMRc,     
                SMRc = SMRc,       
                LMchRc = LMchRc,    
                EMchRc = EMchRc,   
                MMchRc = MMchRc,   
                MchMchRc = MchMchRc,
                RMchRc = RMchRc,   
                RtMchRc = RtMchRc,  
                RcMchRc = RcMchRc,  
                ctMchRc = RctMchRc,
                XMchRc = XMchRc,    
                XtMchRc = XtMchRc,  
                SMchRc = SMchRc,    
                LRRc = LRRc,      
                ERRc = ERRc,        
                MRRc = MRRc,      
                MchRRc = MchRRc,   
                RRRc = RRRc,        
                RtRRc = RtRRc,     
                RcRRc = RcRRc,     
                RctRRc = RctRRc,    
                XRRc = XRRc,       
                XtRRc = XtRRc,   
                SRRc = SRRc,        
                LRtRc = LRtRc,    
                ERtRc = ERtRc,      
                MRtRc = MRtRc,      
                MchRtRc = MchRtRc, 
                RRtRc = RRtRc,     
                RtRtRc = RtRtRc,    
                RcRtRc = RcRtRc,  
                RctRtRc = RctRtRc, 
                XRtRc = XRtRc,      
                XtRtRc = XtRtRc,    
                SRtRc = SRtRc,    
                LRcRc = LRcRc,      
                ERcRc = ERcRc,  
                MRcRc = MRcRc,    
                MchRcRc = MchRcRc,  
                RRcRc = RRcRc,    
                RtRcRc = RtRcRc,   
                RcRcRc = RcRcRc,    
                RctRcRc = RctRcRc,  
               
               #801-900
               XRcRc = XRcRc,    
               XtRcRc = XtRcRc,      
               SRcRc = SRcRc,     
               LRctRc = LRctRc,      
               ERctRc = ERctRc,     
               MRctRc = MRctRc,      
               MchRctRc = MchRctRc,  
               RRctRc = RRctRc,      
               RtRctRc = RtRctRc,   
               RcRctRc = RcRctRc,    
               RctRctRc = RctRctRc, 
               XRctRc = XRctRc,      
               XtRctRc = XtRctRc,   
               SRctRc = SRctRc,      
               LXRc = LXRc,       
               EXRc = EXRc,          
               MXRc = MXRc,       
               MchXRc = MchXRc,      
               RXRc = RXRc,       
               RtXRc = RtXRc,        
               RcXRc = RcXRc,     
               RctXRc = RctXRc,      
               XXRc = XXRc,       
               XtXRc = XtXRc,        
               SXRc = SXRc,        
               LXtRc = LXtRc,        
               EXtRc = EXtRc,      
               MXtRc = MXtRc,        
               MchXtRc = MchXtRc,   
               RXtRc = RXtRc,        
               RtXtRc = RtXtRc,     
               RcXtRc = RcXtRc,      
               RctXtRc = RctXtRc,  
               XXtRc = XXtRc,        
               XtXtRc = XtXtRc,    
               SXtRc = SXtRc,        
               LSRc = LSRc,        
               ESRc = ESRc,          
               MSRc = MSRc,        
               MchSRc = MchSRc,      
               RSRc = RSRc,       
               RtSRc = RtSRc,        
               RcSRc = RcSRc,      
               RctSRc = RctSRc,      
               XSRc = XSRc,       
               XtSRc = XtSRc,        
               SSRc = SSRc,       
               LLRct = LLRct,        
               ELRct = ELRct,      
               MLRct = MLRct,        
               MchLRct = MchLRct,   
               RLRct = RLRct,        
               RtLRct = RtLRct,    
               RcLRct = RcLRct,      
               RctLRct = RctLRct,  
               XLRct = XLRct,        
               XtLRct = XtLRct,   
               SLRct = SLRct,        
               LERct = LERct,      
               EERct = EERct,        
               MERct = MERct,      
               MchERct = MchERct,    
               RERct = RERct,      
               RtERct = RtERct,      
               RcERct = RcERct,    
               RctERct = RctERct,    
               XERct = XERct,      
               XtERct = XtERct,      
               SERct = SERct,      
               LMRct = LMRct,        
               EMRct = EMRct,      
               MMRct = MMRct,        
               MchMRct = MchMRct,  
               RMRct = RMRct,        
               RtMRct = RtMRct,    
               RcMRct = RcMRct,      
               RctMRct = RctMRct,  
               XMRct = XMRct,        
               XtMRct = XtMRct,   
               SMRct = SMRct,        
               LMchRct = LMchRct,  
               EMchRct = EMchRct,    
               MMchRct = MMchRct,  
               MchMchRct = MchMchRct,
               RMchRct = RMchRct,   
               RtMchRct = RtMchRct,  
               RcMchRct = RcMchRct,
               RctMchRct = RctMchRct,
               XMchRct = XMchRct,   
               XtMchRct = XtMchRct,  
               SMchRct = SMchRct, 
               LRRct = LRRct,        
               ERRct = ERRct,       
               MRRct = MRRct,        
               MchRRct = MchRRct,  
               RRRct = RRRct,        
               RtRRct = RtRRct,    
               RcRRct = RcRRct,      
               RctRRct = RctRRct,   
               XRRct = XRRct, 
               
               #901-1000
               XtRRct = XtRRct,  
               SRRct = SRRct,        
               LRtRct = LRtRct,   
               ERtRct = ERtRct,      
               MRtRct = MRtRct,   
               MchRtRct = MchRtRct,  
               RRtRct = RRtRct,   
               RtRtRct = RtRtRct,    
               RcRtRct = RcRtRct,  
               RctRtRct = RctRtRct,  
               XRtRct = XRtRct,  
                XtRtRct = XtRtRct,    
                SRtRct = SRtRct,    
                LRcRct = LRcRct,      
                ERcRct = ERcRct,     
                MRcRct = MRcRct,      
                MchRcRct = MchRcRct,  
                RRcRct = RRcRct,      
                RtRcRct = RtRcRct,   
                RcRcRct = RcRcRct,    
                RctRcRct = RctRcRct, 
                XRcRct = XRcRct,      
                XtRcRct = XtRcRct,   
                SRcRct = SRcRct,      
                LRctRct = LRctRct,   
                ERctRct = ERctRct,    
                MRctRct = MRctRct,   
                MchRctRct = MchRctRct,
                RRctRct = RRctRct,   
                RtRctRct = RtRctRct,  
                RcRctRct = RcRctRct,  
                RctRctRct = RctRctRct,
                XRctRct = XRctRct,   
                XtRctRct = XtRctRct,  
                SRctRct = SRctRct,   
                LXRct = LXRct,        
                EXRct = EXRct,       
                MXRct = MXRct,        
                MchXRct = MchXRct,   
                RXRct = RXRct,        
                RtXRct = RtXRct,     
                RcXRct = RcXRct,      
                RctXRct = RctXRct,   
                XXRct = XXRct,        
                XtXRct = XtXRct,     
                SXRct = SXRct,        
               LXtRct = LXtRct,     
               EXtRct = EXtRct,      
                MXtRct = MXtRct,   
                MchXtRct = MchXtRct,  
                RXtRct = RXtRct,   
                RtXtRct = RtXtRct,    
                RcXtRct = RcXtRct,  
                RctXtRct = RctXtRct,  
                XXtRct = XXtRct,    
                XtXtRct = XtXtRct,    
               SXtRct = SXtRct,      
               LSRct = LSRct,        
                ESRct = ESRct,      
                MSRct = MSRct,        
                MchSRct = MchSRct, 
                RSRct = RSRct,        
                RtSRct = RtSRct,     
                RcSRct = RcSRct,      
               RctSRct = RctSRct,    
               XSRct = XSRct,        
                XtSRct = XtSRct,    
                SSRct = SSRct,        
                LLX = LLX,          
                ELX = ELX,            
                MLX = MLX,          
                MchLX = MchLX,        
                RLX = RLX,          
                RtLX = RtLX,          
                RcLX = RcLX,       
                RctLX = RctLX,        
                XLX = XLX,         
                XtLX = XtLX,          
                SLX = SLX,         
                LEX = LEX,            
                EEX = EEX,          
                MEX = MEX,            
                MchEX = MchEX,     
                REX = REX,            
                RtEX = RtEX,       
                RcEX = RcEX,          
                RctEX = RctEX,    
                XEX = XEX,            
                XtEX = XtEX,      
                SEX = SEX,            
                LMX = LMX,        
                EMX = EMX,            
                MMX = MMX,         
                MchMX = MchMX,        
                RMX = RMX,         
                RtMX = RtMX,          
                RcMX = RcMX,      
                RctMX = RctMX,        
                XMX = XMX,         
                XtMX = XtMX, 
               
               #1001-1100
                SMX = SMX,   
               LMchX = LMchX,   
               EMchX = EMchX,    
               MMchX = MMchX,  
               MchMchX = MchMchX, 
               RMchX = RMchX,    
                RtMchX = RtMchX, 
                RcMchX = RcMchX, 
                RctMchX = RctMchX,
                XMchX = XMchX,   
                XtMchX = XtMchX,  
                SMchX = SMchX,    
                LRX = LRX,       
                ERX = ERX,     
                MRX = MRX,        
                MchRX = MchRX,   
                RRX = RRX,       
                RtRX = RtRX,      
                RcRX = RcRX,   
                RctRX = RctRX, 
                XRX = XRX,        
                XtRX = XtRX,   
                SRX = SRX,      
                LRtX = LRtX,      
                ERtX = ERtX,    
                MRtX = MRtX,     
                MchRtX = MchRtX,  
                RRtX = RRtX,    
                RtRtX = RtRtX,   
                RcRtX = RcRtX,    
                RctRtX = RctRtX, 
                XRtX = XRtX,     
                XtRtX = XtRtX,    
                SRtX = SRtX,    
                LRcX = LRcX,    
                ERcX = ERcX,      
                MRcX = MRcX,    
                MchRcX = MchRcX,
                RRcX = RRcX,      
                RtRcX = RtRcX,   
                RcRcX = RcRcX,   
                RctRcX = RctRcX,  
                XRcX = XRcX,     
                XtRcX = XtRcX,   
                SRcX = SRcX,      
                LRctX = LRctX,   
                ERctX = ERctX,     
                MRctX = MRctX,    
                MchRctX = MchRctX,
                RRctX = RRctX,    
                RtRctX = RtRctX,  
                RcRctX = RcRctX,  
                RctRctX = RctRctX, 
                XRctX = XRctX,    
                XtRctX = XtRctX,  
                SRctX = SRctX,   
                LXX = LXX,        
                EXX = EXX,        
                MXX = MXX,        
                MchXX = MchXX,    
                RXX = RXX,        
                RtXX = RtXX,      
                RcXX = RcXX,      
                RctXX = RctXX,    
                XXX = XXX,       
                XtXX = XtXX,      
                SXX = SXX,        
                LXtX = LXtX,     
                EXtX = EXtX,      
                MXtX = MXtX,     
                MchXtX = MchXtX,  
                RXtX = RXtX,      
                RtXtX = RtXtX,    
                RcXtX = RcXtX,   
                RctXtX = RctXtX,  
                XXtX = XXtX,      
                XtXtX = XtXtX,   
                SXtX = SXtX,      
                LSX = LSX,       
                ESX = ESX,       
                MSX = MSX,        
                MchSX = MchSX,   
                RSX = RSX,       
                RtSX = RtSX,      
                RcSX = RcSX,     
                RctSX = RctSX,   
                XSX = XSX,        
                XtSX = XtSX,     
                SSX = SSX,       
                LLXt = LLXt,      
                ELXt = ELXt,    
                MLXt = MLXt,     
                MchLXt = MchLXt,  
                RLXt = RLXt,     
                RtLXt = RtLXt,   
                RcLXt = RcLXt,    
                RctLXt = RctLXt,  
                XLXt = XLXt,     
                XtLXt = XtLXt,    
                SLXt = SLXt,
                
                #1101 - 1200
                
                 LEXt = LEXt, 
                 EEXt = EEXt,    
                 MEXt = MEXt,        
   MchEXt = MchEXt,   
   REXt = REXt,       
   RtEXt = RtEXt,      
   RcEXt = RcEXt,     
   RctEXt = RctEXt,  
   XEXt = XEXt,        
  XtEXt = XtEXt,    
  SEXt = SEXt,       
  LMXt = LMXt,        
  EMXt = EMXt,      
  MMXt = MMXt,      
  MchMXt = MchMXt,    
  RMXt = RMXt,      
  RtMXt = RtMXt,    
  RcMXt = RcMXt,      
  RctMXt = RctMXt,   
  XMXt = XMXt,      
  XtMXt = XtMXt,      
  SMXt = SMXt,       
  LMchXt = LMchXt,   
  EMchXt = EMchXt,    
  MMchXt = MMchXt,   
  MchMchXt = MchMchXt,
  RMchXt = RMchXt,    
  RtMchXt = RtMchXt,  
  RcMchXt = RcMchXt,  
  RctMchXt = RctMchXt,
  XMchXt = XMchXt,    
  XtMchXt = XtMchXt, 
  SMchXt = SMchXt,    
  LRXt = LRXt,        
  ERXt = ERXt,       
  MRXt = MRXt,        
  MchRXt = MchRXt,    
  RRXt = RRXt,        
  RtRXt = RtRXt,      
  RcRXt = RcRXt,    
  RctRXt = RctRXt,   
  XRXt = XRXt,        
  XtRXt = XtRXt,   
  SRXt = SRXt,       
  LRtXt = LRtXt,      
  ERtXt = ERtXt,     
  MRtXt = MRtXt,      
  MchRtXt = MchRtXt,  
  RRtXt = RRtXt,     
  RtRtXt = RtRtXt,     
  RcRtXt = RcRtXt,    
  RctRtXt = RctRtXt, 
  XRtXt = XRtXt,      
  XtRtXt = XtRtXt,    
  SRtXt = SRtXt,     
  LRcXt = LRcXt,      
  ERcXt = ERcXt,      
  MRcXt = MRcXt,     
  MchRcXt = MchRcXt,  
  RcXt = RRcXt,      
  RtRcXt = RtRcXt,   
  RcRcXt = RcRcXt, 
  RctRcXt = RctRcXt,  
  XRcXt = XRcXt,      
  XtRcXt = XtRcXt, 
  SRcXt = SRcXt,      
  LRctXt = LRctXt,  
  ERctXt = ERctXt,  
  MRctXt = MRctXt,    
 MchRctXt = MchRctXt,
 RRctXt = RRctXt,    
 RtRctXt = RtRctXt,  
  RcRctXt = RcRctXt, 
  RctRctXt = RctRctXt,
  XRctXt = XRctXt,    
  XtRctXt = XtRctXt,  
  SRctXt = SRctXt,   
  LXXt = LXXt,        
  EXXt = EXXt,       
  MXXt = MXXt,       
  MchXXt = MchXXt,    
  RXXt = RXXt,      
  RtXXt = RtXXt,      
  RcXXt = RcXXt,      
  RctXXt = RctXXt,   
  XXXt = XXXt,      
  XtXXt = XtXXt,      
  SXXt = SXXt,      
  LXtXt = LXtXt,   
  EXtXt = EXtXt,      
  MXtXt = MXtXt,    
  MchXtXt = MchXtXt,
  RXtXt = RXtXt,      
  RtXtXt = RtXtXt,   
  RcXtXt = RcXtXt,   
  RctXtXt = RctXtXt,  
  XXtXt = XXtXt,     
  XtXtXt = XtXtXt,   
  SXtXt = SXtXt,      
 LSXt = LSXt,  
                 
                 #1201-13000
                 
                  ESXt = ESXt,   
                  MSXt = MSXt,     
                  MchSXt = MchSXt,  
   RSXt = RSXt,      
   RtSXt = RtSXt,  
   RcSXt = RcSXt,    
   RctSXt = RctSXt, 
   XSXt = XSXt,     
   XtSXt = XtSXt,    
  SSXt = SSXt,     
  LLS = LLS,    
  ELS = ELS,        
  MLS = MLS,     
  MchLS = MchLS,   
  RLS = RLS,        
  RtLS = RtLS,     
  RcLS = RcLS,     
  RctLS = RctLS,    
  XLS = XLS,      
  XtLS = XtLS,    
  SLS = SLS,        
  LES = LES,      
  EES = EES,      
  MES = MES,        
  MchES = MchES,  
  RES = RES,       
  RtES = RtES,      
  RcES = RcES,     
  RctES = RctES,   
  XES = XES,        
  XtES = XtES,     
  SES = SES,     
  LMS = LMS,        
  EMS = EMS,       
  MMS = MMS,      
  MchMS = MchMS,    
  RMS = RMS,       
  RtMS = RtMS,    
  RcMS = RcMS,      
  RctMS = RctMS,   
  XMS = XMS,       
  XtMS = XtMS,      
  SMS = SMS,      
  LMchS = LMchS,    
  EMchS = EMchS,    
  MMchS = MMchS,   
  MchMchS = MchMchS, 
  RMchS = RMchS,    
  RtMchS = RtMchS,   
  RcMchS = RcMchS,  
  RctMchS = RctMchS,
  XMchS = XMchS,    
  XtMchS = XtMchS,  
  SMchS = SMchS,    
  LRS = LRS,        
  ERS = ERS,        
  MRS = MRS,        
  MchRS = MchRS,    
  RRS = RRS,        
  RtRS = RtRS,      
  RcRS = RcRS,     
  RctRS = RctRS,   
  XRS = XRS,        
  XtRS = XtRS,     
  SRS = SRS,       
  LRtS = LRtS,      
  ERtS = ERtS,     
  MRtS = MRtS,    
  MchRtS = MchRtS,  
  RRtS = RRtS,     
  RtRtS = RtRtS,    
  RcRtS = RcRtS,    
  RctRtS = RctRtS,  
  XRtS = XRtS,     
  XtRtS = XtRtS,    
  SRtS = SRtS,     
  LRcS = LRcS,      
  ERcS = ERcS,      
  MRcS = MRcS,      
  MchRcS = MchRcS,  
  RRcS = RRcS,      
  RtRcS = RtRcS,    
  RcRcS = RcRcS,    
  RctRcS = RctRcS,  
  XRcS = XRcS,      
  XtRcS = XtRcS,   
  SRcS = SRcS,      
  LRctS = LRctS,   
  ERctS = ERctS,    
  RctS = MRctS,    
  MchRctS = MchRctS,
  RRctS = RRctS,    
  RtRctS = RtRctS,  
  RcRctS = RcRctS,  
  RctRctS = RctRctS,
  XRctS = XRctS,    
  XtRctS = XtRctS,  
  SRctS = SRctS,    
  LXS = LXS,        
 EXS = EXS,  
 
 #1300-1331
 
  MXS = MXS,    
  MchXS = MchXS,  
  RXS = RXS,      
 RtXS = RtXS,   
 RcXS = RcXS,   
 RctXS = RctXS,  
  XXS = XXS,   
  XtXS = XtXS,   
  SXS = SXS,      
 LXtS = LXtS,    
 EXtS = EXtS,    
 MXtS = MXtS,    
 MchXtS = MchXtS, 
 RXtS = RXtS,    
 RtXtS = RtXtS,  
 RcXtS = RcXtS,   
 RctXtS = RctXtS,
 XXtS = XXtS,    
 XtXtS = XtXtS,  
 SXtS = SXtS,    
 LSS = LSS,     
 ESS = ESS,     
 MSS = MSS,     
 MchSS = MchSS,  
 RSS = RSS,     
 RtSS = RtSS,   
 RcSS = RcSS,    
 RctSS = RctSS, 
 XSS = XSS,     
 XtSS = XtSS,    
 SSS = SSS, 
 
 # IPW
 # 1 - 100
 LL = LL,
 EL = EL,        
 ML = ML, 
 MchL = MchL,    
 RL = RL,     
 RtL = RtL,      
 RcL = RcL,    
 RctL = RctL,    
 XL = XL,   
 XtL = XtL,      
 SL = SL,    
 LE = LE,        
 EE = EE,    
 ME = ME,        
 MchE = MchE,
 RE = RE,        
 RtE = RtE,      
 RcE = RcE,      
 RctE = RctE,     
 XE = XE,        
 XtE = XtE,    
 SE = SE, 
 LM = LM,    
 EM = EM,        
 MM = MM,    
 MchM = MchM,    
 RM = RM,  
 RtM = RtM,      
 RcM = RcM,  
 RctM = RctM,    
 XM = XM,     
 XtM = XtM,      
 SM = SM,      
 LMch = LMch,    
 EMch = EMch,    
 MMch = MMch,    
 MchMch = MchMch,
 RMch = RMch,    
 RtMch = RtMch,  
 RcMch = RcMch,  
 RctMch = RctMch,
 XMch = XMch,    
 XtMch = XtMch,  
 SMch = SMch,    
 LR = LR,      
 ER = ER,        
 MR = MR,      
 MchR = MchR,    
 RR = RR,    
 RtR = RtR,    
 RcR = RcR,  
 RctR = RctR,    
 XR = XR,    
 XtR = XtR,      
 SR = SR,    
 LRt = LRt,      
 ERt = ERt,    
 MRt = MRt,      
 MchRt = MchRt,
 RRt = RRt,      
 RtRt = RtRt,    
 RcRt = RcRt,    
 RctRt = RctRt,  
 XRt = XRt,      
 XtRt = XtRt,
 SRt = SRt,      
 LRc = LRc,
 ERc = ERc,      
 MRc = MRc,    
 MchRc = MchRc,  
 RRc = RRc,  
 RtRc = RtRc,    
 RcRc = RcRc,  
 RctRc = RctRc,  
 XRc = XRc,  
 XtRc = XtRc,    
 SRc = SRc,  
 LRct = LRct,    
 ERct = ERct,  
 MRct = MRct,    
 MchRct = MchRct,
 RRct = RRct,    
 RtRct = RtRct,  
 RcRct = RcRct,  
 RctRct = RctRct,
 XRct = XRct,    
 XtRct = XtRct,
 SRct = SRct,    
 LX = LX,    
 EX = EX,        
 MX = MX,  
 MchX = MchX,    
 RX = RX,    
 RtX = RtX,      
 RcX = RcX,  
 RctX = RctX,    
 XX = XX,  
 XtX = XtX,      
 SX = SX,      
 LXt = LXt,      
 EXt = EXt,  
 MXt = MXt,      
 MchXt = MchXt, 
 RXt = RXt,      
 RtXt = RtXt,  
 RcXt = RcXt,    
 RctXt = RctXt,  
 XXt = XXt,      
 XtXt = XtXt,  
 SXt = SXt,      
 LS = LS,      
 ES = ES,        
 MS = MS,    
 MchS = MchS,    
 RS = RS,    
 RtS = RtS,      
 RcS = RcS,
 RctS = RctS,    
 XS = XS,  
 tS = XtS,      
 SS = SS, 
 ####################
 # Direct estimator #
 ####################
 Dir_tau = Dir_tau, 
 
 S = S, 
 L = L,
 E = E, 
 M = M,
 Mch = Mch, 
 
 R = R,
 Rc = Rc,
 Rt = Rt,
 Rct = Rct,
 X = X, 
 
 Xt = Xt)
stop_time = Sys.time() 
outputName = paste("sae4_", SN, ".RData", sep = "")
outputPath = file.path("/user/work/ze23696/Comp", outputName)
save("Results", file = outputPath)