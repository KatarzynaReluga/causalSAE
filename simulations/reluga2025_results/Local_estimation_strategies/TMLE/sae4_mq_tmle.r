#remotes::install_github("tlverse/sl3")
#remotes::install_github("tlverse/tmle3")
library(Rsolnp)
library(MASS)
library(data.table)
library(dplyr)
library(tmle3)
library(sl3)


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

data_full_i <- data_full_prelim
data_full_i$y <- y_full_imputeMch
########################################################################################
########################################################################################

TM <- NULL
TM_up <- NULL
TM_do <- NULL
TM_se <- NULL


start_time = Sys.time()
for (i in 1:m) {
  data_full_g <- data_full_i[data_full_i$group == group_names[i], -c(11, 13:14)]
  
  
  node_list <- list(
    W = c("X1", "X2", "X3", "X4", "X5",  
          "X6", "X7", "X8", "X9", "X10"),
    A = "A",
    Y = "y"
  )
  
  ate_spec <- tmle_ATE(
    treatment_level = "1",
    control_level = "0"
  )
  
  # choose base learners
  lrnr_mean <- make_learner(Lrnr_mean)
  lrnr_rf <- make_learner(Lrnr_ranger)
  lrnr_xgb <- make_learner(Lrnr_xgboost)
  
  # define metalearners appropriate to data types
  ls_metalearner <- make_learner(Lrnr_nnls)
  mn_metalearner <- make_learner(
    Lrnr_solnp, metalearner_logistic_binomial,
    loss_loglik_binomial
  )
  sl_Y <- Lrnr_sl$new(
    learners = list(lrnr_mean, lrnr_rf, lrnr_xgb),
    metalearner = ls_metalearner
  )
  sl_A <- Lrnr_sl$new(
    learners = list(lrnr_mean, lrnr_rf, lrnr_xgb),
    metalearner = mn_metalearner
  )
  learner_list <- list(A = sl_A, Y = sl_Y)
  
  tmle_fit <- tmle3(ate_spec, data_full_g, node_list, learner_list)
  TM <- c(TM, tmle_fit$summary$tmle_est)
  TM_up <- c(TM_up, tmle_fit$summary$upper)
  TM_do <- c(TM_do, tmle_fit$summary$lower)
  TM_se <- c(TM_se, tmle_fit$summary$se)
  
  
}
#Estimation
stop_time = Sys.time() #Time difference of 12.91595 mins

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
               
               # TM
               TM_Mq = TM,
               TM_up = TM_up,
               TM_do = TM_do,
               TM_se = TM_se)
stop_time = Sys.time()  
outputName = paste("sae4_Mq_TMLE_", SN, ".RData", sep = "")
outputPath = file.path("/user/work/ze23696/Comp", outputName)
save("Results", file = outputPath)