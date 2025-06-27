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

error_csda <- numeric(8000)
TM_mat = matrix(0, nrow = 8000, ncol = 41)

for (SN in 4313:4450) {
  print(paste0("SN = ", SN))

#SN = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
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


impute_Ecsda <-
  impute_y(
    model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A + (1+A|group),
    data_sample,
    data_out_of_sample,
    method = "EBLUP",
    type_model = "gaussian"
  )


smcsda = summary(impute_Ecsda$outcome_fit)

y_full_imputecsda <- impute_Ecsda$y_full_imputed
data_full_i <- data_full_prelim
data_full_i$y <- y_full_imputecsda

noerror <- is.null(smcsda$optinfo$conv$lme4$messages)

if (noerror) {
  error_csda[SN] <- 1
  
  TM <- NULL
  start_time = Sys.time()
  for (i in 1:m) {
    print(paste0("area = ", i))
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
#    TM_up <- c(TM_up, tmle_fit$summary$upper)
#    TM_do <- c(TM_do, tmle_fit$summary$lower)
#    TM_se <- c(TM_se, tmle_fit$summary$se)
    
    
  }
  TM_mat[SN, ] <- TM
  stop_time = Sys.time()
} else {
  error_csda[SN] <- 0
  #stop("Convergence issue")
  TM_mat[SN, ] <- rep(0, 41)

}

########################################################################################
########################################################################################

#Estimation
stop_time = Sys.time() #Time difference of 12.91595 mins

}

TM_sel = TM_mat[which(error_csda == 1), ]
mean(colMeans((TM_sel - matrix(rep(tau_true, 254), nrow = 254, ncol = 41, byrow = T))^2))


write.csv(TM_sel,
          file = 'results_Tcsda2.csv', row.names = FALSE)

tau_csda_tmle <- read.csv("results_Tcsda.csv")

mse = mean(colMeans((tau_csda_tmle[c(1:1000), ] - tau_M)^2))
bias = mean(colMeans(tau_csda_tmle[c(1:1000), ] - tau_M))


y_full_imputecsda <- impute_Ecsda$y_full_imputed
data_full_i <- data_full_prelim
data_full_i$y <- y_full_imputecsda

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
               TM_E = TM,
               TM_up = TM_up,
               TM_do = TM_do,
               TM_se = TM_se)
stop_time = Sys.time() 
outputName = paste("sae4_2E_TMLE_", SN, ".RData", sep = "")
outputPath = file.path("/user/work/ze23696/Comp", outputName)
save("Results", file = outputPath)