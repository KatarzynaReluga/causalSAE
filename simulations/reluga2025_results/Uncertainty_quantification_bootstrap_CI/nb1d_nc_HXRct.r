library("MASS")

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

start_sim = Sys.time()
#for (SN in 296:250) {
SN = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("SN = ", SN))
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


names_sample <- names(data_sample)
data_full_prelim <- rbind(data_sample[, !names_sample %in% "y"],
                          data_out_of_sample[, !names_sample %in% "y"])
data_out_of_sampley <- data_out_of_sample[, !names_sample %in% "y"]
data_sampley <- data_sample[,!names_sample %in% "y"]

data_full <- data_full_prelim[, c(1:11, 13)]

################################################################
#IMpute ------------------------------------------------------------------------------
###############################################################

#       Method       MSE       Bias
#31 AIPW_Ecsda 0.1126749 0.10609411
#62 NIPW_Ecsda 0.1139954 0.10648829
#21    AIPW_E2 0.1194126 0.09578734
#63       OR_E 0.1197558 0.09638847
#52    NIPW_E2 0.1209563 0.09620760


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



############################################################################################
# Propensity score and m_0 and m_1                                                         #
############################################################################################
# E for P(A = 1|X)  #
#########################
ps_hat_X = c(data_sample_ps$ps_hat_X, data_out_of_sample_ps$ps_hat_X)
#which(ps_hat_E >= 0.98)
####################################################################################################################
# Best: AIPW MXMq
# Best: NIPW EMq
#########################
# Controls -------------------------------------
data_sample0 = data_sample[data_sample$A == 0, ]
obj_fit_y0 <- list(data_sample = data_sample0)

# Treated --------------------------------------
data_sample1 = data_sample[data_sample$A == 1, ]
obj_fit_y1 <- list(data_sample = data_sample1)

# Check and format data
model_formula_OR <- y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +X11 + (1|group)


# RF -----------------------------------------------------------------------------
formatted_data <- format_data(model_formula = model_formula_OR,
                              data_sample = data_sample,
                              data_out_of_sample  = data_out_of_sample)
data_full_RX <- rbind(formatted_data$X, formatted_data$X_newdata)
###################################
# Clustering, no tuning
#############################################
class(obj_fit_y0) <- "RF"
mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                               model_formula = model_formula_OR)


OR0_R <- fit_y(mutated_obj0,
               type_model = "continuous", 
               tune_RF = TRUE, 
               clust_RF = TRUE)

mu0_y_R <-  unname(unlist(predict(object = OR0_R$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))


# Treated --------------------------------------
class(obj_fit_y1) <- "RF"
mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                               model_formula = model_formula_OR)


OR1_R <- fit_y(mutated_obj1,
               type_model = "continuous", 
               tune_RF = TRUE,
               clust_RF = TRUE)

mu1_y_R <-  unname(unlist(predict(object = OR1_R$outcome_fit,
                                  newdata = data_full_RX,
                                  allow.new.levels = TRUE)))
# 

# 
# int
AIPW_dataE2 <- data.frame(
  y = impute_Ecsda$y_full_imputed,
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_X,
  mu0_y = mu0_y_R,
  mu1_y = mu1_y_R
)

AIPWfE2 <-
  as.data.frame(calculate_tau(AIPW_dataE2, type_tau = "AIPW"))

tau_hat_AIPW = AIPWfE2$tau



######################################################
# Bootstrapping
#######################################################
#Boot control
boot <- residual_bootstrap(y = data_sample$y,
                           y_hat = impute_Ecsda$y_hat_sample,
                           data_sample = data_sample,
                           n_boot  = 4000,
                           type_boot = c("br1"),
                           boot_seed = 10,
                           cent_resid = FALSE, 
                           c_Win = 3) 


n_boot = 4000
#tau_boot_OR <- matrix(0, nrow = n_boot, ncol = 41)
tau_boot_AIPW <- matrix(0, nrow = n_boot, ncol = 41)
tau_boot2_AIPW <- matrix(0, nrow = n_boot, ncol = 41)


y_full_imputedM <- matrix(0, nrow = n_boot, ncol = 41000) 

error_csda_b = numeric(n_boot)
start_boot <- Sys.time()
for (b in 1:n_boot) {
  
  #print(paste("b =", b))
  
  if(b %% 50 == 0) {
    # Print on the screen some message
    cat(paste0("b: ", b, "\n"))
  }
  
  y_boot <- boot$y_boot[[b]]
  
  data_sampleb = data_sample
  data_sampleb$y <- y_boot
  
  # Impute 
  impute_Ecsdab <-
    impute_y(
      model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A + (1+A|group),
      data_sample = data_sampleb,
      data_out_of_sample,
      method = "EBLUP",
      type_model = "gaussian"
    )
  y_full_imputedM[b, ] <- impute_Ecsdab$y_full_imputed
  smcsda = summary(impute_Ecsdab$outcome_fit)
  
  if (is.null(smcsda$optinfo$conv$lme4$messages)) {
    error_csda_b[b] <- 0
  } else {
    error_csda_b[b] <- 1
  }
  
} 

nw =  which(error_csda_b ==0)[1:1000]

if ( sum(is.na(nw)) !=0 ) {
  stop("Too few bootstraps.")
} 



start_boot <- Sys.time()
for (b in nw) {
  
  print(paste0("b: ", b))
  
  y_boot <- boot$y_boot[[b]]
  
  data_sampleb = data_sample
  data_sampleb$y <- y_boot
  
  # Controls -------------------------------------
  data_sample0b = data_sampleb[data_sampleb$A == 0, ]
  
  obj_fit_y0b <- list(data_sample = data_sample0b)

  
  # Treated --------------------------------------
  data_sample1b = data_sampleb[data_sampleb$A == 1, ]
  
  obj_fit_y1b <- list(data_sample = data_sample1b)

  #########################
  # R for E(Y|X, A = a)  #
  # RF -----------------------------------------------------------------------------
  formatted_datab <- format_data(model_formula = model_formula_OR,
                                 data_sample = data_sampleb,
                                 data_out_of_sample  = data_out_of_sample)
  data_full_RXb <- rbind(formatted_data$X, formatted_datab$X_newdata)
  ###################################
  # No clustering, no tuning
  #############################################
  startR <- Sys.time()
  class(obj_fit_y0b) <- "RF"
  mutated_obj0b <- mutate_obj_fit(obj_fit_y = obj_fit_y0b,
                                  model_formula = model_formula_OR)
  
  
  OR0_Rb <- fit_y(mutated_obj0b,
                  type_model = "continuous", 
                  tune_RF = TRUE, 
                  clust_RF = TRUE)
  
  mu0_y_Rb <-  unname(unlist(predict(object = OR0_Rb$outcome_fit,
                                     newdata = data_full_RXb,
                                     allow.new.levels = TRUE)))
  
  
  # Treated --------------------------------------
  class(obj_fit_y1b) <- "RF"
  mutated_obj1b <- mutate_obj_fit(obj_fit_y = obj_fit_y1b,
                                  model_formula = model_formula_OR)
  
  
  OR1_Rb <- fit_y(mutated_obj1b,
                  type_model = "continuous", 
                  tune_RF = TRUE,
                  clust_RF = TRUE)
  
  mu1_y_Rb <-  unname(unlist(predict(object = OR1_Rb$outcome_fit,
                                     newdata = data_full_RXb,
                                     allow.new.levels = TRUE)))
  stoopR <- Sys.time()
  
  ########################
  # X for P(A = 1|X)  #
  #########################
  #  ps_hat_E = c(data_sample_ps$ps_hat_E, data_out_of_sample_ps$ps_hat_E)
  #  ps_hat_X = c(data_sample_ps$ps_hat_X, data_out_of_sample_ps$ps_hat_X)
  impute_Ecsdab <-
    impute_y(
      model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A + (1+A|group),
      data_sample = data_sampleb,
      data_out_of_sample,
      method = "EBLUP",
      type_model = "gaussian"
    )
  
  smcsda = summary(impute_Ecsdab$outcome_fit)
  

  
  # y_full_imputeEb <- numeric()
  # A0_full <- which(data_full_prelim$A == 0)
  # A1_full <- which(data_full_prelim$A == 1)
  # y_full_imputeEb[A0_full] <- mu0_y_Eb[A0_full]
  # y_full_imputeEb[A1_full] <- mu1_y_Eb[A1_full]
  
  # AIPW MXMq
  
  AIPW_datab <- data.frame(
    y = impute_Ecsdab$y_full_imputed,
    A = c(data_sample$A, data_out_of_sample$A),
    group = c(data_sample$group, data_out_of_sample$group),
    p_score = ps_hat_X, 
    mu0_y = mu0_y_Rb, 
    mu1_y = mu1_y_Rb
  )
  
  AIPWfb <-
    as.data.frame(calculate_tau(AIPW_datab, type_tau = "AIPW"))
  
  tau_boot_AIPW[b, ] <- AIPWfb$tau
  
 
  ##################################
  # Double-bootstrap 
  ######################################

  bootb <- residual_bootstrap(y = data_sampleb$y,
                              y_hat = impute_Ecsdab$y_hat_sample,
                              data_sample = data_sample,
                              n_boot  = 1,
                              type_boot = c("br1"),
                              boot_seed = 11,
                              cent_resid = FALSE, 
                              c_Win = 3) 
  
  y_boot <- bootb$y_boot[[1]]
  #  data_sampleb = data_sample
  data_sampleb$y <- y_boot
  
  
  # Controls -------------------------------------
  data_sample0b = data_sampleb[data_sampleb$A == 0, ]
  
  obj_fit_y0b <- list(data_sample = data_sample0b)
  # class(obj_fit_y0b) <- "EBLUP"
  # mutated_obj0b <- mutate_obj_fit(obj_fit_y = obj_fit_y0b,
  #                                 model_formula = model_formula_OR)
  # 
  # OR0_Eb <- fit_y(mutated_obj0b,
  #                 type_model =  "gaussian")
  # 
  # mu0_y_Eb <-  unname(unlist(predict(object = OR0_Eb$outcome_fit,
  #                                    newdata = data_full_prelim,
  #                                    allow.new.levels = TRUE)))
  # 
  
  # Treated --------------------------------------
  data_sample1b = data_sampleb[data_sampleb$A == 1, ]
  
  obj_fit_y1b <- list(data_sample = data_sample1b)
  # class(obj_fit_y1b) <- "EBLUP"
  # mutated_obj1b <- mutate_obj_fit(obj_fit_y = obj_fit_y1b,
  #                                 model_formula = model_formula_OR)
  # 
  # OR1_Eb <- fit_y(mutated_obj1b,
  #                 type_model =  "gaussian")
  # 
  # mu1_y_Eb <-  unname(unlist(predict(object = OR1_Eb$outcome_fit,
  #                                    newdata = data_full_prelim,
  #                                    allow.new.levels = TRUE)))
  # 
  # 
  # # OR E Boot
  # df_Eb <- data.frame(mu1_y = mu1_y_Eb,
  #                     mu0_y = mu0_y_Eb,
  #                     group = c(data_sample$group, 
  #                               data_out_of_sample$group))
  # 
  # tau_treatEb = aggregate(df_Eb$mu1_y, list(df_Eb$group), FUN = mean)$x
  # tau_untreatEb = aggregate(df_Eb$mu0_y, list(df_Eb$group), FUN = mean)$x
  #  tau_Eb = tau_treatEb - tau_untreatEb
  
  # tau_boot_OR[b, ] <- tau_treatEb - tau_untreatEb
  
  # R for E(Y|X, A = a)  #
  # RF -----------------------------------------------------------------------------
  formatted_datab <- format_data(model_formula = model_formula_OR,
                                 data_sample = data_sampleb,
                                 data_out_of_sample  = data_out_of_sample)
  data_full_RXb <- rbind(formatted_data$X, formatted_datab$X_newdata)
  ###################################
  # No clustering, no tuning
  #############################################
  startR <- Sys.time()
  class(obj_fit_y0b) <- "RF"
  mutated_obj0b <- mutate_obj_fit(obj_fit_y = obj_fit_y0b,
                                  model_formula = model_formula_OR)
  
  
  OR0_Rb <- fit_y(mutated_obj0b,
                  type_model = "continuous", 
                  tune_RF = TRUE, 
                  clust_RF = TRUE)
  
  mu0_y_Rb <-  unname(unlist(predict(object = OR0_Rb$outcome_fit,
                                     newdata = data_full_RXb,
                                     allow.new.levels = TRUE)))
  
  
  # Treated --------------------------------------
  class(obj_fit_y1b) <- "RF"
  mutated_obj1b <- mutate_obj_fit(obj_fit_y = obj_fit_y1b,
                                  model_formula = model_formula_OR)
  
  
  OR1_Rb <- fit_y(mutated_obj1b,
                  type_model = "continuous", 
                  tune_RF = TRUE,
                  clust_RF = TRUE)
  
  mu1_y_Rb <-  unname(unlist(predict(object = OR1_Rb$outcome_fit,
                                     newdata = data_full_RXb,
                                     allow.new.levels = TRUE)))
  stoopR <- Sys.time()
  
  #########################
  # X for P(A = 1|X)  #
  #########################
  impute_Ecsdab <-
    impute_y(
      model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A + (1+A|group),
      data_sample = data_sampleb,
      data_out_of_sample,
      method = "EBLUP",
      type_model = "gaussian"
    )
  
  smcsda = summary(impute_Ecsdab$outcome_fit)
  

  
  # y_full_imputeEb <- numeric()
  # A0_full <- which(data_full_prelim$A == 0)
  # A1_full <- which(data_full_prelim$A == 1)
  # y_full_imputeEb[A0_full] <- mu0_y_Eb[A0_full]
  # y_full_imputeEb[A1_full] <- mu1_y_Eb[A1_full]
  
  # AIPW MXMq
  
  AIPW_datab <- data.frame(
    y = impute_Ecsdab$y_full_imputed,
    A = c(data_sample$A, data_out_of_sample$A),
    group = c(data_sample$group, data_out_of_sample$group),
    p_score = ps_hat_X, 
    mu0_y = mu0_y_Rb, 
    mu1_y = mu1_y_Rb
  )
  
  AIPWfb <-
    as.data.frame(calculate_tau(AIPW_datab, type_tau = "AIPW"))
  
  tau_boot2_AIPW[b, ] <- AIPWfb$tau

}
stop_boot <- Sys.time()
stop_boot - start_boot


#### OR ######
#  MSE
# tau_boot_dif_OR <- tau_boot_OR - matrix(rep(tau_hat_OR, n_boot), nrow = n_boot, ncol = 41, byrow = TRUE)
# tau_boot_s_OR <- (tau_boot_dif_OR)^2
# MSE_boot_OR <- colMeans(tau_boot_s_OR)
# 
# # CI with MSE
# upInt_OR <- tau_hat_OR + sqrt(MSE_boot_OR) * 1.96
# loInt_OR <- tau_hat_OR - sqrt(MSE_boot_OR) * 1.96
# covB_OR <- (loInt_OR <= tau_true ) & (upInt_OR >= tau_true )
# lenB_OR <- upInt_OR - loInt_OR
# 
# # CI with percentiles
# tau_boot_025_OR <- apply(tau_boot_OR, 2,  quantile, 0.025)
# tau_boot_975_OR <- apply(tau_boot_OR, 2,  quantile, 0.975)
# covBper_OR <- (tau_boot_025_OR <= tau_true ) & (tau_boot_975_OR >= tau_true )
# lenBper_OR <- tau_boot_975_OR - tau_boot_025_OR
# 

#### AIPW ######
#  MSE
tau_boot_dif_AIPW <- tau_boot_AIPW[nw, ] - matrix(rep(tau_hat_AIPW, 1000), nrow = 1000, ncol = 41, byrow = TRUE)
tau_boot_s_AIPW <- (tau_boot_dif_AIPW)^2
MSE_boot_AIPW <- colMeans(tau_boot_s_AIPW)

tau_boot2_dif_AIPW <- tau_boot2_AIPW[nw, ]  - tau_boot_AIPW[nw, ] 
tau_boot2_s_AIPW <- (tau_boot2_dif_AIPW)^2
MSE_boot2_AIPW <- colMeans(tau_boot2_s_AIPW)

#MSE_boot_AIPWbc <- 2 * MSE_boot_AIPW - MSE_boot2_AIPW

Bias_boot_AIPWbc <- 2 * colMeans(tau_boot_dif_AIPW) - colMeans(tau_boot2_dif_AIPW)
Bias_boot_AIPW <-  colMeans(tau_boot_dif_AIPW)

MSE_boot_AIPWbc = numeric(m) 
for (i in 1:m){
  if (MSE_boot_AIPW[i] >= MSE_boot2_AIPW[i] ) {
    MSE_boot_AIPWbc[i] = 2 * MSE_boot_AIPW[i] - MSE_boot2_AIPW[i]
  } else {
    MSE_boot_AIPWbc[i] = MSE_boot_AIPW[i]*exp((MSE_boot_AIPW[i]-MSE_boot2_AIPW[i])/MSE_boot2_AIPW[i])}
}

# CI with MSE
upInt_AIPW <- tau_hat_AIPW + sqrt(MSE_boot_AIPW) * 1.96
loInt_AIPW <- tau_hat_AIPW - sqrt(MSE_boot_AIPW) * 1.96
covB_AIPW <- (loInt_AIPW <= tau_true ) & (upInt_AIPW >= tau_true )
lenB_AIPW <- upInt_AIPW - loInt_AIPW

#bias corrected MSE
upInt_AIPWbc <- tau_hat_AIPW + sqrt(MSE_boot_AIPWbc) * 1.96
loInt_AIPWbc <- tau_hat_AIPW - sqrt(MSE_boot_AIPWbc) * 1.96
covB_AIPWbc <- (loInt_AIPWbc <= tau_true ) & (upInt_AIPWbc >= tau_true )
lenB_AIPWbc <- upInt_AIPWbc - loInt_AIPWbc

# bias corrected estimator 1
upInt_AIPWb <- tau_hat_AIPW + Bias_boot_AIPW + sqrt(MSE_boot_AIPW) * 1.96
loInt_AIPWb <- tau_hat_AIPW + Bias_boot_AIPW - sqrt(MSE_boot_AIPW) * 1.96
covB_AIPWb <- (loInt_AIPWb <= tau_true ) & (upInt_AIPWb >= tau_true )
lenB_AIPWb <- upInt_AIPWb - loInt_AIPWb

# bias corrected estimator 2
upInt_AIPWb2 <- tau_hat_AIPW + Bias_boot_AIPWbc + sqrt(MSE_boot_AIPW) * 1.96
loInt_AIPWb2 <- tau_hat_AIPW + Bias_boot_AIPWbc - sqrt(MSE_boot_AIPW) * 1.96
covB_AIPWb2 <- (loInt_AIPWb2 <= tau_true ) & (upInt_AIPWb2 >= tau_true )
lenB_AIPWb2 <- upInt_AIPWb2 - loInt_AIPWb2

# bias corrected estimator 1 + MSE
upInt_AIPWbbc <- tau_hat_AIPW + Bias_boot_AIPW + sqrt(MSE_boot_AIPWbc) * 1.96
loInt_AIPWbbc <- tau_hat_AIPW + Bias_boot_AIPW - sqrt(MSE_boot_AIPWbc) * 1.96
covB_AIPWbbc <- (loInt_AIPWbbc <= tau_true ) & (upInt_AIPWbbc >= tau_true )
lenB_AIPWbbc <- upInt_AIPWbbc - loInt_AIPWbbc


# bias corrected estimator 2 + MSE
upInt_AIPWb2bc <- tau_hat_AIPW + Bias_boot_AIPWbc + sqrt(MSE_boot_AIPWbc) * 1.96
loInt_AIPWb2bc <- tau_hat_AIPW + Bias_boot_AIPWbc - sqrt(MSE_boot_AIPWbc) * 1.96
covB_AIPWb2bc <- (loInt_AIPWbbc <= tau_true ) & (upInt_AIPWbbc >= tau_true )
lenB_AIPWb2bc <- upInt_AIPWb2bc - loInt_AIPWb2bc


# CI with percentiles
tau_boot_025_AIPW <- apply(tau_boot_AIPW, 2,  quantile, 0.025)
tau_boot_975_AIPW <- apply(tau_boot_AIPW, 2,  quantile, 0.975)
covBper_AIPW <- (tau_boot_025_AIPW <= tau_true ) & (tau_boot_975_AIPW >= tau_true )
lenBper_AIPW <- tau_boot_975_AIPW - tau_boot_025_AIPW


tau_boot_025_AIPWb <- apply(tau_boot_AIPW + Bias_boot_AIPW, 2,  quantile, 0.025)
tau_boot_975_AIPWb <- apply(tau_boot_AIPW + Bias_boot_AIPW, 2,  quantile, 0.975)
covBper_AIPWb <- (tau_boot_025_AIPWb <= tau_true ) & (tau_boot_975_AIPWb >= tau_true )
lenBper_AIPWb <- tau_boot_975_AIPWb - tau_boot_025_AIPWb

tau_boot_025_AIPWb2 <- apply(tau_boot_AIPW + Bias_boot_AIPWbc, 2,  quantile, 0.025)
tau_boot_975_AIPWb2 <- apply(tau_boot_AIPW + Bias_boot_AIPWbc, 2,  quantile, 0.975)
covBper_AIPWb2 <- (tau_boot_025_AIPWb2 <= tau_true ) & (tau_boot_975_AIPWb2 >= tau_true )
lenBper_AIPWb2 <- tau_boot_975_AIPWb2 - tau_boot_025_AIPWb2


Results = list(
  ni = ni,
  n = sum(ni),
  
  Ni = Ni,
  N = sum(Ni),
  frac_nN = frac_nN,
  
  
  #Tau
  tau_true = tau_true,
  tau_hat_AIPW = tau_hat_AIPW,
  #               tau_hat_OR = tau_hat_OR,
  
  # Boot AIPW
  #  MSE
  tau_boot_dif_AIPW  = tau_boot_dif_AIPW,
  tau_boot_s_AIPW  = tau_boot_s_AIPW,
  MSE_boot_AIPW = MSE_boot_AIPW,
  
  tau_boot2_dif_AIPW = tau_boot2_dif_AIPW,
  tau_boot2_s_AIPW = tau_boot2_s_AIPW,
  MSE_boot2_AIPW = MSE_boot2_AIPW,
  
  Bias_boot_AIPWbc = Bias_boot_AIPWbc,
  Bias_boot_AIPW = Bias_boot_AIPW,
  
  # CI with MSE
  upInt_AIPW  = upInt_AIPW,
  loInt_AIPW = loInt_AIPW,
  covB_AIPW = covB_AIPW,
  lenB_AIPW = lenB_AIPW,
  
  #bias corrected MSE
  upInt_AIPWbc = upInt_AIPWbc,
  loInt_AIPWbc = loInt_AIPWbc,
  covB_AIPWbc  = covB_AIPWbc,
  lenB_AIPWbc = lenB_AIPWbc,
  
  # bias corrected estimator 1
  upInt_AIPWb  = upInt_AIPWb,
  loInt_AIPWb = loInt_AIPWb,
  covB_AIPWb  = covB_AIPWb,
  lenB_AIPWb  = lenB_AIPWb,
  
  # bias corrected estimator 2
  upInt_AIPWb2 = upInt_AIPWb2,
  loInt_AIPWb2 = loInt_AIPWb2,
  covB_AIPWb2 = covB_AIPWb2,
  lenB_AIPWb2 = lenB_AIPWb2,
  
  # bias corrected estimator 1 + MSE
  upInt_AIPWbbc = upInt_AIPWbbc,
  loInt_AIPWbbc = loInt_AIPWbbc,
  covB_AIPWbbc = covB_AIPWbbc,
  lenB_AIPWbbc = lenB_AIPWbbc,
  
  # bias corrected estimator 2 + MSE
  upInt_AIPWb2bc = upInt_AIPWb2bc,
  loInt_AIPWb2bc = loInt_AIPWb2bc,
  covB_AIPWb2bc = covB_AIPWb2bc,
  lenB_AIPWb2bc = lenB_AIPWb2bc,
  
  # Percentiles
  # CI with percentiles
  tau_boot_025_AIPW = tau_boot_025_AIPW,
  tau_boot_975_AIPW = tau_boot_975_AIPW,
  covBper_AIPW = covBper_AIPW,
  lenBper_AIPW = lenBper_AIPW,
  
  tau_boot_025_AIPWb = tau_boot_025_AIPWb,
  tau_boot_975_AIPWb  = tau_boot_975_AIPWb,
  covBper_AIPWb = covBper_AIPWb,
  lenBper_AIPWb = lenBper_AIPWb,
  
  tau_boot_025_AIPWb2  = tau_boot_025_AIPWb2,
  tau_boot_975_AIPWb2 = tau_boot_975_AIPWb2,
  covBper_AIPWb2  = covBper_AIPWb2,
  lenBper_AIPWb2  = lenBper_AIPWb2,
  

  error_csda = error_csda
)

outputName = paste("nb1d_nc_HXRct_", SN, ".RData", sep = "")
outputPath = file.path("/user/work/ze23696/Comp", outputName)
save("Results", file = outputPath)