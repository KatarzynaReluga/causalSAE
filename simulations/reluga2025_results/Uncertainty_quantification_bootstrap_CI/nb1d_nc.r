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


############################################################################################
# Propensity score and m_0 and m_1                                                         #
############################################################################################
ps_hat_E = c(data_sample_ps$ps_hat_E, data_out_of_sample_ps$ps_hat_E)

####################################################################################################################
# Best: AIPW MXMq
# Best: NIPW EMq
#########################
# Controls -------------------------------------
data_sample0 = data_sample[data_sample$A == 0, ]
obj_fit_y0 <- list(data_sample = data_sample0)

# Check and format data
model_formula_OR <- y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +X11 + (1|group)


class(obj_fit_y0) <- "MQ"
mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                               model_formula = model_formula_OR)


OR0_M <- fit_y(mutated_obj0,
               type_model =  "continuous")

mu0_hat_M <-  unname(unlist(predict(object = OR0_M$outcome_fit,
                                    newdata = data_full_prelim,
                                    allow.new.levels = TRUE)))


# Treated --------------------------------------
data_sample1 = data_sample[data_sample$A == 1, ]
obj_fit_y1 <- list(data_sample = data_sample1)

class(obj_fit_y1) <- "MQ"
mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                               model_formula = model_formula_OR)


OR1_M <- fit_y(mutated_obj1,
               type_model = "continuous")

mu1_hat_M <-  unname(unlist(predict(object = OR1_M$outcome_fit,
                                    newdata = data_full_prelim,
                                    allow.new.levels = TRUE)))


################################################################
# EBLUP ------------------------------------------------------------------------------
###############################################################

# Controls -------------------------------------
class(obj_fit_y0) <- "EBLUP"
mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                               model_formula = model_formula_OR)

OR0_E <- fit_y(mutated_obj0,
               type_model =  "gaussian")

mu0_y_E <-  unname(unlist(predict(object = OR0_E$outcome_fit,
                                  newdata = data_full_prelim,
                                  allow.new.levels = TRUE)))


# Treated --------------------------------------

obj_fit_y1 <- list(data_sample = data_sample1)
class(obj_fit_y1) <- "EBLUP"
mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                               model_formula = model_formula_OR)


OR1_E <- fit_y(mutated_obj1,
               type_model = "gaussian")

mu1_y_E <-  unname(unlist(predict(object = OR1_E$outcome_fit,
                                  newdata = data_full_prelim,
                                  allow.new.levels = TRUE)))


# OR E 
df_E <- data.frame(mu1_y = mu1_y_E,
                   mu0_y = mu0_y_E,
                   group = c(data_sample$group, 
                             data_out_of_sample$group))

tau_treatE = aggregate(df_E$mu1_y, list(df_E$group), FUN = mean)$x
tau_untreatE = aggregate(df_E$mu0_y, list(df_E$group), FUN = mean)$x
tau_hat_OR = tau_treatE - tau_untreatE

#       Method       MSE       Bias
#31 AIPW_Ecsda 0.1126749 0.10609411
#62 NIPW_Ecsda 0.1139954 0.10648829
#21    AIPW_E2 0.1194126 0.09578734
#63       OR_E 0.1197558 0.09638847
#52    NIPW_E2 0.1209563 0.09620760

sm = summary(OR0_E$outcome_fit)

if (is.null(sm$optinfo$conv$lme4$messages)) {
  error_E <- 0
} else {
  error_E <- 1 #1 means error
}


y_full_imputeE2 <- numeric()
A0_full <- which(data_full_prelim$A == 0)
A1_full <- which(data_full_prelim$A == 1)
y_full_imputeE2[A0_full] <- mu0_y_E[A0_full]
y_full_imputeE2[A1_full] <- mu1_y_E[A1_full]

# int
AIPW_dataE2 <- data.frame(
  y = y_full_imputeE2,
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_E,
  mu0_y = mu0_hat_M,
  mu1_y = mu1_hat_M
)

AIPWfE2 <-
  as.data.frame(calculate_tau(AIPW_dataE2, type_tau = "AIPW"))

tau_hat_AIPW = AIPWfE2$tau

IPW_dataE2 <- data.frame(
  y = y_full_imputeE2,
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_E
)

IPWfE2 <-
  as.data.frame(calculate_tau(IPW_dataE2, type_tau = "H"))

tau_hat_NIPW = IPWfE2$tau

#######################################################
# Bootstrapping
#######################################################
# Bootstrapping using EBLUP, the same y_hat for all methods
model_formula_OR_f <- y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A+ (1|group)

obj_fit_y <- list(data_sample = data_sample)

class(obj_fit_y) <- "EBLUP"
mutated_obj <- mutate_obj_fit(obj_fit_y = obj_fit_y,
                              model_formula = model_formula_OR_f)

OR_E <- fit_y(mutated_obj, type_model =  "gaussian")


y_hat <- unname(unlist(predict(object = OR_E$outcome_fit,
                               newdata = data_sample,
                               allow.new.levels = TRUE)))


#Boot control
boot <- residual_bootstrap(y = data_sample$y,
                           y_hat = y_hat,
                           data_sample = data_sample,
                           n_boot  = 1000,
                           type_boot = c("br1"),
                           boot_seed = 10,
                           cent_resid = FALSE, 
                           c_Win = 3) 

n_boot = 100
tau_boot_OR <- matrix(0, nrow = n_boot, ncol = 41)
tau_boot_AIPW <- matrix(0, nrow = n_boot, ncol = 41)
tau_boot_NIPW <- matrix(0, nrow = n_boot, ncol = 41)

tau_boot2_OR <- matrix(0, nrow = n_boot, ncol = 41)
tau_boot2_AIPW <- matrix(0, nrow = n_boot, ncol = 41)
tau_boot2_NIPW <- matrix(0, nrow = n_boot, ncol = 41)

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
  
  # Controls -------------------------------------
  data_sample0b = data_sampleb[data_sampleb$A == 0, ]
  
  obj_fit_y0b <- list(data_sample = data_sample0b)
  class(obj_fit_y0b) <- "EBLUP"
  mutated_obj0b <- mutate_obj_fit(obj_fit_y = obj_fit_y0b,
                                  model_formula = model_formula_OR)
  
  OR0_Eb <- fit_y(mutated_obj0b,
                  type_model =  "gaussian")
  
  mu0_y_Eb <-  unname(unlist(predict(object = OR0_Eb$outcome_fit,
                                     newdata = data_full_prelim,
                                     allow.new.levels = TRUE)))
  
  
  # Treated --------------------------------------
  data_sample1b = data_sampleb[data_sampleb$A == 1, ]
  
  obj_fit_y1b <- list(data_sample = data_sample1b)
  class(obj_fit_y1b) <- "EBLUP"
  mutated_obj1b <- mutate_obj_fit(obj_fit_y = obj_fit_y1b,
                                  model_formula = model_formula_OR)
  
  OR1_Eb <- fit_y(mutated_obj1b,
                  type_model =  "gaussian")
  
  mu1_y_Eb <-  unname(unlist(predict(object = OR1_Eb$outcome_fit,
                                     newdata = data_full_prelim,
                                     allow.new.levels = TRUE)))
  
  
  # OR E Boot
  df_Eb <- data.frame(mu1_y = mu1_y_Eb,
                      mu0_y = mu0_y_Eb,
                      group = c(data_sample$group, 
                                data_out_of_sample$group))
  
  tau_treatEb = aggregate(df_Eb$mu1_y, list(df_Eb$group), FUN = mean)$x
  tau_untreatEb = aggregate(df_Eb$mu0_y, list(df_Eb$group), FUN = mean)$x
  #  tau_Eb = tau_treatEb - tau_untreatEb
  
  tau_boot_OR[b, ] <- tau_treatEb - tau_untreatEb

  #########################
  # MQ for E(Y|X, A = a)  #
  #########################
  # Controls -------------------------------------

  class(obj_fit_y0b) <- "MQ"
  mutated_obj0b <- mutate_obj_fit(obj_fit_y = obj_fit_y0b,
                                  model_formula = model_formula_OR)
  
  
  OR0_Mb <- fit_y(mutated_obj0b,
                  type_model =  "continuous")
  
  mu0_hat_Mb <-  unname(unlist(predict(object = OR0_Mb$outcome_fit,
                                       newdata = data_full_prelim,
                                       allow.new.levels = TRUE)))
  
  
  # Treated --------------------------------------
  #data_sample1b = data_sampleb[data_sampleb$A == 1, ]
  #obj_fit_y1b <- list(data_sample = data_sample1b)
  
  
  class(obj_fit_y1b) <- "MQ"
  mutated_obj1b <- mutate_obj_fit(obj_fit_y = obj_fit_y1b,
                                  model_formula = model_formula_OR)
  
  
  OR1_Mb <- fit_y(mutated_obj1b,
                  type_model =  "continuous")
  
  mu1_hat_Mb <-  unname(unlist(predict(object = OR1_Mb$outcome_fit,
                                       newdata = data_full_prelim,
                                       allow.new.levels = TRUE)))
  
  #########################
  # X for P(A = 1|X)  #
  #########################
  ps_hat_E = c(data_sample_ps$ps_hat_E, data_out_of_sample_ps$ps_hat_E)
  
 
   y_full_imputeEb <- numeric()
   A0_full <- which(data_full_prelim$A == 0)
   A1_full <- which(data_full_prelim$A == 1)
   y_full_imputeEb[A0_full] <- mu0_y_Eb[A0_full]
   y_full_imputeEb[A1_full] <- mu1_y_Eb[A1_full]

  # AIPW MXMq
  
  AIPW_datab <- data.frame(
    y = y_full_imputeEb,
    A = c(data_sample$A, data_out_of_sample$A),
    group = c(data_sample$group, data_out_of_sample$group),
    p_score = ps_hat_E, 
    mu0_y = mu0_hat_Mb, 
    mu1_y = mu1_hat_Mb
  )
  
  AIPWfb <-
    as.data.frame(calculate_tau(AIPW_datab, type_tau = "AIPW"))
  
  tau_boot_AIPW[b, ] <- AIPWfb$tau
  
  
  # NIPW EMq
  
  IPW_datab <- data.frame(
    y = y_full_imputeEb,
    A = c(data_sample$A, data_out_of_sample$A),
    group = c(data_sample$group, data_out_of_sample$group),
    p_score = ps_hat_E
  )
  
  IPWfb <-
    as.data.frame(calculate_tau(IPW_datab, type_tau = "H"))
  
  tau_boot_NIPW[b, ] <- IPWfb$tau
  ##################################
  # Double-bootstrap 
  ######################################
#  model_formula_OR_f <- y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A+ (1|group)
  
  obj_fit_yb <- list(data_sample = data_sampleb)
  
  class(obj_fit_yb) <- "EBLUP"
  mutated_objb <- mutate_obj_fit(obj_fit_y = obj_fit_yb,
                                model_formula = model_formula_OR_f)
  
  OR_Eb <- fit_y(mutated_objb, type_model =  "gaussian")
  
  
  y_hatb <- unname(unlist(predict(object = OR_Eb$outcome_fit,
                                 newdata = data_sampleb,
                                 allow.new.levels = TRUE)))
  
  
  bootb <- residual_bootstrap(y = data_sampleb$y,
                             y_hat = y_hatb,
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
  class(obj_fit_y0b) <- "EBLUP"
  mutated_obj0b <- mutate_obj_fit(obj_fit_y = obj_fit_y0b,
                                  model_formula = model_formula_OR)
  
  OR0_Eb <- fit_y(mutated_obj0b,
                  type_model =  "gaussian")
  
  mu0_y_Eb <-  unname(unlist(predict(object = OR0_Eb$outcome_fit,
                                     newdata = data_full_prelim,
                                     allow.new.levels = TRUE)))
  
  
  # Treated --------------------------------------
  data_sample1b = data_sampleb[data_sampleb$A == 1, ]
  
  obj_fit_y1b <- list(data_sample = data_sample1b)
  class(obj_fit_y1b) <- "EBLUP"
  mutated_obj1b <- mutate_obj_fit(obj_fit_y = obj_fit_y1b,
                                  model_formula = model_formula_OR)
  
  OR1_Eb <- fit_y(mutated_obj1b,
                  type_model =  "gaussian")
  
  mu1_y_Eb <-  unname(unlist(predict(object = OR1_Eb$outcome_fit,
                                     newdata = data_full_prelim,
                                     allow.new.levels = TRUE)))
  
  
  # OR E Boot
  df_Eb <- data.frame(mu1_y = mu1_y_Eb,
                      mu0_y = mu0_y_Eb,
                      group = c(data_sample$group, 
                                data_out_of_sample$group))
  
  tau_treatEb = aggregate(df_Eb$mu1_y, list(df_Eb$group), FUN = mean)$x
  tau_untreatEb = aggregate(df_Eb$mu0_y, list(df_Eb$group), FUN = mean)$x
  #  tau_Eb = tau_treatEb - tau_untreatEb
  
  tau_boot2_OR[b, ] <- tau_treatEb - tau_untreatEb
  
  #########################
  # MQ for E(Y|X, A = a)  #
  #########################
  # Controls -------------------------------------
  
  class(obj_fit_y0b) <- "MQ"
  mutated_obj0b <- mutate_obj_fit(obj_fit_y = obj_fit_y0b,
                                  model_formula = model_formula_OR)
  
  
  OR0_Mb <- fit_y(mutated_obj0b,
                  type_model =  "continuous")
  
  mu0_hat_Mb <-  unname(unlist(predict(object = OR0_Mb$outcome_fit,
                                       newdata = data_full_prelim,
                                       allow.new.levels = TRUE)))
  
  
  # Treated --------------------------------------
  #data_sample1b = data_sampleb[data_sampleb$A == 1, ]
  #obj_fit_y1b <- list(data_sample = data_sample1b)
  
  
  class(obj_fit_y1b) <- "MQ"
  mutated_obj1b <- mutate_obj_fit(obj_fit_y = obj_fit_y1b,
                                  model_formula = model_formula_OR)
  
  
  OR1_Mb <- fit_y(mutated_obj1b,
                  type_model =  "continuous")
  
  mu1_hat_Mb <-  unname(unlist(predict(object = OR1_Mb$outcome_fit,
                                       newdata = data_full_prelim,
                                       allow.new.levels = TRUE)))
  
  #########################
  # X for P(A = 1|X)  #
  #########################
  ps_hat_E = c(data_sample_ps$ps_hat_E, data_out_of_sample_ps$ps_hat_E)
  
  
  y_full_imputeEb <- numeric()
  A0_full <- which(data_full_prelim$A == 0)
  A1_full <- which(data_full_prelim$A == 1)
  y_full_imputeEb[A0_full] <- mu0_y_Eb[A0_full]
  y_full_imputeEb[A1_full] <- mu1_y_Eb[A1_full]
  
  # AIPW MXMq
  
  AIPW_datab <- data.frame(
    y = y_full_imputeEb,
    A = c(data_sample$A, data_out_of_sample$A),
    group = c(data_sample$group, data_out_of_sample$group),
    p_score = ps_hat_E, 
    mu0_y = mu0_hat_Mb, 
    mu1_y = mu1_hat_Mb
  )
  
  AIPWfb <-
    as.data.frame(calculate_tau(AIPW_datab, type_tau = "AIPW"))
  
  tau_boot2_AIPW[b, ] <- AIPWfb$tau
  
  
  # NIPW EMq
  
  IPW_datab <- data.frame(
    y = y_full_imputeEb,
    A = c(data_sample$A, data_out_of_sample$A),
    group = c(data_sample$group, data_out_of_sample$group),
    p_score = ps_hat_E
  )
  
  IPWfb <-
    as.data.frame(calculate_tau(IPW_datab, type_tau = "H"))
  
  tau_boot2_NIPW[b, ] <- IPWfb$tau
  
  
}
stop_boot <- Sys.time()
stop_boot - start_boot

#### OR ######
#  MSE
tau_boot_dif_OR <- tau_boot_OR - matrix(rep(tau_hat_OR, n_boot), nrow = n_boot, ncol = 41, byrow = TRUE)
tau_boot_s_OR <- (tau_boot_dif_OR)^2
MSE_boot_OR <- colMeans(tau_boot_s_OR)

tau_boot2_dif_OR <- tau_boot2_OR - tau_boot_OR
tau_boot2_s_OR <- (tau_boot2_dif_OR)^2
MSE_boot2_OR <- colMeans(tau_boot2_s_OR)

#MSE_boot_ORbc <- 2 * MSE_boot_OR - MSE_boot2_OR

Bias_boot_ORbc <- 2 * colMeans(tau_boot_dif_OR) - colMeans(tau_boot2_dif_OR)
Bias_boot_OR <-  colMeans(tau_boot_dif_OR)

MSE_boot_ORbc = numeric(m) 
for (i in 1:m){
  if (MSE_boot_OR[i] >= MSE_boot2_OR[i] ) {
    MSE_boot_ORbc[i] = 2 * MSE_boot_OR[i] - MSE_boot2_OR[i]
  } else {
    MSE_boot_ORbc[i] = MSE_boot_OR[i]*exp((MSE_boot_OR[i]-MSE_boot2_OR[i])/MSE_boot2_OR[i])}
  }

############################################################
# CI with MSE
###############################################################
#plain
upInt_OR <- tau_hat_OR + sqrt(MSE_boot_OR) * 1.96
loInt_OR <- tau_hat_OR - sqrt(MSE_boot_OR) * 1.96
covB_OR <- (loInt_OR <= tau_true ) & (upInt_OR >= tau_true )
lenB_OR <- upInt_OR - loInt_OR
#bias corrected MSE
upInt_ORbc <- tau_hat_OR + sqrt(MSE_boot_ORbc) * 1.96
loInt_ORbc <- tau_hat_OR - sqrt(MSE_boot_ORbc) * 1.96
covB_ORbc <- (loInt_ORbc <= tau_true ) & (upInt_ORbc >= tau_true )
lenB_ORbc <- upInt_ORbc - loInt_ORbc

# bias corrected estimator 1
upInt_ORb <- tau_hat_OR + Bias_boot_OR + sqrt(MSE_boot_OR) * 1.96
loInt_ORb <- tau_hat_OR + Bias_boot_OR - sqrt(MSE_boot_OR) * 1.96
covB_ORb <- (loInt_ORb <= tau_true ) & (upInt_ORb >= tau_true )
lenB_ORb <- upInt_ORb - loInt_ORb

# bias corrected estimator 2
upInt_ORb2 <- tau_hat_OR + Bias_boot_ORbc + sqrt(MSE_boot_OR) * 1.96
loInt_ORb2 <- tau_hat_OR + Bias_boot_ORbc - sqrt(MSE_boot_OR) * 1.96
covB_ORb2 <- (loInt_ORb2 <= tau_true ) & (upInt_ORb2 >= tau_true )
lenB_ORb2 <- upInt_ORb2 - loInt_ORb2

# bias corrected estimator 1 + MSE
upInt_ORbbc <- tau_hat_OR + Bias_boot_OR + sqrt(MSE_boot_ORbc) * 1.96
loInt_ORbbc <- tau_hat_OR + Bias_boot_OR - sqrt(MSE_boot_ORbc) * 1.96
covB_ORbbc <- (loInt_ORbbc <= tau_true ) & (upInt_ORbbc >= tau_true )
lenB_ORbbc <- upInt_ORbbc - loInt_ORbbc


# bias corrected estimator 2 + MSE
upInt_ORb2bc <- tau_hat_OR + Bias_boot_ORbc + sqrt(MSE_boot_ORbc) * 1.96
loInt_ORb2bc <- tau_hat_OR + Bias_boot_ORbc - sqrt(MSE_boot_ORbc) * 1.96
covB_ORb2bc <- (loInt_ORbbc <= tau_true ) & (upInt_ORbbc >= tau_true )
lenB_ORb2bc <- upInt_ORb2bc - loInt_ORb2bc

##########################################
# CI with percentiles
##############################################
tau_boot_025_OR <- apply(tau_boot_OR, 2,  quantile, 0.025)
tau_boot_975_OR <- apply(tau_boot_OR, 2,  quantile, 0.975)
covBper_OR <- (tau_boot_025_OR <= tau_true ) & (tau_boot_975_OR >= tau_true )
lenBper_OR <- tau_boot_975_OR - tau_boot_025_OR


tau_boot_025_ORb <- apply(tau_boot_OR + Bias_boot_OR, 2,  quantile, 0.025)
tau_boot_975_ORb <- apply(tau_boot_OR + Bias_boot_OR, 2,  quantile, 0.975)
covBper_ORb <- (tau_boot_025_ORb <= tau_true ) & (tau_boot_975_ORb >= tau_true )
lenBper_ORb <- tau_boot_975_ORb - tau_boot_025_ORb

tau_boot_025_ORb2 <- apply(tau_boot_OR + Bias_boot_ORbc, 2,  quantile, 0.025)
tau_boot_975_ORb2 <- apply(tau_boot_OR + Bias_boot_ORbc, 2,  quantile, 0.975)
covBper_ORb2 <- (tau_boot_025_ORb2 <= tau_true ) & (tau_boot_975_ORb2 >= tau_true )
lenBper_ORb2 <- tau_boot_975_ORb2 - tau_boot_025_ORb2

#################################################################
#### AIPW ######
####################################################################

#  MSE
tau_boot_dif_AIPW <- tau_boot_AIPW - matrix(rep(tau_hat_AIPW, n_boot), nrow = n_boot, ncol = 41, byrow = TRUE)
tau_boot_s_AIPW <- (tau_boot_dif_AIPW)^2
MSE_boot_AIPW <- colMeans(tau_boot_s_AIPW)

tau_boot2_dif_AIPW <- tau_boot2_AIPW - tau_boot_AIPW
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

######################################################
# CI with MSE
#######################################################
#plain
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


##########################################
# CI with percentiles
##############################################
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

###############################################
#### NIPW ######
#  MSE
tau_boot_dif_NIPW <- tau_boot_NIPW - matrix(rep(tau_hat_NIPW, n_boot), nrow = n_boot, ncol = 41, byrow = TRUE)
tau_boot_s_NIPW <- (tau_boot_dif_NIPW)^2
MSE_boot_NIPW <- colMeans(tau_boot_s_NIPW)

tau_boot2_dif_NIPW <- tau_boot2_NIPW - tau_boot_NIPW
tau_boot2_s_NIPW <- (tau_boot2_dif_NIPW)^2
MSE_boot2_NIPW <- colMeans(tau_boot2_s_NIPW)


#MSE_boot_NIPWbc <- 2 * MSE_boot_NIPW - MSE_boot2_NIPW

Bias_boot_NIPWbc <- 2 * colMeans(tau_boot_dif_NIPW) - colMeans(tau_boot2_dif_NIPW)
Bias_boot_NIPW <-  colMeans(tau_boot_dif_NIPW)

MSE_boot_NIPWbc = numeric(m) 
for (i in 1:m){
  if (MSE_boot_NIPW[i] >= MSE_boot2_NIPW[i] ) {
    MSE_boot_NIPWbc[i] = 2 * MSE_boot_NIPW[i] - MSE_boot2_NIPW[i]
  } else {
    MSE_boot_NIPWbc[i] = MSE_boot_NIPW[i]*exp((MSE_boot_NIPW[i]-MSE_boot2_NIPW[i])/MSE_boot2_NIPW[i])}
}


# CI with MSE
upInt_NIPW <- tau_hat_NIPW + sqrt(MSE_boot_NIPW) * 1.96
loInt_NIPW <- tau_hat_NIPW - sqrt(MSE_boot_NIPW) * 1.96
covB_NIPW <- (loInt_NIPW <= tau_true ) & (upInt_NIPW >= tau_true )
lenB_NIPW <- upInt_NIPW - loInt_NIPW

#bias corrected MSE
upInt_NIPWbc <- tau_hat_NIPW + sqrt(MSE_boot_NIPWbc) * 1.96
loInt_NIPWbc <- tau_hat_NIPW - sqrt(MSE_boot_NIPWbc) * 1.96
covB_NIPWbc <- (loInt_NIPWbc <= tau_true ) & (upInt_NIPWbc >= tau_true )
lenB_NIPWbc <- upInt_NIPWbc - loInt_NIPWbc

# bias corrected estimator 1
upInt_NIPWb <- tau_hat_NIPW + Bias_boot_NIPW + sqrt(MSE_boot_NIPW) * 1.96
loInt_NIPWb <- tau_hat_NIPW + Bias_boot_NIPW - sqrt(MSE_boot_NIPW) * 1.96
covB_NIPWb <- (loInt_NIPWb <= tau_true ) & (upInt_NIPWb >= tau_true )
lenB_NIPWb <- upInt_NIPWb - loInt_NIPWb

# bias corrected estimator 2
upInt_NIPWb2 <- tau_hat_NIPW + Bias_boot_NIPWbc + sqrt(MSE_boot_NIPW) * 1.96
loInt_NIPWb2 <- tau_hat_NIPW + Bias_boot_NIPWbc - sqrt(MSE_boot_NIPW) * 1.96
covB_NIPWb2 <- (loInt_NIPWb2 <= tau_true ) & (upInt_NIPWb2 >= tau_true )
lenB_NIPWb2 <- upInt_NIPWb2 - loInt_NIPWb2

# bias corrected estimator 1 + MSE
upInt_NIPWbbc <- tau_hat_NIPW + Bias_boot_NIPW + sqrt(MSE_boot_NIPWbc) * 1.96
loInt_NIPWbbc <- tau_hat_NIPW + Bias_boot_NIPW - sqrt(MSE_boot_NIPWbc) * 1.96
covB_NIPWbbc <- (loInt_NIPWbbc <= tau_true ) & (upInt_NIPWbbc >= tau_true )
lenB_NIPWbbc <- upInt_NIPWbbc - loInt_NIPWbbc


# bias corrected estimator 2 + MSE
upInt_NIPWb2bc <- tau_hat_NIPW + Bias_boot_NIPWbc + sqrt(MSE_boot_NIPWbc) * 1.96
loInt_NIPWb2bc <- tau_hat_NIPW + Bias_boot_NIPWbc - sqrt(MSE_boot_NIPWbc) * 1.96
covB_NIPWb2bc <- (loInt_NIPWbbc <= tau_true ) & (upInt_NIPWbbc >= tau_true )
lenB_NIPWb2bc <- upInt_NIPWb2bc - loInt_NIPWb2bc


# CI with percentiles
tau_boot_025_NIPW <- apply(tau_boot_NIPW, 2,  quantile, 0.025)
tau_boot_975_NIPW <- apply(tau_boot_NIPW, 2,  quantile, 0.975)
covBper_NIPW <- (tau_boot_025_NIPW <= tau_true ) & (tau_boot_975_NIPW >= tau_true )
lenBper_NIPW <- tau_boot_975_NIPW - tau_boot_025_NIPW

tau_boot_025_NIPWb <- apply(tau_boot_NIPW + Bias_boot_NIPW, 2,  quantile, 0.025)
tau_boot_975_NIPWb <- apply(tau_boot_NIPW + Bias_boot_NIPW, 2,  quantile, 0.975)
covBper_NIPWb <- (tau_boot_025_NIPWb <= tau_true ) & (tau_boot_975_NIPWb >= tau_true )
lenBper_NIPWb <- tau_boot_975_NIPWb - tau_boot_025_NIPWb

tau_boot_025_NIPWb2 <- apply(tau_boot_NIPW + Bias_boot_NIPWbc, 2,  quantile, 0.025)
tau_boot_975_NIPWb2 <- apply(tau_boot_NIPW + Bias_boot_NIPWbc, 2,  quantile, 0.975)
covBper_NIPWb2 <- (tau_boot_025_NIPWb2 <= tau_true ) & (tau_boot_975_NIPWb2 >= tau_true )
lenBper_NIPWb2 <- tau_boot_975_NIPWb2 - tau_boot_025_NIPWb2



Results = list(ni = ni,
                  n = sum(ni),
                  
                  Ni = Ni,
                  N = sum(Ni),
                  frac_nN = frac_nN,
                  
                  
               #Tau
               tau_true = tau_true,
               tau_hat_AIPW = tau_hat_AIPW, 
               tau_hat_NIPW = tau_hat_NIPW, 
               tau_hat_OR = tau_hat_OR, 
               
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
               
               # Boot NIPW
               #  MSE
               tau_boot_dif_NIPW  = tau_boot_dif_NIPW, 
               tau_boot_s_NIPW  = tau_boot_s_NIPW, 
               MSE_boot_NIPW = MSE_boot_NIPW,
               
               tau_boot2_dif_NIPW = tau_boot2_dif_NIPW, 
               tau_boot2_s_NIPW = tau_boot2_s_NIPW, 
               MSE_boot2_NIPW = MSE_boot2_NIPW, 
               
               Bias_boot_NIPWbc = Bias_boot_NIPWbc, 
               Bias_boot_NIPW = Bias_boot_NIPW, 
               
               # CI with MSE
               upInt_NIPW  = upInt_NIPW, 
               loInt_NIPW = loInt_NIPW, 
               covB_NIPW = covB_NIPW, 
               lenB_NIPW = lenB_NIPW,

               #bias corrected MSE
               upInt_NIPWbc = upInt_NIPWbc, 
               loInt_NIPWbc = loInt_NIPWbc, 
               covB_NIPWbc  = covB_NIPWbc,  
               lenB_NIPWbc = lenB_NIPWbc, 
               
               # bias corrected estimator 1
               upInt_NIPWb  = upInt_NIPWb, 
               loInt_NIPWb = loInt_NIPWb, 
               covB_NIPWb  = covB_NIPWb, 
               lenB_NIPWb  = lenB_NIPWb, 
               
               # bias corrected estimator 2
               upInt_NIPWb2 = upInt_NIPWb2, 
               loInt_NIPWb2 = loInt_NIPWb2, 
               covB_NIPWb2 = covB_NIPWb2, 
               lenB_NIPWb2 = lenB_NIPWb2,
               
               # bias corrected estimator 1 + MSE
               upInt_NIPWbbc = upInt_NIPWbbc, 
               loInt_NIPWbbc = loInt_NIPWbbc, 
               covB_NIPWbbc = covB_NIPWbbc,
               lenB_NIPWbbc = lenB_NIPWbbc, 
               
               # bias corrected estimator 2 + MSE
               upInt_NIPWb2bc = upInt_NIPWb2bc, 
               loInt_NIPWb2bc = loInt_NIPWb2bc, 
               covB_NIPWb2bc = covB_NIPWb2bc, 
               lenB_NIPWb2bc = lenB_NIPWb2bc, 
               
               # CI with percentiles
               tau_boot_025_NIPW = tau_boot_025_NIPW, 
               tau_boot_975_NIPW = tau_boot_975_NIPW, 
               covBper_NIPW = covBper_NIPW,
               lenBper_NIPW = lenBper_NIPW,
               
               tau_boot_025_NIPWb = tau_boot_025_NIPWb, 
               tau_boot_975_NIPWb  = tau_boot_975_NIPWb,  
               covBper_NIPWb = covBper_NIPWb, 
               lenBper_NIPWb = lenBper_NIPWb, 
               
               tau_boot_025_NIPWb2  = tau_boot_025_NIPWb2, 
               tau_boot_975_NIPWb2 = tau_boot_975_NIPWb2, 
               covBper_NIPWb2  = covBper_NIPWb2, 
               lenBper_NIPWb2  = lenBper_NIPWb2, 
               
               
               # Boot OR
               # MSE
               #  MSE
               tau_boot_dif_OR  = tau_boot_dif_OR, 
               tau_boot_s_OR  = tau_boot_s_OR, 
               MSE_boot_OR = MSE_boot_OR,
               
               tau_boot2_dif_OR = tau_boot2_dif_OR, 
               tau_boot2_s_OR = tau_boot2_s_OR, 
               MSE_boot2_OR = MSE_boot2_OR, 
               
               Bias_boot_ORbc = Bias_boot_ORbc, 
               Bias_boot_OR = Bias_boot_OR, 
               
               # CI with MSE
               upInt_OR  = upInt_OR, 
               loInt_OR = loInt_OR, 
               covB_OR = covB_OR, 
               lenB_OR = lenB_OR,
               
               #bias corrected MSE
               upInt_ORbc = upInt_ORbc, 
               loInt_ORbc = loInt_ORbc, 
               covB_ORbc  = covB_ORbc,  
               lenB_ORbc = lenB_ORbc, 
               
               # bias corrected estimator 1
               upInt_ORb  = upInt_ORb, 
               loInt_ORb = loInt_ORb, 
               covB_ORb  = covB_ORb, 
               lenB_ORb  = lenB_ORb, 
               
               # bias corrected estimator 2
               upInt_ORb2 = upInt_ORb2, 
               loInt_ORb2 = loInt_ORb2, 
               covB_ORb2 = covB_ORb2, 
               lenB_ORb2 = lenB_ORb2,
               
               # bias corrected estimator 1 + MSE
               upInt_ORbbc = upInt_ORbbc, 
               loInt_ORbbc = loInt_ORbbc, 
               covB_ORbbc = covB_ORbbc,
               lenB_ORbbc = lenB_ORbbc, 
               
               # bias corrected estimator 2 + MSE
               upInt_ORb2bc = upInt_ORb2bc, 
               loInt_ORb2bc = loInt_ORb2bc, 
               covB_ORb2bc = covB_ORb2bc, 
               lenB_ORb2bc = lenB_ORb2bc, 
               
               # Percentiles
               # CI with percentiles
               tau_boot_025_OR = tau_boot_025_OR, 
               tau_boot_975_OR = tau_boot_975_OR, 
               covBper_OR = covBper_OR,
               lenBper_OR = lenBper_OR, 
               
               
               tau_boot_025_ORb = tau_boot_025_ORb, 
               tau_boot_975_ORb  = tau_boot_975_ORb,  
               covBper_ORb = covBper_ORb, 
               lenBper_ORb = lenBper_ORb, 
               
               tau_boot_025_ORb2  = tau_boot_025_ORb2, 
               tau_boot_975_ORb2 = tau_boot_975_ORb2, 
               covBper_ORb2  = covBper_ORb2, 
               lenBper_ORb2  = lenBper_ORb2, 
               

               error_E = error_E)


outputName = paste("nb1d_nc_", SN, ".RData", sep = "")
outputPath = file.path("/user/work/ze23696/Comp", outputName)
save("Results", file = outputPath)