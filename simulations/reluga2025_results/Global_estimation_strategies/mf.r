library("MASS")
library("LongituRF")

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
############################################################################################
### Other estimators 
############################################################################################
############################################################################################
# Propensity score and m_0 and m_1                                                         #
############################################################################################
# E for P(A = 1|X)  #
#########################
ps_hat_E = c(data_sample_ps$ps_hat_E, data_out_of_sample_ps$ps_hat_E)
#which(ps_hat_E >= 0.98)
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



time_al <- NULL
time_al_oos <- NULL

for (i in 1:41) {
  time_al <- c(time_al, 1:ni[i])
  time_al_oos <- c(time_al_oos, 1:(1000 - ni[i]))
}



###################################
# MERF
#############################################
# Regular one model 
smerf <- MERF(X = data_sample[, c(1:12)],
              Y = data_sample$y,
              Z = matrix(1, nrow = nrow(data_sample), ncol = 1),
              id = data_sample$group,
              time = time_al,
              mtry=2, ntree=500, sto="BM")

y_full_impute_MRs <- predict(smerf, 
                             X = as.matrix(data_sample[, c(1:12)]),
                             Z = matrix(1, nrow = nrow(data_sample), ncol = 1),
                             id = data_sample$group,
                             time = time_al)

y_full_impute_MRoos <- predict(smerf, 
                               X = as.matrix(data_out_of_sample[, c(1:12)]),
                               Z = matrix(1, nrow = nrow(data_out_of_sample), ncol = 1),
                               id = data_out_of_sample$group,
                               time = time_al_oos)

# AIPW 1M

AIPW_dataMR <- data.frame(
  y = c(y_full_impute_MRs, y_full_impute_MRoos),
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_E,
  mu0_y = mu0_hat_M,
  mu1_y = mu1_hat_M
)

AIPWfMR <-
  as.data.frame(calculate_tau(AIPW_dataMR, type_tau = "AIPW"))

tau_hat_AIPW_MR = AIPWfMR$tau
#tau_hat_AIPW_MR[SN,] <- tau_hat_AIPWMR

# IPW 1M
IPW_dataMR <- data.frame(
  y =  c(y_full_impute_MRs, y_full_impute_MRoos),
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_E
)

IPWfMR <-
  as.data.frame(calculate_tau(IPW_dataMR, type_tau = "H"))

tau_hat_NIPW_MR = IPWfMR$tau


#tau_hat_NIPW_MR[SN,] <- tau_hat_NIPWMR

time_al0 <- NULL
time_al1 <- NULL

for (i in 1:41) {
  time_al0 <- c(time_al0, 1:nc[i])
  time_al1 <- c(time_al1, 1:nt[i])
}


# OR -- controls  --------------------------------

smerf0 <- MERF(X = data_sample0[, c(1:11)],
               Y = data_sample0$y,
               Z = matrix(1, nrow = nrow(data_sample0), ncol = 1),
               id = data_sample0$group,
               time = time_al0,
               mtry=2, ntree=500, sto="BM")


mu0_y_MRs <- predict(smerf0, 
                     X = as.matrix(data_sample[, c(1:11)]),
                     Z = matrix(1, nrow = nrow(data_sample), ncol = 1),
                     id = data_sample$group,
                     time = time_al)

mu0_y_MRoos <- predict(smerf0, 
                       X = as.matrix(data_out_of_sample[, c(1:11)]),
                       Z = matrix(1, nrow = nrow(data_out_of_sample), ncol = 1),
                       id = data_out_of_sample$group,
                       time = time_al_oos)

mu0_y_MR = c(mu0_y_MRs, mu0_y_MRoos)

# OR Treated --------------------


smerf1 <- MERF(X = data_sample1[, c(1:11)],
               Y = data_sample1$y,
               Z = matrix(1, nrow = nrow(data_sample1), ncol = 1),
               id = data_sample1$group,
               time = time_al1,
               mtry=2, ntree=500, sto="BM")

mu1_y_MRs <- predict(smerf1, 
                     X = as.matrix(data_sample[, c(1:11)]),
                     Z = matrix(1, nrow = nrow(data_sample), ncol = 1),
                     id = data_sample$group,
                     time = time_al)

mu1_y_MRoos <- predict(smerf1, 
                       X = as.matrix(data_out_of_sample[, c(1:11)]),
                       Z = matrix(1, nrow = nrow(data_out_of_sample), ncol = 1),
                       id = data_out_of_sample$group,
                       time = time_al_oos)

mu1_y_MR = c(mu1_y_MRs, mu1_y_MRoos)


# OR R 
df_MR <- data.frame(mu1_y = mu1_y_MR,
                    mu0_y = mu0_y_MR,
                    group =  c(data_sample$group, data_out_of_sample$group))

tau_treatMR = aggregate(df_MR$mu1_y, list(df_MR$group), FUN = mean)$x
tau_untreatMR = aggregate(df_MR$mu0_y, list(df_MR$group), FUN = mean)$x
tau_hat_OR_MR = tau_treatMR - tau_untreatMR


y_full_imputeMR2 <- numeric()
A0_full <- which(data_full_prelim$A == 0)
A1_full <- which(data_full_prelim$A == 1)
y_full_imputeMR2[A0_full] <- mu0_y_MR[A0_full]
y_full_imputeMR2[A1_full] <- mu1_y_MR[A1_full]


# AIPW 1M

AIPW_dataMR2 <- data.frame(
  y = y_full_imputeMR2,
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_E,
  mu0_y = mu0_hat_M,
  mu1_y = mu1_hat_M
)

AIPWfMR2 <-
  as.data.frame(calculate_tau(AIPW_dataMR2, type_tau = "AIPW"))

tau_hat_AIPW_MR2 = AIPWfMR2$tau

# IPW 1M

IPW_dataMR2 <- data.frame(
  y = y_full_imputeMR2,
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_E
)

IPWfMR2 <-
  as.data.frame(calculate_tau(IPW_dataMR2, type_tau = "H"))

tau_hat_NIPW_MR2 = IPWfMR2$tau



# M interaction
data_sampleEX <- data_sample

data_sampleEX$AX1 <- data_sample$A * data_sample$X1
data_sampleEX$AX2 <- data_sample$A * data_sample$X2

data_sampleEX$AX3 <- data_sample$A * data_sample$X3
data_sampleEX$AX4 <- data_sample$A * data_sample$X4

data_sampleEX$AX5 <- data_sample$A * data_sample$X5
data_sampleEX$AX6 <- data_sample$A * data_sample$X6

data_sampleEX$AX7 <- data_sample$A * data_sample$X7
data_sampleEX$AX8 <- data_sample$A * data_sample$X8

data_sampleEX$AX9 <- data_sample$A * data_sample$X9
data_sampleEX$AX10 <- data_sample$A * data_sample$X10

data_sampleEX$AX11 <- data_sample$A * data_sample$X11


data_out_of_sampleEX <- data_out_of_sample

data_out_of_sampleEX$AX1 <- data_out_of_sample$A * data_out_of_sample$X1
data_out_of_sampleEX$AX2 <- data_out_of_sample$A * data_out_of_sample$X2

data_out_of_sampleEX$AX3 <- data_out_of_sample$A * data_out_of_sample$X3
data_out_of_sampleEX$AX4 <- data_out_of_sample$A * data_out_of_sample$X4

data_out_of_sampleEX$AX5 <- data_out_of_sample$A * data_out_of_sample$X5
data_out_of_sampleEX$AX6 <- data_out_of_sample$A * data_out_of_sample$X6

data_out_of_sampleEX$AX7 <- data_out_of_sample$A * data_out_of_sample$X7
data_out_of_sampleEX$AX8 <- data_out_of_sample$A * data_out_of_sample$X8

data_out_of_sampleEX$AX9 <- data_out_of_sample$A * data_out_of_sample$X9
data_out_of_sampleEX$AX10 <- data_out_of_sample$A * data_out_of_sample$X10

data_out_of_sampleEX$AX11 <- data_out_of_sample$A * data_out_of_sample$X11



smerfi <- MERF(X = data_sampleEX[, c(1:12, 16:26)],
               Y = data_sampleEX$y,
               Z = matrix(1, nrow = nrow(data_sample), ncol = 1),
               id = data_sampleEX$group,
               time = time_al,
               mtry=2, ntree=500, sto="BM")

y_full_impute_MRis <- predict(smerfi, 
                              X = as.matrix(data_sampleEX[, c(1:12, 16:26)]),
                              Z = matrix(1, nrow = nrow(data_sample), ncol = 1),
                              id = data_sample$group,
                              time = time_al)

y_full_impute_MRioos <- predict(smerfi, 
                                X = as.matrix(data_out_of_sampleEX[,  c(1:12, 16:26)]),
                                Z = matrix(1, nrow = nrow(data_out_of_sample), ncol = 1),
                                id = data_out_of_sample$group,
                                time = time_al_oos)


# AIPW 1M Ri

AIPW_dataMRi <- data.frame(
  y = c(y_full_impute_MRis, y_full_impute_MRioos),
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_E,
  mu0_y = mu0_hat_M,
  mu1_y = mu1_hat_M
)

AIPWfMRi <-
  as.data.frame(calculate_tau(AIPW_dataMRi, type_tau = "AIPW"))

tau_hat_AIPW_MRi = AIPWfMRi$tau

# IPW 1M

IPW_dataMRi <- data.frame(
  y =  c(y_full_impute_MRis, y_full_impute_MRioos),
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_E
)

IPWfMRi <-
  as.data.frame(calculate_tau(IPW_dataMRi, type_tau = "H"))

tau_hat_NIPW_MRi = IPWfMRi$tau


# Regular one model 
smerf_csda <- MERF(X = data_sample[, c(1:12)],
              Y = data_sample$y,
              Z = cbind(1, data_sample$A),
              id = data_sample$group,
              time = time_al,
              mtry=2, ntree=500, sto="BM")

y_full_impute_MR_csda_s <- predict(smerf_csda, 
                             X = as.matrix(data_sample[, c(1:12)]),
                             Z = cbind(1, data_sample$A),
                             id = data_sample$group,
                             time = time_al)

y_full_impute_MR_csda_oos <- predict(smerf_csda, 
                               X = as.matrix(data_out_of_sample[, c(1:12)]),
                               Z = cbind(1, data_out_of_sample$A),
                               id = data_out_of_sample$group,
                               time = time_al_oos)

# AIPW 1M

AIPW_dataMR_csda <- data.frame(
  y = c(y_full_impute_MR_csda_s, y_full_impute_MR_csda_oos),
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_E,
  mu0_y = mu0_hat_M,
  mu1_y = mu1_hat_M
)

AIPWfMR_csda <-
  as.data.frame(calculate_tau(AIPW_dataMR_csda, type_tau = "AIPW"))

tau_hat_AIPW_MR_csda = AIPWfMR_csda$tau
#tau_hat_AIPW_MR[SN,] <- tau_hat_AIPWMR

# IPW 1M
IPW_dataMR_csda <- data.frame(
  y =  c(y_full_impute_MRs, y_full_impute_MRoos),
  A = c(data_sample$A, data_out_of_sample$A),
  group = c(data_sample$group, data_out_of_sample$group),
  p_score  = ps_hat_E
)

IPWfMR_csda <-
  as.data.frame(calculate_tau(IPW_dataMR_csda, type_tau = "H"))

tau_hat_NIPW_MR_csda = IPWfMR_csda$tau


Results = list(ni = ni,
               n = sum(ni),
               
               Ni = Ni,
               N = sum(Ni),
               frac_nN = frac_nN,
               
               
               #Tau
               tau_true = tau_true,
               tau_hat_AIPW_MR = tau_hat_AIPW_MR, 
               tau_hat_NIPW_MR = tau_hat_NIPW_MR,
               
               tau_hat_OR_MR = tau_hat_OR_MR, 
               
               tau_hat_AIPW_MR2 = tau_hat_AIPW_MR2,
               tau_hat_NIPW_MR2 = tau_hat_NIPW_MR2, 
               
               tau_hat_AIPW_MRi = tau_hat_AIPW_MRi, 
               tau_hat_NIPW_MRi = tau_hat_NIPW_MRi, 
               
               tau_hat_NIPW_MR_csda = tau_hat_NIPW_MR_csda, 
               tau_hat_AIPW_MR_csda = tau_hat_AIPW_MR_csda)


outputName = paste("mf_", SN, ".RData", sep = "")
outputPath = file.path("/user/work/ze23696/Comp", outputName)
save("Results", file = outputPath)
