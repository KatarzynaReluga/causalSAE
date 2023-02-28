#setwd("C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/Mplustune25_2")
#setwd("C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/PositiveEffects/Tuned/Mplustune100_2")
setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/PositiveEffects/Mplus100_2s")
library(data.table)
#setwd("C:/Users/katar/Documents/GitHub/causalSAE")


####################
### Load results ###
####################
SimNum = 1000
res <- list()

file_list <- list.files()
for (i in 1:SimNum) {
  load(file_list[i])

#  resM <- data.frame(do.call(rbind, Results))[-c(1:3), ]
  resM <- data.frame(do.call(rbind, Results))
  colnames(resM) <- paste("C", 1:dim(resM)[2], sep = "")
  nMethotds <- dim(resM)[1]
  nClusters <- dim(resM)[2]
  tau_trueM <- matrix(rep(Results$tau_true, nMethotds),
                      ncol = nClusters,
                      nrow = nMethotds, byrow  = T)
  Method <- rownames(resM)
  rownames(resM) <- NULL

  resMC <- resM - tau_trueM
  colnames(resMC) <- paste("CC", 1:dim(resMC)[2], sep = "")

  resMCS <- (resM - tau_trueM) ^ 2
  colnames(resMCS) <- paste("CCS", 1:dim(resMCS)[2], sep = "")


  res[[i]] <- cbind(resM, resMC, resMCS, tau_trueM)
  res[[i]]$Method <- Method

#  colnames(all_results[[i]]) <- paste("C", 1:dim(all_results[[i]])[2], sep = "")

}

 # Join results and compute
dt_res <- data.table(do.call(rbind, res))
dt_Mean <- dt_res[, lapply(.SD, mean), by = "Method"]

# 25
tau_trueM <- dt_Mean[c(4:72), c(77:101)]
# 50
tau_trueM <- dt_Mean[c(4:72), c(152:201)]
# 100
tau_trueM <- dt_Mean[c(4:72), c(302:401)]

# design based
tau_trueM <- dt_Mean[c(4:72), c(125:165)]

Methods <- dt_res$Method[4:72]


# 25
tau_trueM <- dt_Mean[c(4:9), c(77:101)]
# 50
tau_trueM <- dt_Mean[c(4:9), c(152:201)]
# 100
tau_trueM <- dt_Mean[c(4:9), c(302:401)]

# design based
tau_trueM <- dt_Mean[c(4:9), c(125:165)]

Methods <- dt_res$Method[4:9]

# Bias ----------------------------------------------
RB <- dt_Mean[c(4:72), c(27:51)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(27:51)]/abs(tau_trueM) * 100

RB <- dt_Mean[c(4:72), c(52:101)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(52:101)]/abs(tau_trueM) * 100

RB <- dt_Mean[c(4:72), c(102:201)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(102:201)]/abs(tau_trueM) * 100

# Design based
RB <- dt_Mean[c(4:9), c(43:83)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:9), c(43:83)]/abs(tau_trueM) * 100



RB <- dt_Mean[c(4:9), c(27:51)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:9), c(27:51)]/abs(tau_trueM) * 100

RB <- dt_Mean[c(4:9), c(52:101)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:9), c(52:101)]/abs(tau_trueM) * 100

RB <- dt_Mean[c(4:9), c(102:201)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:9), c(102:201)]/abs(tau_trueM) * 100

# Design based
RB <- dt_Mean[c(4:9), c(43:83)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:9), c(43:83)]/abs(tau_trueM) * 100



# Relative Mean -------------------------------------------
RBPMean <- rowMeans(RBP)
# Relative bias ------------------------------------------
RBPMedian <- apply(RBP, 1, quantile, probs = 0.5)

# Absolute Bias ----------------------------------------------
ARB <- abs(dt_Mean[c(4:72), c(27:51)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:72), c(27:51)]/abs(tau_trueM) * 100)

ARB <- abs(dt_Mean[c(4:72), c(52:101)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:72), c(52:101)]/abs(tau_trueM) * 100)

ARB <- abs(dt_Mean[c(4:72), c(102:201)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:72), c(102:201)]/abs(tau_trueM) * 100)

# Design bias
ARB <- abs(dt_Mean[c(4:72), c(43:83)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:72),c(43:83)]/abs(tau_trueM) * 100)


# Absolute Bias ----------------------------------------------
ARB <- abs(dt_Mean[c(4:9), c(27:51)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:9), c(27:51)]/abs(tau_trueM) * 100)

ARB <- abs(dt_Mean[c(4:9), c(52:101)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:9), c(52:101)]/abs(tau_trueM) * 100)

ARB <- abs(dt_Mean[c(4:9), c(102:201)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:9), c(102:201)]/abs(tau_trueM) * 100)

# Design bias
ARB <- abs(dt_Mean[c(4:9), c(43:83)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:9),c(43:83)]/abs(tau_trueM) * 100)

# Means --------------------------------------------------
ARBPMean <- rowMeans(ARBP)
# Medians ------------------------------------------------
ARBPMedian <- apply(ARBP, 1, quantile, probs = 0.5)

# RRMSE ----------------------------------------------
RMSE <- sqrt(dt_Mean[c(4:72), c(52:76)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:72), c(52:76)])/abs(tau_trueM) * 100

RMSE <- sqrt(dt_Mean[c(4:72), c(102:151)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:72), c(102:151)])/abs(tau_trueM) * 100

RMSE <- sqrt(dt_Mean[c(4:72), c(202:301)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:72), c(202:301)])/abs(tau_trueM) * 100

# Design based --------------------------------------------
RMSE <- sqrt(dt_Mean[c(4:72), c(84:124)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:72), c(84:124)])/abs(tau_trueM) * 100


# RRMSE ----------------------------------------------
RMSE <- sqrt(dt_Mean[c(4:9), c(52:76)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:9), c(52:76)])/abs(tau_trueM) * 100

RMSE <- sqrt(dt_Mean[c(4:9), c(102:151)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:9), c(102:151)])/abs(tau_trueM) * 100

RMSE <- sqrt(dt_Mean[c(4:9), c(202:301)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:9), c(202:301)])/abs(tau_trueM) * 100

# Design based --------------------------------------------
RMSE <- sqrt(dt_Mean[c(4:9), c(84:124)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:9), c(84:124)])/abs(tau_trueM) * 100


RMSEMean <- rowMeans(RMSEP)
RMSEPMedian <- apply(RMSEP, 1, quantile, probs = 0.5)


resultsMPlus100tune <- data.frame(Method = Methods,

           RBPMean = RBPMean,
           RBPMedian = RBPMedian,

           ARBPMean = ARBPMean,
           ARBPMedian = ARBPMedian,

           RMSEMean = RMSEMean,
           RMSEPMedian = RMSEPMedian)

write.csv(resultsMPlus100tune,
          file = 'resultsMPlus100tune.csv', row.names = FALSE)

# Not tuned

Methods <- resultsMPlus50$Method

MPlus25 <- read.csv('C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/PositiveEffects/resultsMPlus25.csv')
ARBPMedian25 <- MPlus25$ARBPMedian
ARBPMediansort25 <- sort(ARBPMedian25, index.return = TRUE)$x
ARBPMediansorti25 <- sort(ARBPMedian25, index.return = TRUE)$ix
Methods[ARBPMediansorti25]


RMSEPMedian25 <- MPlus25$RMSEPMedian
RMSEPMediansort25 <- sort(RMSEPMedian25, index.return = TRUE)$x
RMSEPMediansorti25 <- sort(RMSEPMedian25, index.return = TRUE)$ix
Methods[RMSEPMediansorti25]

MPlus50 <- read.csv('C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/PositiveEffects/resultsMPlus50.csv')
ARBPMedian50 <- MPlus50$ARBPMedian
ARBPMediansort50 <- sort(ARBPMedian50, index.return = TRUE)$x
ARBPMediansorti50 <- sort(ARBPMedian50, index.return = TRUE)$ix
Methods[ARBPMediansorti50]


RMSEPMedian50 <- MPlus50$RMSEPMedian
RMSEPMediansort50 <- sort(RMSEPMedian50, index.return = TRUE)$x
RMSEPMediansorti50 <- sort(RMSEPMedian50, index.return = TRUE)$ix
Methods[RMSEPMediansorti50]
Methods <- MPlus50$Method


MPlus100 <- read.csv('C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/PositiveEffects/resultsMPlus100.csv')
ARBPMedian100 <- MPlus100$ARBPMedian
ARBPMediansort100 <- sort(ARBPMedian100, index.return = TRUE)$x
ARBPMediansorti100 <- sort(ARBPMedian100, index.return = TRUE)$ix
Methods[ARBPMediansorti100]


RMSEPMedian100 <- MPlus100$RMSEPMedian
RMSEPMediansort100 <- sort(RMSEPMedian100, index.return = TRUE)$x
RMSEPMediansorti100 <- sort(RMSEPMedian100, index.return = TRUE)$ix
Methods[RMSEPMediansorti100]

# Tuned
MPlus25t <- read.csv('C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/PositiveEffects/resultsMPlus25tune.csv')
ARBPMedian25t <- MPlus25t$ARBPMedian
ARBPMediansort25t <- sort(ARBPMedian25t, index.return = TRUE)$x
ARBPMediansorti25t <- sort(ARBPMedian25t, index.return = TRUE)$ix
Methods[ARBPMedian25t]


RMSEPMedian25t <- MPlus25t$RMSEPMedian
RMSEPMediansort25t <- sort(RMSEPMedian25t, index.return = TRUE)$x
RMSEPMediansorti25t <- sort(RMSEPMedian25t, index.return = TRUE)$ix
Methods[RMSEPMediansorti25t]

MPlus50t <- read.csv('C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/PositiveEffects/resultsMPlus50tune.csv')
Methods <- MPlus50$Method

ARBPMedian50t <- MPlus50t$ARBPMedian
ARBPMediansort50t <- sort(ARBPMedian50t, index.return = TRUE)$x
ARBPMediansorti50t <- sort(ARBPMedian50t, index.return = TRUE)$ix
Methods[ARBPMediansorti50t]


RMSEPMedian50t <- MPlus50$RMSEPMedian
RMSEPMediansort50t <- sort(RMSEPMedian50, index.return = TRUE)$x
RMSEPMediansorti50t <- sort(RMSEPMedian50, index.return = TRUE)$ix
Methods[RMSEPMediansorti50t]
Methods <- MPlus50t$Method


MPlus100t <- read.csv('C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/PositiveEffects/resultsMPlus100tune.csv')
ARBPMedian100t <- MPlus100t$ARBPMedian
ARBPMediansort100t <- sort(ARBPMedian100t, index.return = TRUE)$x
ARBPMediansorti100t <- sort(ARBPMedian100t, index.return = TRUE)$ix
Methods[ARBPMediansorti100t]


RMSEPMedian100t <- MPlus100t$RMSEPMedian
RMSEPMediansort100t <- sort(RMSEPMedian100t, index.return = TRUE)$x
RMSEPMediansorti100t <- sort(RMSEPMedian100t, index.return = TRUE)$ix
Methods[RMSEPMediansorti100t]

# All together


MPlus25F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/PositiveEffects/resultsMPlus25F.csv')
Methods <- MPlus25F$Method
ARBPMedian25F <- MPlus25F$ARBPMedian
ARBPMediansort25F <- sort(ARBPMedian25F, index.return = TRUE)$x
ARBPMediansorti25F <- sort(ARBPMedian25F, index.return = TRUE)$ix
Methods[ARBPMediansorti25F][1:10]
plot(density(ARBPMediansort25F))
plot(1:127, ARBPMedian25F)

RMSEPMedian25F <- MPlus25F$RMSEPMedian
RMSEPMediansort25F <- sort(RMSEPMedian25F, index.return = TRUE)$x
RMSEPMediansorti25F <- sort(RMSEPMedian25F, index.return = TRUE)$ix
Methods[RMSEPMediansorti25F][1:10]
plot(density(RMSEPMediansort25F))
plot(1:127, RMSEPMediansort25F)


MPlus50F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/PositiveEffects/resultsMPlus50F.csv')
ARBPMedian50F <- MPlus50F$ARBPMedian
ARBPMediansort50F <- sort(ARBPMedian50F, index.return = TRUE)$x
ARBPMediansorti50F <- sort(ARBPMedian50F, index.return = TRUE)$ix
Methods[ARBPMediansorti50F]
Methods[ARBPMediansorti50F][1:10]
plot(density(ARBPMediansort50F))
plot(1:127, ARBPMedian50F)


RMSEPMedian50F <- MPlus50F$RMSEPMedian
RMSEPMediansort50F <- sort(RMSEPMedian50F, index.return = TRUE)$x
RMSEPMediansorti50F <- sort(RMSEPMedian50F, index.return = TRUE)$ix
Methods[RMSEPMediansorti50F]
Methods[ARBPMediansorti50F][1:10]
plot(density(ARBPMediansort50F))
plot(1:127, ARBPMedian50F)


MPlus100 <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/PositiveEffects/resultsMPlus100.csv')
ARBPMedian100 <- MPlus100$ARBPMedian
ARBPMediansort100 <- sort(ARBPMedian100, index.return = TRUE)$x
ARBPMediansorti100 <- sort(ARBPMedian100, index.return = TRUE)$ix
Methods[ARBPMediansorti100]


RMSEPMedian100 <- MPlus100$RMSEPMedian
RMSEPMediansort100 <- sort(RMSEPMedian100, index.return = TRUE)$x
RMSEPMediansorti100 <- sort(RMSEPMedian100, index.return = TRUE)$ix
Methods[RMSEPMediansorti100]


# Means processing -----------------------------------------
# Very skewed results, some areas not estimated precisely, we look at medians instead
# ARBPMean <- rowMeans(ARBP)
#############################################################
ARBPMeansorti <- sort(ARBPMean, index.return = T)$ix
ARBPMeansort <- sort(ARBPMean)
ARBPMeanMin <- which.min(ARBPMean)
Methods[ARBPMeansorti]
plot(1:length(Methods), ARBPMean, type = "l")
##############################################################
# Medians processing ----------------------------------
# Safe medians -----------------------------
#ARBPMedian <- apply(ARBP, 1, quantile, probs = 0.5)
##############################################################
ARBPMediansorti <- sort(ARBPMedian, index.return = T)$ix
ARBPMediansort <- sort(ARBPMedian)
ARBPMedianMin <- which.min(ARBPMedian)
Methods[ARBPMediansorti]
##############################################################

RBPQ <- apply(RBP, 1, quantile, probs = c(0.25, 0.75))
ARBPQ <- abs(RBPQ)
ARBPQMin <- apply(ARBPQ, 1, which.min)
plot(1:ncol(RBP), ARBP[ARBPMeanMin, ], type = "l")
lines(1:ncol(RBP), ARBP[ARBPMedianMin, ], type = "l", col = 2)


plot(density(unlist(ARBP[ARBPMeanMin, ])))
lines(density(unlist(ARBP[ARBPMedianMin, ])), col = 2)
sort(unlist(ARBP[ARBPMedianMin, ]), index.return = T)
sort(unlist(ARBP[ARBPMeanMin, ]), index.return = T)
Q9 <- apply(ARBP, 1, quantile, probs = 0.9)
Q84 <- apply(ARBP, 1, quantile, probs = 0.84)
plot(density(unlist(ARBP[1,])))
sort(unname(unlist(ARBP[1,])))

x <- unname(unlist(ARBP[ARBPMeanMin, ]))
funsort <- function(x) {
  sortx <- sort(x, index.return = T)
  lengthx <- length(x)
  indexs <- sortx$ix[c(0.9 * lengthx + 1):lengthx]
  indexs <- sortx$ix[c(0.84 * lengthx + 1):lengthx]

}
plot(density(unlist(ARBP[60, ])))
sort(unname(unlist(ARBP[60, ])))

indapp <- apply(ARBP, 1, funsort)
table(indapp[1, ])
table(indapp[2, ])
table(indapp[3, ])
table(indapp[4, ])
table(indapp[5, ])
table(indapp)
# 25
# 7 14 18 19


# 50
# 6, 10, 11, 12, 25, 30, 34
# 10, 11, 12, 25, 30
#10, 11, 12, 25, 26, 30, 44, 49

# 100
# 30 46 48 64 65 71 86 93
tau_true[c(6, 10, 11, 12, 25, 30, 34)]
sort(tau_true[c(10, 11, 12, 25, 26, 30, 44, 49)])
Methods[c(ARBPMeanMin, ARBPMedianMin)]

#pp <- sqrt(colMeans((E_OR - tau_trueM)^2))/abs(tau_true) * 100

RMSEMeanMin <- which.min(RMSEMean)

RMSEPMedianMin <- which.min(RMSEPMedian)
RMSEsort <- sort(RMSEPMedian, index.return = TRUE)$ix
Methods[RMSEsort]
Methods[RMSEPMedianMin]
Methods[RMSEMeanMin]

indapp <- apply(RMSE, 1, funsort)
table(indapp[1, ])
table(indapp[2, ])
table(indapp[3, ])
table(indapp[4, ])
table(indapp[5, ])
table(indapp)
# 7 14 18 19

RMSEPQ <- apply(RMSEP, 1, quantile, probs = c(0.25, 0.75))
RMSEPQMin <- apply(RMSEPQ, 1, which.min)

plot(1:ncol(RMSE), RMSEP[RMSEMeanMin,], type = "l")
lines(1:ncol(RMSE), RMSEP[RMSEPMedianMin,], type = "l", col = 2)

plot(density(as.numeric(RMSEP[3,])))



ARBPMean <- abs(RBPMean)
ARBPMeanMin <- which.min(ARBPMean)

RBPMedian <- apply(RBP, 1, quantile, probs = 0.5)
ARBPMedian <- abs(RBPMedian)
ARBPMedianMin <- which.min(ARBPMedian)

Methods[c(ARBPMeanMin,ARBPMedianMin)]

#########################################################################
#########################################################################
###########################
### Comparison criteria ###
###########################
#Relative Bias percent
RBias_percent <- function(estimator, true_value) {

  value <- (apply(estimator, 2, mean) - true_value) / abs(true_value) * 100
  value

}

#Relative Bias
RBias <- function(estimator, true_value) {

  value <- (apply(estimator, 2, mean) - true_value) / abs(true_value)
  value

}

# Relative RMSE
RRMSE_percent <- function(estimator, true_value) {

  value <- (sqrt(apply(apply(estimator, 1, function(x) {
    (x - true_value) ^ 2
  }), 1, mean, na.rm = TRUE
  )) / abs(true_value)) * 100
  value
}

#RRMSE
RRMSE <- function(estimator, true_value) {

  value <- (sqrt(apply(apply(estimator, 1, function(x) {
    (x - true_value) ^ 2
  }), 1, mean, na.rm = TRUE
  )) / abs(true_value))
  value
}

#RRMSE percent
RRMSE_percent <- function(estimator, true_value) {

  value <- (sqrt(apply(apply(estimator, 1, function(x) {
    (x - true_value) ^ 2
  }), 1, mean, na.rm = TRUE
  )) / abs(true_value)) * 100
  value
}

#MSE
MSE <- function(estimator, true_value) {

  value <- apply(apply(estimator, 1, function(x) {
    (x - true_value) ^ 2
  }), 1, mean, na.rm = TRUE
  )
  value
}




D = 100

tau_trueM <- matrix(0, nrow = SimNum, ncol = D)

# OR ---------------------------
E_OR <- matrix(0, nrow = SimNum, ncol = D)
M_OR <- matrix(0, nrow = SimNum, ncol = D)
R_OR <- matrix(0, nrow = SimNum, ncol = D)
X_OR <- matrix(0, nrow = SimNum, ncol = D)

# NIPW ----------------------
EE_NIPW <- matrix(0, nrow = SimNum, ncol = D)
EM_NIPW <- matrix(0, nrow = SimNum, ncol = D)
ER_NIPW <- matrix(0, nrow = SimNum, ncol = D)
EX_NIPW <- matrix(0, nrow = SimNum, ncol = D)

MM_NIPW <- matrix(0, nrow = SimNum, ncol = D)
ME_NIPW <- matrix(0, nrow = SimNum, ncol = D)
MR_NIPW <- matrix(0, nrow = SimNum, ncol = D)
MX_NIPW <- matrix(0, nrow = SimNum, ncol = D)

RR_NIPW <- matrix(0, nrow = SimNum, ncol = D)
RE_NIPW <- matrix(0, nrow = SimNum, ncol = D)
RM_NIPW <- matrix(0, nrow = SimNum, ncol = D)
RX_NIPW <- matrix(0, nrow = SimNum, ncol = D)

XX_NIPW <- matrix(0, nrow = SimNum, ncol = D)
XE_NIPW <- matrix(0, nrow = SimNum, ncol = D)
XM_NIPW <- matrix(0, nrow = SimNum, ncol = D)
XR_NIPW <- matrix(0, nrow = SimNum, ncol = D)

# AIPW
# EBLUP -------
EEM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
EER_AIPW <- matrix(0, nrow = SimNum, ncol = D)
EEX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

EMM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
EMR_AIPW <- matrix(0, nrow = SimNum, ncol = D)
EMX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

ERM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
ERR_AIPW <- matrix(0, nrow = SimNum, ncol = D)
ERX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

EXM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
EXR_AIPW <- matrix(0, nrow = SimNum, ncol = D)
EXX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

# MQ ----------------
MEE_AIPW <- matrix(0, nrow = SimNum, ncol = D)
MER_AIPW <- matrix(0, nrow = SimNum, ncol = D)
MEX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

MME_AIPW <- matrix(0, nrow = SimNum, ncol = D)
MMR_AIPW <- matrix(0, nrow = SimNum, ncol = D)
MMX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

MRE_AIPW <- matrix(0, nrow = SimNum, ncol = D)
MRR_AIPW <- matrix(0, nrow = SimNum, ncol = D)
MRX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

MXE_AIPW <- matrix(0, nrow = SimNum, ncol = D)
MXR_AIPW <- matrix(0, nrow = SimNum, ncol = D)
MXX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

# R --------------
REE_AIPW <- matrix(0, nrow = SimNum, ncol = D)
REM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
REX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

RME_AIPW <- matrix(0, nrow = SimNum, ncol = D)
RMM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
RMX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

RRE_AIPW <- matrix(0, nrow = SimNum, ncol = D)
RRM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
RRX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

RXE_AIPW <- matrix(0, nrow = SimNum, ncol = D)
RXM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
RXX_AIPW <- matrix(0, nrow = SimNum, ncol = D)

# X --------------
XEE_AIPW <- matrix(0, nrow = SimNum, ncol = D)
XEM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
XER_AIPW <- matrix(0, nrow = SimNum, ncol = D)

XME_AIPW <- matrix(0, nrow = SimNum, ncol = D)
XMM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
XMR_AIPW <- matrix(0, nrow = SimNum, ncol = D)

XRE_AIPW <- matrix(0, nrow = SimNum, ncol = D)
XRM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
XRR_AIPW <- matrix(0, nrow = SimNum, ncol = D)

XXE_AIPW <- matrix(0, nrow = SimNum, ncol = D)
XXM_AIPW <- matrix(0, nrow = SimNum, ncol = D)
XXR_AIPW <- matrix(0, nrow = SimNum, ncol = D)

tauDirect <- matrix(0, nrow = SimNum, ncol = D)

file_list <- list.files()
for (i in 1:SimNum) {

  load(file_list[i])

  tau_trueM[i, ] <- Results$tau_true

  # OR ---------------------------
  E_OR[i, ] <- Results$E_OR
  M_OR[i, ] <- Results$M_OR
  R_OR[i, ] <- Results$R_OR
  X_OR[i, ] <- Results$X_OR

  # NIPW ----------------------
  EE_NIPW[i, ] <- Results$EE_NIPW
  EM_NIPW[i, ] <- Results$EM_NIPW
  ER_NIPW[i, ] <- Results$ER_NIPW
  EX_NIPW[i, ] <- Results$EX_NIPW

  MM_NIPW[i, ] <- Results$MM_NIPW
  ME_NIPW[i, ] <- Results$ME_NIPW
  MR_NIPW[i, ] <- Results$MR_NIPW
  MX_NIPW[i, ] <- Results$MX_NIPW

  RR_NIPW[i, ] <- Results$RR_NIPW
  RE_NIPW[i, ] <- Results$RE_NIPW
  RM_NIPW[i, ] <- Results$RM_NIPW
  RX_NIPW[i, ] <- Results$RX_NIPW

  XX_NIPW[i, ] <- Results$XX_NIPW
  XE_NIPW[i, ] <- Results$XE_NIPW
  XM_NIPW[i, ] <- Results$XM_NIPW
  XR_NIPW[i, ] <- Results$XR_NIPW


  # EBLUP -------
  EEM_AIPW[i, ] <- Results$EEM_AIPW
  EER_AIPW[i, ] <- Results$EER_AIPW
  EEX_AIPW[i, ] <- Results$EEX_AIPW

  EMM_AIPW[i, ] <- Results$EMM_AIPW
  EMR_AIPW[i, ] <- Results$EMR_AIPW
  EMX_AIPW[i, ] <- Results$EMX_AIPW

  ERM_AIPW[i, ] <- Results$ERM_AIPW
  ERR_AIPW[i, ] <- Results$ERR_AIPW
  ERX_AIPW[i, ] <- Results$ERX_AIPW

  EXM_AIPW[i, ] <- Results$EXM_AIPW
  EXR_AIPW[i, ] <- Results$EXR_AIPW
  EXX_AIPW[i, ] <- Results$EXX_AIPW

  # MQ ----------------
  MEE_AIPW[i, ] <- Results$MEE_AIPW
  MER_AIPW[i, ] <- Results$MER_AIPW
  MEX_AIPW[i, ] <- Results$MEX_AIPW

  MME_AIPW[i, ] <- Results$MME_AIPW
  MMR_AIPW[i, ] <- Results$MMR_AIPW
  MMX_AIPW[i, ] <- Results$MMX_AIPW

  MRE_AIPW[i, ] <- Results$MRE_AIPW
  MRR_AIPW[i, ] <- Results$MRR_AIPW
  MRX_AIPW[i, ] <- Results$MRX_AIPW

  MXE_AIPW[i, ] <- Results$MXE_AIPW
  MXR_AIPW[i, ] <- Results$MXR_AIPW
  MXX_AIPW[i, ] <- Results$MXX_AIPW

  # R --------------
  REE_AIPW[i, ] <- Results$REE_AIPW
  REM_AIPW[i, ] <- Results$REM_AIPW
  REX_AIPW[i, ] <- Results$REX_AIPW

  RME_AIPW[i, ] <- Results$RME_AIPW
  RMM_AIPW[i, ] <- Results$RMM_AIPW
  RMX_AIPW[i, ] <- Results$RMX_AIPW

  RRE_AIPW[i, ] <- Results$RRE_AIPW
  RRM_AIPW[i, ] <- Results$RRM_AIPW
  RRX_AIPW[i, ] <- Results$RRX_AIPW

  RXE_AIPW[i, ] <- Results$RXE_AIPW
  RXM_AIPW[i, ] <- Results$RXM_AIPW
  RXX_AIPW[i, ] <- Results$RXX_AIPW

  # X --------------
  XEE_AIPW[i, ] <- Results$XEE_AIPW
  XEM_AIPW[i, ] <- Results$XEM_AIPW
  XER_AIPW[i, ] <- Results$XER_AIPW

  XME_AIPW[i, ] <- Results$XME_AIPW
  XMM_AIPW[i, ] <- Results$XMM_AIPW
  XMR_AIPW[i, ] <- Results$XMR_AIPW

  XRE_AIPW[i, ] <- Results$XRE_AIPW
  XRM_AIPW[i, ] <- Results$XRM_AIPW
  XRR_AIPW[i, ] <- Results$XRR_AIPW

  XXE_AIPW[i, ] <- Results$XXE_AIPW
  XXM_AIPW[i, ] <- Results$XXM_AIPW
  XXR_AIPW[i, ] <- Results$XXR_AIPW

  tauDirect[i, ] <- Results$Dir_tau

}

#########################
# Relative Bias percent #
#########################

tau_true <- tau_trueM[1, ]

RBperDirect <- RBias_percent(tauDirect, tau_true)

# OR ----------------------------------------
RBperE <- RBias_percent(E_OR, tau_true)
RBperM <- RBias_percent(M_OR, tau_true)
RBperR <- RBias_percent(R_OR, tau_true)
RBperX <- RBias_percent(X_OR, tau_true)

mean(RBperE)
mean(RBperM)
mean(RBperR)
mean(RBperX)
mean(RBperDirect)

quantile(RBperE, 0.5, type = 1)
quantile(RBperM, 0.5, type = 1)
quantile(RBperR, 0.5, type = 1)
quantile(RBperX, 0.5, type = 1)
quantile(RBperDirect, 0.5, type = 1)

#
boxplot(
  cbind(e = RBperE,
        m = RBperM,
        r = RBperR,
        x = RBperX,

        d = RBperDirect),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# IPW ---------------------------------------
# NIPW ----------------------
RBperEE <- RBias_percent(EE_NIPW, tau_true)
RBperEM <- RBias_percent(EM_NIPW, tau_true)
RBperER <- RBias_percent(ER_NIPW, tau_true)
RBperEX <- RBias_percent(EX_NIPW, tau_true)

RBperMM <- RBias_percent(MM_NIPW, tau_true)
RBperME <- RBias_percent(ME_NIPW, tau_true)
RBperMR <- RBias_percent(MR_NIPW, tau_true)
RBperMX <- RBias_percent(MX_NIPW, tau_true)

RBperRR <- RBias_percent(RR_NIPW, tau_true)
RBperRE <- RBias_percent(RE_NIPW, tau_true)
RBperRM <- RBias_percent(RM_NIPW, tau_true)
RBperRX <- RBias_percent(RX_NIPW, tau_true)

RBperXX <- RBias_percent(XX_NIPW, tau_true)
RBperXE <- RBias_percent(XE_NIPW, tau_true)
RBperXM <- RBias_percent(XM_NIPW, tau_true)
RBperXR <- RBias_percent(XR_NIPW, tau_true)


mean(RBperEE)
mean(RBperEM)
mean(RBperER)
mean(RBperEX)

mean(RBperMM)
mean(RBperME)
mean(RBperMR)
mean(RBperMX)

mean(RBperRR)
mean(RBperRE)
mean(RBperRM)
mean(RBperRX)

mean(RBperXX)
mean(RBperXE)
mean(RBperXM)
mean(RBperXR)

boxplot(
  cbind(ee = RBperEE,
        em = RBperEM,
        er = RBperER,
        ex = RBperEX,

        me = RBperME,
        mm = RBperMM,
        mr = RBperMR,
        mx = RBperMX,

        rr = RBperRR,
        re = RBperRE,
        rm = RBperRM,
        rx = RBperRX,

        xe = RBperXE,
        xm = RBperXM,
        xr = RBperXR,
        xx = RBperXX,

        d = RBperDirect),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# AIPW
# EBLUP --------------------------------------
RBperEEM <- RBias_percent(EEM_AIPW, tau_true)
RBperEER <- RBias_percent(EER_AIPW, tau_true)
RBperEEX <- RBias_percent(EEX_AIPW, tau_true)

RBperEMM <- RBias_percent(EMM_AIPW, tau_true)
RBperEMR <- RBias_percent(EMR_AIPW, tau_true)
RBperEMX <- RBias_percent(EMX_AIPW, tau_true)

RBperERM <- RBias_percent(ERM_AIPW, tau_true)
RBperERR <- RBias_percent(ERR_AIPW, tau_true)
RBperERX <- RBias_percent(ERX_AIPW, tau_true)

RBperEXM <- RBias_percent(EXM_AIPW, tau_true)
RBperEXR <- RBias_percent(EXR_AIPW, tau_true)
RBperEXX <- RBias_percent(EXX_AIPW, tau_true)

# MQ ----------------------------------------
RBperMEE = RBias_percent(MEE_AIPW, tau_true)
RBperMER = RBias_percent(MER_AIPW, tau_true)
RBperMEX = RBias_percent(MEX_AIPW, tau_true)

RBperMME = RBias_percent(MME_AIPW, tau_true)
RBperMMR = RBias_percent(MMR_AIPW, tau_true)
RBperMMX = RBias_percent(MMX_AIPW, tau_true)

RBperMRE = RBias_percent(MRE_AIPW, tau_true)
RBperMRR = RBias_percent(MRR_AIPW, tau_true)
RBperMRX = RBias_percent(MRX_AIPW, tau_true)

RBperMXE = RBias_percent(MXE_AIPW, tau_true)
RBperMXR = RBias_percent(MXR_AIPW, tau_true)
RBperMXX = RBias_percent(MXX_AIPW, tau_true)

# R ------------------------------------------
RBperREE = RBias_percent(REE_AIPW, tau_true)
RBperREM = RBias_percent(REM_AIPW, tau_true)
RBperREX = RBias_percent(REX_AIPW, tau_true)

RBperRME = RBias_percent(RME_AIPW, tau_true)
RBperRMM = RBias_percent(RMM_AIPW, tau_true)
RBperRMX = RBias_percent(RMX_AIPW, tau_true)

RBperRRE = RBias_percent(RRE_AIPW, tau_true)
RBperRRM = RBias_percent(RRM_AIPW, tau_true)
RBperRRX = RBias_percent(RRX_AIPW, tau_true)

RBperRXE = RBias_percent(RXE_AIPW, tau_true)
RBperRXM = RBias_percent(RXM_AIPW, tau_true)
RBperRXX = RBias_percent(RXX_AIPW, tau_true)

# X ------------------------------------------
RBperXEE = RBias_percent(XEE_AIPW, tau_true)
RBperXEM = RBias_percent(XEM_AIPW, tau_true)
RBperXER = RBias_percent(XER_AIPW, tau_true)

RBperXME = RBias_percent(XME_AIPW, tau_true)
RBperXMM = RBias_percent(XMM_AIPW, tau_true)
RBperXMR = RBias_percent(XMR_AIPW, tau_true)

RBperXRE = RBias_percent(XRE_AIPW, tau_true)
RBperXRM = RBias_percent(XRM_AIPW, tau_true)
RBperXRR = RBias_percent(XRR_AIPW, tau_true)

RBperXXE = RBias_percent(XXE_AIPW, tau_true)
RBperXXM = RBias_percent(XXM_AIPW, tau_true)
RBperXXR = RBias_percent(XXR_AIPW, tau_true)

#
boxplot(
  cbind(eem = RBperEEM,
        eer = RBperEER,
        eex = RBperEEX,

        emm = RBperEMM,
        emr = RBperEMR,
        emx = RBperEMX,

        erm = RBperERM,
        err = RBperERR,
        erx = RBperERX,

        exm = RBperEXM,
        exr = RBperEXR,
        eex = RBperEXX,

        # MQ ----------------------------------------
        mee = RBperMEE,
        mer = RBperMER,
        mex = RBperMEX,

        mme = RBperMME,
        mmr = RBperMMR,
        mmx = RBperMMX,

        mre = RBperMRE,
        mrr = RBperMRR,
        mrx = RBperMRX,

        mxe = RBperMXE,
        mxr = RBperMXR,
        mxx = RBperMXX,

        # R ------------------------------------------
        ree = RBperREE,
        rem = RBperREM,
        rex = RBperREX,

        rme = RBperRME,
        rmm = RBperRMM,
        rmx = RBperRMX,

        rre = RBperRRE,
        rrm = RBperRRM,
        rrx = RBperRRX,

        rxe = RBperRXE,
        rxm = RBperRXM,
        rxx = RBperRXX,

        # X ------------------------------------------
        xee = RBperXEE,
        xem = RBperXEM,
        xer = RBperXER,

        xme = RBperXME,
        xmm = RBperXMM,
        xmr = RBperXMR,

        xre = RBperXRE,
        xrm = RBperXRM,
        xrr = RBperXRR,

        xxe = RBperXXE,
        xxm = RBperXXM,
        xxr = RBperXXR,

        d = RBperDirect),
  main = "Relative Bias in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

########
# Bias #
########
##########
# EBLUP -------------------------------------------
#########
boxplot(
  cbind(eem = RBperEEM,
        eer = RBperEER,
        eex = RBperEEX,

        emm = RBperEMM,
        emr = RBperEMR,
        emx = RBperEMX,

        erm = RBperERM,
        err = RBperERR,
        erx = RBperERX,

        exm = RBperEXM,
        exr = RBperEXR,
        eex = RBperEXX,

        d = RBperDirect),
  main = "Relative Bias in %",
  cex.axis = 0.8,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
##########
# MQ ---------------------------------------
#########
boxplot(
  cbind(# MQ ----------------------------------------
        mee = RBperMEE,
        mer = RBperMER,
        mex = RBperMEX,

        mme = RBperMME,
        mmr = RBperMMR,
        mmx = RBperMMX,

        mre = RBperMRE,
        mrr = RBperMRR,
        mrx = RBperMRX,

        mxe = RBperMXE,
        mxr = RBperMXR,
        mxx = RBperMXX,

        d = RBperDirect),
  main = "Relative Bias in %",
  cex.axis = 0.8,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

##########
# RF ---------------------------------------
#########
boxplot(
  cbind( # R ------------------------------------------
         ree = RBperREE,
         rem = RBperREM,
         rex = RBperREX,

         rme = RBperRME,
         rmm = RBperRMM,
         rmx = RBperRMX,

         rre = RBperRRE,
         rrm = RBperRRM,
         rrx = RBperRRX,

         rxe = RBperRXE,
         rxm = RBperRXM,
         rxx = RBperRXX,

         d = RBperDirect),
  main = "Relative Bias in %",
  cex.axis = 0.8,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)


##########
# XGB ---------------------------------------
#########
boxplot(
  cbind( # X ------------------------------------------
         xee = RBperXEE,
         xem = RBperXEM,
         xer = RBperXER,

         xme = RBperXME,
         xmm = RBperXMM,
         xmr = RBperXMR,

         xre = RBperXRE,
         xrm = RBperXRM,
         xrr = RBperXRR,

         xxe = RBperXXE,
         xxm = RBperXXM,
         xxr = RBperXXR,

         d = RBperDirect),
  main = "Relative Bias in %",
  cex.axis = 0.8,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

########
# RMSE #
########
#Dir
RMSEperDirect <- RRMSE_percent(tauDirect, tau_true)

# OR ----------------------------------------
RMSEperE <- RRMSE_percent(E_OR, tau_true)
RMSEperM <- RRMSE_percent(M_OR, tau_true)
RMSEperR <- RRMSE_percent(R_OR, tau_true)
RMSEperX <- RRMSE_percent(X_OR, tau_true)

tau_trueBig <- matrix(rep(tau_true, 1000),
                      nrow = 1000,
                      ncol = 100, byrow = TRUE)
pp <- sqrt(colMeans((E_OR - tau_trueBig)^2))/abs(tau_true) * 100



RRMSE_percent <- function(estimator, true_value) {

  value <- (sqrt(apply(apply(estimator, 1, function(x) {
    (x - true_value) ^ 2
  }), 1, mean, na.rm = TRUE
  )) / abs(true_value)) * 100
  value
}

estimator = E_OR
true_value = tau_true
kk <- apply(estimator, 1, function(x) {
  (x - true_value) ^ 2
})


par(mfrow = c(1, 1))
#
boxplot(
  cbind(e = RMSEperE,
        m = RMSEperM,
        r = RMSEperR,
        x = RMSEperX,

        d = RMSEperDirect),
  main = "Relative MSE in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)


# IPW ---------------------------------------
# NIPW ----------------------
RMSEperEE <- RRMSE_percent(EE_NIPW, tau_true)
RMSEperEM <- RRMSE_percent(EM_NIPW, tau_true)
RMSEperER <- RRMSE_percent(ER_NIPW, tau_true)
RMSEperEX <- RRMSE_percent(EX_NIPW, tau_true)

RMSEperMM <- RRMSE_percent(MM_NIPW, tau_true)
RMSEperME <- RRMSE_percent(ME_NIPW, tau_true)
RMSEperMR <- RRMSE_percent(MR_NIPW, tau_true)
RMSEperMX <- RRMSE_percent(MX_NIPW, tau_true)

RMSEperRR <- RRMSE_percent(RR_NIPW, tau_true)
RMSEperRE <- RRMSE_percent(RE_NIPW, tau_true)
RMSEperRM <- RRMSE_percent(RM_NIPW, tau_true)
RMSEperRX <- RRMSE_percent(RX_NIPW, tau_true)

RMSEperXX <- RRMSE_percent(XX_NIPW, tau_true)
RMSEperXE <- RRMSE_percent(XE_NIPW, tau_true)
RMSEperXM <- RRMSE_percent(XM_NIPW, tau_true)
RMSEperXR <- RRMSE_percent(XR_NIPW, tau_true)

#
boxplot(
  cbind(ee = RMSEperEE,
        em = RMSEperEM,
        er = RMSEperER,
        ex = RMSEperEX,

        me = RMSEperME,
        mm = RMSEperMM,
        mr = RMSEperMR,
        mx = RMSEperMX,

        rr = RMSEperRR,
        re = RMSEperRE,
        rm = RMSEperRM,
        rx = RMSEperRX,

        xe = RMSEperXE,
        xm = RMSEperXM,
        xr = RMSEperXR,
        xx = RMSEperXX,

        d = RMSEperDirect),
  main = "Relative RMSE in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5,
  ylime = c(-1, 200)
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

# AIPW
# EBLUP --------------------------------------
RMSEperEEM <- RRMSE_percent(EEM_AIPW, tau_true)
RMSEperEER <- RRMSE_percent(EER_AIPW, tau_true)
RMSEperEEX <- RRMSE_percent(EEX_AIPW, tau_true)

RMSEperEMM <- RRMSE_percent(EMM_AIPW, tau_true)
RMSEperEMR <- RRMSE_percent(EMR_AIPW, tau_true)
RMSEperEMX <- RRMSE_percent(EMX_AIPW, tau_true)

RMSEperERM <- RRMSE_percent(ERM_AIPW, tau_true)
RMSEperERR <- RRMSE_percent(ERR_AIPW, tau_true)
RMSEperERX <- RRMSE_percent(ERX_AIPW, tau_true)

RMSEperEXM <- RRMSE_percent(EXM_AIPW, tau_true)
RMSEperEXR <- RRMSE_percent(EXR_AIPW, tau_true)
RMSEperEXX <- RRMSE_percent(EXX_AIPW, tau_true)

#Means

mRMSEperEEM <- mean(RMSEperEEM)
mRMSEperEER <- mean(RMSEperEER)
mRMSEperEEX <- mean(RMSEperEER)

mRMSEperEMM <- mean(RMSEperEMM)
mRMSEperEMR <- mean(RMSEperEMR)
mRMSEperEMX <- mean(RMSEperEMX)

mRMSEperERM <- mean(RMSEperERM)
mRMSEperERR <- mean(RMSEperERR)
mRMSEperERX <- mean(RMSEperERX)

mRMSEperEXM <- mean(RMSEperEXM)
mRMSEperEXR <- mean(RMSEperEXR)
mRMSEperEXX <- mean(RMSEperEXX)

#Median
mRMSEperEEM <- quantile(RMSEperEEM, 0.5)
mRMSEperEER <- quantile(RMSEperEER, 0.5)
mRMSEperEEX <- quantile(RMSEperEER, 0.5)

mRMSEperEMM <- quantile(RMSEperEMM, 0.5)
mRMSEperEMR <- quantile(RMSEperEMR, 0.5)
mRMSEperEMX <- quantile(RMSEperEMX, 0.5)

mRMSEperERM <- quantile(RMSEperERM, 0.5)
mRMSEperERR <- quantile(RMSEperERR, 0.5)
mRMSEperERX <- quantile(RMSEperERX, 0.5)

mRMSEperEXM <- quantile(RMSEperEXM, 0.5)
mRMSEperEXR <- quantile(RMSEperEXR, 0.5)
mRMSEperEXX <- quantile(RMSEperEXX, 0.5)

cE <- c(mRMSEperEEM, mRMSEperEER, mRMSEperEEX,
        mRMSEperEMM, mRMSEperEMR, mRMSEperEMX,
        mRMSEperERM, mRMSEperERR, mRMSEperERX,
        mRMSEperEXM, mRMSEperEXR, mRMSEperEXX)

plot(1:12, cE)
# MQ ----------------------------------------
RMSEperMEE = RRMSE_percent(MEE_AIPW, tau_true)
RMSEperMER = RRMSE_percent(MER_AIPW, tau_true)
RMSEperMEX = RRMSE_percent(MEX_AIPW, tau_true)

RMSEperMME = RRMSE_percent(MME_AIPW, tau_true)
RMSEperMMR = RRMSE_percent(MMR_AIPW, tau_true)
RMSEperMMX = RRMSE_percent(MMX_AIPW, tau_true)

RMSEperMRE = RRMSE_percent(MRE_AIPW, tau_true)
RMSEperMRR = RRMSE_percent(MRR_AIPW, tau_true)
RMSEperMRX = RRMSE_percent(MRX_AIPW, tau_true)

RMSEperMXE = RRMSE_percent(MXE_AIPW, tau_true)
RMSEperMXR = RRMSE_percent(MXR_AIPW, tau_true)
RMSEperMXX = RRMSE_percent(MXX_AIPW, tau_true)

# R ------------------------------------------
RMSEperREE = RRMSE_percent(REE_AIPW, tau_true)
RMSEperREM = RRMSE_percent(REM_AIPW, tau_true)
RMSEperREX = RRMSE_percent(REX_AIPW, tau_true)

RMSEperRME = RRMSE_percent(RME_AIPW, tau_true)
RMSEperRMM = RRMSE_percent(RMM_AIPW, tau_true)
RMSEperRMX = RRMSE_percent(RMX_AIPW, tau_true)

RMSEperRRE = RRMSE_percent(RRE_AIPW, tau_true)
RMSEperRRM = RRMSE_percent(RRM_AIPW, tau_true)
RMSEperRRX = RRMSE_percent(RRX_AIPW, tau_true)

RMSEperRXE = RRMSE_percent(RXE_AIPW, tau_true)
RMSEperRXM = RRMSE_percent(RXM_AIPW, tau_true)
RMSEperRXX = RRMSE_percent(RXX_AIPW, tau_true)

# X ------------------------------------------
RMSEperXEE = RRMSE_percent(XEE_AIPW, tau_true)
RMSEperXEM = RRMSE_percent(XEM_AIPW, tau_true)
RMSEperXER = RRMSE_percent(XER_AIPW, tau_true)

RMSEperXME = RRMSE_percent(XME_AIPW, tau_true)
RMSEperXMM = RRMSE_percent(XMM_AIPW, tau_true)
RMSEperXMR = RRMSE_percent(XMR_AIPW, tau_true)

RMSEperXRE = RRMSE_percent(XRE_AIPW, tau_true)
RMSEperXRM = RRMSE_percent(XRM_AIPW, tau_true)
RMSEperXRR = RRMSE_percent(XRR_AIPW, tau_true)

RMSEperXXE = RRMSE_percent(XXE_AIPW, tau_true)
RMSEperXXM = RRMSE_percent(XXM_AIPW, tau_true)
RMSEperXXR = RRMSE_percent(XXR_AIPW, tau_true)





par(mfrow = c(2,2))
par(mfrow = c(1,1))
##########
# All #
#######
boxplot(
  cbind(# EBLUP--------------------
        eem = RMSEperEEM,
        eer = RMSEperEER,
        eex = RMSEperEEX,

        emm = RMSEperEMM,
        emr = RMSEperEMR,
        emx = RMSEperEMX,

        erm = RMSEperERM,
        err = RMSEperERR,
        erx = RMSEperERX,

        exm = RMSEperEXM,
        exr = RMSEperEXR,
        eex = RMSEperEXX,

        # MQ-----------------
        mee = RMSEperMEE,
        mer = RMSEperMER,
        mex = RMSEperMEX,

        mme = RMSEperMME,
        mmr = RMSEperMMR,
        mmx = RMSEperMMX,

        mre = RMSEperMRE,
        mrr = RMSEperMRR,
        mrx = RMSEperMRX,

        mxe = RMSEperMXE,
        mxr = RMSEperMXR,
        mxx = RMSEperMXX,

        #RF-------------------
        ree = RMSEperREE,
        rem = RMSEperREM,
        rex = RMSEperREX,

        rme = RMSEperRME,
        rmm = RMSEperRMM,
        rmx = RMSEperRMX,

        rre = RMSEperRRE,
        rrm = RMSEperRRM,
        rrx = RMSEperRRX,

        rxe = RMSEperRXE,
        rxm = RMSEperRXM,
        rxx = RMSEperRXX,

        # Xgb -----------------
        xee = RMSEperXEE,
        xem = RMSEperXEM,
        xer = RMSEperXER,

        xme = RMSEperXME,
        xmm = RMSEperXMM,
        xmr = RMSEperXMR,

        xre = RMSEperXRE,
        xrm = RMSEperXRM,
        xrr = RMSEperXRR,

        xxe = RMSEperXXE,
        xxm = RMSEperXXM,
        xxr = RMSEperXXR,

        d = RMSEperDirect),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 0.8,
  font.lab = 2,
  lwd = 1.5,
  ylim = c(0, 200)
)

abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)


##########
# EBLUP -------------------------------------------
#########
boxplot(
  cbind(eem = RMSEperEEM,
        eer = RMSEperEER,
        eex = RMSEperEEX,

        emm = RMSEperEMM,
        emr = RMSEperEMR,
        emx = RMSEperEMX,

        erm = RMSEperERM,
        err = RMSEperERR,
        erx = RMSEperERX,

        exm = RMSEperEXM,
        exr = RMSEperEXR,
        eex = RMSEperEXX,

        d = RMSEperDirect),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 0.8,
  font.lab = 2,
  lwd = 1.5,
  ylim = c(0, 300)
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
##########
# MQ ---------------------------------------
#########
boxplot(
  cbind(# MQ ----------------------------------------
        mee = RMSEperMEE,
        mer = RMSEperMER,
        mex = RMSEperMEX,

        mme = RMSEperMME,
        mmr = RMSEperMMR,
        mmx = RMSEperMMX,

        mre = RMSEperMRE,
        mrr = RMSEperMRR,
        mrx = RMSEperMRX,

        mxe = RMSEperMXE,
        mxr = RMSEperMXR,
        mxx = RMSEperMXX,

        d = RMSEperDirect),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 0.8,
  font.lab = 2,
  lwd = 1.5,
  ylim = c(0, 250)
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

##########
# RF ---------------------------------------
#########
boxplot(
  cbind( # R ------------------------------------------
         ree = RMSEperREE,
         rem = RMSEperREM,
         rex = RMSEperREX,

         rme = RMSEperRME,
         rmm = RMSEperRMM,
         rmx = RMSEperRMX,

         rre = RMSEperRRE,
         rrm = RMSEperRRM,
         rrx = RMSEperRRX,

         rxe = RMSEperRXE,
         rxm = RMSEperRXM,
         rxx = RMSEperRXX,

         d = RMSEperDirect),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 0.8,
  font.lab = 2,
  lwd = 1.5,
  ylim = c(0, 355)
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)


##########
# XGB ---------------------------------------
#########
boxplot(
  cbind( # X ------------------------------------------
         xee = RMSEperXEE,
         xem = RMSEperXEM,
         xer = RMSEperXER,

         xme = RMSEperXME,
         xmm = RMSEperXMM,
         xmr = RMSEperXMR,

         xre = RMSEperXRE,
         xrm = RMSEperXRM,
         xrr = RMSEperXRR,

         xxe = RMSEperXXE,
         xxm = RMSEperXXM,
         xxr = RMSEperXXR,

         d = RMSEperDirect),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 0.8,
  font.lab = 2,
  lwd = 1.5,
  ylim = c(0, 355)
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)

