#setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/OR")
setwd("C:/Users/katar/Documents/GitHub/causalSAE/simulations/results_test")
#setwd("C:/Users/katar/Documents/GitHub/causalSAE")

# Set seed
#set.seed(100)

#m = 50
#ni = rep(10, m)
#Ni = rep(100, m)
#N = sum(Ni)
#n = sum(ni)

#n_boot = 500

# Generate covariates
#X <- generate_X(
#  n = N,
#  p = 1,
#  covariance_norm = NULL,
#  cov_type = "unif",
#  seed = 1
#)

#X_outcome <- generate_X(
#  n = N,
#  p = 1,
#  covariance_norm = NULL,
#  cov_type = "lognorm",
#  seed = 1
#)

# Generate populations
#populations <- generate_pop(X, X_outcome,
#                            coeffs = get_default_coeffs(),
#                            errors_outcome = get_default_errors_outcome(),
#                            rand_eff_outcome = get_default_rand_eff_outcome(),
#                            rand_eff_p_score = get_default_rand_eff_p_score(),
#                            regression_type = "continuous",
#                            Ni_size = 100,
#                            m = 50,
#                            no_sim = 1,
#                            seed = 10)

#tau_true <- calculate_tau(list(populations), type_tau = "H")[[1]]$tau


RBias_percent <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RB.EBLUP <- (apply(EBLUP_tau, 2, mean) - True) / abs(True) * 100
  value <- (apply(estimator, 2, mean) - true_value) / abs(true_value) * 100
  value

}

RBias <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RB.EBLUP <- (apply(EBLUP_tau, 2, mean) - True) / abs(True) * 100
  value <- (apply(estimator, 2, mean) - true_value) / abs(true_value)
  value

}

RRMSE_percent <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RRMSE.Direct <-
  #  (sqrt(apply(
  #    apply(Direct, 1, function(x) {
  #      (x - True) ^ 2
  #    }), 1, mean, na.rm = TRUE
  #  )) / abs(True)) * 100
  value <- (sqrt(apply(apply(estimator, 1, function(x) {
    (x - true_value) ^ 2
  }), 1, mean, na.rm = TRUE
  )) / abs(true_value)) * 100
  value
}

RRMSE <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RRMSE.Direct <-
  #  (sqrt(apply(
  #    apply(Direct, 1, function(x) {
  #      (x - True) ^ 2
  #    }), 1, mean, na.rm = TRUE
  #  )) / abs(True)) * 100
  value <- (sqrt(apply(apply(estimator, 1, function(x) {
    (x - true_value) ^ 2
  }), 1, mean, na.rm = TRUE
  )) / abs(true_value))
  value
}

RRMSE_percent <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RRMSE.Direct <-
  #  (sqrt(apply(
  #    apply(Direct, 1, function(x) {
  #      (x - True) ^ 2
  #    }), 1, mean, na.rm = TRUE
  #  )) / abs(True)) * 100
  value <- (sqrt(apply(apply(estimator, 1, function(x) {
    (x - true_value) ^ 2
  }), 1, mean, na.rm = TRUE
  )) / abs(true_value)) * 100
  value
}

MSE <- function(estimator, true_value) {

  #True = tau_EB[[1]]$tau
  #RRMSE.Direct <-
  #  (sqrt(apply(
  #    apply(Direct, 1, function(x) {
  #      (x - True) ^ 2
  #    }), 1, mean, na.rm = TRUE
  #  )) / abs(True)) * 100
  value <- apply(apply(estimator, 1, function(x) {
    (x - true_value) ^ 2
  }), 1, mean, na.rm = TRUE
  )
  value
}

# Retrive vectors ------------------------------------------------------
SimNum = 500
D = 50

tau_true <- matrix(0, nrow = SimNum, ncol = D)

tauEBLUP <- matrix(0, nrow = SimNum, ncol = D)
tauMQ <- matrix(0, nrow = SimNum, ncol = D)
tauRF <- matrix(0, nrow = SimNum, ncol = D)
tauXGB <- matrix(0, nrow = SimNum, ncol = D)

tauEEIPW <- matrix(0, nrow = SimNum, ncol = D)
tauEMIPW <- matrix(0, nrow = SimNum, ncol = D)
tauERIPW <- matrix(0, nrow = SimNum, ncol = D)
tauEXIPW <- matrix(0, nrow = SimNum, ncol = D)

tauMMIPW <- matrix(0, nrow = SimNum, ncol = D)
tauMEIPW <- matrix(0, nrow = SimNum, ncol = D)
tauMRIPW <- matrix(0, nrow = SimNum, ncol = D)
tauMXIPW <- matrix(0, nrow = SimNum, ncol = D)

tauRRIPW <- matrix(0, nrow = SimNum, ncol = D)
tauREIPW <- matrix(0, nrow = SimNum, ncol = D)
tauRMIPW <- matrix(0, nrow = SimNum, ncol = D)
tauRXIPW <- matrix(0, nrow = SimNum, ncol = D)

tauXXIPW <- matrix(0, nrow = SimNum, ncol = D)
tauXEIPW <- matrix(0, nrow = SimNum, ncol = D)
tauXMIPW <- matrix(0, nrow = SimNum, ncol = D)
tauXRIPW <- matrix(0, nrow = SimNum, ncol = D)


tauEER <- matrix(0, nrow = SimNum, ncol = D)
tauEEX <- matrix(0, nrow = SimNum, ncol = D)
tauMMR <- matrix(0, nrow = SimNum, ncol = D)
tauMMX <- matrix(0, nrow = SimNum, ncol = D)

tauDirect <- matrix(0, nrow = SimNum, ncol = D)

file_list <- list.files()
for (i in 1:SimNum) {

  load(file_list[i])

  tau_true[i, ] <- Results$tau_true$tau

  tauEEIPW[i, ]  <- Results$EE_NIPW
  tauEMIPW[i, ]  <- Results$EM_NIPW
  tauERIPW[i, ]  <- Results$ER_NIPW
  tauEXIPW[i, ]  <- Results$EX_NIPW

  tauMMIPW[i, ] <- Results$MM_NIPW
  tauMEIPW[i, ] <- Results$ME_NIPW
  tauMRIPW[i, ] <- Results$MR_NIPW
  tauMXIPW[i, ] <- Results$MX_NIPW

  tauRRIPW[i, ] <- Results$RR_NIPW
  tauREIPW[i, ] <- Results$RE_NIPW
  tauRMIPW[i, ] <- Results$RM_NIPW
  tauRXIPW[i, ] <- Results$RX_NIPW

  tauXXIPW[i, ] <- Results$RR_NIPW
  tauXXIPW[i, ] <- Results$RR_NIPW
  tauXXIPW[i, ] <- Results$RR_NIPW
  tauXXIPW[i, ] <- Results$RR_NIPW


#  tauEER[i, ] <- Results$EER_AIPW
#  tauEEX[i, ] <- Results$EEX_AIPW
#  tauMMR[i, ] <- Results$MMR_AIPW
#  tauMMX[i, ] <- Results$MMX_AIPW

  tauDirect[i, ] <- Results$Dir_tau

}
# Check tests --------------------
covE <- mean(colMeans(TestE))
covM <- mean(colMeans(TestM))
covR <- mean(colMeans(TestR))
covX <- mean(colMeans(TestX))
covE
covM
covR
covX
#######################################################
tau_true <- tau_true[1, ]

#RBperEB <- RBias_percent(tauEBLUP, tau_true)
#RBperMQ <- RBias_percent(tauMQ, tau_true)
#RBperRF <- RBias_percent(tauRF, tau_true)
#RBperXGB <- RBias_percent(tauXGB, tau_true)

RBperEE <- RBias_percent(tauEEIPW, tau_true)
RBperEM <- RBias_percent(tauEMIPW, tau_true)
RBperER <- RBias_percent(tauERIPW, tau_true)
RBperEX <- RBias_percent(tauEXIPW, tau_true)

RBperMM <- RBias_percent(tauMMIPW, tau_true)
RBperME <- RBias_percent(tauMEIPW, tau_true)
RBperMR <- RBias_percent(tauMRIPW, tau_true)
RBperMX <- RBias_percent(tauMXIPW, tau_true)

RBperRR <- RBias_percent(tauRRIPW, tau_true)
RBperRE <- RBias_percent(tauREIPW, tau_true)
RBperRM <- RBias_percent(tauRMIPW, tau_true)
RBperRX <- RBias_percent(tauRXIPW, tau_true)

RBperXX <- RBias_percent(tauXXIPW, tau_true)
RBperXE <- RBias_percent(tauXEIPW, tau_true)
RBperXM <- RBias_percent(tauXMIPW, tau_true)
RBperXR <- RBias_percent(tauXRIPW, tau_true)

#RBperEER <- RBias_percent(tauEER, tau_true)
#RBperEEX <- RBias_percent(tauEEX, tau_true)
#RBperMMR <- RBias_percent(tauMMR, tau_true)
#RBperMMX <- RBias_percent(tauMMX, tau_true)

RBperDirect <- RBias_percent(tauDirect, tau_true)


RBEBLUP <- RBias(tauEBLUP, tau_true)
RBMQ <- RBias(tauMQ, tau_true)
RBRF <- RBias(tauRF, tau_true)
RBXGB <- RBias(tauXGB, tau_true)
#######################################################
#RMSEperEBLUP <- RRMSE_percent(tauEBLUP, tau_true)
#RMSEperMQ <- RRMSE_percent(tauMQ, tau_true)
#RMSEperRF <- RRMSE_percent(tauRF, tau_true)
#RMSEperXGB <- RRMSE_percent(tauXGB, tau_true)

RMSEperEE <- RRMSE_percent(tauEEIPW, tau_true)
RMSEperEM <- RRMSE_percent(tauEMIPW, tau_true)
RMSEperER <- RRMSE_percent(tauERIPW, tau_true)
RMSEperEX <- RRMSE_percent(tauEXIPW, tau_true)

RMSEperMM <- RRMSE_percent(tauMMIPW, tau_true)
RMSEperME <- RRMSE_percent(tauMEIPW, tau_true)
RMSEperMR <- RRMSE_percent(tauMRIPW, tau_true)
RMSEperMX <- RRMSE_percent(tauMXIPW, tau_true)

RMSEperRR <- RRMSE_percent(tauRRIPW, tau_true)
RMSEperRE <- RRMSE_percent(tauREIPW, tau_true)
RMSEperRM <- RRMSE_percent(tauRMIPW, tau_true)
RMSEperRX <- RRMSE_percent(tauRXIPW, tau_true)

RMSEperXX <- RRMSE_percent(tauXXIPW, tau_true)
RMSEperXE <- RRMSE_percent(tauXEIPW, tau_true)
RMSEperXM <- RRMSE_percent(tauXMIPW, tau_true)
RMSEperXR <- RRMSE_percent(tauXRIPW, tau_true)
#RMSEperEER <- RRMSE_percent(tauEER, tau_true)
#RMSEperEEX <- RRMSE_percent(tauEEX, tau_true)
#RMSEperMMR <- RRMSE_percent(tauMMR, tau_true)
#RMSEperMMX <- RRMSE_percent(tauMMX, tau_true)

RMSEperDirect <- RRMSE_percent(tauDirect, tau_true)

RMSEEBLUP <- RRMSE(tauEBLUP, tau_true)
RMSEMQ <- RRMSE(tauMQ, tau_true)
RMSERF <- RRMSE(tauRF, tau_true)
RMSEXGB <- RRMSE(tauXGB, tau_true)

####################################################
MSEEBLUP <- MSE(tauEBLUP, tau_true)
index_sortE =  sort(MSEEBLUP, decreasing = T, index.return = T)$ix

MSEMQ <- MSE(tauMQ, tau_true)
index_sortM =  sort(MSEMQ, decreasing = T, index.return = T)$ix

MSERF <- MSE(tauRF, tau_true)
index_sortR =  sort(MSERF, decreasing = T, index.return = T)$ix

MSEXGB <- MSE(tauXGB, tau_true)
index_sortX =  sort(MSEXGB, decreasing = T, index.return = T)$ix

EBLUPvar <- colMeans(tauEBLUP_var)
plot(1:D, MSEEBLUP[index_sortE], col = 2, lwd = 2, pch = 0, type = "b")
lines(EBLUPvar[index_sortE], col = 1, lwd = 2, pch = 1, type = "b")


MQvar <- colMeans(tauMQ_var)
plot(1:D, MSEMQ[index_sortM], col = 2, lwd = 2, pch = 0, type = "b")
lines(MQvar[index_sortM], col = 1, lwd = 2, pch = 1, type = "b")


RFvar <- colMeans(tauRF_var)
plot(1:D, MSERF[index_sortR], col = 2, lwd = 2, pch = 0, type = "b")
lines(RFvar[index_sortR], col = 1, lwd = 2, pch = 1, type = "b")

XGBvar <- colMeans(tauXGB_var)
plot(1:D, MSEXGB[index_sortX], col = 2, lwd = 2, pch = 0, type = "b")
lines(XGBvar[index_sortX], col = 1, lwd = 2, pch = 1, type = "b")


#
boxplot(
  cbind(ee = RBperEE,
        em = RBperEM,
        er = RBperER,
        ex = RBperEX,

        mm = RBperMM,
        me = RBperME,
        mr = RBperMR,
        mx = RBperMX,

        rr = RBperRR,
        re = RBperRE,
        rm = RBperRM,
        rx = RBperRX,

        xx = RBperXX,
        xe = RBperXE,
        xr = RBperXR,
        xm = RBperXM,

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


boxplot(
  cbind(e = RMSEperEBLUP,
        m =  RMSEperMQ,
        r = RMSEperRF,
        x = RMSEperXGB,

        ee = RMSEperEE,
        mm = RMSEperMM,
        rr = RMSEperRR,
        xx = RMSEperXX,

        eex = RMSEperEEX,
        eer = RMSEperEER,
        mmx = RMSEperMMX,
        mmr = RMSEperMMR,

        d = RMSEperDirect),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5,
  ylim = c(0, 60)
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)
#

# Test for the same means ---------------------------------------------------
# setwd("C:/Users/katar/Documents/GitHub/causalSAE/simulations/results_test")

# SimNum = 50
# D = 50
# true_tau <- matrix(0, nrow = SimNum, ncol = D)

# file_list <- list.files()
# for (i in 1:SimNum) {

#  load(file_list[i])
#  true_tau[i, ] = Results$tau[[1]]$tau

#}
#means_tau <- colMeans(true_tau)
#MatrixMeansTau <- matrix(rep(means_tau, D), nrow = SimNum, ncol = D, byrow = T)
#colMeans(true_tau - MatrixMeansTau)^2

