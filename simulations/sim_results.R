#setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/OR")
#setwd("C:/Users/katar/Documents/GitHub/causalSAE")
#setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/OR_boot_correct")
setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/OR_sample_boot_correct")
# Set seed
set.seed(100)

m = 50
ni = rep(10, m)
Ni = rep(100, m)
N = sum(Ni)
n = sum(ni)

n_boot = 500

# Generate covariates
X <- generate_X(
  n = N,
  p = 1,
  covariance_norm = NULL,
  cov_type = "unif",
  seed = 1
)

X_outcome <- generate_X(
  n = N,
  p = 1,
  covariance_norm = NULL,
  cov_type = "lognorm",
  seed = 1
)

# Generate populations
populations <- generate_pop(X, X_outcome,
                            coeffs = get_default_coeffs(),
                            errors_outcome = get_default_errors_outcome(),
                            rand_eff_outcome = get_default_rand_eff_outcome(),
                            rand_eff_p_score = get_default_rand_eff_p_score(),
                            regression_type = "continuous",
                            Ni_size = 100,
                            m = 50,
                            no_sim = 1,
                            seed = 10)

tau_true <- calculate_tau(list(populations), type_tau = "H")[[1]]$tau


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
SimNum = 200
D = 50

tau_trueM <- matrix(0, nrow = SimNum, ncol = D)

tauEBLUP <- matrix(0, nrow = SimNum, ncol = D)
tauMQ <- matrix(0, nrow = SimNum, ncol = D)
tauRF <- matrix(0, nrow = SimNum, ncol = D)
tauXGB <- matrix(0, nrow = SimNum, ncol = D)

tauEBLUP_var <- matrix(0, nrow = SimNum, ncol = D)
tauMQ_var <- matrix(0, nrow = SimNum, ncol = D)
tauRF_var <- matrix(0, nrow = SimNum, ncol = D)
tauXGB_var <- matrix(0, nrow = SimNum, ncol = D)

TestE <- matrix(0, nrow = SimNum, ncol = D)
TestM <- matrix(0, nrow = SimNum, ncol = D)
TestR <- matrix(0, nrow = SimNum, ncol = D)
TestX <- matrix(0, nrow = SimNum, ncol = D)

widthE <- matrix(0, nrow = SimNum, ncol = D)
widthM <- matrix(0, nrow = SimNum, ncol = D)
widthR <- matrix(0, nrow = SimNum, ncol = D)
widthX <- matrix(0, nrow = SimNum, ncol = D)


file_list <- list.files()
for (i in 1:SimNum) {

  load(file_list[i])

  tau_trueM[i, ] <- Results$tau_true[[1]]$tau

  tauEBLUP[i, ] <- Results$EBLUP_OR
  tauMQ[i, ] <- Results$MQ_OR
  tauRF[i, ] <- Results$RF_OR
  tauXGB[i, ] <- Results$XGB_OR

  tauEBLUP_var[i, ] <- Results$EBLUP_OR_var
  tauMQ_var[i, ] <- Results$MQ_OR_var
  tauRF_var[i, ] <- Results$RF_OR_var
  tauXGB_var[i, ] <- Results$XGB_OR_var

  tau_true <- Results$tau_true[[1]]$tau

  TestE[i, ] <- ((Results$EBLUP_OR - 1.96 * sqrt(Results$EBLUP_OR_var)) <= tau_true) & (tau_true <= (1.96 * sqrt(Results$EBLUP_OR_var) + Results$EBLUP_OR))
  TestM[i, ] <- ((Results$MQ_OR - 1.96 * sqrt(Results$MQ_OR_var)) <= tau_true) & (tau_true <= (1.96 * sqrt(Results$MQ_OR_var) + Results$MQ_OR))
  TestR[i, ] <- ((Results$RF_OR - 1.96 * sqrt(Results$RF_OR_var)) <= tau_true) & (tau_true <= (1.96 * sqrt(Results$RF_OR_var) + Results$RF_OR))
  TestX[i, ] <- ((Results$XGB_OR - 1.96 * sqrt(Results$XGB_OR_var)) <= tau_true) & (tau_true <= (1.96 * sqrt(Results$XGB_OR_var) + Results$XGB_OR))

  widthE[i, ] <- 2 * 1.96 * sqrt(Results$EBLUP_OR_var)
  widthM[i, ] <- 2 * 1.96 * sqrt(Results$MQ_OR_var)
  widthR[i, ] <- 2 * 1.96 * sqrt(Results$RF_OR_var)
  widthX[i, ] <- 2 * 1.96 * sqrt(Results$XGB_OR_var)

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
RBperEB <- RBias_percent(tauEBLUP, tau_true)
RBperMQ <- RBias_percent(tauMQ, tau_true)
RBperRF <- RBias_percent(tauRF, tau_true)
RBperXGB <- RBias_percent(tauXGB, tau_true)

RBEBLUP <- RBias(tauEBLUP, tau_true)
RBMQ <- RBias(tauMQ, tau_true)
RBRF <- RBias(tauRF, tau_true)
RBXGB <- RBias(tauXGB, tau_true)


#######################################################
RMSEperEBLUP <- RRMSE_percent(tauEBLUP, tau_true)
RMSEperMQ <- RRMSE_percent(tauMQ, tau_true)
RMSEperRF <- RRMSE_percent(tauRF, tau_true)
RMSEperXGB <- RRMSE_percent(tauXGB, tau_true)

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






boxplot(
  cbind(EBLUP = RBperEB,
        MQ =  RBperMQ,
        RF = RBperRF,
        XGB = RBperXGB),
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
  cbind(EBLUP = RMSEperEBLUP,
        MQ =  RMSEperMQ,
        RF = RMSEperRF,
        XGB = RMSEperXGB),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5,
  ylim = c(0, 30)
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

