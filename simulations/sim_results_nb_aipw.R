#setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/OR")
setwd("C:/Users/katar/Documents/GitHub/causalSAE/simulations/results")
#setwd("C:/Users/katar/Documents/GitHub/causalSAE")


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

tau_trueM <- matrix(0, nrow = SimNum, ncol = D)


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

  tau_trueM[i, ] <- Results$tau_true$tau


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
tau_true <- tau_trueM[1, ]


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

RBperDirect <- RBias_percent(tauDirect, tau_true)


RBEBLUP <- RBias(tauEBLUP, tau_true)
RBMQ <- RBias(tauMQ, tau_true)
RBRF <- RBias(tauRF, tau_true)
RBXGB <- RBias(tauXGB, tau_true)
#######################################################
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
par(mfrow = c(2,2))
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
  ylim = c(0, 55)
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
  ylim = c(0, 55)
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
  ylim = c(0, 55)
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
  ylim = c(0, 55)
)
abline(h = 0,
       col = "red",
       lty = "dashed",
       lwd = 2)


boxplot(
  cbind(ee = RMSEperEE,
        em = RMSEperEM,
        er = RMSEperER,
        ex = RMSEperEX,

        mm = RMSEperMM,
        me = RMSEperME,
        mr = RMSEperMR,
        mx = RMSEperMX,

        rr = RMSEperRR,
        re = RMSEperRE,
        rm = RMSEperRM,
        rx = RMSEperRX,

        xx = RMSEperXX,
        xe = RMSEperXE,
        xm = RMSEperXM,
        xr = RMSEperXR,

        d = RMSEperDirect),
  main = "Relative Root Mean Square Error in %",
  cex.axis = 1.5,
  font.lab = 2,
  lwd = 1.5,
  ylim = c(0, 20)
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

