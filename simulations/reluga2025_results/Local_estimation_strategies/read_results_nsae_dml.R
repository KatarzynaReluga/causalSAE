setwd(".DML/csda")

library("RColorBrewer")
library(data.table)

####################
### Load results ###
####################
file_list <- list.files()
SimNum = 1000
res <- list()
for (i in 1:SimNum) {
  load(file_list[i])
  
  ResultsM <- Results[c(10:13, 17, 21, 25, 29, 33, 37, 41, 45)]
  resM <- data.frame(do.call(rbind, ResultsM))
  colnames(resM) <- paste("C", 1:dim(resM)[2], sep = "")
  nMethotds <- dim(resM)[1]
  nClusters <- dim(resM)[2]
  tau_trueM <- matrix(rep(ResultsM$tau_true, nMethotds),
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
  
}


# Join results and compute
dt_res <- data.table(do.call(rbind, res))
dt_Mean <- dt_res[, lapply(.SD, mean), by = "Method"]

# 41
tau_trueM <- dt_Mean[c(4:12), c(125:165)]
Methods <- dt_res$Method[4:12]

TMLE_df <- as.data.frame(dt_res[Method == "XXcsda_DML", .SD, .SDcols = 1:41])


####################################################
# D = 25
# Bias -----------------------------------------------
B <- dt_Mean[c(4:12), c(43:83)]
BMean <- rowMeans(B)
BMedian <- apply(B, 1, quantile, probs = 0.5)
# Absolute Bias ----------------------------------------------
AB <- abs(B)
ABMean <- rowMeans(AB)
ABMedian <- apply(AB, 1, quantile, probs = 0.5)
# Bias Relative ----------------------------------------------
RB <- AB/abs(tau_trueM)
RBMean <- rowMeans(RB)
RBMedian <- apply(RB, 1, quantile, probs = 0.5)
# MSE --------------------------------
MSE <- dt_Mean[c(4:12), c(84:124)]
MSEMean <- rowMeans(MSE)
MSEMedian <- apply(MSE, 1, quantile, probs = 0.5)

RRMSEi <- sqrt(MSE)/abs(tau_trueM)
RRMSEMean <- rowMeans(RRMSEi)

# csda
results_csda <- data.frame(Method = Methods,
                         mse = MSEMean,
                         Bias = BMean,
                         ABias = ABMean,
                         RABias = RBMean,
                         rrmse = RRMSEMean)

# E
results_E2 <- data.frame(Method = Methods,
                        mse = MSEMean,
                        Bias = BMean,
                        ABias = ABMean,
                        RABias = RBMean,
                        rrmse = RRMSEMean)

# E
results_E <- data.frame(Method = Methods,
                          mse = MSEMean,
                          Bias = BMean,
                          ABias = ABMean,
                          RABias = RBMean,
                          rrmse = RRMSEMean)

# L 
results_L <- data.frame(Method = Methods,
                        mse = MSEMean,
                        Bias = BMean,
                        ABias = ABMean,
                        RABias = RBMean, 
                        rrmse = RRMSEMean)

# M
results_M <- data.frame(Method = Methods,
                        mse = MSEMean,
                        Bias = BMean,
                        ABias = ABMean,
                        RABias = RBMean, 
                        rrmse = RRMSEMean)

# Mq
results_Mq <- data.frame(Method = Methods,
                        mse = MSEMean,
                        Bias = BMean,
                        ABias = ABMean,
                        RABias = RBMean,
                        rrmse = RRMSEMean)
# R
results_R <- data.frame(Method = Methods,
                         mse = MSEMean,
                         Bias = BMean,
                         ABias = ABMean,
                         RABias = RBMean,
                        rrmse = RRMSEMean)


# Rc
results_Rc <- data.frame(Method = Methods,
                           mse = MSEMean,
                           Bias = BMean,
                           ABias = ABMean,
                           RABias = RBMean,
                         rrmse = RRMSEMean)


# Rct
results_Rct <- data.frame(Method = Methods,
                            mse = MSEMean,
                            Bias = BMean,
                            ABias = ABMean,
                            RABias = RBMean,
                          rrmse = RRMSEMean)


# Rt
results_Rt <- data.frame(Method = Methods,
                         mse = MSEMean,
                         Bias = BMean,
                         ABias = ABMean,
                         RABias = RBMean,
                         rrmse = RRMSEMean)


# X
results_X <- data.frame(Method = Methods,
                          mse = MSEMean,
                          Bias = BMean,
                          ABias = ABMean,
                          RABias = RBMean,
                        rrmse = RRMSEMean)

#Xt

results_Xt <- data.frame(Method = Methods,
                           mse = MSEMean,
                           Bias = BMean,
                           ABias = ABMean,
                           RABias = RBMean,
                         rrmse = RRMSEMean)

#S

results_S <- data.frame(Method = Methods,
                         mse = MSEMean,
                         Bias = BMean,
                         ABias = ABMean,
                         RABias = RBMean,
                        rrmse = RRMSEMean)


resultsF <- rbind(results_E, results_L, results_M, results_Mq, 
                    results_R, results_Rc, results_Rct, results_Rt, 
                  results_S, results_X, results_Xt)


resultsF <- rbind(results_E2, results_csda)

resultsF$Impute <- c(rep("E", 9),
                     rep("L", 9),
                     rep("M", 9),
                     rep("Mq", 9),
                     
                     rep("R", 9),
                     rep("Rc", 9),
                     rep("Rct", 9),
                     rep("Rt", 9),
                     
                     rep("S", 9),
                     rep("X", 9),
                     rep("Xt", 9))

resultsF$Impute <- c(rep("2E", 9),
                     rep("csda", 9))

