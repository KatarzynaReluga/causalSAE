setwd("./results/ModelBased/DBscenarios/DB25t")
library("RColorBrewer")
library(data.table)

####################
### Load results ###
####################
SimNum = 1000
res <- list()

file_list <- list.files()
for (i in 1:SimNum) {
  load(file_list[i])

  ResultsM <- Results[- c(1:10)]
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

# 25
tau_trueM <- dt_Mean[c(4:72), c(77:101)]
# 50
tau_trueM <- dt_Mean[c(4:72), c(152:201)]
# 75
tau_trueM <- dt_Mean[c(4:72), c(227:301)]

Methods <- dt_res$Method[4:72]

# Bias ----------------------------------------------
# 25
RB <- dt_Mean[c(4:72), c(27:51)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(27:51)]/abs(tau_trueM) * 100

# 50
RB <- dt_Mean[c(4:72), c(52:101)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(52:101)]/abs(tau_trueM) * 100

# 75
RB <- dt_Mean[c(4:72), c(77:151)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(77:151)]/abs(tau_trueM) * 100


# Relative Mean -------------------------------------------
RBPMean <- rowMeans(RBP)
# Relative bias ------------------------------------------
RBPMedian <- apply(RBP, 1, quantile, probs = 0.5)

# Absolute Bias ----------------------------------------------
# 25
ARB <- abs(dt_Mean[c(4:72), c(27:51)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:72), c(27:51)]/abs(tau_trueM) * 100)

# 50
ARB <- abs(dt_Mean[c(4:72), c(52:101)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:72), c(52:101)]/abs(tau_trueM) * 100)

#75
ARB <- abs(dt_Mean[c(4:72), c(77:151)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:72), c(77:151)]/abs(tau_trueM) * 100)

# Means --------------------------------------------------
ARBPMean <- rowMeans(ARBP)
# Medians ------------------------------------------------
ARBPMedian <- apply(ARBP, 1, quantile, probs = 0.5)

# RRMSE ----------------------------------------------
# 25
RMSE <- sqrt(dt_Mean[c(4:72), c(52:76)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:72), c(52:76)])/abs(tau_trueM) * 100

# 50
RMSE <- sqrt(dt_Mean[c(4:72), c(102:151)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:72), c(102:151)])/abs(tau_trueM) * 100

# 75
RMSE <- sqrt(dt_Mean[c(4:72), c(152:226)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:72), c(152:226)])/abs(tau_trueM) * 100

RMSEMean <- rowMeans(RMSEP)
RMSEPMedian <- apply(RMSEP, 1, quantile, probs = 0.5)


resultsDB25t <- data.frame(Method = Methods,

           RBPMean = RBPMean,
           RBPMedian = RBPMedian,

           ARBPMean = ARBPMean,
           ARBPMedian = ARBPMedian,

           RMSEMean = RMSEMean,
           RMSEPMedian = RMSEPMedian)

write.csv(resultsDB25t,
          file = 'resultsDB25t.csv', row.names = FALSE)
