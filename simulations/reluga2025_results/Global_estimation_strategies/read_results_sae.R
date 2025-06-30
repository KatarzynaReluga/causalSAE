setwd("./Model-based/res4/Global/csda_2e")

library("RColorBrewer")
library(data.table)
library(ggplot2)
library(cowplot)
####################
### Load results ###
####################
file_list <- list.files()

SimNum = 1000
res <- list()
for (i in 1:SimNum) {
  load(file_list[i])

  ResultsM <- Results[- c(1:9)]
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
#tau_trueM <- dt_Mean[c(4:1467), c(125:165)]
#Methods <- dt_res$Method[4:1467]

tau_trueM <- dt_Mean[c(4:267), c(125:165)]
Methods <- dt_res$Method[4:267]

AIPW_df <- as.data.frame(dt_res[Method == "MXcsda", .SD, .SDcols = 1:41])

########################
####################################################
# D = 25
# Bias -----------------------------------------------
B <- dt_Mean[c(4:267), c(43:83)]
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
MSE <- dt_Mean[c(4:267), c(84:124)]
MSEMean <- rowMeans(MSE)
MSEMedian <- apply(MSE, 1, quantile, probs = 0.5)

RRMSEi <- sqrt(MSE)/abs(tau_trueM)
RRMSEMean <- rowMeans(RRMSEi)

results_sae <- data.frame(Method = Methods,
                                mse = MSEMean,
                                Bias = BMean,
                                ABias = ABMean,
                                RABias = RBMean,
                          rrmse = RRMSEMean)

results_saeord <- results_sae[order(results_sae$mse, decreasing = FALSE), ]
head(results_saeord)

write.csv(results_sae,
          file = 'results_global_model_csda_2e.csv', row.names = FALSE)


############
# Plot #####
################################################################################
# Attention!!!
# DML_df can be obtained after launching simulations under the local estimation strategy
# and reading the results, see line 47 in R file read_results_nsae_dml.R
# DML_df <- as.data.frame(dt_res[Method == "XXcsda_DML", .SD, .SDcols = 1:41])
#####################################################################################

par(mar = c(5, 5, 1, 1))
plot(as.numeric(tau_trueM[1,]), type = "b", ylim = c(0, 4), lwd = 4,
     xlab = "Subpopulation", ylab = "MSE",
     cex.lab = 3,
     yaxt="n", xaxt="n",
     font.axis = 2,
     font.lab  = 2)
axis(1, at = seq(1, 41,
                 by = 2), cex.axis=2, font.axis = 2)
axis(2, at = seq(0, 4, by = 1), cex.axis=2, font.axis = 2)


tau_trueV <- unlist(as.vector(tau_trueM[1,]))
tau_trueV_sort <- sort(tau_trueV, index.return = TRUE)


par(mgp = c(2.5, 0.7, 0))
plot(1:41, tau_trueV, lwd = 4, type = "b",
     ylim = c(0, 2.5),
     xlab = "Subpopulation",
     ylab = expression(tau[j]),
     cex.lab = 2,    # Enlarges the axis labels
     cex.axis = 1.5,
     xlim = c(1,41))
for (i in  seq(1, 500, 100)) {
  #  lines(tau_hat_AIPW[i, ], lwd = 2, col = "lightgreen") # 13, 16, 17, 21, 25
  #lines(up_AIPW[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25
  #lines(do_AIPW[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25
  lines(unlist(AIPW_df[i, ]), lwd = 2, col = "blue", lty = "dashed") # 13, 16, 17, 21, 25
  lines(unlist(DML_df[i, ]), lwd = 2, col = "red", lty = "dashed") # 13, 16, 17, 21, 25
#  lines(tau_boot_975_AIPWb2M[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25
#  lines(tau_boot_025_NIPWb2M[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25
#  lines(tau_boot_975_NIPWb2M[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25

}
legend("topleft", legend = c(expression(tau[j]),
                         expression(hat(tau)[AIPWj]),
                         expression(hat(tau)[DMLj])),
       fill = c(1, "blue", "red"), bty = "n", cex = 2, horiz = "TRUE")




par(mgp = c(2.5, 0.7, 0))
plot(1:41, tau_trueV[tau_trueV_sort$ix], lwd = 4, type = "b",
     ylim = c(0, 2.5),
     xlab = "Subpopulation",
     ylab = expression(tau[j]),
     cex.lab = 2,    # Enlarges the axis labels
     cex.axis = 1.5,
     xlim = c(1,41))
for (i in  seq(1, 300, 100)) {

  AIPW_i <- unlist(AIPW_df[i, ])
  TMLE_i <- unlist(DML_df[i, ])

  lines(AIPW_i[tau_trueV_sort$ix], lwd = 4, col = "blue",  lty = 4) # 13, 16, 17, 21, 25
  lines(TMLE_i[tau_trueV_sort$ix], lwd = 4, col = "red",  lty = 4) # 13, 16, 17, 21, 25


}
legend("top", legend = c(expression(tau[j]),
                             expression(hat(tau)[AIPWj]),
                             expression(hat(tau)[DMLj])),
       fill = c(1, "blue", "red"), bty = "n", cex = 2, horiz = "TRUE")



## Sort index
ix <- tau_trueV_sort$ix

# Combine all 6 estimators: 6 × 41 matrix
bar_matrix <- rbind(unlist(AIPW_df[1, ix]),
                    unlist(AIPW_df[101, ix]),
                    unlist(AIPW_df[201, ix]),
                    unlist(DML_df[1, ix]),
                    unlist(DML_df[101, ix]),
                    unlist(DML_df[201, ix]))  # dimensions: 6 x 41

# Define light color palette
bar_cols <- c("steelblue1", "steelblue1", "steelblue1",
              "indianred1", "indianred1", "indianred1")



bar_cols <- c(rep("#0072B2", 3),  # 3 AIPW bars
              rep("#E69F00", 3))

bar_cols <- c(rep("#20B2AA", 3),  # AIPW bars (teal)
              rep("#FFB347", 3))

bar_cols <- c(rep("#66C2A5", 3),  # AIPW bars (aqua green)
              rep("#FC8D62", 3))  # TMLE bars (soft coral)
# TMLE bars (gold)

bar_cols <- c(rep("#1CA9C9", 3),  # AIPW bars (cyan)
              rep("#FDAE61", 3))  # TMLE bars (orange)



# Plot
par(mar = c(5, 5, 2, 2), mgp = c(2.5, 0.7, 0), xaxs = "i")
bp <- barplot(bar_matrix, beside = TRUE,
              col = bar_cols,
              border = NA,  # <- removes black borders
              ylim = c(0, 2.5),
              xlab = "Subpopulation",
              ylab = expression(tau[j]),
              cex.lab = 2,
              cex.axis = 1.5,
              las = 1,
              names.arg = rep("", ncol(bar_matrix)),
              space = c(0.2, 1))  # adjust spacing

# Compute center of each group of 6 bars
J <- ncol(bar_matrix)
x_centers <- numeric(J)
for (j in 1:J) {
  idx <- (j - 1) * 6 + 1
  x_centers[j] <- mean(bp[idx:(idx + 5)])
}

# True values
true_vals <- tau_trueV[ix]
points(x_centers, true_vals, pch = 16, col = "black", cex = 1.3)
axis(side = 1, at = x_centers[c(1, 10, 20, 30, 41)],
     labels = c(1, 10, 20, 30, 41),
     cex.axis = 1.5, tick = TRUE)

# Add legend

legend("top", legend = c(expression(tau[j]),
                         expression(hat(tau)[AIPWj]),
                         expression(hat(tau)[DMLj])),
       fill = c(NA, "#1CA9C9","#FDAE61"),
       #fill = c(NA, "#66C2A5","#FC8D62"),
       #fill = c(NA, "#20B2AA","#FFB347"),
       #fill = c(NA, "#0072B2", "#E69F00"),
       #fill = c(NA, "steelblue1", "indianred1"),
       pch = c(16, rep(NA, 3)),
       #col = c("black", "steelblue1", "indianred1"),
       #col = c("black", "#0072B2", "#E69F00"),
       #col = c("black", "#20B2AA","#FFB347"),
       #col = c("black", "#66C2A5","#FC8D62"),
       col = c("black", "#1CA9C9","#FDAE61"),
       border = NA, bty = "n", cex = 2, horiz = "TRUE")


