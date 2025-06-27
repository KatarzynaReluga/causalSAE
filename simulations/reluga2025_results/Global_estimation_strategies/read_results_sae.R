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


##################################################################################
setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/ModelBased/DBscenarios/Nonparam_PS/Old simulations/OR_NC_25")
####################
### Load results ###
####################
file_list <- list.files()
#SimNum = length(file_list)
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
dim(dt_Mean)
tau_trueM <- dt_Mean[c(4:15), c(77:101)]
Methods <- dt_res$Method[4:15]


# Bias -----------------------------------------------
B <- dt_Mean[c(4:15), c(27:51)]
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
MSE <- dt_Mean[c(4:15), c(52:76)]
MSEMean <- rowMeans(MSE)
MSEMedian <- apply(MSE, 1, quantile, probs = 0.5)




resultsOR_25 <- data.frame(Method = Methods,
                        
                        mse = MSEMean,
                        Bias = BMean,
                        ABias = ABMean,
                        RABias = RBMean)

write.csv(resultsOR_25,
          file = 'resultsOR_25.csv', row.names = FALSE)


resultsF_25 <- rbind(resultsOR_25, results_25)
resultsF_25ord <- resultsF_25[order(resultsF_25$mse, decreasing = FALSE), ]  
write.csv(resultsF_25ord,
          file = 'resultsF_25ord.csv', row.names = FALSE)


#
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
  lines(unlist(TMLE_df[i, ]), lwd = 2, col = "red", lty = "dashed") # 13, 16, 17, 21, 25
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
  TMLE_i <- unlist(TMLE_df[i, ])
  
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
                    unlist(TMLE_df[1, ix]), 
                    unlist(TMLE_df[101, ix]),
                    unlist(TMLE_df[201, ix]))  # dimensions: 6 x 41

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








lines(Results$MXcsda,  col = 2, lwd = 3, type = "b")
lines(Results$MXMch, col = 3, lwd = 3, type = "b")
lines(Results$EXMch, col = 4, lwd = 3, type = "b")
lines(Results$RcXMch, col = 5, lwd = 3, type = "b")
lines(Results6$MRctM, col = 6, lwd = 3, type = "b")

lines(Results15$XtXtM, col = 2, lwd = 3, type = "b")
lines(Results15$MSM, col = 3, lwd = 3, type = "b")
lines(Results15$MRtM, col = 4, lwd = 3, type = "b")
lines(Results15$RctXM, col = 5, lwd = 3, type = "b")
lines(Results15$MRctM, col = 6, lwd = 3, type = "b")

lines(Results2$XtXtM, col = 2, lwd = 3, type = "b")
lines(Results2$MSM, col = 3, lwd = 3, type = "b")
lines(Results2$MRtM, col = 4, lwd = 3, type = "b")
lines(Results2$RctXM, col = 5, lwd = 3, type = "b")
lines(Results2$MRctM, col = 6, lwd = 3, type = "b")


legend("bottomleft", legend = c("tau", "XtXtM", "MSM", "MRtM", "RctXM", "MRctM"),
       lwd = 3, col = c(1, 2, 3, 4, 5, 6), 
       bty = "n", cex = 1.9, horiz = TRUE)

plot(x, y, main = "My title", sub = "Subtitle",
     cex.main = 2,   # Title size
     cex.sub = 1.5,  # Subtitle size
     cex.lab = 3,    # X-axis and Y-axis labels size
     cex.axis = 0.5) # Axis labels size


i = 6
load(file_list[i])
ResultsM <- Results[- c(1:9)]
Results6 <- ResultsM
i = 15
load(file_list[i])
ResultsM <- Results[- c(1:9)]
Results15 <- ResultsM
i = 2
load(file_list[i])
ResultsM <- Results[- c(1:9)]
Results2 <- ResultsM

lines(ResultsM$XtXtM, col = 2)
lines(ResultsM$MSM, col = 3)
lines(ResultsM$MRtM, col = 4)
lines(ResultsM$RctXM, col = 5)
lines(ResultsM$MRctM, col = 6)

MSE <- c(unlist(tau_trueM[1,]), unlist(Results6$XtXtM), unlist(Results6$MSM), unlist(Results6$MRtM), 
         unlist(Results6$RctXM), unlist(Results6$MRctM),
         
         unlist(Results15$XtXtM), unlist(Results15$MSM), unlist(Results15$MRtM), 
          unlist(Results15$RctXM), unlist(Results15$MRctM), 
         
         unlist(Results2$XtXtM), unlist(Results2$MSM), unlist(Results2$MRtM), 
         unlist(Results2$RctXM), unlist(Results2$MRctM))

Method <- c(rep("tau", 25), 
            rep("XtXtM", 25), rep("MSM", 25), rep("MRtM", 25),
            rep("RctXM", 25),  rep("MRctM", 25), 
            
            rep("XtXtM", 25), rep("MSM", 25), rep("MRtM", 25),
            rep("RctXM", 25),  rep("MRctM", 25), 
            
            rep("XtXtM", 25), rep("MSM", 25), rep("MRtM", 25),
            rep("RctXM", 25),  rep("MRctM", 25))

Subpopulation <- c(1:25, 
                   1:25, 1:25, 1:25, 1:25, 1:25, 
                   1:25, 1:25, 1:25, 1:25, 1:25, 
                   1:25, 1:25, 1:25, 1:25, 1:25)

plot25 <- data.frame(MSE = MSE, 
                     Method = Method, 
                     Subpopulation = Subpopulation)


mse25 <- ggplot(plot25, aes(x = Subpopulation,
                          y = MSE,
                          col = Method, palette = "Spectral")) +
  geom_point(size = 3) +
#  geom_line() +
  labs(y = "MSE", x = "Subpopulation")+ 
  theme_bw() +
#  theme(axis.text.x=element_blank(), #remove x axis labels
#        axis.ticks.x=element_blank(), 
#        panel.grid.major.x = element_blank()#remove x axis ticks
#  ) + 
  coord_cartesian(ylim=c(0, 25)) +
  scale_y_continuous(breaks = seq(0, 19, 1)) +
  scale_x_discrete(limits = NULL) +
  theme(legend.direction = "horizontal",
        legend.position = c(0.5, 0.15),
        legend.box = "horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold")) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=20,face="bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 22)) +
  theme(axis.title.x = element_text(size = 24)) +
  guides(col = guide_legend(nrow = 3, byrow = T, 
                            title.position = "top",
                            title.hjust = 0.5, order = 1)) 
mse25



# Analyse results

MPlus25F <- read.csv("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/ModelBased/DBscenarios/Nonparam_PS/results_25F.csv")
#MPlus25F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/Simulations/ModelBased/DBscenarios/resultsDB25F.csv')
Methods <- MPlus25F$Method
ARBPMedian25F <- MPlus25F$ABias
#ARBPMedian25F <- MPlus25F$ARBPMean
ARBPMediansort25F <- sort(ARBPMedian25F, index.return = TRUE)$x
ARBPMediansorti25F <- sort(ARBPMedian25F, index.return = TRUE)$ix
Methods[ARBPMediansorti25F][1:10]
plot(density(ARBPMediansort25F[1:1000]))
plot(1:1464, ARBPMedian25F,  ylim = c(1.2,1.70))
#plot(1:127, ARBPMedian25F, ylim = c(7,35))
abline(a = 10, b = 0)
which(ARBPMedian25F <= 11)
which(ARBPMedian25F <= 7.5)
Methods[which(ARBPMedian25F <=11)]
#DB
# 1   2   3   4   6   9  15  19  21  24  27  30  46  49  52  55  58  61  64  67  69  77
# 81  87  90 104 107 110 113 116

which(ARBPMedian25F <= 7.5)
# 2   6   9  52  55  77  81 110 113
# 2   6  55  77 110
RMSEPMedian25F <- MPlus25F$mse
#RMSEPMedian25F <- MPlus25F$RMSEMean
RMSEPMediansort25F <- sort(RMSEPMedian25F, index.return = TRUE)$x
RMSEPMediansorti25F <- sort(RMSEPMedian25F, index.return = TRUE)$ix
Methods[RMSEPMediansorti25F][1:20]
plot(density(RMSEPMedian25F))
plot(1:1464, RMSEPMedian25F, ylim = c(16, 17))
plot(1:126, RMSEPMedian25F[-127])
plot(1:126, RMSEPMedian25F[-127], ylim = c(30, 90))
abline(35.5,0)
which(RMSEPMedian25F <= 17 & RMSEPMedian25F >= 16)
Methods[which(RMSEPMedian25F <= 16.4 & RMSEPMedian25F >= 16.1)]
#DB
# 2   3   6   9  15  19  21  24  27  30  46  49  52  55  58  61  64  67  69  77  81  87
# 90 104 107 110 113 116 119 122 125
