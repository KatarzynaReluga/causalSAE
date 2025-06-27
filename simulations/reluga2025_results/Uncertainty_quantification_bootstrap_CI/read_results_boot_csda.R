
setwd("./nb2_c_csda")
setwd("./nb1d_nc_csda_bg")


library("RColorBrewer")
library(data.table)

####################
### Load results ###
####################
file_list <- list.files()
SimNum = 1000

tau_true <- matrix(0, nrow = SimNum, ncol = 41)

tau_boot_025_AIPWM <- matrix(0, nrow = SimNum, ncol = 41) 
tau_boot_975_AIPWM <- matrix(0, nrow = SimNum, ncol = 41) 

tau_boot_025_AIPWbM <- matrix(0, nrow = SimNum, ncol = 41) 
tau_boot_975_AIPWbM <- matrix(0, nrow = SimNum, ncol = 41) 

tau_boot_025_AIPWb2M <- matrix(0, nrow = SimNum, ncol = 41) 
tau_boot_975_AIPWb2M <- matrix(0, nrow = SimNum, ncol = 41) 

covBper_AIPWM <- matrix(0, nrow = SimNum, ncol = 41)
covBper_AIPWbM <- matrix(0, nrow = SimNum, ncol = 41)
covBper_AIPWb2M <- matrix(0, nrow = SimNum, ncol = 41)

for (i in 1:SimNum) {
  
  load(file_list[i])
  tau_true[i, ] <- Results$tau_true
  
  tau_boot_AIPW <- Results$tau_boot_dif_AIPW + matrix(rep(Results$tau_hat_AIPW, 1000), nrow = 1000, ncol = 41, byrow = TRUE)
  

  # CI with percentiles
  tau_boot_025_AIPWM[i, ] <- apply(tau_boot_AIPW, 2,  quantile, 0.025)
  tau_boot_975_AIPWM[i, ] <- apply(tau_boot_AIPW, 2,  quantile, 0.975)
  covBper_AIPWM[i, ] <- (tau_boot_025_AIPWM[i, ] <= tau_true[i, ] ) & (tau_boot_975_AIPWM[i,] >= tau_true[i, ] )
#  lenBper_AIPW <- tau_boot_975_AIPW - tau_boot_025_AIPW
  
  
  tau_boot_025_AIPWbM[i, ]  <- apply(tau_boot_AIPW - Results$Bias_boot_AIPW, 2,  quantile, 0.025)
  tau_boot_975_AIPWbM[i, ]  <- apply(tau_boot_AIPW - Results$Bias_boot_AIPW, 2,  quantile, 0.975)
  covBper_AIPWbM[i, ] <- (tau_boot_025_AIPWbM[i, ] <= tau_true[i, ] ) & (tau_boot_975_AIPWbM[i,] >= tau_true[i, ] )
  #  lenBper_AIPWb <- tau_boot_975_AIPWb - tau_boot_025_AIPWb
  
  tau_boot_025_AIPWb2M[i, ] <- apply(tau_boot_AIPW - Results$Bias_boot_AIPWbc, 2,  quantile, 0.025)
  tau_boot_975_AIPWb2M[i, ] <- apply(tau_boot_AIPW - Results$Bias_boot_AIPWbc, 2,  quantile, 0.975)
  covBper_AIPWb2M <- (tau_boot_025_AIPWb2M[i, ]  <= tau_true ) & ( tau_boot_975_AIPWb2M[i, ] >= tau_true )
#  lenBper_AIPWb2 <- tau_boot_975_AIPWb2 - tau_boot_025_AIPWb2
  
  
}

mean(colMeans(covBper_AIPWM))
mean(colMeans(covBper_AIPWbM))
mean(colMeans(covBper_AIPWb2M))




covBM_NIPW <- matrix(0, nrow = SimNum, ncol = 41)
covBperM_NIPW <- matrix(0, nrow = SimNum, ncol = 41)
covBM_NIPW2 <- matrix(0, nrow = SimNum, ncol = 41)
covB_NIPWbcM <- matrix(0, nrow = SimNum, ncol = 41)
covB_NIPWbM <- matrix(0, nrow = SimNum, ncol = 41)
covB_NIPWb2M <- matrix(0, nrow = SimNum, ncol = 41)
covB_NIPWbbcM <- matrix(0, nrow = SimNum, ncol = 41)
covB_NIPWb2bcM <- matrix(0, nrow = SimNum, ncol = 41)
covBper_NIPWbM <- matrix(0, nrow = SimNum, ncol = 41)
covBper_NIPWb2M <- matrix(0, nrow = SimNum, ncol = 41)

tau_boot_025_NIPWb2M <- matrix(0, nrow = SimNum, ncol = 41) 
tau_boot_975_NIPWb2M <- matrix(0, nrow = SimNum, ncol = 41) 

covBM_AIPW <- matrix(0, nrow = SimNum, ncol = 41)
covBperM_AIPW <- matrix(0, nrow = SimNum, ncol = 41)
covBM_AIPW2 <- matrix(0, nrow = SimNum, ncol = 41)
covB_AIPWbcM <- matrix(0, nrow = SimNum, ncol = 41)
covB_AIPWbM <- matrix(0, nrow = SimNum, ncol = 41)
covB_AIPWb2M <- matrix(0, nrow = SimNum, ncol = 41)
covB_AIPWbbcM <- matrix(0, nrow = SimNum, ncol = 41)
covB_AIPWb2bcM <- matrix(0, nrow = SimNum, ncol = 41)
covBper_AIPWbM <- matrix(0, nrow = SimNum, ncol = 41)
covBper_AIPWb2M <- matrix(0, nrow = SimNum, ncol = 41)

tau_boot_025_AIPWb2M <- matrix(0, nrow = SimNum, ncol = 41) 
tau_boot_975_AIPWb2M <- matrix(0, nrow = SimNum, ncol = 41) 

tau_boot_025_NIPWbM <- matrix(0, nrow = SimNum, ncol = 41) 
tau_boot_975_NIPWbM <- matrix(0, nrow = SimNum, ncol = 41) 

tau_boot_025_AIPWbM <- matrix(0, nrow = SimNum, ncol = 41) 
tau_boot_975_AIPWbM <- matrix(0, nrow = SimNum, ncol = 41) 

tau_true <- matrix(0, nrow = SimNum, ncol = 41)

tau_hat_AIPW <- matrix(0, nrow = SimNum, ncol = 41)
tau_hat_OR <- matrix(0, nrow = SimNum, ncol = 41)
tau_hat_NIPW <- matrix(0, nrow = SimNum, ncol = 41)

MSE_OR <- matrix(0, nrow = SimNum, ncol = 41)
MSE_AIPW <- matrix(0, nrow = SimNum, ncol = 41)
MSE_NIPW <- matrix(0, nrow = SimNum, ncol = 41)


tau_hat_OR_nb <- matrix(0, nrow = SimNum, ncol = 41)
tau_hat_AIPW_nb <- matrix(0, nrow = SimNum, ncol = 41)
tau_hat_NIPW_nb <- matrix(0, nrow = SimNum, ncol = 41)

Bias_OR <- matrix(0, nrow = SimNum, ncol = 41)
Bias_AIPW <- matrix(0, nrow = SimNum, ncol = 41)
Bias_NIPW <- matrix(0, nrow = SimNum, ncol = 41)



b_OR5 <- matrix(0, nrow = length(file_list), ncol = 512)
b_AIPW5 <- matrix(0, nrow = length(file_list), ncol = 512)
b_NIPW5 <- matrix(0, nrow = length(file_list), ncol = 512)

error_csda <- numeric(SimNum)

#b_OR30 <- matrix(0, nrow = length(file_list), ncol = 512)
#b_AIPW30 <- matrix(0, nrow = length(file_list), ncol = 512)
#b_NIPW30 <- matrix(0, nrow = length(file_list), ncol = 512)

qlAIPW <- matrix(0, nrow = SimNum, ncol = 41)
quAIPW <- matrix(0, nrow = SimNum, ncol = 41)
perAIPW <- matrix(0, nrow = SimNum, ncol = 41)

qlNIPW <- matrix(0, nrow = SimNum, ncol = 41)
quNIPW <- matrix(0, nrow = SimNum, ncol = 41)
perNIPW <- matrix(0, nrow = SimNum, ncol = 41)

for (i in 1:SimNum) {
  
  load(file_list[i])
  
  error_csda[i] <- Results$error_csda
  #tau_true[i, ] <- taus[[1]]$tau
  tau_true[i, ] <- Results$tau_true
  
  tau_hat_AIPW[i, ] <- Results$tau_hat_AIPW
#  tau_hat_OR[i, ] <- Results$tau_hat_OR
  tau_hat_NIPW[i, ] <- Results$tau_hat_NIPW
  
  covBM_AIPW[i, ] <- Results$covB_AIPW
  covBM_NIPW[i, ] <- Results$covB_NIPW
#  covBM_OR[i, ] <- Results$covB_OR
  
  covBperM_AIPW[i, ] <- Results$covBper_AIPW
  covBperM_NIPW[i, ] <- Results$covBper_NIPW
#  covBperM_OR[i, ] <- Results$covBper_OR
  
 # MSE_OR[i, ] <- Results$MSE_boot_OR
  MSE_NIPW[i, ] <- Results$MSE_boot_NIPW
  MSE_AIPW[i, ] <- Results$MSE_boot_AIPW
  
#  Bias_OR[i, ]  <- colMeans(Results$tau_boot_dif_OR)
  Bias_AIPW[i, ]  <- colMeans(Results$tau_boot_dif_AIPW)
  Bias_NIPW[i, ]  <- colMeans(Results$tau_boot_dif_NIPW)
  
#  tau_hat_OR_nb[i, ] <- - Bias_OR[i, ] + Results$tau_hat_OR
  tau_hat_AIPW_nb[i, ] <- - Bias_AIPW[i, ] + Results$tau_hat_AIPW
  tau_hat_NIPW_nb[i, ] <- - Bias_NIPW[i, ] + Results$tau_hat_NIPW
  
#  upInt_NIPW <- tau_hat_NIPW + sqrt(MSE_boot_NIPW) * 1.96
#  loInt_NIPW <- tau_hat_NIPW - sqrt(MSE_boot_NIPW) * 1.96
#  covB_NIPW <- (loInt_NIPW <= tau_true ) & (upInt_NIPW >= tau_true )
  
  
  covBM_AIPW2[i, ] <- ((tau_hat_AIPW_nb[i, ] - sqrt(MSE_AIPW[i, ]) * 1.96) <= tau_true[i, ]) & ((tau_hat_AIPW_nb[i, ] + sqrt(MSE_AIPW[i, ]) * 1.96) >= tau_true[i, ]) 
  covBM_NIPW2[i, ] <- ((tau_hat_NIPW_nb[i, ] - sqrt(MSE_NIPW[i, ]) * 1.96) <= tau_true[i, ]) & ((tau_hat_NIPW_nb[i, ] + sqrt(MSE_NIPW[i, ]) * 1.96) >= tau_true[i, ]) 
#  covBM_OR2[i, ] <-   ((tau_hat_OR_nb[i, ] - sqrt(MSE_OR[i, ]) * 1.96) <= tau_true[i, ]) & ((tau_hat_OR_nb[i, ] + sqrt(MSE_OR[i, ]) * 1.96) >= tau_true[i, ]) 
  
  
  covB_AIPWbcM[i, ] <- Results$covB_AIPWbc
  covB_AIPWbM[i, ] <- Results$covB_AIPWb
  covB_AIPWb2M[i, ] <- Results$covB_AIPWb2
  covB_AIPWbbcM[i, ] <- Results$covB_AIPWbbc
  covB_AIPWb2bcM[i, ] <- Results$covB_AIPWb2bc
  covBper_AIPWbM[i, ] <- Results$covBper_AIPWb
  covBper_AIPWb2M[i, ] <- Results$covBper_AIPWb2
  
  tau_boot_025_NIPWb2M[i, ] <- Results$tau_boot_025_NIPWb2
  tau_boot_975_NIPWb2M[i, ] <- Results$tau_boot_975_NIPWb2
 
   
  covB_NIPWbcM[i, ] <- Results$covB_NIPWbc
  covB_NIPWbM[i, ] <- Results$covB_NIPWb
  covB_NIPWb2M[i, ] <- Results$covB_NIPWb2
  covB_NIPWbbcM[i, ] <- Results$covB_NIPWbbc
  covB_NIPWb2bcM[i, ] <- Results$covB_NIPWb2bc
  covBper_NIPWbM[i, ] <- Results$covBper_NIPWb
  covBper_NIPWb2M[i, ] <- Results$covBper_NIPWb2
  
  tau_boot_025_AIPWb2M[i, ] <- Results$tau_boot_025_AIPWb2
  tau_boot_975_AIPWb2M[i, ] <- Results$tau_boot_975_AIPWb2
  
  
  tau_boot_025_NIPWbM[i, ] <- Results$tau_boot_025_NIPWb
  tau_boot_975_NIPWbM[i, ] <- Results$tau_boot_975_NIPWb
  
  tau_boot_025_AIPWbM[i, ] <- Results$tau_boot_025_AIPWb
  tau_boot_975_AIPWbM[i, ] <- Results$tau_boot_975_AIPWb
  # Density
#  b_OR0 <- density(Results$tau_boot_dif_OR[, 5])
#  b_OR5[i, ] <- b_OR0$y
  
  b_AIPW0 <- density(Results$tau_boot_dif_AIPW[, 5])
  b_AIPW5[i, ] <- b_AIPW0$y
  
  b_NIPW0 <- density(Results$tau_boot_dif_NIPW[, 5])
  b_NIPW5[i, ] <- b_NIPW0$y
  
  
#  b_OR0 <- density(Results$tau_boot_dif_OR[, 30])
#  b_OR30[i, ] <- b_OR0$y
  
#  b_AIPW0 <- density(Results$tau_boot_dif_AIPW[, 30])
#  b_AIPW30[i, ] <- b_AIPW0$y
  
#  b_NIPW0 <- density(Results$tau_boot_dif_NIPW[, 30])
#  b_NIPW30[i, ] <- b_NIPW0$y
  
  qlAIPW[i, ] <- apply(Results$tau_boot_dif_AIPW, 2,  quantile, 0.025)
  quAIPW[i, ] <- apply(Results$tau_boot_dif_AIPW, 2,  quantile, 0.975) 
  
  perAIPW[i, ] <- (tau_true[i, ] >= qlAIPW[i, ]) & (tau_true[i, ] <= quAIPW[i, ])
  
  qlNIPW[i, ] <- apply(Results$tau_boot_dif_NIPW, 2,  quantile, 0.025)
  quNIPW[i, ] <- apply(Results$tau_boot_dif_NIPW, 2,  quantile, 0.975) 
  
  perNIPW[i, ] <- (tau_true[i, ] >= qlNIPW[i, ]) & (tau_true[i, ] <= quNIPW[i, ])
}

nw <- which(error_csda == 0)

mean(colMeans(perAIPW))
mean(colMeans(perNIPW))

###################
# Plots density OR
###################
par(mfrow=c(1, 1))

par(mar = c(4, 4, 1, 1) + 0.5)
plot(density(tau_hat_OR[, 5] - tau_true[, 5]), lwd = 5, col = 2, ylim = c(0, 2), 
     xlim = c(-1, 1),
     lty = 1, main = "", xlab = "", cex.axis = 1.3, cex.lab = 1.5)

for (i in c(21:50)) {
  lines(b_OR0$x, b_OR5[i, ], lwd = 2, col = 3, lty = "dashed") # 13, 16, 17, 21, 25
  
}

legend("topleft", legend = c(expression(hat(tau)-tau),expression(hat(tau)^b-hat(tau))),  
       fill = c(2, 3), bty = "n", cex = 1.4)
##############################
BiasOR = colMeans(tau_hat_OR) - colMeans(tau_true)

par(mfrow=c(2,3))

plot(1:41, tau_true[1, ], lwd = 2, type = "l", ylim = c(0, 3), ylab = "tau")
for (i in seq(1, 1000, 100)) {
  lines(tau_hat_OR[i, ], lwd = 2, col = 3, lty = "dashed") # 13, 16, 17, 21, 25
  
}

legend("topright", legend = c("true", "OR"),  fill = c(1, 3), bty = "n", cex = 1.4)

boots300 <- Results$tau_boot_dif_OR +  matrix(rep( tau_hat_OR[300, ], 1000), nrow = 1000, ncol = 41, byrow = TRUE)
plot(1:41, tau_hat_OR[300, ], lwd = 2, type = "l", ylim = c(0, 3))
for (i in seq(1, 300, 20)) {
  lines(boots300[i, ], lwd = 2, col = 3, lty = "dashed") # 13, 16, 17, 21, 25
  
}

legend("topright", legend = c("OR", "OR boot"),  fill = c(1, 3), bty = "n", cex = 1.4)
###############################
###############################
# Plots denisty AIPW
###############################
##############################

par(mar = c(4, 4, 1, 1) + 0.5)
plot(density(tau_hat_AIPW[, 5] - tau_true[, 5]), lwd = 5, col = 2, ylim = c(0, 4), 
     xlim = c(-1, 1),
     lty = 1, main = "", xlab = "", cex.axis = 1.3, cex.lab = 1.5)

for (i in seq(1, 100, 10)) {
  lines(b_AIPW0$x, b_AIPW5[i, ], lwd = 2, col = 3, lty = "dashed") # 13, 16, 17, 21, 25
  
}

legend("topleft", legend = c(expression(hat(tau)-tau),expression(hat(tau)^b-hat(tau))),  
       fill = c(2, 3), bty = "n", cex = 1.4)

###################
plot(1:41, tau_true[1, ], lwd = 2, type = "l", ylim = c(0, 3), ylab = "tau")
for (i in nw[20:40]) {
  lines(tau_hat_AIPW[i, ], lwd = 2, col = 3, lty = "dashed") # 13, 16, 17, 21, 25
  
}
legend("topright", legend = c("true", "AIPW"),  fill = c(1, 3), bty = "n", cex = 1.4)


BiasAIPW = colMeans(tau_hat_AIPW) - colMeans(tau_true)
  
boots300AIPW <- Results$tau_boot_dif_AIPW +  matrix(rep( tau_hat_AIPW[300, ], 1000), nrow = 1000, ncol = 41, byrow = TRUE)
plot(1:41, tau_hat_AIPW[300, ], lwd = 2, type = "l", ylim = c(0, 3))
for (i in seq(1, 300, 20)) {
  lines(boots300AIPW[i, ], lwd = 2, col = 3, lty = "dashed") # 13, 16, 17, 21, 25
  
}
legend("topright", legend = c("AIPW", "AIPW boot"),  fill = c(1, 3), bty = "n", cex = 1.4)


# Plots NIPW
####################
par(mar = c(4, 4, 1, 1) + 0.5)
plot(density(tau_hat_NIPW[, 5] - tau_true[, 5]), lwd = 5, col = 2, ylim = c(0, 4), 
     xlim = c(-1, 1),
     lty = 1, main = "", xlab = "", cex.axis = 1.3, cex.lab = 1.5)

for (i in seq(1, 100, 10)) {
  lines(b_NIPW0$x, b_NIPW5[i, ], lwd = 2, col = 3, lty = "dashed") # 13, 16, 17, 21, 25
  
}


legend("topleft", legend = c(expression(hat(tau)-tau),expression(hat(tau)^b-hat(tau))),  
       fill = c(2, 3), bty = "n", cex = 1.4)

##########################
plot(1:41, tau_true[1, ], lwd = 2, type = "l", ylim = c(0, 3), ylab = "tau")
for (i in seq(1, 1000, 50)) {
  lines(tau_hat_NIPW[i, ], lwd = 2, col = 3, lty = "dashed") # 13, 16, 17, 21, 25
  
}
legend("topright", legend = c("true", "NIPW"),  fill = c(1, 3), bty = "n", cex = 1.4)


BiasNIPW = colMeans(tau_hat_NIPW) - colMeans(tau_true)

boots300NIPW <- Results$tau_boot_dif_NIPW +  matrix(rep( tau_hat_NIPW[300, ], 1000), nrow = 1000, ncol = 41, byrow = TRUE)
plot(1:41, tau_hat_NIPW[300, ], lwd = 2, type = "l", ylim = c(0, 3))
for (i in seq(1, 300, 20)) {
  lines(boots300NIPW[i, ], lwd = 2, col = 3, lty = "dashed") # 13, 16, 17, 21, 25
  
}
legend("topright", legend = c("NIPW", "NIPW boot"),  fill = c(1, 3), bty = "n", cex = 1.4)


#Percenitle does not work -- bias 

true_MSE_OR <- colMeans((tau_hat_OR - tau_true)^2)
true_MSE_AIPW <- colMeans((tau_hat_AIPW - tau_true)^2)
true_MSE_NIPW <- colMeans((tau_hat_NIPW - tau_true)^2)
plot(true_MSE_AIPW, type = "l", ylim = c(0, 1.2), ylab = "MSE")
lines(true_MSE_AIPW, type = "l", col = 2, lwd = 2)
#lines(MSE_boot_AIPW, col = "orange", lwd = 2)
lines(colMeans(MSE_AIPW), col = "orange", lwd = 2)

par(mfrow=c(1,1))
plot(true_MSE_OR, type = "l", ylim = c(0, 1.2), ylab = "MSE", col = 1, lwd = 2)
lines(MSE_OR[1, ], col = "grey", lwd = 2)
lines(colMeans(MSE_OR), col = "grey", lwd = 2)
plot(true_MSE_AIPW, type = "l", col = 2, lwd = 2)
lines(MSE_AIPW[1, ], col = "orange", lwd = 2)
lines(colMeans(MSE_AIPW), col = "orange", lwd = 2)
lines(true_MSE_NIPW, type = "l", col = "blue", lwd = 2)
lines(MSE_NIPW[1, ], col = "lightblue", lwd = 2)
lines(colMeans(MSE_NIPW), col = "lightblue", lwd = 2)
legend("top", legend = c(expression("OR MSE"[true]),expression("MSE"[OR]),
                              expression("AIPW MSE"[true]), expression("MSE"[AIPW]), 
                              expression("NIPW MSE"[true]), expression("MSE"[NIPW])),  
       fill = c(1, "grey", 2, "orange", "blue", "lightblue"), bty = "n", cex = 1.4)


lines(colMeans(MSE_AIPW), col = 4)

plot(density(tau_hat_NIPW[, 4] - tau_true[, 4]), col = 4)
lines(density(Results$tau_boot_dif_NIPW[, 4]), lty = 2)

true_MSE_NIPW <- colMeans((tau_hat_NIPW - tau_true)^2)
lines(true_MSE_NIPW, type = "l", col = 3)
lines(colMeans(MSE_NIPW), col = 4)

legend("top", legend = c(expression("MSE"[true]),expression("MSE"[OR]), expression("MSE"[AIPW]), expression("MSE"[NIPW])),  
       fill = c(1, 2, 3, 4), bty = "n", cex = 1.4)


mean(colMeans(covBM_AIPW[,]))
mean(colMeans(covBM_AIPW2[, ]))
mean(colMeans(covBperM_AIPW))
mean(colMeans(covB_AIPWbcM))
mean(colMeans(covB_AIPWbM))
mean(colMeans(covB_AIPWb2M))
mean(colMeans(covB_AIPWbbcM))
mean(colMeans(covB_AIPWb2bcM))
mean(colMeans(covBper_AIPWbM))




mean(colMeans(covBM_NIPW[, ]))
mean(colMeans(covBM_NIPW2[, ]))
mean(colMeans(covBperM_NIPW[, ]))
mean(colMeans(covB_NIPWbcM))
mean(colMeans(covB_NIPWbM))
mean(colMeans(covB_NIPWb2M))
mean(colMeans(covB_NIPWbbcM))
mean(colMeans(covB_NIPWb2bcM))
mean(colMeans(covBper_NIPWbM))

mean(colMeans(covBper_NIPWb2M))
mean(colMeans(covBper_AIPWb2M))

tau_boot_025_AIPWb2M[i, ] <- Results$tau_boot_025_AIPWb2
tau_boot_975_AIPWb2M[i, ] <- Results$tau_boot_975_AIPWb2
par(mgp = c(2, 1, 0)) 
plot(1:41, tau_true[2, ], lwd = 2, type = "l", 
     ylim = c(-0.5, 3.5), 
     xlab = "Subpopulation", 
     ylab = expression(tau[j]), 
     cex.lab = 1.5,    # Enlarges the axis labels
     cex.axis = 1.2, 
     xlim = c(1,41)) 
for (i in  seq(1, 200, 20)) {
  #  lines(tau_hat_AIPW[i, ], lwd = 2, col = "lightgreen") # 13, 16, 17, 21, 25
  #lines(up_AIPW[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25
  #lines(do_AIPW[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25
  lines(tau_boot_025_AIPWb2M[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25
  lines(tau_boot_975_AIPWb2M[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25
  lines(tau_boot_025_NIPWb2M[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25
  lines(tau_boot_975_NIPWb2M[i, ], lwd = 2, col = "lightblue", lty = "dashed") # 13, 16, 17, 21, 25
  
}
legend("top", legend = c(expression(tau[j]), 
                               expression(hat(tau)[j])),  
       fill = c(1, "lightblue"), bty = "n", cex = 1.4, horiz = "TRUE")





tau_trueM <- colMeans(tau_true)
tau_hatM <- colMeans(tau_hat_AIPW)
tau_hatMB <- - colMeans(Results$tau_boot_dif_AIPW) + tau_hat_AIPW[i,]

plot(1:41, tau_trueM, type = "l")
lines(tau_hatM, col = 2)
lines(tau_hatMB, col = 3)
# Join results and compute
dt_res <- data.table(do.call(rbind, res))
dt_Mean <- dt_res[, lapply(.SD, mean), by = "Method"]

