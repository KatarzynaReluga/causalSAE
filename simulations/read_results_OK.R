#setwd("C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/Mplustune25_2")
#setwd("C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/PositiveEffects/Tuned/Mplustune100_2")
setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/dbXtune")
setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/DBscenarios/DB50")


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

#  colnames(all_results[[i]]) <- paste("C", 1:dim(all_results[[i]])[2], sep = "")

}

 # Join results and compute
dt_res <- data.table(do.call(rbind, res))
#E_OR <- dt_res[dt_res$Method == "E_OR"]
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

# Bias ----------------------------------------------
RB <- dt_Mean[c(4:72), c(27:51)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(27:51)]/abs(tau_trueM) * 100

RB <- dt_Mean[c(4:72), c(52:101)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(52:101)]/abs(tau_trueM) * 100

RB <- dt_Mean[c(4:72), c(102:201)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(102:201)]/abs(tau_trueM) * 100

# Design based
RB <- dt_Mean[c(4:72), c(43:83)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(43:83)]/abs(tau_trueM) * 100


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


RMSEMean <- rowMeans(RMSEP)
RMSEPMedian <- apply(RMSEP, 1, quantile, probs = 0.5)


resultsDB50 <- data.frame(Method = Methods,

           RBPMean = RBPMean,
           RBPMedian = RBPMedian,

           ARBPMean = ARBPMean,
           ARBPMedian = ARBPMedian,

           RMSEMean = RMSEMean,
           RMSEPMedian = RMSEPMedian)

write.csv(resultsDB50,
          file = 'resultsDB50.csv', row.names = FALSE)

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
which(ARBPMedian25F <=21.57)
Methods[which(ARBPMedian25F <=21.57)]
# 1   3   6  21  27  30  32  46  49  52  55  58  90 104 127


# 2   3   6   9  15  19  21  24  27  30  46  49  52  55  58  64
# 67  69  77  81  87  90 104 107 110 113 116
RMSEPMedian25F <- MPlus25F$RMSEPMedian
RMSEPMediansort25F <- sort(RMSEPMedian25F, index.return = TRUE)$x
RMSEPMediansorti25F <- sort(RMSEPMedian25F, index.return = TRUE)$ix
Methods[RMSEPMediansorti25F][1:14]
plot(density(RMSEPMedian25F))
plot(1:127, RMSEPMedian25F)
plot(1:126, RMSEPMedian25F[-127])
plot(1:126, RMSEPMedian25F[-127], ylim = c(0, 36))
Methods[which(RMSEPMedian25F <= 36)]
#2   3   6   9  15  19  21  24  27  30  31  43  44  46  49  52
#55  56  58  64  67  69  77  81  87  90 104 107 110 113 116

#[1]   2   3   6   9  15  19  21  24  27  30  46  49  52  55
# 58  69  77  81  87  90 104 107 110 113 116

#2   3   6   9  15  19  21  24  27  30  46  49  52  55
# 58  69 77  87  90 104 107 110 113 116

#2  3  6  9 15 19 21 24 27 30 46 49 52 55
# 58 64 67

MPlus50F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/PositiveEffects/resultsMPlus50F.csv')
MPlus50F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/DBscenarios/resultsDB50.csv')
Methods <- MPlus50F$Method
ARBPMedian50F <- MPlus50F$ARBPMedian
ARBPMediansort50F <- sort(ARBPMedian50F, index.return = TRUE)$x
ARBPMediansorti50F <- sort(ARBPMedian50F, index.return = TRUE)$ix
Methods[ARBPMediansorti50F]
Methods[ARBPMediansorti50F][1:10]
plot(density(ARBPMediansort50F))
plot(1:69, ARBPMedian50F)
which(ARBPMedian50F <=19)
Methods[which(ARBPMedian50F <=22)]

# 5  10  14  18  33  36  39  42  45  48  51  54  57  60  63  66
# 76  80  97 100 103 106 109 112 115

# 3 6 21 27 30 46 49 52 55 58 90 104

# 2   3   6   9  15  19  21  24  27  30  46  49  52  55  58  64
# 67  69  77  81  87  90 104 107 110 113 116

RMSEPMedian50F <- MPlus50F$RMSEPMedian
Methods <- MPlus50F$Method
RMSEPMediansort50F <- sort(RMSEPMedian50F, index.return = TRUE)$x
RMSEPMediansorti50F <- sort(RMSEPMedian50F, index.return = TRUE)$ix
Methods[RMSEPMediansorti50F][1:14]
Methods[ARBPMediansorti50F][1:10]
plot(density(RMSEPMedian50F))
plot(1:69, RMSEPMedian50F)
plot(1:126, RMSEPMedian50F[-127], ylim = c(0,30))
plot(1:126, RMSEPMedian50F[-127], ylim = c(25,35))
Methods[which(RMSEPMedian50F <= 30)]
#[1]   2   3   6   9  15  19  21  24  27  30  46  49  52  55
# 69  77 81  87  90 104 107 110 113


#2  3  6  9 15 19 21 24 27 30 46 49 52 55
# 58 64 67
MPlus100F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/PositiveEffects/resultsMPlus100F.csv')
ARBPMedian100F <- MPlus100F$ARBPMedian
ARBPMediansort100F <- sort(ARBPMedian100F, index.return = TRUE)$x
ARBPMediansorti100F <- sort(ARBPMedian100F, index.return = TRUE)$ix
Methods[ARBPMediansorti100F]
plot(density(ARBPMediansort100F))
plot(1:127, ARBPMedian100F)
which(ARBPMedian100F <= 18.63702)
Methods[which(ARBPMedian100F <= 18.63702)]
# 2   3   6   9  15  19  21  24  27  30  46  49  52  55  58  64
# 67  69  77  81  87  90 104 107 110 113 116


v100 <- which(ARBPMedian100F <= 18.63702)

RMSEPMedian100F <- MPlus100F$RMSEPMedian
RMSEPMediansort100F <- sort(RMSEPMedian100F, index.return = TRUE)$x
RMSEPMediansorti100F <- sort(RMSEPMedian100F, index.return = TRUE)$ix
Methods[RMSEPMediansorti100F]
plot(density(RMSEPMedian100F))
plot(1:127, RMSEPMedian100F)
plot(1:126, RMSEPMedian100F[-127], ylim = c(22,40))
Methods[which(RMSEPMedian100F <= 30)]
#2  3  6  9 15 19 21 24 27 30 46 49 52 55
#58 64 67
ind = c(2,3,6,9,15,   19,21, 24, 27, 30,    46, 49, 52, 55)
ind100  = c(58, 64, 67)
ind50  = c(69, 77, 81, 87, 90, 104, 107, 110, 113)
ind25  = c(58, 69, 77, 87, 90, 104, 107, 110, 113, 116)

vec <- numeric(127)
vec100 <- numeric(127)
vec50 <- numeric(127)
vec25 <- numeric(127)

v100 <- numeric(127)
v50 <- numeric(127)
v25 <- numeric(127)


vec[ind] <- RMSEPMedian100F[ind]
vec[ind] <- RMSEPMedian50F[ind]
vec[ind] <- RMSEPMedian25F[ind]

Methods[ind]

vec100[ind100] <- RMSEPMedian100F[ind100]
vec50[ind50] <- RMSEPMedian50F[ind50]
vec25[ind25] <- RMSEPMedian25F[ind25]

v100[which(ARBPMedian100F <= 18.63702)] <- ARBPMedian100F[which(ARBPMedian100F <= 18.63702)]
v50[which(ARBPMedian50F <= 19)] <- ARBPMedian50F[which(ARBPMedian50F <=19)]
v25[which(ARBPMedian25F <= 21.57)] <- ARBPMedian25F[which(ARBPMedian25F <=21.57)]

# Plots
plot(1:127, RMSEPMedian100F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "RRMSE", ylim = c(20, 40),
     xaxt="n", yaxt="n", cex.main = 2)
lines(1:127, vec, col = "blue", lwd = 4, type = "p", lty = 0)
lines(1:127, vec100, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, vec50, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, vec25, col = "darkgreen", lwd = 4, type = "p", lty = 0)

#title(ylab = "RRMSE", line = 2, cex.lab = 1.2)
title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b=c(10, 20, 25,  30, 35,  40, 60, 80, 100, 120, 140)
axis(2, at = b, labels = b, cex.axis = 1.5)

title(xlab = "Estimation methods", line = 2.1, cex.lab = 1.2)

plot(1:127, ARBPMedian100F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB", ylim = c(10, 40),
     xaxt="n", yaxt="n", cex.main = 2)
plot(1:127, ARBPMedian50F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB", ylim = c(10, 40),
     xaxt="n", yaxt="n", cex.main = 2)
plot(1:127, ARBPMedian25F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB", ylim = c(10, 40),
     xaxt="n", yaxt="n", cex.main = 2)
lines(1:127, vec, col = "blue", lwd = 4, type = "p", lty = 0)
lines(1:127, v100, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, v50, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, v25, col = "darkgreen", lwd = 4, type = "p", lty = 0)

#title(ylab = "RRMSE", line = 2, cex.lab = 1.2)
title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b = c(10, 20, 25,  30, 35,  40, 60, 80, 100, 120, 140)
b = c(10, 15, 20, 25,  30, 35, 40, 50, 60, 70, 80, 100, 120, 140)
axis(2, at=b,labels=b,cex.axis=1.5)

###############################
dbE <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsEF.csv')
Methods <- dbE$Method

ARBPMediandbE <- dbE$ARBPMedian
ARBPMediansortdbE <- sort(ARBPMediandbE, index.return = TRUE)$x
ARBPMediansortidbE <- sort(ARBPMediandbE, index.return = TRUE)$ix
Methods[ARBPMediansortidbE]
plot(density(ARBPMediansortdbE))
plot(1:127, ARBPMediandbE)
#127  70   3  69   4   1   2

ind = c(127, 70, 3, 69, 4, 1, 2)
vec <- numeric(127)

vec[ind] <- ARBPMediandbE[ind]

indE = ind
vecE <- numeric(127)
vecE[indE] <- ARBPMediandbE[indE]

plot(1:127, ARBPMediandbE,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB",
     xaxt="n", yaxt="n", cex.main = 2)
title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b=c(10, 20, 30, 40, 50,  60, 70, 80, 90, 100, 120, 140)
axis(2, at = b, labels = b, cex.axis = 1.5)


#lines(1:127, vecE, col = "darkgreen", lwd = 4, type = "p", lty = 0)
#lines(1:127, vecM, col = "darkgreen", lwd = 4, type = "p", lty = 0)
#lines(1:127, vecR, col = "darkgreen", lwd = 4, type = "p", lty = 0)
#lines(1:127, vecX, col = "darkgreen", lwd = 4, type = "p", lty = 0)
#lines(1:127, vecRtune, col = "darkgreen", lwd = 4, type = "p", lty = 0)
#lines(1:127, vecXtune, col = "darkgreen", lwd = 4, type = "p", lty = 0)

plot(1:127, ARBPMediandbE,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB", ylim = c(10, 90),
     xaxt="n", yaxt="n", cex.main = 2)

vec[ind] <- ARBPMediandbE[ind]
lines(1:127, vec, col = "blue", lwd = 4, type = "p", lty = 0)
#title(ylab = "RRMSE", line = 2, cex.lab = 1.2)
title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b=c(10, 20, 30, 40, 50,  60, 70, 80, 90, 100, 120, 140)
axis(2, at = b, labels = b, cex.axis = 1.5)

RMSEPMediandbE <- dbE$RMSEPMedian
RMSEPMediansortdbE <- sort(RMSEPMediandbE, index.return = TRUE)$x
RMSEPMediansortidbE <- sort(RMSEPMediandbE, index.return = TRUE)$ix
Methods[RMSEPMediansortidbE]
plot(density(RMSEPMediandbE))
plot(1:127, RMSEPMediandbE)
#  69   3   4   2  55    54  56   1  70  52    53
# 69   3   4   2 1 53 54 55 56 70

ind = c(3, 69)
vec <- numeric(127)

vec[ind] <- RMSEPMediandbE[ind]

indE = c(4, 2, 55, 54, 56, 1, 70, 52, 53)
vecE <- numeric(127)
vecE[indE] <- RMSEPMediandbE[indE]


plot(1:127, RMSEPMediandbXtune,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "RRMSE", ylim = c(35, 100),
     xaxt="n", yaxt="n", cex.main = 2)

lines(1:127, vecE, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, vecM, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, vecR, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, vecX, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, vecRtune, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, vecXtune, col = "darkgreen", lwd = 4, type = "p", lty = 0)

lines(1:127, vec, col = "blue", lwd = 4, type = "p", lty = 0)
#title(ylab = "RRMSE", line = 2, cex.lab = 1.2)
title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b=c(10, 20, 25,  30, 40, 50,  60, 70, 80, 90, 100, 120, 140)
axis(2, at = b, labels = b, cex.axis = 1.5)

#title(xlab = "Estimation methods", line = 2.1, cex.lab = 1.2)



dbM <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsMF.csv')

ARBPMediandbM <- dbM$ARBPMedian
ARBPMediansortdbM <- sort(ARBPMediandbM, index.return = TRUE)$x
ARBPMediansortidbM <- sort(ARBPMediandbM, index.return = TRUE)$ix
Methods[ARBPMediansortidbM]
plot(density(ARBPMediansortdbM))
plot(1:127, ARBPMediandbM)
#127  70  69   3   4   2   1

vec[ind] <- ARBPMediandbM[ind]

indM = ind
vecM <- numeric(127)
vecM[indM] <- ARBPMediandbM[indM]

RMSEPMediandbM <- dbM$RMSEPMedian
RMSEPMediansortdbM <- sort(RMSEPMediandbM, index.return = TRUE)$x
RMSEPMediansortidbM <- sort(RMSEPMediandbM, index.return = TRUE)$ix
Methods[RMSEPMediansortidbM]
plot(density(RMSEPMediandbM))
plot(1:127, RMSEPMediandbM)
#DR TMLE (double robust TMLE)
# 3  69   4   2     55 54  70  56    1  52    53 110
# 69   3   4   2 1 53 54 55 56 70

ind = c(3, 69)
vec <- numeric(127)

vec[ind] <- RMSEPMediandbM[ind]

indM = c(3, 69, 4, 2, 55, 54, 70, 56, 1, 52, 53, 110)
vecM <- numeric(127)
vecM[indM] <- RMSEPMediandbM[indM]


dbX <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsXF.csv')
ARBPMediandbX <- dbX$ARBPMedian
ARBPMediansortdbX <- sort(ARBPMediandbX, index.return = TRUE)$x
ARBPMediansortidbX <- sort(ARBPMediandbX, index.return = TRUE)$ix
Methods[ARBPMediansortidbX]
plot(density(ARBPMediansortdbX))
plot(1:127, ARBPMediandbX)
# 127  70   3  69   4   1   2  53  55  54  51  56 112  66  67
#127  70   3  69   4   1   2

vec[ind] <- ARBPMediandbX[ind]

indX = c(127, 70, 3, 69, 4, 1, 2, 53, 55, 54, 51, 56, 112, 66, 67)
vecX <- numeric(127)
vecX[indX] <- ARBPMediandbX[indX]


RMSEPMediandbX <- dbX$RMSEPMedian
RMSEPMediansortdbX <- sort(RMSEPMediandbX, index.return = TRUE)$x
RMSEPMediansortidbX <- sort(RMSEPMediandbX, index.return = TRUE)$ix
Methods[RMSEPMediansortidbX]
plot(density(RMSEPMediandbX))
plot(1:127, RMSEPMediandbX)
# 3  69   4  70   2   1


vec[ind] <- RMSEPMediandbX[ind]

indX = c(3, 69, 4, 70, 2, 1)
vecX <- numeric(127)
vecX[indX] <- RMSEPMediandbX[indX]


dbR <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsRF.csv')
Methods <- dbR$Method
ARBPMediandbR <- dbR$ARBPMedian
ARBPMediansortdbR <- sort(ARBPMediandbR, index.return = TRUE)$x
ARBPMediansortidbR <- sort(ARBPMediandbR, index.return = TRUE)$ix
Methods[ARBPMediansortidbR]
plot(density(ARBPMediansortdbR))
plot(1:127, ARBPMediandbR, ylim = c())
#127  70   3  69   4   1   2
#127  70   3  69   4   1   2


vec[ind] <- ARBPMediandbR[ind]

indR = ind
vecR <- c(127, 70, 3, 69, 4, 1, 2)
vecR[indR] <- ARBPMediandbR[indR]


RMSEPMediandbR <- dbR$RMSEPMedian
RMSEPMediansortdbR <- sort(RMSEPMediandbR, index.return = TRUE)$x
RMSEPMediansortidbR <- sort(RMSEPMediandbR, index.return = TRUE)$ix
Methods[RMSEPMediansortidbR]
plot(density(RMSEPMediandbR))
plot(1:127, RMSEPMediandbR)
# 3  69  70   4      2  54  55   1    52  56  51  53 110
# 69   3   4   2 1           53 54 55 56 70


vec[ind] <- ARBPMediandbR[ind]

indR = c(3, 69, 70, 4, 2, 54, 55, 1, 52, 56, 51, 53, 110)
vecR <- numeric(127)
vecR[indR] <- RMSEPMediandbR[indR]

#ind = c(3, 69)
#vec <- numeric(127)
vec[ind] <- RMSEPMediandbR[ind]


dbRtune <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsRtuneF.csv')
ARBPMediandbRtune <- dbRtune$ARBPMedian
ARBPMediansortdbRtune <- sort(ARBPMediandbRtune, index.return = TRUE)$x
ARBPMediansortidbRtune <- sort(ARBPMediandbRtune, index.return = TRUE)$ix
Methods[ARBPMediansortidbRtune]
plot(density(ARBPMediansortdbRtune))
plot(1:127, ARBPMediandbRtune, ylim = c())
# 127  70   3   4  69   1   2


RMSEPMediandbRtune <- dbRtune$RMSEPMedian
RMSEPMediansortdbRtune <- sort(RMSEPMediandbRtune, index.return = TRUE)$x
RMSEPMediansortidbRtune <- sort(RMSEPMediandbRtune, index.return = TRUE)$ix
Methods[RMSEPMediansortidbRtune]
plot(density(RMSEPMediandbRtune))
plot(1:127, RMSEPMediandbRtune)
#3  69  70   4   2  55  56   1  54  53  52  51
# 69   3   4   2 1           53 54 55 56 70


#ind = c(3, 69)
#vec <- numeric(127)
vec[ind] <- RMSEPMediandbRtune[ind]

indRtune = c(3, 69, 70, 4, 2, 55, 56, 1, 54, 53, 52, 51)
vecRtune <- numeric(127)
vecRtune[indRtune] <- RMSEPMediandbRtune[indRtune]



dbXtune <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsXtuneF.csv')
ARBPMediandbXtune <- dbXtune$ARBPMedian
ARBPMediansortdbXtune <- sort(ARBPMediandbXtune, index.return = TRUE)$x
ARBPMediansortidbXtune <- sort(ARBPMediandbXtune, index.return = TRUE)$ix
Methods[ARBPMediansortidbXtune]
plot(density(ARBPMediansortdbXtune))
plot(1:127, ARBPMediandbXtune)
#127  69   3   4   1   2  70  56  54  53
#127  70   3  69   4   1   2
RMSEPMediandbXtune <- dbXtune$RMSEPMedian
RMSEPMediansortdbXtune <- sort(RMSEPMediandbXtune, index.return = TRUE)$x
RMSEPMediansortidbXtune <- sort(RMSEPMediandbXtune, index.return = TRUE)$ix
Methods[RMSEPMediansortidbXtune]
plot(density(RMSEPMediandbXtune))
plot(1:127, RMSEPMediandbXtune)
# 3  69   4   2   1

indXtune = c(3, 69, 70, 4, 2, 1)
vecXtune <- numeric(127)
vecXtune[indXtune] <- RMSEPMediandbXtune[indXtune]

vec[ind] <- RMSEPMediandbXtune[ind]
