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
E_OR <- dt_res[dt_res$Method == "E_OR"]
dt_Mean <- dt_res[, lapply(.SD, mean), by = "Method"]

# 25
tau_trueM <- dt_Mean[c(4:72), c(77:101)]
# 50
tau_trueM <- dt_Mean[c(4:72), c(152:201)]
# 75
tau_trueM <- dt_Mean[c(4:72), c(227:301)]
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

RB <- dt_Mean[c(4:72), c(77:151)]/abs(tau_trueM)
RBP <- dt_Mean[c(4:72), c(77:151)]/abs(tau_trueM) * 100

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

ARB <- abs(dt_Mean[c(4:72), c(77:151)]/abs(tau_trueM))
ARBP <- abs(dt_Mean[c(4:72), c(77:151)]/abs(tau_trueM) * 100)

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

RMSE <- sqrt(dt_Mean[c(4:72), c(152:226)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:72), c(152:226)])/abs(tau_trueM) * 100

RMSE <- sqrt(dt_Mean[c(4:72), c(202:301)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:72), c(202:301)])/abs(tau_trueM) * 100

# Design based --------------------------------------------
RMSE <- sqrt(dt_Mean[c(4:72), c(84:124)])/abs(tau_trueM)
RMSEP <- sqrt(dt_Mean[c(4:72), c(84:124)])/abs(tau_trueM) * 100


RMSEMean <- rowMeans(RMSEP)
RMSEPMedian <- apply(RMSEP, 1, quantile, probs = 0.5)


resultsDB75t <- data.frame(Method = Methods,

           RBPMean = RBPMean,
           RBPMedian = RBPMedian,

           ARBPMean = ARBPMean,
           ARBPMedian = ARBPMedian,

           RMSEMean = RMSEMean,
           RMSEPMedian = RMSEPMedian)

write.csv(resultsDB75t,
          file = 'resultsDB75t.csv', row.names = FALSE)

# Analyse results

MPlus25F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/PositiveEffects/resultsMPlus25F.csv')
MPlus25F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/DBscenarios/resultsDB25F.csv')
Methods <- MPlus25F$Method
ARBPMedian25F <- MPlus25F$ARBPMedian
ARBPMediansort25F <- sort(ARBPMedian25F, index.return = TRUE)$x
ARBPMediansorti25F <- sort(ARBPMedian25F, index.return = TRUE)$ix
Methods[ARBPMediansorti25F][1:5]
plot(density(ARBPMediansort25F))
plot(1:127, ARBPMedian25F)
plot(1:127, ARBPMedian25F, ylim = c(7,15))
abline(a = 7.5, b = 0)
which(ARBPMedian25F <= 11)
which(ARBPMedian25F <= 7.5)
Methods[which(ARBPMedian25F <=11)]
#DB
# 1   2   3   4   6   9  15  19  21  24  27  30  46  49  52  55  58  61  64  67  69  77
# 81  87  90 104 107 110 113 116

which(ARBPMedian25F <= 7.5)
# 2   6   9  52  55  77  81 110 113
# 2   6  55  77 110
RMSEPMedian25F <- MPlus25F$RMSEPMedian
RMSEPMediansort25F <- sort(RMSEPMedian25F, index.return = TRUE)$x
RMSEPMediansorti25F <- sort(RMSEPMedian25F, index.return = TRUE)$ix
Methods[RMSEPMediansorti25F][1:10]
plot(density(RMSEPMedian25F))
plot(1:127, RMSEPMedian25F)
plot(1:126, RMSEPMedian25F[-127])
plot(1:126, RMSEPMedian25F[-127], ylim = c(30, 40))
abline(35.5,0)
which(RMSEPMedian25F <= 35.5)
Methods[which(RMSEPMedian25F <= 41)]
#DB
# 2   3   6   9  15  19  21  24  27  30  46  49  52  55  58  61  64  67  69  77  81  87
# 90 104 107 110 113 116 119 122 125


MPlus50F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/PositiveEffects/resultsMPlus50F.csv')
MPlus50F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/DBscenarios/resultsDB50F.csv')
ARBPMedian50F <- MPlus50F$ARBPMedian
ARBPMediansort50F <- sort(ARBPMedian50F, index.return = TRUE)$x
ARBPMediansorti50F <- sort(ARBPMedian50F, index.return = TRUE)$ix
Methods[ARBPMediansorti50F]
Methods[ARBPMediansorti50F][1:10]
plot(density(ARBPMediansort50F))
plot(1:127, ARBPMedian50F, ylim = c(10,13))
abline(10.5,0)
which(ARBPMedian50F <= 13)
Methods[which(ARBPMedian50F <= 14)]
#DB 25
# 1   2   3   4   6   9  15  19  21  24  27  30  46  49  52  55  58  61  64  67  69  77
# 81  87  90 104 107 110 113 116

#DB 50
#   1   2   4   6   8   9  15  19  21  23  24  26  27  29  30  32  35  38  41  44  46  47  49  50  52  53  55
#  56  58  61  64  67  77  81  87  90 104 107 110 113 116 127

which(ARBPMedian50F <= 11)
#2  21  24  27  30  46  49  58  61  64  67  77  87  90 127
which(ARBPMedian50F <= 10.5)
# 27 30 46 58 61 64 67 87 90

RMSEPMedian50F <- MPlus50F$RMSEPMedian
RMSEPMediansort50F <- sort(RMSEPMedian50F, index.return = TRUE)$x
RMSEPMediansorti50F <- sort(RMSEPMedian50F, index.return = TRUE)$ix
Methods[RMSEPMediansorti50F]
Methods[RMSEPMediansorti50F][1:10]
plot(density(RMSEPMedian50F))
plot(1:127, RMSEPMedian50F)
plot(1:126, RMSEPMedian50F[-127], ylim = c(20,40))
abline(30,0)
abline(26.5,0)
which(RMSEPMedian50F <= 30)
Methods[which(RMSEPMedian50F <= 30)]
# DB 36
#  2   3   4   6   8   9  12  15  16  17  19  21  23  24  26  27  29  30  32  35  38  41
#  44  46  47  49  50  52  53  55  56  58  61  64  67  69  77  81  87  90 104 107 110 113
# 116 119 122 125

#   2   3   6   9  15  19  21  24  27  30  46  49  52  55  58  61  64  67  69  77  81  87  90 104 107 110 113
#116
MPlus100F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/PositiveEffects/resultsMPlus100F.csv')
MPlus100F <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/DBscenarios/resultsDB75F.csv')
ARBPMedian100F <- MPlus100F$ARBPMedian
ARBPMediansort100F <- sort(ARBPMedian100F, index.return = TRUE)$x
ARBPMediansorti100F <- sort(ARBPMedian100F, index.return = TRUE)$ix
Methods[ARBPMediansorti100F][1:10]
plot(density(ARBPMediansort100F))
plot(1:127, ARBPMedian100F, ylim = c(7,14))
which(ARBPMedian100F <= 12.5)
abline(12.5,0)
which(ARBPMedian100F <= 11)
abline(11, 0)
which(ARBPMedian100F <= 7.9)
abline(7.9,0)
Methods[which(ARBPMedian100F <= 12.5)]
#   1   2   3   4   6   8   9  12  15  19  21  23  24  26  27  29  30  32  35  38  41  44  46  47  49  50  52
#  53  55  56  58  61  64  67  69  77  81  87  90 104 107 110 113 116


which(ARBPMedian100F <= 11)
# 1   2   3   4   6   9  15  19  21  24  27  30  46  49  52  55  58  61  64  67  69  77  81  87  90
# 104 107 110 113 116


which(ARBPMedian100F <= 8)
#6   9  19  21  77  81 104 116

v100 <- which(ARBPMedian100F <= 18.63702)

RMSEPMedian100F <- MPlus100F$RMSEPMedian
RMSEPMediansort100F <- sort(RMSEPMedian100F, index.return = TRUE)$x
RMSEPMediansorti100F <- sort(RMSEPMedian100F, index.return = TRUE)$ix
Methods[RMSEPMediansorti100F][1:10]
RMSEPMediansort100F[1:10]
plot(density(RMSEPMedian100F))
plot(1:127, RMSEPMedian100F)
plot(1:126, RMSEPMedian100F[-127], ylim = c(18,40))
abline(21.2,0)
Methods[which(RMSEPMedian100F <= 21.2)]
#DB
#2   3   6   9  15  19  21  24  27  30  46  49  52  55  58  61  64  67  69  77  81  87
#  90 104 107 110 113 116

# 2   6   9  15  19  21  46  49  52  55  77  81 104 107 110 113 116
ind = c(2, 3, 6, 9, 15, 19, 21, 24, 27, 30, 46, 49, 52,
        55, 58, 61, 64, 67, 69, 77, 81, 87, 90, 104, 107,
        110, 113, 116)

Methods[ind100]
ind100  = c(2, 6, 9,  15,  19,  21,  46,  49,  52,  55,  77,  81, 104, 107, 110, 113, 116)
#which(RMSEPMedian50F <= 26.3)
#          2   6   9  15  19  21  46  49  52  55 104 107 110 113 116
Methods[ind100]
Methods[ind50]
Methods[ind25]
ind50  = c(2, 6, 9,  15, 19, 21,  46,  49,  52,  55,  104, 107, 110, 113, 116)
ind25  = c(2,   6,   9,  15,  19,  21,  24,  27,  30,  46,  49,  52,  55,  69,  77,  81,  87,  90, 104, 107, 110, 113, 116)

vec <- numeric(127)
vec100 <- numeric(127)
vec50 <- numeric(127)
vec25 <- numeric(127)


vec[ind] <- RMSEPMedian100F[ind]
vec[ind] <- RMSEPMedian50F[ind]
vec[ind] <- RMSEPMedian25F[ind]

Methods[ind]

vec100[ind100] <- RMSEPMedian100F[ind100]
vec50[ind50] <- RMSEPMedian50F[ind50]
vec25[ind25] <- RMSEPMedian25F[ind25]




ARBP <- c(Methods[ARBPMediansorti100F][1:15],
          Methods[ARBPMediansorti50F][1:15],
          Methods[ARBPMediansorti25F][1:15])
namesARB <- names(table(ARBP))
occARB <- as.vector(table(ARBP))
indsort <- sort(occARB,  decreasing = TRUE, index.return = TRUE)$ix
bestmethods <- namesARB[indsort][1:3]
# Bias "EEM_AIPW"  "tRM_NIPW"  "XEM_AIPW"

RMSEP <- c(Methods[RMSEPMediansorti100F][1:15],
           Methods[RMSEPMediansorti50F][1:15],
           Methods[RMSEPMediansorti25F][1:15])
namesRMS <- names(table(RMSEP))
occRMS <- as.vector(table(RMSEP))
indsortRMS <- sort(occRMS,  decreasing = TRUE, index.return = TRUE)$ix
bestmethodsR <- namesRMS[indsortRMS][1:10]
# "EM_NIPW"   "MM_NIPW"   "REM_AIPW"  "RM_NIPW"   "RMM_AIPW"  "RRM_AIPW"
# "tREM_AIPW" "tRRM_AIPW" "tRXM_AIPW" "XM_NIPW"


#"EEM_AIPW" "EMM_AIPW" "M_OR"     "MM_NIPW"  "R_OR"     "RRM_AIPW" "XM_NIPW"

# Plots
plot(1:127, RMSEPMedian25F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "RRMSE", ylim = c(20, 90),
     xaxt="n", yaxt="n", cex.main = 2)
plot(1:127, RMSEPMedian50F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "RRMSE", ylim = c(20, 90),
     xaxt="n", yaxt="n", cex.main = 2)
plot(1:127, RMSEPMedian100F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "RRMSE", ylim = c(20, 90),
     xaxt="n", yaxt="n", cex.main = 2)
title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b=c(10, 20, 30, 40, 50, 60, 70,  80, 90,  100, 120, 140)
axis(2, at=b,labels=b,cex.axis=1.5)

plot(1:127, RMSEPMedian25F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "RRMSE", ylim = c(20, 40),
     xaxt="n", yaxt="n", cex.main = 2)
plot(1:127, RMSEPMedian50F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "RRMSE", ylim = c(20, 40),
     xaxt="n", yaxt="n", cex.main = 2)
plot(1:127, RMSEPMedian100F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "RRMSE", ylim = c(20, 40),
     xaxt="n", yaxt="n", cex.main = 2)
lines(1:127, vec, col = "blue", lwd = 4, type = "p", lty = 0)
lines(1:127, vec100, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, vec50, col = "darkgreen", lwd = 4, type = "p", lty = 0)
lines(1:127, vec25, col = "darkgreen", lwd = 4, type = "p", lty = 0)

title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b=c(10, 20, 25,  30, 35,  40, 60, 80, 100, 120, 140)
axis(2, at=b,labels=b,cex.axis=1.5)
title(xlab = "Estimation methods", line = 2.1, cex.lab = 1.2)

ind = c(2, 3, 6, 9, 15,19,  21,  24,
        27,  30,  46,  49,  52,  55,
        58,  61,  64,  67,  69,  77,
        81,  87,  90, 104, 107, 110,
        113, 116)

ind = c(6, 77)

ind25  = c(1,   2,   3,   4,   6,   9,  15,  19,  21,  24,
            27,  30,  46,  49,  52,  55,  58,  61,  64,  67,
            69,  77, 81,  87,  90, 104, 107, 110, 113, 116)

ind50  = c(1,   2,   4,   6,   8,   9,  15,  19,  21,  23,
           24,  26,  27,  29,  30,  32,  35,  38,  41,  44,
           46,  47,  49,  50,  52,  53,  55, 56,  58,  61,  64,
           67,  77,  81,  87,  90, 104, 107, 110, 113, 116, 127)

#ind100  = c(1,   2,   3,   4,   6,   8,   9,  12,  15,  19,  21,
#            23,  24,  26,  27,  29,  30,  32,  35,  38,  41,  44,
#            46,  47,  49,  50,  52,  53,  55,  56,  58,  61,  64,
#            67,  69,  77,  81,  87,  90, 104, 107, 110, 113, 116)

ind100  = c(1,   2,   3,   4,   6,   9,  15,  19,  21,  24,  27,
            30, 46,  49,  52,  55,  58,  61,  64,  67,  69,  77,
            81,  87,  90, 104, 107, 110, 113, 116)


vec <- numeric(127)
vec100 <- numeric(127)
vec50 <- numeric(127)
vec25 <- numeric(127)


vec100[ind100] <- ARBPMedian100F[ind100]
vec50[ind50] <- ARBPMedian50F[ind50]
vec25[ind25] <- ARBPMedian25F[ind25]


vec1002 <- numeric(127)
vec502 <- numeric(127)
vec252 <- numeric(127)

ind252 <- c(2,   6,  55,  77, 110)
Methods[c(2,   6,  55,  77, 110)]
ind502 <- c(27, 30, 46, 58, 61, 64, 67, 87, 90)
Methods[ind502]
ind1002 <- c(6, 9,  19,  21,  77,  81, 104, 116)
Methods[ind1002]

vec1002[ind1002] <- ARBPMedian100F[ind1002]
vec502[ind502] <- ARBPMedian50F[ind502]
vec252[ind252] <- ARBPMedian25F[ind252]


plot(1:127, ARBPMedian100F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB", ylim = c(7, 85),
     xaxt="n", yaxt="n", cex.main = 2)
plot(1:127, ARBPMedian50F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB", ylim = c(7, 85),
     xaxt="n", yaxt="n", cex.main = 2)
plot(1:127, ARBPMedian25F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB", ylim = c(7, 85),
     xaxt="n", yaxt="n", cex.main = 2)
title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b = c(10, 20,  30,  40, 50,  60, 70,  80, 90,  100, 120, 140)
axis(2, at=b,labels=b,cex.axis=1.5)

Meth <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/DBscenarios/Methods.csv')
Labels <- Meth$X

brewer.pal(7, "Set2")
#"#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
# OR        NIPW       AIPW       tOR     tNIPW      tAIPW      dir
which(RMSEPMedian25F <= 35.5)
Methods[which(RMSEPMedian25F <= 41)]

ARBPMedian25Fs <- ARBPMedian25F[which(ARBPMedian25F <= 12)]
Methods25 <- Methods[which(ARBPMedian25F <= 12)]
Labels25 <- Labels[which(ARBPMedian25F <= 12)]
Type25 <- c(rep("OR", 4),
          rep("NIPW", 4),
          rep("AIPW", 12),
          rep("tOR", 1),
          rep("tNIPW", 2),
          rep("tAIPW", 7))

df25 <- data.frame(method = c(1:30),
                   RRMSE = ARBPMedian25Fs,
                   Label = Labels25,
                   Type = Type25)

RMSEPMedian25Fs <- RMSEPMedian25F[which(RMSEPMedian25F <= 35.5)]
Methods25 <- Methods[which(RMSEPMedian25F <= 35.5)]
Labels25 <- Labels[which(RMSEPMedian25F <= 35.5)]
Type25 <- c(rep("OR", 1),
            rep("NIPW", 4),
            rep("AIPW", 8),
            rep("tOR", 1),
            rep("tNIPW", 2),
            rep("tAIPW", 7))

df25 <- data.frame(method = c(1:23),
                   RRMSE = RMSEPMedian25Fs,
                   Label = Labels25,
                   Type = Type25)


ARBPMedian50Fs <- ARBPMedian50F[which(ARBPMedian50F <= 13)]
Methods50 <- Methods[(which(ARBPMedian50F <= 13))]
Labels50 <-  Labels[which(ARBPMedian50F <= 13)]


Type50 <- c(rep("OR", 3),
          rep("NIPW", 4),
          rep("AIPW", 19),
          rep("tNIPW", 2),
          rep("tAIPW", 7),
          rep("dir", 1))

df50 <- data.frame(method = c(1:36),
                   ARB = ARBPMedian50Fs,
                   Label = Labels50,
                   Type = Type50)


which(RMSEPMedian50F <= 30)

RMSEPMedian50Fs <- RMSEPMedian50F[which(RMSEPMedian50F <= 30)]
Methods50 <- Methods[which(RMSEPMedian50F <= 30)]
Labels50 <- Labels[which(RMSEPMedian50F <= 30)]

Type50 <- c(rep("OR", 2),
            rep("NIPW", 4),
            rep("AIPW", 12),
            rep("tOR", 1),
            rep("tNIPW", 2),
            rep("tAIPW", 7))

df50 <- data.frame(method = c(1:28),
                   RRMSE = RMSEPMedian50Fs,
                   Label = Labels50,
                   Type = Type50)



ARBPMedian100Fs <- ARBPMedian100F[which(ARBPMedian100F <= 11)]
Methods100 <- Methods[(which(ARBPMedian100F <= 11))]
Labels100 <-  Labels[which(ARBPMedian100F <= 11)]


Type100 <- c(rep("OR", 4),
            rep("NIPW", 4),
            rep("AIPW", 12),
            rep("tOR", 1),
            rep("tNIPW", 2),
            rep("tAIPW", 7))

df100 <- data.frame(method = c(1:30),
                   ARB = ARBPMedian100Fs,
                   Label = Labels100,
                   Type = Type100)


RMSEPMedian100Fs <- RMSEPMedian100F[which(RMSEPMedian100F <= 25)]
Methods100 <- Methods[(which(RMSEPMedian100F <= 25))]
Labels100 <-  Labels[which(RMSEPMedian100F <= 25)]


Type100 <- c(rep("OR", 2),
             rep("NIPW", 4),
             rep("AIPW", 12),
             rep("tOR", 1),
             rep("tNIPW", 2),
             rep("tAIPW", 7))

df100 <- data.frame(method = c(1:28),
                    RRMSE = RMSEPMedian100Fs,
                    Label = Labels100,
                    Type = Type100)




g25 <- ggplot(df25, aes(x = method,
#                 y = RMSEPMedian25Fs,
                 y = ARBPMedian25Fs,
                 col = Type)) +
#  labs(y = "RRMSE", x = "Method")+theme_bw()+
  labs(y = "ARB", x = "Method")+theme_bw()+
#  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(19, 41)) +
#  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(19,41,2)) +
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(5, 14)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(5,20,2)) +

  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "NIPW"="#FC8D62",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tNIPW"="#A6D854",
                                "tAIPW"="#FFD92F" ))+
#  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = Labels25), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

g25 +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

#library(ggrepel)
#"#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"

g50 <- ggplot(df50, aes(x = method,
#                        y = RMSEPMedian50Fs,
                        y = ARBPMedian50Fs,
                        col = Type)) +
#  labs(y = "RRMSE", x = "Method")+theme_bw()+
  labs(y = "ARB", x = "Method")+theme_bw()+
#  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(19, 41)) +
#  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(19,41,2))+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(5, 14)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(5,20,2)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "NIPW"="#FC8D62",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tNIPW"="#A6D854",
                                "tAIPW"="#FFD92F"))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = Labels50), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

g50 +   ggtitle("m = 50") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )



g100 <- ggplot(df100, aes(x = method,
                          y = RMSEPMedian100Fs,
                          col = Type)) +
  labs(y = "RRMSE", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(19, 41)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(19,41,2)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "NIPW"="#FC8D62",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tNIPW"="#A6D854",
                                "tAIPW"="#FFD92F"))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = Labels100), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

g100 +   ggtitle("m = 100") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )


#install.packages("cowplot")
#library(cowplot)
plot_grid(g25, g50, g100,
          labels = c("(a)", "(b)", "(c)"), nrow = 1)




plot(1:30, ARBPMedian25Fs,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB", ylim = c(7, 14),
     xaxt="n", yaxt="n", cex.main = 2)

lines(1:30, vecor, col = "1", lwd = 1, pch = 16, type = "p", lty = 0)
lines(1:30, vecnipw, col = "2", lwd = 1, pch = 16, type = "p", lty = 0)
lines(1:30, vecaipw, col = "3", lwd = 1, pch = 16, type = "p", lty = 0)

lines(1:30, vecort, col = "4", lwd = 1, pch = 16, type = "p", lty = 0)
lines(1:30, vecnipwt, col = "5", lwd = 1, pch = 16, type = "p", lty = 0)
lines(1:30, vecaipwt, col = "6", lwd = 1, pch = 16, type = "p", lty = 0)

title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b = seq(4,20,2)
axis(2, at = b, labels = b, cex.axis = 1.5)
#axis(1, at = 1:30, labels = 1:30, cex.axis=1)




plot(1:127, ARBPMedian100F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB", ylim = c(7, 20),
     xaxt="n", yaxt="n", cex.main = 2)
lines(1:127, vec100, col = "blue", lwd = 4, type = "p", lty = 0)
lines(1:127, vec1002, col = "darkgreen", lwd = 4, type = "p", lty = 0)
title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b = seq(4,20,2)
axis(2, at=b,labels=b,cex.axis=1.5)


plot(1:127, ARBPMedian50F,  col = 1, lwd = 2, pch = 20,
     xlab = " ", ylab = " ", main = "ARB", ylim = c(7, 20),
     xaxt="n", yaxt="n", cex.main = 2)
lines(1:127, vec50, col = "blue", lwd = 4, type = "p", lty = 0)
lines(1:127, vec502, col = "darkgreen", lwd = 4, type = "p", lty = 0)
title(xlab = "Estimation methods", line = 1, cex.lab = 2)
b = seq(4,20,2)
axis(2, at=b,labels=b,cex.axis=1.5)
legend("top",c("OR","NIPW","AIPW",
               "ORt", "NIPWt", "AIPWt"),
       bty = "n",
       lty=1,
       col=c(1:6),
       cex=1,  pt.cex = 2, y.intersp=0.8, ncol = 6)


#title(ylab = "RRMSE", line = 2, cex.lab = 1.2)

dbE <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsEF.csv')
Methods <- dbE$Method

ARBPMediandbE <- dbE$ARBPMedian
ARBPMediansortdbE <- sort(ARBPMediandbE, index.return = TRUE)$x
ARBPMediansortidbE <- sort(ARBPMediandbE, index.return = TRUE)$ix
Methods[ARBPMediansortidbE]
plot(density(ARBPMediansortdbE))
plot(1:127, ARBPMediandbE)

# "Dir_tau"   "tX_OR"     "R_OR"      "tR_OR"     "X_OR"      "E_OR"      "M_OR"
# "RRM_AIPW"  "RRE_AIPW"  "RXM_AIPW"  "RXE_AIPW"  "RRX_AIPW"  "RXX_AIPW"  "tRXE_AIPW"

RMSEPMediandbE <- dbE$RMSEPMedian
RMSEPMediansortdbE <- sort(RMSEPMediandbE, index.return = TRUE)$x
RMSEPMediansortidbE <- sort(RMSEPMediandbE, index.return = TRUE)$ix
Methods[RMSEPMediansortidbE]
plot(density(RMSEPMediandbE))
plot(1:127, RMSEPMediandbE, ylim = c(35, 100))

# [1] "tR_OR"     "R_OR"      "X_OR"      "M_OR"      "RXM_AIPW"  "RXE_AIPW"  "RXX_AIPW"
#[8] "E_OR"      "tX_OR"     "RRM_AIPW"  "RRX_AIPW"  "tRXM_AIPW" "tRXE_AIPW" "tRRX_AIPW"

dbM <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsM.csv')
dbX <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsX.csv')
dbR <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsR.csv')
dbRtune <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsRtune.csv')
dbXtune <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsXtune.csv')


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

