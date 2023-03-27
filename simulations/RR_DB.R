#setwd("C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/Mplustune25_2")
#setwd("C:/Users/katar/Documents/GitHub/causalSAE/simulations/results/ModelBased/PositiveEffects/Tuned/Mplustune100_2")
setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/dbXtune")
setwd("C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/DBscenarios/DB50")
#library(ggrepel)
library("RColorBrewer")
library(data.table)
library(cowplot)
library(ggplot2)
library(ggrepel)
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

###############################
dbE <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsEF.csv')
#dbE <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsEoutF.csv')

Methods <- dbE$Method

ARBPMediandbE <- dbE$ARBPMedian
ARBPMediansortdbE <- sort(ARBPMediandbE, index.return = TRUE)$x
ARBPMediansortidbE <- sort(ARBPMediandbE, index.return = TRUE)$ix
Methods[ARBPMediansortidbE]
plot(density(ARBPMediansortdbE))
plot(1:127, ARBPMediandbE, ylim = c(10,80))
abline(58, 0)
#127  70   3  69   4   1   2

RMSEPMediandbE <- dbE$RMSEPMedian
RMSEPMediansortdbE <- sort(RMSEPMediandbE, index.return = TRUE)$x
RMSEPMediansortidbE <- sort(RMSEPMediandbE, index.return = TRUE)$ix
Methods[RMSEPMediansortidbE]
plot(density(RMSEPMediandbE))
plot(1:127, RMSEPMediandbE, ylim = c(35, 80))
abline(73.5, 0)
plot(1:127, RMSEPMediandbM, ylim = c(35, 80))
abline(72.5, 0)
plot(1:127, RMSEPMediandbR, ylim = c(35, 80))
abline(70.5, 0)
plot(1:127, RMSEPMediandbX, ylim = c(35, 80))
abline(74, 0)
plot(1:127, RMSEPMediandbRtune, ylim = c(35, 80))
abline(71.5, 0)
plot(1:127, RMSEPMediandbXtune, ylim = c(35, 80))
abline(73, 0)
#  69   3   4   2  55    54  56   1  70  52    53
# 69   3   4   2 1 53 54 55 56 70


dbM <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsMF.csv')
#dbM <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsMoutF.csv')

ARBPMediandbM <- dbM$ARBPMedian
ARBPMediansortdbM <- sort(ARBPMediandbM, index.return = TRUE)$x
ARBPMediansortidbM <- sort(ARBPMediandbM, index.return = TRUE)$ix
Methods[ARBPMediansortidbM][1:12]
plot(density(ARBPMediansortdbM))
plot(1:127, ARBPMediandbM, ylim = c(10,80))
abline(58,0)
#127  70  69   3   4   2   1

RMSEPMediandbM <- dbM$RMSEPMedian
RMSEPMediansortdbM <- sort(RMSEPMediandbM, index.return = TRUE)$x
RMSEPMediansortidbM <- sort(RMSEPMediandbM, index.return = TRUE)$ix
Methods[RMSEPMediansortidbM]
plot(density(RMSEPMediandbM))
plot(1:127, RMSEPMediandbM)
#DR TMLE (double robust TMLE)
# 3  69   4   2     55 54  70  56    1  52    53 110
# 69   3   4   2 1 53 54 55 56 70

dbX <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsXF.csv')
dbX <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsXoutF.csv')

ARBPMediandbX <- dbX$ARBPMedian
ARBPMediansortdbX <- sort(ARBPMediandbX, index.return = TRUE)$x
ARBPMediansortidbX <- sort(ARBPMediandbX, index.return = TRUE)$ix
Methods[ARBPMediansortidbX]
plot(density(ARBPMediansortdbX))
plot(1:127, ARBPMediandbX)
abline(58,0)
# 127  70   3  69   4   1   2  53  55  54  51  56 112  66  67
#127  70   3  69   4   1   2

RMSEPMediandbX <- dbX$RMSEPMedian
RMSEPMediansortdbX <- sort(RMSEPMediandbX, index.return = TRUE)$x
RMSEPMediansortidbX <- sort(RMSEPMediandbX, index.return = TRUE)$ix
Methods[RMSEPMediansortidbX]
plot(density(RMSEPMediandbX))
plot(1:127, RMSEPMediandbX)
# 3  69   4  70   2   1

dbR <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsRF.csv')
dbR <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsRoutF.csv')

Methods <- dbR$Method
ARBPMediandbR <- dbR$ARBPMedian
ARBPMediansortdbR <- sort(ARBPMediandbR, index.return = TRUE)$x
ARBPMediansortidbR <- sort(ARBPMediandbR, index.return = TRUE)$ix
Methods[ARBPMediansortidbR]
plot(density(ARBPMediansortdbR))
plot(1:127, ARBPMediandbR, ylim = c(10, 80))
abline(58, 0)
#127  70   3  69   4   1   2
#127  70   3  69   4   1   2


RMSEPMediandbR <- dbR$RMSEPMedian
RMSEPMediansortdbR <- sort(RMSEPMediandbR, index.return = TRUE)$x
RMSEPMediansortidbR <- sort(RMSEPMediandbR, index.return = TRUE)$ix
Methods[RMSEPMediansortidbR]
plot(density(RMSEPMediandbR))
plot(1:127, RMSEPMediandbR)
# 3  69  70   4      2  54  55   1    52  56  51  53 110
# 69   3   4   2 1           53 54 55 56 70

dbRtune <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsRtuneF.csv')
dbRtune <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsRtuneoutF.csv')

ARBPMediandbRtune <- dbRtune$ARBPMedian
ARBPMediansortdbRtune <- sort(ARBPMediandbRtune, index.return = TRUE)$x
ARBPMediansortidbRtune <- sort(ARBPMediandbRtune, index.return = TRUE)$ix
Methods[ARBPMediansortidbRtune]
plot(density(ARBPMediansortdbRtune))
plot(1:127, ARBPMediandbRtune, ylim = c(10, 80))
abline(58,0)
# 127  70   3   4  69   1   2


RMSEPMediandbRtune <- dbRtune$RMSEPMedian
RMSEPMediansortdbRtune <- sort(RMSEPMediandbRtune, index.return = TRUE)$x
RMSEPMediansortidbRtune <- sort(RMSEPMediandbRtune, index.return = TRUE)$ix
Methods[RMSEPMediansortidbRtune]
plot(density(RMSEPMediandbRtune))
plot(1:127, RMSEPMediandbRtune)
#3  69  70   4   2  55  56   1  54  53  52  51
# 69   3   4   2 1           53 54 55 56 70


dbXtune <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsXtuneF.csv')
dbXtune <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/DesignBased/resultsXtuneoutF.csv')
ARBPMediandbXtune <- dbXtune$ARBPMedian
ARBPMediansortdbXtune <- sort(ARBPMediandbXtune, index.return = TRUE)$x
ARBPMediansortidbXtune <- sort(ARBPMediandbXtune, index.return = TRUE)$ix
Methods[ARBPMediansortidbXtune]
plot(density(ARBPMediansortdbXtune))
plot(1:127, ARBPMediandbXtune, ylim = c(10,80))
abline(58,0)
#127  69   3   4   1   2  70  56  54  53
#127  70   3  69   4   1   2
RMSEPMediandbXtune <- dbXtune$RMSEPMedian
RMSEPMediansortdbXtune <- sort(RMSEPMediandbXtune, index.return = TRUE)$x
RMSEPMediansortidbXtune <- sort(RMSEPMediandbXtune, index.return = TRUE)$ix
Methods[RMSEPMediansortidbXtune]
plot(density(RMSEPMediandbXtune))
plot(1:127, RMSEPMediandbXtune)
# 3  69   4   2   1




ARBP <- c(Methods[ARBPMediansortidbE][1:15],
          Methods[ARBPMediansortidbM][1:15],
          Methods[ARBPMediansortidbXtune][1:15],
          Methods[ARBPMediansortidbX][1:15],
          Methods[ARBPMediansortidbRtune][1:15],
          Methods[ARBPMediansortidbR][1:15])
namesARB <- names(table(ARBP))
occARB <- as.vector(table(ARBP))
indsort <- sort(occARB,  decreasing = TRUE, index.return = TRUE)$ix
bestmethods <- namesARB[indsort][1:9]
bestmethods <- namesARB[indsort][1:10]
# Bias "Dir_tau"   "E_OR"      "M_OR"      "R_OR"      "RXM_AIPW"  "tR_OR"
# "tX_OR"     "tXXR_AIPW" "X_OR"

# "Dir_tau"  "E_OR"     "M_OR"     "R_OR"     "RXE_AIPW" "RXM_AIPW"
# "RXX_AIPW" "tR_OR"    "tX_OR" "X_OR"
RMSEP <- c(Methods[RMSEPMediansortidbE][1:15],
           Methods[RMSEPMediansortidbM][1:15],
           Methods[RMSEPMediansortidbR][1:15],
           Methods[RMSEPMediansortidbRtune][1:15],
           Methods[RMSEPMediansortidbX][1:15],
           Methods[RMSEPMediansortidbXtune][1:15])
namesRMS <- names(table(RMSEP))
occRMS <- as.vector(table(RMSEP))
indsortRMS <- sort(occRMS,  decreasing = TRUE, index.return = TRUE)$ix
bestmethodsR <- namesRMS[indsortRMS][1:11]
bestmethodsR <- namesRMS[indsortRMS][1:10]
#RMSE  "E_OR"      "M_OR"      "R_OR"      "RRE_AIPW"  "RRM_AIPW"  "RXE_AIPW"
# "RXM_AIPW"  "tR_OR"     "tRXM_AIPW" "tX_OR"     "X_OR"

# "E_OR"     "M_OR"     "R_OR"     "RRX_AIPW" "RXE_AIPW" "RXM_AIPW" "RXX_AIPW"
# "tR_OR"    "tX_OR"    "X_OR"

# Plots
###################################################################
#ARB
###################################################################
#E
ARBPMediandbEs <- ARBPMediandbE[which(ARBPMediandbE <= 58)]
MethodsE <- Methods[(which(ARBPMediandbE <= 58))]
LabelsE <-  Labels[which(ARBPMediandbE <= 58)]
brewer.pal(7, "Set2")
#"#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
# OR        NIPW       AIPW       tOR     tNIPW      tAIPW      dir

TypeE <- c(rep("OR", 4),
             rep("AIPW", 6),
             rep("tOR", 2),
             rep("tAIPW", 1),
             rep("dir", 1))

TypeE <- c(rep("OR", 4),
           rep("AIPW", 8),
           rep("tOR", 2),
           rep("tAIPW", 12),
           rep("dir", 1))

dfE <- data.frame(method = c(1:27),
                    ARB = ARBPMediandbEs,
                    Label = LabelsE,
                    Type = TypeE)


gE <- ggplot(dfE, aes(x = method,
                      y = ARBPMediandbEs,
                      col = TypeE)) +
  labs(y = "ARB", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(10, 60)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(9,61,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F",
                                "dir"="#E5C494" ))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.1),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsE), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gE +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gE

#M
ARBPMediandbMs <- ARBPMediandbM[which(ARBPMediandbM <= 58)]
MethodsM <- Methods[(which(ARBPMediandbM <= 58))]
LabelsM <-  Labels[which(ARBPMediandbM <= 58)]

TypeM <- c(rep("OR", 4),
           rep("AIPW", 6),
           rep("tOR", 2),
           rep("tAIPW", 2),
           rep("dir", 1))

TypeM <- c(rep("OR", 4),
           rep("AIPW", 12),
           rep("tOR", 2),
           rep("tAIPW", 12),
           rep("dir", 1))

dfM <- data.frame(method = c(1:31),
                  ARB = ARBPMediandbMs,
                  Label = LabelsM,
                  Type = TypeM)


gM <- ggplot(dfM, aes(x = method,
                      y = ARBPMediandbMs,
                      col = TypeM)) +
  labs(y = "ARB", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(10, 60)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(9,61,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F",
                                "dir"="#E5C494" ))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.1),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsM), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gM +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gM


#R
ARBPMediandbRs <- ARBPMediandbR[which(ARBPMediandbR <= 58)]
MethodsR <- Methods[(which(ARBPMediandbR <= 58))]
LabelsR <-  Labels[which(ARBPMediandbR <= 58)]


TypeR <- c(rep("OR", 4),
           rep("AIPW", 6),
           rep("tOR", 2),
           rep("tAIPW", 5),
           rep("dir", 1))


TypeR <- c(rep("OR", 4),
           rep("AIPW", 8),
           rep("tOR", 2),
           rep("tAIPW", 9),
           rep("dir", 1))

dfR <- data.frame(method = c(1:24),
                  ARB = ARBPMediandbRs,
                  Label = LabelsR,
                  Type = TypeR)



gR <- ggplot(dfR, aes(x = method,
                      y = ARBPMediandbRs,
                      col = TypeR)) +
  labs(y = "ARB", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(10, 60)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(9,61,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F",
                                "dir"="#E5C494" ))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.1),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsR), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gR +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gR

#Rtune
ARBPMediandbRtunes <- ARBPMediandbRtune[which(ARBPMediandbRtune <= 58)]
MethodsRtune <- Methods[(which(ARBPMediandbRtune <= 58))]
LabelsRtune <-  Labels[which(ARBPMediandbRtune <= 58)]


TypeRtune <- c(rep("OR", 4),
           rep("AIPW", 6),
           rep("tOR", 2),
           rep("tAIPW", 5),
           rep("dir", 1))

TypeRtune <- c(rep("OR", 4),
               rep("AIPW", 8),
               rep("tOR", 2),
               rep("tAIPW", 9),
               rep("dir", 1))

dfRtune <- data.frame(method = c(1:24),
                  ARB = ARBPMediandbRtunes,
                  Label = LabelsRtune,
                  Type = TypeRtune)


gRtune <- ggplot(dfRtune, aes(x = method,
                      y = ARBPMediandbRtunes,
                      col = TypeRtune)) +
  labs(y = "ARB", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(10, 60)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(9,61,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F",
                                "dir"="#E5C494" ))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.1),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsRtune), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gRtune +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gRtune


# X
ARBPMediandbXs <- ARBPMediandbX[which(ARBPMediandbX <= 58)]
MethodsX <- Methods[(which(ARBPMediandbX <= 58))]
LabelsX <-  Labels[which(ARBPMediandbX <= 58)]


TypeX <- c(rep("OR", 4),
           rep("AIPW", 8),
           rep("tOR", 2),
           rep("tAIPW", 2),
           rep("dir", 1))

TypeX <- c(rep("OR", 4),
           rep("AIPW", 10),
           rep("tOR", 2),
           rep("tAIPW", 10),
           rep("dir", 1))

dfX <- data.frame(method = c(1:27),
                  ARB = ARBPMediandbXs,
                  Label = LabelsX,
                  Type = TypeX)



gX <- ggplot(dfX, aes(x = method,
                      y = ARBPMediandbXs,
                      col = TypeX)) +
  labs(y = "ARB", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(10, 60)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(9,61,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F",
                                "dir"="#E5C494" ))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.1),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsX), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gX +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gX

# X tune
ARBPMediandbXtunes <- ARBPMediandbXtune[which(ARBPMediandbXtune <= 58)]
MethodsXtune <- Methods[(which(ARBPMediandbXtune <= 58))]
LabelsXtune <-  Labels[which(ARBPMediandbXtune <= 58)]


TypeXtune <- c(rep("OR", 4),
           rep("AIPW", 2),
           rep("tOR", 2),
           rep("dir", 1))

TypeXtune <- c(rep("OR", 4),
               rep("AIPW", 8),
               rep("tOR", 2),
               rep("tAIPW", 10),
               rep("dir", 1))

dfXtune <- data.frame(method = c(1:25),
                  ARB = ARBPMediandbXtunes,
                  Label = LabelsXtune,
                  Type = TypeXtune)

gXtune <- ggplot(dfXtune, aes(x = method,
                              y = ARBPMediandbXtunes,
                              col = TypeXtune)) +
  labs(y = "ARB", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(10, 60)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(9,61,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F",
                                "dir"="#E5C494" ))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.1),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsXtune), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gXtune +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gXtune


plot_grid(gE, gM, gR, gRtune, gX, gXtune,
          labels = c("(a)", "(b)", "(c)",
                     "(d)", "(e)", "(f)"), nrow = 2)

Meth <- read.csv('C:/Users/katar/Documents/Kasia/4_PostDoc/rok_2022_2023/simultaions_causalSAE/results/ModelBased/DBscenarios/Methods.csv')
Labels <- Meth$X

###########################################################################
#RRMSE
###########################################################################
#plot(1:127, RMSEPMediandbE, ylim = c(35, 80))
#abline(73.5, 0)
#E
RMSEPMediandbEs <- RMSEPMediandbE[which(RMSEPMediandbE <= 73.5)]
MethodsE <- Methods[(which(RMSEPMediandbE <= 73.5))]
LabelsE <-  Labels[which(RMSEPMediandbE <= 73.5)]
#brewer.pal(7, "Set2")
#"#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
# OR        NIPW       AIPW       tOR     tNIPW      tAIPW      dir

TypeE <- c(rep("OR", 4),
           rep("AIPW", 14),
           rep("tOR", 2),
           rep("tAIPW", 12))

dfE <- data.frame(method = c(1:32),
                  ARB = RMSEPMediandbEs,
                  Label = LabelsE,
                  Type = TypeE)


gE <- ggplot(dfE, aes(x = method,
                      y = RMSEPMediandbEs,
                      col = TypeE)) +
  labs(y = "RRMSE", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(36, 84)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(10,88,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F"))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsE), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gE +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gE


#plot(1:127, RMSEPMediandbM, ylim = c(35, 80))
#abline(72.5, 0)

RMSEPMediandbMs <- RMSEPMediandbM[which(RMSEPMediandbM <= 72.5)]
MethodsM <- Methods[(which(RMSEPMediandbM <= 72.5))]
LabelsM <-  Labels[which(RMSEPMediandbM <= 72.5)]



TypeM <- c(rep("OR", 4),
           rep("AIPW", 12),
           rep("tOR", 2),
           rep("tAIPW", 12))

dfM <- data.frame(method = c(1:30),
                  ARB = RMSEPMediandbMs,
                  Label = LabelsM,
                  Type = TypeM)


gM <- ggplot(dfM, aes(x = method,
                      y = RMSEPMediandbMs,
                      col = TypeM)) +
  labs(y = "RRMSE", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(36, 84)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(10,88,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F"))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsM), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gM +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gM


#plot(1:127, RMSEPMediandbR, ylim = c(35, 80))
#abline(70.5, 0)


#R
RMSEPMediandbRs <- RMSEPMediandbR[which(RMSEPMediandbR <= 70.5)]
MethodsR <- Methods[(which(RMSEPMediandbR <= 70.5))]
LabelsR <-  Labels[which(RMSEPMediandbR <= 70.5)]


TypeR <- c(rep("OR", 4),
           rep("AIPW", 9),
           rep("tOR", 2),
           rep("tAIPW", 12))

dfR <- data.frame(method = c(1:27),
                  ARB = RMSEPMediandbRs,
                  Label = LabelsR,
                  Type = TypeR)



gR <- ggplot(dfR, aes(x = method,
                      y = RMSEPMediandbRs,
                      col = TypeR)) +
  labs(y = "RRMSE", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(36, 84)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(10,88,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F" ))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsR), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gR +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gR



#plot(1:127, RMSEPMediandbRtune, ylim = c(35, 80))
#abline(71.5, 0)


#Rtune
RMSEPMediandbRtunes <- RMSEPMediandbRtune[which(RMSEPMediandbRtune <= 71.5)]
MethodsRtune <- Methods[(which(RMSEPMediandbRtune <= 71.5))]
LabelsRtune <-  Labels[which(RMSEPMediandbRtune <= 71.5)]


TypeRtune <- c(rep("OR", 4),
               rep("AIPW", 13),
               rep("tOR", 2),
               rep("tAIPW", 12))

dfRtune <- data.frame(method = c(1:31),
                      ARB = RMSEPMediandbRtunes,
                      Label = LabelsRtune,
                      Type = TypeRtune)


gRtune <- ggplot(dfRtune, aes(x = method,
                              y = RMSEPMediandbRtunes,
                              col = TypeRtune)) +
  labs(y = "RRMSE", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(36, 84)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(10,88,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F" ))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsRtune), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gRtune +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gRtune

#plot(1:127, RMSEPMediandbX, ylim = c(35, 80))
#abline(74, 0)


# X
RMSEPMediandbXs <- RMSEPMediandbX[which(RMSEPMediandbX <= 74)]
MethodsX <- Methods[(which(RMSEPMediandbX <= 74))]
LabelsX <-  Labels[which(RMSEPMediandbX <= 74)]


TypeX <-  c(rep("OR", 4),
            rep("AIPW", 20),
            rep("tOR", 2),
            rep("tAIPW", 12))

dfX <- data.frame(method = c(1:38),
                  ARB = RMSEPMediandbXs,
                  Label = LabelsX,
                  Type = TypeX)



gX <- ggplot(dfX, aes(x = method,
                      y = RMSEPMediandbXs,
                      col = TypeX)) +
  labs(y = "RRMSE", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(36, 84)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(10,88,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F" ))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsX), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gX +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gX

#plot(1:127, RMSEPMediandbXtune, ylim = c(35, 80))
#abline(73, 0)

# X tune
RMSEPMediandbXtunes <- RMSEPMediandbXtune[which(RMSEPMediandbXtune <= 73)]
MethodsXtune <- Methods[(which(RMSEPMediandbXtune <= 73))]
LabelsXtune <-  Labels[which(RMSEPMediandbXtune <= 73)]


TypeXtune <-  c(rep("OR", 4),
                rep("AIPW", 12),
                rep("tOR", 2),
                rep("tAIPW", 12))

dfXtune <- data.frame(method = c(1:30),
                      ARB = RMSEPMediandbXtunes,
                      Label = LabelsXtune,
                      Type = TypeXtune)

gXtune <- ggplot(dfXtune, aes(x = method,
                              y = RMSEPMediandbXtunes,
                              col = TypeXtune)) +
  labs(y = "RRMSE", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(36, 84)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(10,89,5)) +
  scale_color_manual(name = " ",
                     values = c("OR" ="#66C2A5",
                                "AIPW"="#8DA0CB",
                                "tOR"="#E78AC3",
                                "tAIPW" = "#FFD92F" ))+
  #  annotate('text', x = 5, y = 7.145041, label = 'C') +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  geom_text_repel(aes(label = LabelsXtune), show.legend=F, size = 5)+
  guides(color = guide_legend(nrow = 2, byrow = T))

gXtune +   ggtitle("m = 25") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(size=14, face="bold")
  )

gXtune


plot_grid(gE, gM, gR, gRtune, gX, gXtune,
          labels = c("(a)", "(b)", "(c)",
                     "(d)", "(e)", "(f)"), nrow = 2)

