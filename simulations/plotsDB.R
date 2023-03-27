#############################################################################
# Model based --------------------------------------------------------------
#############################################################################
df25 <- data.frame(subpop = c(1:25),
                   tau_true = tau_true)

df50 <- data.frame(subpop = c(1:50),
                   tau_true = tau_true)

df75 <- data.frame(subpop = c(1:75),
                   tau_true = tau_true)



g25 <- ggplot(df25, aes(x = subpop,
                        y = tau_true)) +
  labs(y = "ATE", x = "Subpopulation")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(0, 17)) +
  scale_y_continuous(breaks=seq(0,88,5)) +
  scale_x_continuous(breaks=seq(0,75,5)) +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))

g50 <- ggplot(df50, aes(x = subpop,
                        y = tau_true)) +
  labs(y = "ATE", x = "Subpopulation")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(0, 17)) +
  scale_y_continuous(breaks=seq(0,88,5)) +
  scale_x_continuous(breaks=seq(0,75,5)) +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))

g75 <- ggplot(df75, aes(x = subpop,
                        y = tau_true)) +
  labs(y = "ATE", x = "Subpopulation")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(0, 17)) +
  scale_y_continuous(breaks=seq(0,88,5)) +
  scale_x_continuous(breaks=seq(0,75,5)) +
  theme(legend.position = c(0.5, 0.9),
        legend.direction="horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))

plot_grid(g25, g50, g75,
          labels = c("(a)", "(b)", "(c)"), nrow = 1)

#################################################################################
# Design based ------------------------------------------------------------------
#################################################################################
data_survey <- read.csv2("data_filtered_survey.csv")

# Read and pre_process data

data_survey_inc <- data_survey[, -c(17)]
data_survey_inc <- data_survey_inc %>% filter(eq_house_disp_inc > 0)
y = log(data_survey_inc$eq_house_disp_inc)
data_survey_inc$y <- y

# Standardize count variables
data_survey_inc$age <- scale(data_survey_inc$age,
                             center = T, scale = T)

data_survey_inc$house_size <- scale(data_survey_inc$house_size,
                                    center = T, scale = T)

data_pop <- data_survey_inc %>%
  rename(A  = type_contract,
         group  = province) %>%
  select(-eq_house_disp_inc)

####################################################################################
# Estimate the propensity score
####################################################################################
# Drop: single, edu0

#formula_p_score = A ~ X1 + (1|group)
formula_p_score = A ~ sex + nationality + age + house_size + married +
  separated + widowed + divorced  + edu1 + edu2 + edu3 + (1|group)

# EBLUP

obj_p_score_EBLUP <- list(data_p_score = data_pop)
class(obj_p_score_EBLUP) <- "EBLUP"

ps_hat_EBLUP <-  p_score(obj_p_score = obj_p_score_EBLUP,
                         model_formula = formula_p_score)
data_popE <- data_pop
data_popE$p_score <- ps_hat_EBLUP

tauE <- calculate_tau(list(data_popE), type_tau = "H")
tau_trueE <- tauE[[1]]$tau

# MQ

obj_p_score_MQ <- list(data_p_score = data_pop)
class(obj_p_score_MQ) <- "MQ"

ps_hat_MQ <-  p_score(obj_p_score = obj_p_score_MQ,
                      model_formula = formula_p_score)

data_popM <- data_pop
data_popM$p_score <- ps_hat_MQ

tauM <- calculate_tau(list(data_popM), type_tau = "H")
tau_trueM <- tauM[[1]]$tau

# RF

obj_p_score_RF <- list(data_p_score = data_pop)
class(obj_p_score_RF) <- "RF"

# No tune
ps_hat_RF <-  p_score(obj_p_score = obj_p_score_RF,
                      model_formula = formula_p_score,
                      tune_RF = FALSE)

data_popR <- data_pop
data_popR$p_score <- ps_hat_RF

tauR <- calculate_tau(list(data_popR), type_tau = "H")
tau_trueR <- tauR[[1]]$tau

# Tune
ps_hat_RFt <-  p_score(obj_p_score = obj_p_score_RF,
                      model_formula = formula_p_score,
                      tune_RF = TRUE)

data_popRt <- data_pop
data_popRt$p_score <- ps_hat_RFt

tauRt <- calculate_tau(list(data_popRt), type_tau = "H")
tau_trueRt <- tauRt[[1]]$tau

# XGB
obj_p_score_XGB <- list(data_p_score = data_pop)
class(obj_p_score_XGB) <- "XGB"

# No tune
ps_hat_XGB <-  p_score(obj_p_score = obj_p_score_XGB,
                       model_formula = formula_p_score,
                       xgboost_params = list(CV_XGB = FALSE,
                                             nfolds = 5,
                                             nrounds = 50))
data_popX <- data_pop
data_popX$p_score <- ps_hat_XGB

tauX <- calculate_tau(list(data_popX), type_tau = "H")
tau_trueX <- tauX[[1]]$tau

# Tune
ps_hat_XGBt <-  p_score(obj_p_score = obj_p_score_XGB,
                       model_formula = formula_p_score,
                       xgboost_params = list(CV_XGB = TRUE,
                                             nfolds = 5,
                                             nrounds = 50))
data_popXt <- data_pop
data_popXt$p_score <- ps_hat_XGBt

tauXt <- calculate_tau(list(data_popXt), type_tau = "H")
tau_trueXt <- tauXt[[1]]$tau

Method <- c(rep("E", 41),
            rep("M", 41),
            rep("R", 41),
            rep("X", 41),
            rep("Rt", 41),
            rep("Xt", 41))

tau <- c(tau_trueE, tau_trueM,
         tau_trueR, tau_trueX,
         tau_trueRt, tau_trueXt)

dftr <- data.frame(tau_trueE, Region, Province)

dfL <- dftr[dftr$Region == "L", ]
dfLsort <- dfL[with(dfL, order(tau_trueE)), ]

dfT <- dftr[dftr$Region == "T", ]
dfTsort <- dfT[with(dfT, order(tau_trueE)), ]

dfU <- dftr[dftr$Region == "U", ]
dfUsort <- dfU[with(dfU, order(tau_trueE)), ]

dfM <- dftr[dftr$Region == "M", ]
dfMsort <- dfM[with(dfM, order(tau_trueE)), ]

dfC <- dftr[dftr$Region == "C", ]
dfCsort <- dfC[with(dfC, order(tau_trueE)), ]

dfS <- dftr[dftr$Region == "S", ]
dfSsort <- dfS[with(dfS, order(tau_trueE)), ]


ProvSort <- c(dfLsort$Province, dfTsort$Province,
              dfUsort$Province, dfMsort$Province,
              dfCsort$Province, dfSsort$Province)




#C  L  M  S  T  U
#5 11  4  9 10  2
regionGeo <- numeric(41)
regionGeo[which(Region == "L")] <- "1L"
regionGeo[which(Region == "T")] <- "2T"
regionGeo[which(Region == "U")] <- "3U"
regionGeo[which(Region == "M")] <- "4M"
regionGeo[which(Region == "C")] <- "5C"
regionGeo[which(Region == "S")] <- "6S"

dfM <- data.frame(subpopulation  = rep(Province, 6),
                  tau = tau,
                  Method = Method,
                  region = rep(Region, 6),
                  regiongeo = rep(regionGeo, 6))


#dfMs <- dfM[with(dfM, order(regiongeo, tau,subpopulation)), ]
#brewer.pal(7, "Set2")
#"#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
# OR        NIPW       AIPW       tOR     tNIPW      tAIPW      dir

values = c("L" ="#66C2A5",
           "T"="#FC8D62",
           "U"="#8DA0CB",
           "M" = "#E78AC3",
           "C" = "#A6D854",
           "S" = "#FFD92F")


                    values = c("E" ="#66C2A5",
                                "M"="#FC8D62",
                                "R"="#8DA0CB",
                                "X" = "#E78AC3",
                                "Rt" = "#A6D854",
                                "Xt" = "#FFD92F"))+

gDB <- ggplot(dfM, aes(x = subpopulation,
                       y = tau,
                       col = region)) +
  labs(y = "ATE", x = "Province")+
  geom_point(aes(shape = Method), size = 3) + theme_bw() + coord_cartesian(ylim=c(-0.75, 1.8)) +
  scale_y_continuous(breaks = seq(-1,1.5,0.5)) +
  scale_x_discrete(limits = ProvSort)+
  scale_color_manual(name = " Region",
                     values = c("L" ="#66C2A5",
                                "T"="#FC8D62",
                                "U"="#8DA0CB",
                                "M" = "#E78AC3",
                                "C" = "#A6D854",
                                "S" = "#FFD92F"))+
  scale_shape_manual(name = "Method",
                     values = c("E" = 8,
                                "M"= 16,
                                "R"= 17,
                                "X" = 15,
                                "Rt" = 2,
                                "Xt" = 0))+
  theme(legend.direction = "horizontal",
        legend.position = c(0.5, 0.9),
        legend.box = "horizontal",
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12, face="bold"))+
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90, vjust =  0.35, hjust = 0.7)) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=24,face="bold"),
        legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 24)) +
  theme(axis.title.x = element_text(size = 24)) +
  theme(axis.title.x = element_text(vjust = 1)) +
  guides(color = guide_legend(nrow = 1, byrow = T, title.position = "top")) +
  guides(shape = guide_legend(nrow = 1, byrow = T, title.position = "top"))

gDB




gM <- ggplot(dfM, aes(x = method,
                      y = ARBPMediandbMs,
                      col = TypeM)) +
  labs(y = "ARB", x = "Method")+theme_bw()+
  geom_point(size = 2) + theme_bw() + coord_cartesian(ylim=c(10, 60)) +
  scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=seq(10,61,5)) +
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
