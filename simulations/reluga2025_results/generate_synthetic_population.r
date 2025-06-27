# rm(list=ls())
# .rs.restartR()

library("MASS")

setwd("/user/home/ze23696/Comp/causalSAE")
devtools::load_all()

setwd("/user/home/ze23696/Comp/SuperLearner")
devtools::load_all()

setwd("/user/home/ze23696/Comp/BinaryMQ")
devtools::load_all()

source("/user/home/ze23696/Comp/QRLM.r")

set.seed(10)

m = 41
Nii = 1000
Ni = rep(Nii, m)
N = sum(Ni)

var_norm <- matrix(c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 0.5,
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5,
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5, 0.5,
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 0.5,
                     0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1),
                   nrow = 10, ncol = 10, byrow = T)

## Retrieve group indicator ---------------------------------------------------
group = rep(1:m, times = Ni)

X <- generate_X(
  n = N,
  p = 10,
  covariance_norm = var_norm,
  cov_type = "norm",
  seed = 3
)

untreat1 <- list()
untreat2 <- list()

treat1 <- list()
treat2 <- list()

for (i in 1:m) {
  untreat1[[i]] <- matrix(rep(rnorm(10, 0, 3), Nii), nrow = Nii, ncol = 10, byrow = T)
  untreat2[[i]] <- 0.1 *  matrix(rep(rnorm(10, 0, 3), Nii), nrow = Nii, ncol = 10, byrow = T)

  treat1[[i]] <- untreat1[[i]] + 2 * matrix(rep(rnorm(10, 4, 3), Nii), nrow = Nii, ncol = 10, byrow = T)
  treat2[[i]] <- untreat2[[i]] + 0.1 * matrix(rep(rnorm(10, 4, 3), Nii), nrow = Nii, ncol = 10, byrow = T)
}

coef_01 <- do.call("rbind", untreat1)
coef_02 <- do.call("rbind", untreat2)

coef_11 <- do.call("rbind", treat1)
coef_12 <- do.call("rbind", treat2)

ref0 <- runif(m, 1, 2)
ref1 <- runif(m, 2, 3)

c0 <- rep(ref0, each = Nii)
c1 <- rep(ref1, each = Nii)

set.seed(1)
X11 <- rep(rnorm(m, 0, 2.6), each = Nii)

set.seed(39)
y0ne <- rowMeans(coef_01 * X + exp(coef_02 * X)) + X11 + c0
y0 <- log(abs(y0ne + rnorm(Nii * m, 0, mean(y0ne))))

y1ne <- rowMeans(coef_11 * X + exp(coef_12 * X)) + X11 + c1
y1 <- log(abs(y1ne + rnorm(Nii * m, 0, mean(y1ne))))

var(y0ne)/var(rnorm(Nii * m, 0, mean(y0ne)))
var(y1ne)/var(rnorm(Nii * m, 0, mean(y1ne)))

tau_treat <- aggregate(y1, list(group), FUN = mean)$x
tau_untreat <- aggregate(y0, list(group), FUN = mean)$x
tau_true = tau_treat - tau_untreat
range(tau_true);max(tau_true)-min(tau_true)

var(y0);var(y1)
#tau_true
plot(1:41, tau_true)

X <- cbind(X, X11)

# Propensity score ---------------------------------------
set.seed(31)
X_g <- data.frame(X, group)
p_scoreM <- matrix(0, ncol = m, nrow = Nii)
group_names <- as.data.frame(table(group))$group
for (i in 1:m) {
  X_gg <- X_g[X_g$group == group_names[i], ]
  coef_p_score = c(rnorm(9, 2, 6), rnorm(2, 0, 0.25))
  X_group <- X_gg[, -12]
  Xreg_p_score = abs(rowMeans(coef_p_score * X_group))
  exp_p_score = exp(Xreg_p_score)
  p_scoreM[, i] = exp_p_score * (1 + exp_p_score)^(-1)
}

p_score <- c(p_scoreM)
#
A <- rbinom(N, 1, p_score)
A_group <- aggregate(A, list(group), FUN = mean)$x
range(A_group) #0.636 0.878
mean(A_group) #0.7727073

#range(A_group)
#0.6056504 0.8991213
#mean(A_group)
#0.8030633
#A_census <-  c(0.8915738 0.8950139 0.8315111 0.8886850 0.8837953 0.8675712 0.8756629
#0.8771440 0.8642530 0.8340261 0.8410248 0.8393271 0.8229929 0.8486160
#0.8483193 0.8427683 0.8603392 0.8045844 0.8522335 0.8487627 0.8165947
#0.7806220 0.8298834 0.8467180 0.7693839 0.8055868 0.7813291 0.7247728
#0.7289499 0.6994361 0.7330159 0.7087267 0.6441270 0.6514668 0.6487688
#0.7055741 0.6056504 0.6886817 0.8991213 0.8901192 0.8488613)

colnames(X) <-  paste0("X", 1:ncol(X))
y = A * y1 + (1 - A) * y0

populations <- data.frame(X, A, group, p_score, y)

#Super learner 

# Fit propensity score
y_train_ps <- populations$A
x_train_ps <- populations[, -c(12, 14:15)]
D = Sys.time()
SLps = SuperLearner(Y = y_train_ps,
                    X = x_train_ps,
                    id = populations$group,
                    family = binomial(),
                    SL.library = c("SL.glm",
                                   "SL.glmmb", 
                                   "SL.mqb","SL.mqchb",
                                   "SL.grf", "SL.grf_nc", "SL.grf_t", "SL.grf_nct",
                                   "SL.xgboost", "SL.xgboost_t"),
                    cvControl = list(V = 5))
B = Sys.time() #33.41462 mins
ps_hat_S <- SLps$SL.predict[, 1]


C = Sys.time()
# PS methods ------------------------------------------------------------
# GLM


ps_fit_L <- glm(A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +X11,
                data = populations, 
                family = binomial)
ps_hat_L <- ps_fit_L$fitted.values



# EBLUP -----------------------------------------------------------------

obj_p_score_EBLUP <- list(data_p_score = populations[, -c(14:15)])
class(obj_p_score_EBLUP) <- "EBLUP"

ps_hat_E <-  p_score(obj_p_score = obj_p_score_EBLUP,
                     model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +X11+ (1|group))


# MQ -------------------------------------------------------------------

obj_p_score_MQ <- list(data_p_score = populations[, -c(14:15)])
class(obj_p_score_MQ) <- "MQ"

ps_hat_M <-  p_score(obj_p_score = obj_p_score_MQ,
                     model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10+X11 + (1|group))

# Mch ----------------------------------------------------------------
#Q = sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))  
Q = c(0.25, 0.30, 0.4, seq(from = 0.45, to = 0.55, by = 0.005), 0.60, 0.65, 0.75)

M_p0 <- glm.mq.binom(y = populations$A,
                     x = cbind(1, as.matrix(x_train_ps[, -c(12)])),
                     q = Q,
                     k = 1.6,
                     maxit = 100) 

tmp.scores <- QSCORE(populations$A, 
                     M_p0$fitted, 
                     Q)
scores = (cbind(populations$group, 
                tmp.scores))

MQ_p = rep(0, m)
for (i in 1:m){
  MQ_p[i] = (mean(scores[, 2][scores[, 1]==group_names[i]]))
}


M_p <- glm.mq.binom(y = populations$A,
                    x = cbind(1, as.matrix(x_train_ps[, -c(12)])),
                    q = MQ_p,
                    k = 1.6,
                    maxit = 100) 

ps_hat_Mch <- NULL

# Check this against the old example
for (i in 1:m){
  ps_hat_Mch <- c(ps_hat_Mch, exp(cbind(1, as.matrix(populations[populations$group == group_names[i], c(1:11)])) %*% M_p$coefficients[, i])/
                    (1 + exp(cbind(1, as.matrix(populations[populations$group == group_names[i], c(1:11)])) %*% M_p$coefficients[, i])))  
}



# RF ----------------------------------------------------------------
# No cluster no tuning 
obj_p_score_RF <- list(data_p_score = populations[, -c(14:15)])
class(obj_p_score_RF) <- "RF"

ps_hat_R <-  p_score(obj_p_score = obj_p_score_RF,
                     model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + (1|group),
                     tune_RF = FALSE, 
                     clust_RF = FALSE)
time_e = numeric(10)


# Cluster, no tuning
ps_hat_Rc <-  p_score(obj_p_score = obj_p_score_RF,
                      model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ (1|group),
                      tune_RF = FALSE, 
                      clust_RF = TRUE)


# Cluster, no tuning
ps_hat_Rt <-  p_score(obj_p_score = obj_p_score_RF,
                      model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ (1|group),
                      tune_RF = TRUE, 
                      clust_RF = FALSE)

# Cluster, tuning
ps_hat_Rct <-  p_score(obj_p_score = obj_p_score_RF,
                       model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ (1|group),
                       tune_RF = TRUE, 
                       clust_RF = TRUE)


# XGB ----------------------------------------------------------------------
for (i in 1: 10) {
  print(i)
  start_time =  Sys.time()
obj_p_score_XGB <- list(data_p_score = populations[, -c(14:15)])
class(obj_p_score_XGB) <- "XGB"

ps_hat_X <-  p_score(obj_p_score = obj_p_score_XGB,
                     model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ (1|group),
                     xgboost_params = list(CV_XGB = FALSE,
                                           nfolds = 5,
                                           nrounds = 100))
stop_time =  Sys.time()
time_e[i] <- stop_time - start_time

}
mean(time_e)
# XGBt ----------------------------------------------------------------------

ps_hat_Xt <-  p_score(obj_p_score = obj_p_score_XGB,
                      model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11+ (1|group),
                      xgboost_params = list(CV_XGB = TRUE,
                                            nfolds = 5,
                                            nrounds = 100))


populations$ps_hat_S <- ps_hat_S

populations$ps_hat_L <- ps_hat_L
populations$ps_hat_E <- ps_hat_E
populations$ps_hat_M <- ps_hat_M
populations$ps_hat_Mch <- ps_hat_Mch

populations$ps_hat_R <- ps_hat_R
populations$ps_hat_Rc <- ps_hat_Rc
populations$ps_hat_Rt <- ps_hat_Rt
populations$ps_hat_Rct <- ps_hat_Rct

populations$ps_hat_X <- ps_hat_X
populations$ps_hat_Xt <- ps_hat_Xt


write.csv(populations,
          file = 'population_sae4.csv', row.names = FALSE)

D = Sys.time() #Time difference of 59.26068 mins

# S ----------------------------------------------------
IPWS <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_S
)

tauS <-
  as.data.frame(calculate_tau(IPWS, type_tau = "H"))

Sm <- mean((tauS$tau - tau_true)^2)



# E ----------------------------------------------------------
IPWE <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_E
)

tauE <-
  as.data.frame(calculate_tau(IPWE, type_tau = "H"))

Em <- mean((tauE$tau - tau_true)^2)


# L ----------------------------------------------------------
IPWL <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_L
)

tauL <-
  as.data.frame(calculate_tau(IPWL, type_tau = "H"))

Lm <- mean((tauL$tau - tau_true)^2)


# M ----------------------------------------------------------
IPWM <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_M
)

tauM <-
  as.data.frame(calculate_tau(IPWM, type_tau = "H"))

Mm <- mean((tauM$tau - tau_true)^2)


# Mch ----------------------------------------------------------
IPWMch <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_Mch
)

tauMch <-
  as.data.frame(calculate_tau(IPWMch, type_tau = "H"))

Mchm <- mean((tauMch$tau - tau_true)^2)



IPWMchS <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = (populations$ps_hat_Mch + populations$ps_hat_S + populations$ps_hat_E)/3
)

tauMchS <-
  as.data.frame(calculate_tau(IPWMchS, type_tau = "H"))

Mchm <- mean((tauMchS$tau - tau_true)^2)


# R ----------------------------------------------------------
IPWR <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_R
)

tauR <-
  as.data.frame(calculate_tau(IPWR, type_tau = "H"))

Rm <- mean((tauR$tau - tau_true)^2)


# Rc ----------------------------------------------------------
IPWRc <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_Rc
)

tauRc <-
  as.data.frame(calculate_tau(IPWRc, type_tau = "H"))

Rcm <- mean((tauRc$tau - tau_true)^2)

# Rt ----------------------------------------------------------
IPWRt <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_Rt
)

tauRt <-
  as.data.frame(calculate_tau(IPWRt, type_tau = "H"))

Rtm <- mean((tauRt$tau - tau_true)^2)

# Rct ----------------------------------------------------------
IPWRct <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_Rct
)

tauRct <-
  as.data.frame(calculate_tau(IPWRct, type_tau = "H"))
Rctm <- mean((tauRct$tau - tau_true)^2)

# X ----------------------------------------------------------
IPWX <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_X
)

tauX <-
  as.data.frame(calculate_tau(IPWX, type_tau = "H"))

Xm <- mean((tauX$tau - tau_true)^2)

# Xt ----------------------------------------------------------
IPWXt <- data.frame(
  y = populations$y,
  A = populations$A,
  group = populations$group,
  p_score  = populations$ps_hat_Xt
)

tauXt <-
  as.data.frame(calculate_tau(IPWXt, type_tau = "H"))

Xtm <- mean((tauXt$tau - tau_true)^2)


tau_df <- data.frame(tau = tau_true, 
                     tau_treat = tau_treat, 
                     tau_untreat = tau_untreat, 
                     
                     tau_S  = tauS$tau, 
                     
                     tau_M  = tauM$tau, 
                     tau_Mch  = tauMch$tau, 
                     tau_E  = tauE$tau, 
                     tau_L  = tauL$tau, 
                     
                     tau_R  = tauR$tau, 
                     tau_Rc  = tauRc$tau, 
                     tau_Rt  = tauRt$tau, 
                     tau_Rct  = tauRct$tau, 
                     
                     tau_X  = tauX$tau, 
                     tau_Xt  = tauXt$tau, 
                     
                     group = group_names)

write.csv(tau_df,
          file = 'tau_sae4.csv', row.names = FALSE)


