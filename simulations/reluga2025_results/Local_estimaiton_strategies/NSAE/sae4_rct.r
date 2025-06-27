library("MASS")


setwd("/user/work/ze23696/Comp/causalSAE")
devtools::load_all()


setwd("/user/work/ze23696/Comp/SuperLearner")
devtools::load_all()

setwd("/user/work/ze23696/Comp/BinaryMQ")
devtools::load_all()


source("/user/work/ze23696/Comp/QRLM.r")



populations <- read.csv('/user/work/ze23696/Comp/population_sae4.csv')
taus <- read.csv('/user/work/ze23696/Comp/tau_sae4.csv')
tau_true <- taus$tau
tau_treat <- taus$tau_treat
tau_untreat <- taus$tau_untreat



m = 41

Nii = 1000
Ni = rep(Nii, m)
N = sum(Ni)
group_names <- as.data.frame(table(populations$group))$Var1
Nc = as.numeric(table(populations$group[populations$A == 0]))
Nt = as.numeric(table(populations$group[populations$A == 1]))

frac_nt = 0.02
nt <- round(frac_nt * Nt)
frac_nc0 <- c(runif(25, 0.01, 0.5), runif(10, 0.51, 1), runif(6, 0.9, 1))
nc <- ceiling(frac_nc0 * nt)
frac_nc <- nc/Nc

SN = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(SN * 2022)
subpopulation <- sample_subpopulations(populations,
                                       frac_nc = frac_nc, frac_nt = frac_nt,
                                       seed = set.seed(SN * 2022))
data_sample <- data.frame(populations[subpopulation, c(1:15)])
data_out_of_sample <- populations[-subpopulation, c(1:15)]

data_sample_ps <-  data.frame(populations[subpopulation, -c(1:15)])
data_out_of_sample_ps <-  data.frame(populations[-subpopulation, -c(1:15)])

ni = as.numeric(table(data_sample$group))
nc = as.numeric(table(data_sample$group[data_sample$A == 0]))
nt = as.numeric(table(data_sample$group[data_sample$A == 1]))
frac_nN <- dim(data_sample)[1]/dim(populations)[1]


########################################################################################
########################################################################################
# Super Learner ------------------------------------------------------------------

#y_train <- data_sample$y
#x_train <- data_sample[, c(1:11, 13)]

names_sample <- names(data_sample)
#names_out_of_sample <- names(data_out_of_sample)
data_full_prelim <- rbind(data_sample[, !names_sample %in% "y"],
                          data_out_of_sample[, !names_sample %in% "y"])
data_out_of_sampley <- data_out_of_sample[, !names_sample %in% "y"]

data_full <- data_full_prelim[, c(1:11, 13)]

# Imputation --------------------------------------------------------------------

################################################################
# Linear regression -------------------------------------------
################################################################
impute_Rct <- impute_y(model_formula = y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + A + (1|group),
                       data_sample,
                       data_out_of_sample,
                       method = "RF",
                       tune_RF = TRUE, 
                       clust_RF  = TRUE)

y_full_imputeRct <- impute_Rct$y_full_imputed

data_full_i <- data_full_prelim
data_full_i$y <- y_full_imputeRct

########################################################################################
########################################################################################

# OR
mu1_y_L = NULL
mu1_y_M = NULL

mu1_y_R = NULL
mu1_y_Rt = NULL

mu1_y_X = NULL
mu1_y_Xt = NULL

mu1_y_S = NULL

mu0_y_L = NULL 
mu0_y_M = NULL

mu0_y_R = NULL
mu0_y_Rt = NULL

mu0_y_X = NULL
mu0_y_Xt = NULL

mu0_y_S = NULL


# PS
ps_hat_L = NULL
ps_hat_M = NULL 

ps_hat_R = NULL
ps_hat_Rt = NULL

ps_hat_X = NULL
ps_hat_Xt = NULL

ps_hat_S = NULL

A_g = NULL
y_g = NULL
group_g <- NULL


start_time = Sys.time()
for (i in 1:m) {
  
  print(i)
  data_full_g <- data_full_i[data_full_i$group == group_names[i], ]
  
  y_g <- c(y_g, data_full_g$y)
  A_g <- c(A_g, data_full_g$A)
  group_g <- c(group_g, data_full_g$group)
  # OR --------------------------------------------------------  
  # Super Learner ------------------------------------------------------------------
  y_train <- data_full_g$y
  x_train <- data_full_g[, c(1:10)]
  
  
  # Controls -------------------------------------
  data_sample0 = data_full_g[data_full_g$A == 0, ]
  obj_fit_y0 <- list(data_sample = data_sample0)
  
  y_train0 <- data_sample0$y
  x_train0 <- data_sample0[, c(1:10)]
  
  SL0 = SuperLearner(Y = y_train0,
                     X = x_train0,
                     family = gaussian(),
                     newX = x_train,
                     SL.library = c("SL.lm",
                                    "SL.mq",  
                                    "SL.grf_nc", "SL.grf_nct",
                                    "SL.xgboost", "SL.xgboost_t"),
                     cvControl = list(V = 5))
  mu0_y_S <- c(mu0_y_S, SL0$SL.predict[, 1])
  
  
  # Treat -----------------------------------------
  data_sample1 = data_full_g[data_full_g$A == 1, ]
  obj_fit_y1 <- list(data_sample = data_sample1)
  
  
  y_train1 <- data_sample1$y
  x_train1 <- data_sample1[, c(1:10)]
  
  SL1 = SuperLearner(Y = y_train1,
                     X = x_train1,
                     family = gaussian(),
                     newX = x_train,
                     SL.library = c("SL.lm",
                                    "SL.mq",  
                                    "SL.grf_nc", "SL.grf_nct",
                                    "SL.xgboost", "SL.xgboost_t"),
                     cvControl = list(V = 5))
  
  mu1_y_S <- c(mu1_y_S, SL1$SL.predict[,1])
  
  ################################################################
  # Linear regression -------------------------------------------
  ################################################################
  
  OR1_L <- lm(y ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = data_sample1)
  mu1_y_L <- c(mu1_y_L, predict(OR1_L, newdata = x_train))
  
  OR0_L <- lm(y ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data = data_sample0)
  mu0_y_L <- c(mu0_y_L, predict(OR0_L, newdata = x_train))
  
  ###########################################################
  # MQ -----------------------------------------------------------------------------
  #############################################################
  model_formula_OR <- y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group)
  
  
  class(obj_fit_y0) <- "MQ"
  mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                                 model_formula = model_formula_OR)
  
  
  OR0_M <- fit_y(mutated_obj0,
                 type_model =  "continuous")
  
  mu0_y_M <- c(mu0_y_M, unname(unlist(predict(object = OR0_M$outcome_fit,
                                              newdata =  x_train,
                                              allow.new.levels = TRUE))))
  
  
  # Treated --------------------------------------
  class(obj_fit_y1) <- "MQ"
  mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                                 model_formula = model_formula_OR)
  
  
  OR1_M <- fit_y(mutated_obj1,
                 type_model = "continuous")
  
  mu1_y_M <-  c(mu1_y_M, unname(unlist(predict(object = OR1_M$outcome_fit,
                                               newdata =  x_train,
                                               allow.new.levels = TRUE))))
  
  ###################################
  # No clustering, no tuning
  #############################################
  obj_fit_y0 <- list(data_sample = data_sample0)
  class(obj_fit_y0) <- "RF"
  mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                                 model_formula = model_formula_OR)
  
  
  OR0_R <- fit_y(mutated_obj0,
                 type_model = "continuous", 
                 tune_RF = FALSE, 
                 clust_RF = FALSE)
  
  mu0_y_R <-  c(mu0_y_R, unname(unlist(predict(object = OR0_R$outcome_fit,
                                               newdata = x_train,
                                               allow.new.levels = TRUE))))
  
  
  # Treated --------------------------------------
  class(obj_fit_y1) <- "RF"
  mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                                 model_formula = model_formula_OR)
  
  
  OR1_R <- fit_y(mutated_obj1,
                 type_model = "continuous", 
                 tune_RF = FALSE,
                 clust_RF = FALSE)
  
  mu1_y_R <-  c(mu1_y_R, unname(unlist(predict(object = OR1_R$outcome_fit,
                                               newdata = x_train,
                                               allow.new.levels = TRUE))))
  
  #############################################
  # RF no clustering, tunning --------------------------------------------------------------------------
  #######################################
  OR0_Rt <- fit_y(mutated_obj0,
                  type_model = "continuous", 
                  tune_RF = TRUE, 
                  clust_RF = FALSE)
  
  
  
  mu0_y_Rt <-  c(mu0_y_Rt, unname(unlist(predict(object = OR0_Rt$outcome_fit,
                                                 newdata = x_train,
                                                 allow.new.levels = TRUE))))
  
  
  # Treated --------------------------------------
  OR1_Rt <- fit_y(mutated_obj1,
                  type_model = "continuous", 
                  tune_RF = TRUE, 
                  clust_RF = FALSE)
  
  mu1_y_Rt <-  c(mu1_y_Rt, unname(unlist(predict(object = OR1_Rt$outcome_fit,
                                                 newdata = x_train,
                                                 allow.new.levels = TRUE))))
  
  
  
  
  ######################################################
  # XGboost
  #####################################################
  # X --------------------------------------------------------------------------------
  class(obj_fit_y0) <- "XGB"
  mutated_obj0 <- mutate_obj_fit(obj_fit_y = obj_fit_y0,
                                 model_formula = model_formula_OR)
  
  
  OR0_X <- fit_y(mutated_obj0,
                 type_model = "continuous", 
                 xgboost_params = list(CV_XGB = FALSE,
                                       nfolds = 5,
                                       nrounds = 100))
  
  
  
  mu0_y_X <-  c(mu0_y_X, unname(unlist(predict(object = OR0_X$outcome_fit,
                                               newdata = as.matrix(x_train),
                                               allow.new.levels = TRUE))))
  
  
  # Treated --------------------------------------
  class(obj_fit_y1) <- "XGB"
  mutated_obj1 <- mutate_obj_fit(obj_fit_y = obj_fit_y1,
                                 model_formula = model_formula_OR)
  
  
  OR1_X <- fit_y(mutated_obj1,
                 type_model = "continuous", 
                 xgboost_params = list(CV_XGB = FALSE,
                                       nfolds = 5,
                                       nrounds = 100))
  
  mu1_y_X <-  c(mu1_y_X, unname(unlist(predict(object = OR1_X$outcome_fit,
                                               newdata = as.matrix(x_train),
                                               allow.new.levels = TRUE))))
  
  # X t --------------------------------------------------------------------------------
  OR0_Xt <- fit_y(mutated_obj0,
                  type_model = "continuous", 
                  xgboost_params = list(CV_XGB = TRUE,
                                        nfolds = 5,
                                        nrounds = 100))
  
  
  
  mu0_y_Xt <-  c(mu0_y_Xt, unname(unlist(predict(object = OR0_Xt$outcome_fit,
                                                 newdata = as.matrix(x_train),
                                                 allow.new.levels = TRUE))))
  
  
  # Treated --------------------------------------
  OR1_Xt <- fit_y(mutated_obj1,
                  type_model = "continuous", 
                  xgboost_params = list(CV_XGB = TRUE,
                                        nfolds = 5,
                                        nrounds = 100))
  
  mu1_y_Xt <- c(mu1_y_Xt, unname(unlist(predict(object = OR1_Xt$outcome_fit,
                                                newdata = as.matrix(x_train),
                                                allow.new.levels = TRUE))))
  
  
  
  ############################################################
  # Propensity scores 
  ############################################################
  #SL -----------------------------------------------------------------------
  # Fit propensity score
  y_train_ps <- data_full_g$A
  x_train_ps <- x_train
  A = Sys.time()
  SLps = SuperLearner(Y = y_train_ps,
                      X = x_train_ps,
                      family = binomial(),
                      SL.library = c("SL.glm", 
                                     "SL.mqb",
                                     "SL.grf_nc", "SL.grf_nct",
                                     "SL.xgboost", "SL.xgboost_t"),
                      cvControl = list(V = 5))
  B = Sys.time()
  ps_hat_S <- c(ps_hat_S, SLps$SL.predict[, 1])
  # B - A   Time difference of 2.075268 hours]
  
  C = Sys.time()
  # PS methods ------------------------------------------------------------
  # GLM
  
  ps_fit_L <- glm(A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                  data = data_full_g, 
                  family = binomial)
  ps_hat_L <- c(ps_hat_L, ps_fit_L$fitted.values)
  
  # MQ -------------------------------------------------------------------
  
  obj_p_score_MQ <- list(data_p_score = data_full_g)
  class(obj_p_score_MQ) <- "MQ"
  
  ps_hat_M <- c(ps_hat_M, p_score(obj_p_score = obj_p_score_MQ,
                                  model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group)))
  
  
  
  # RF ----------------------------------------------------------------
  # No cluster no tuning 
  obj_p_score_RF <- list(data_p_score = data_full_g[, -15])
  class(obj_p_score_RF) <- "RF"
  
  ps_hat_R <-  c(ps_hat_R, p_score(obj_p_score = obj_p_score_RF,
                                   model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                   tune_RF = FALSE, 
                                   clust_RF = FALSE))
  
  
  # Cluster, no tuning
  ps_hat_Rt <-  c(ps_hat_Rt, p_score(obj_p_score = obj_p_score_RF,
                                     model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                     tune_RF = TRUE, 
                                     clust_RF = FALSE))
  
  
  
  # XGB ----------------------------------------------------------------------
  
  obj_p_score_XGB <- list(data_p_score = data_full_g[, -15])
  class(obj_p_score_XGB) <- "XGB"
  
  ps_hat_X <- c(ps_hat_X, p_score(obj_p_score = obj_p_score_XGB,
                                  model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                  xgboost_params = list(CV_XGB = FALSE,
                                                        nfolds = 5,
                                                        nrounds = 100)))
  
  
  # XGBt ----------------------------------------------------------------------
  
  
  ps_hat_Xt <-  c(ps_hat_Xt, p_score(obj_p_score = obj_p_score_XGB,
                                     model_formula = A ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + (1|group),
                                     xgboost_params = list(CV_XGB = TRUE,
                                                           nfolds = 5,
                                                           nrounds = 100)))
  
  
  
  
  
  
  
}
stop_time = Sys.time() #Time difference of 16.94238 mins


df_mu1 <- data.frame(L = mu1_y_L, 
                     
                     M = mu1_y_M, 
                     
                     R = mu1_y_R,
                     Rt = mu1_y_Rt,
                     
                     
                     X = mu1_y_X,
                     Xt = mu1_y_Xt, 
                     
                     S = mu1_y_S)

df_mu0 <- data.frame(L = mu0_y_L, 
                     
                     M = mu0_y_M, 
                     
                     R = mu0_y_R,
                     Rt = mu0_y_Rt,
                     
                     X = mu0_y_X,
                     Xt = mu0_y_Xt, 
                     
                     S = mu0_y_S)



df_ps <- data.frame(L = ps_hat_L, 
                    
                    M = ps_hat_M, 
                    
                    R = ps_hat_R,
                    Rt = ps_hat_Rt,
                    
                    X = ps_hat_X,
                    Xt = ps_hat_Xt,
                    
                    S = ps_hat_S)


# Compute AIPW -----------------------------------------------
names_method <- c("L", "M", "R", "Rt", "X", "Xt", "S")

AIPW_names <- expand.grid(names_method, names_method)
IPW_names <- expand.grid(names_method)

C = Sys.time()
for (k in 1:dim(AIPW_names)[1]) {
  
  mu1_hat <- df_mu1[, names(df_mu1) %in% unlist(unname(AIPW_names[k, 1]))]
  mu0_hat <- df_mu0[, names(df_mu0) %in% unlist(unname(AIPW_names[k, 1]))]
  
  ps_hat <- df_ps[, names(df_ps) %in% unlist(unname(AIPW_names[k, 2]))]
  
  
  AIPW_data <- data.frame(
    y = y_g,
    A = A_g,
    group = group_g,
    p_score  = ps_hat, 
    mu0_y = mu0_hat, 
    mu1_y = mu1_hat
  )
  
  AIPWf <-
    as.data.frame(calculate_tau(AIPW_data, type_tau = "AIPW"))
  
  assign(paste0(paste0(unlist(unname(AIPW_names[k, ])), collapse = ""), "Rct"), 
         AIPWf$tau)
  
  
}

for (k in 1:dim(IPW_names)[1]) {
  #for (k in 1:2) {
  
  #  group <- X[ ,names(X) %in% c("group")]
  
  ps_hat <- df_ps[, names(df_ps) %in% unlist(unname(IPW_names[k, ]))]
  
  
  IPW_data <- data.frame(
    y = y_g,
    A = A_g,
    group = group_g,
    p_score  = ps_hat
  )
  IPWf <-
    as.data.frame(calculate_tau(IPW_data, type_tau = "H"))
  
  assign(paste0(paste0(unlist(unname(IPW_names[k, ])), collapse = ""), "Rct"), IPWf$tau)
}


D = Sys.time() #Time difference of 16.61508 mins


# OR S
df_S <- data.frame(mu1_y = mu1_y_S,
                   mu0_y = mu0_y_S,
                   group = group_g)

tau_treatS = aggregate(df_S$mu1_y, list(df_S$group), FUN = mean)$x
tau_untreatS = aggregate(df_S$mu0_y, list(df_S$group), FUN = mean)$x
S = tau_treatS - tau_untreatS


# OR L 
df_L <- data.frame(mu1_y = mu1_y_L,
                   mu0_y = mu0_y_L,
                   group = group_g)

tau_treatL = aggregate(df_L$mu1_y, list(df_L$group), FUN = mean)$x
tau_untreatL = aggregate(df_L$mu0_y, list(df_L$group), FUN = mean)$x
L = tau_treatL - tau_untreatL


# OR M 
df_M <- data.frame(mu1_y = mu1_y_M,
                   mu0_y = mu0_y_M,
                   group = group_g)

tau_treatM = aggregate(df_M$mu1_y, list(df_M$group), FUN = mean)$x
tau_untreatM = aggregate(df_M$mu0_y, list(df_M$group), FUN = mean)$x
M = tau_treatM - tau_untreatM


# OR R
df_R <- data.frame(mu1_y = mu1_y_R,
                   mu0_y = mu0_y_R,
                   group = group_g)

tau_treatR = aggregate(df_R$mu1_y, list(df_R$group), FUN = mean)$x
tau_untreatR = aggregate(df_R$mu0_y, list(df_R$group), FUN = mean)$x
R = tau_treatR - tau_untreatR

# OR Rt
df_Rt <- data.frame(mu1_y = mu1_y_Rt,
                    mu0_y = mu0_y_Rt,
                    group = group_g)

tau_treatRt = aggregate(df_Rt$mu1_y, list(df_Rt$group), FUN = mean)$x
tau_untreatRt = aggregate(df_Rt$mu0_y, list(df_Rt$group), FUN = mean)$x
Rt = tau_treatRt - tau_untreatRt

# OR X
df_X <- data.frame(mu1_y = mu1_y_X,
                   mu0_y = mu0_y_X,
                   group = group_g)

tau_treatX = aggregate(df_X$mu1_y, list(df_X$group), FUN = mean)$x
tau_untreatX = aggregate(df_X$mu0_y, list(df_X$group), FUN = mean)$x
X = tau_treatX - tau_untreatX

# OR Xt
df_Xt <- data.frame(mu1_y = mu1_y_Xt,
                    mu0_y = mu0_y_Xt,
                    group = group_g)

tau_treatXt = aggregate(df_Xt$mu1_y, list(df_Xt$group), FUN = mean)$x
tau_untreatXt = aggregate(df_Xt$mu0_y, list(df_Xt$group), FUN = mean)$x
Xt = tau_treatXt - tau_untreatXt




###########################################################################
# Store results in the list - standard for baobab.                         #
############################################################################
Results = list(ni = ni,
               nc = nc,
               nt = nt,
               n = sum(ni),
               
               Ni = Ni,
               Nc = Nc,
               Nt = Nt,
               N = sum(Ni),
               frac_nN = frac_nN,
               
               
               #Tau
               tau_true = tau_true,
               tau_treat = tau_treat,
               tau_untreat = tau_untreat,
               
               # 1 - 100
               LLRct =  LLRct,
               MLRct =  MLRct,        
               RLRct =  RLRct,
               RtLRct =  RtLRct,      
               XLRct =  XLRct,      
               XtLRct =  XtLRct,
               SLRct =  SLRct,
               
               LMRct =  LMRct,      
               MMRct =  MMRct,       
               RMRct =  RMRct,        
               RtMRct =  RtMRct,      
               XMRct =  XMRct,       
               XtMRct =  XtMRct,    
               SMRct =  SMRct,  
               
               LRRct =  LRRct, 
               MRRct =  MRRct, 
               RRRct =  RRRct,       
               RtRRct =  RtRRct, 
               XRRct =  XRRct, 
               XtRRct =  XtRRct, 
               SRRct =  SRRct,   
               
               LRtRct =  LRtRct,
               MRtRct =  MRtRct,
               RRtRct =  RRtRct,      
               RtRtRct =  RtRtRct,  
               XRtRct =  XRtRct,   
               XtRtRct =  XtRtRct, 
               SRtRct =  SRtRct,  
               
               LXRct =  LXRct,        
               MXRct =  MXRct,      
               RXRct =  RXRct,        
               RtXRct =  RtXRct,   
               XXRct =  XXRct,     
               XtXRct =  XtXRct,    
               SXRct =  SXRct,  
               
               LXtRct =  LXtRct,
               MXtRct =  MXtRct, 
               RXtRct =  RXtRct, 
               RtXtRct =  RtXtRct,
               XXtRct =  XXtRct, 
               XtXtRct =  XtXtRct,    
               SXtRct =  SXtRct,  
               
               LSRct =  LSRct,        
               MSRct =  MSRct, 
               RSRct =  RSRct,        
               RtSRct =  RtSRct,    
               XSRct =  XSRct,     
               XtSRct =  XtSRct, 
               SSRct =  SSRct,  
               
               
               # IPW
               # 1 - 100
               LRct =  LRct,
               MRct =  MRct, 
               RRct =  RRct,     
               RtRct =  RtRct,      
               XRct =  XRct,   
               XtRct =  XtRct,      
               SRct =  SRct, 
               
               S_Rct = S, 
               L_Rct = L,
               M_Rct = M,
               R_Rct = R,
               Rt_Rct = Rt,
               X_Rct = X, 
               Xt_Rct = Xt)
stop_time = Sys.time() 
outputName = paste("sae4_Rct_", SN, ".RData", sep = "")
outputPath = file.path("/user/work/ze23696/Comp", outputName)
save("Results", file = outputPath)