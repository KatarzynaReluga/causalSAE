# Model based simulations
#setwd("./causalSAE")
#devtools::load_all()
library(foreach)
library(doParallel)


# Set seed
set.seed(100)

m = 50
ni = rep(10, m)
Ni = rep(100, m)
N = sum(Ni)
n = sum(ni)

n_boot = 500

# Generate covariates
X <- generate_X(
  n = N,
  p = 1,
  covariance_norm = NULL,
  cov_type = "unif",
  seed = 1
)

X_outcome <- generate_X(
  n = N,
  p = 1,
  covariance_norm = NULL,
  cov_type = "lognorm",
  seed = 1
)

# Generate populations
populations <- generate_pop(X, X_outcome,
                            coeffs = get_default_coeffs(),
                            errors_outcome = get_default_errors_outcome(),
                            rand_eff_outcome = get_default_rand_eff_outcome(),
                            rand_eff_p_score = get_default_rand_eff_p_score(),
                            regression_type = "continuous",
                            Ni_size = 100,
                            m = 50,
                            no_sim = 1,
                            seed = 10)

NoSim <- 100

EBLUP_boot = matrix(0, nrow = n_boot, ncol = m)
MQ_boot = matrix(0, nrow = n_boot, ncol = m)
RF_boot = matrix(0, nrow = n_boot, ncol = m)
XGB_boot = matrix(0, nrow = n_boot, ncol = m)

a  = Sys.time()
registerDoParallel(14)
resEBLUP <- foreach (i = 15:28) %dopar% {
  set.seed(i * 2022)

  devtools::load_all()

  subpopulation <- sample_subpopulations(populations, frac_nc = 0.05, frac_nt = 0.05)
  data_sample <- data.frame(populations[subpopulation, ])
  data_out_of_sample <- populations[-subpopulation, ]


  # EBLUP OR --------------------------------------------------------------------------------
  EBLUP_ORf <- hte(type_hte = "OR",
                data_sample,
                data_out_of_sample,
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                  method = "EBLUP",
                                  tune_RF = FALSE,
                                  xgboost_params = list(CV_XGB = TRUE,
                                                        nfolds = 5,
                                                        nrounds = 50),
                                                        type_model = "gaussian"))
  EBLUP_OR <- EBLUP_ORf$tau

  # MQ OR ------------------------------------------------------------------------------------
  MQ_ORf <- hte(type_hte = "OR",
                   data_sample,
                   data_out_of_sample,
                   params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                    method = "MQ",
                                    tune_RF = FALSE,
                                    xgboost_params = list(CV_XGB = TRUE,
                                                          nfolds = 5,
                                                          nrounds = 50),
                                    type_model = "continuous"))
  MQ_OR <- MQ_ORf$tau

  # RF OR -------------------------------------------------------------------------------------
  RF_ORf <- hte(type_hte = "OR",
                data_sample,
                data_out_of_sample,
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                 method = "RF",
                                 tune_RF = FALSE,
                                 xgboost_params = list(CV_XGB = TRUE,
                                                       nfolds = 5,
                                                       nrounds = 50),
                                 type_model = "continuous"))
  RF_OR <- RF_ORf$tau

  # EBLUP XGB ---------------------------------------------------------------------------------
  XGB_ORf <- hte(type_hte = "OR",
                data_sample,
                data_out_of_sample,
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                 method = "XGB",
                                 tune_RF = FALSE,
                                 xgboost_params = list(CV_XGB = FALSE,
                                                       nfolds = 5,
                                                       nrounds = 50),
                                 type_model = "continuous"))
  XGB_OR <- XGB_ORf$tau

  # Bootstrap samples ----------------------------------------------------------------------
  bootstrap_indices <- sample_bootstrap_indices(sample_sizes = as.data.frame(table(data_sample$group))$Freq,
                                                out_of_sample_sizes = as.data.frame(table(data_out_of_sample$group))$Freq,
                                                n_boot = n_boot,
                                                seed = 2 * i)
  for (b in 1:n_boot) {
  # Bootstrap samples ---------------------------------------------------------

  index_sample <- bootstrap_indices[[i]]$ind_sample
  data_sample_boot <- data_sample[index_sample, ]
  row.names(data_sample_boot) <- 1 : dim(data_sample_boot)[1]

  # Bootstrap out of sample ----------------------------------------------------

  index_out_of_sample <- bootstrap_indices[[i]]$ind_population
  data_out_of_sample_boot <- data_out_of_sample[index_out_of_sample, ]
  row.names(data_out_of_sample_boot) <- 1 : dim(data_out_of_sample_boot)[1]


  # EBLUP --------------------------------------
  EBLUP_boot[b, ] <-  hte(type_hte = "OR",
                         data_sample = data_sample_boot,
                         data_out_of_sample = data_out_of_sample_boot,
                         params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                          method = "EBLUP",
                                          tune_RF = FALSE,
                                          xgboost_params = list(CV_XGB = TRUE,
                                                                nfolds = 5,
                                                                nrounds = 50),
                                          type_model = "gaussian"))$tau

  # MQ ------------------------------------------------------
  MQ_boot[b, ] <-  hte(type_hte = "OR",
                      data_sample = data_sample_boot,
                      data_out_of_sample = data_out_of_sample_boot,
                      params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                       method = "MQ",
                                       tune_RF = FALSE,
                                       xgboost_params = list(CV_XGB = TRUE,
                                                             nfolds = 5,
                                                             nrounds = 50),
                                       type_model = "continuous"))$tau

  # RF ----------------------------------------------------------
  RF_boot[b, ] <- hte(type_hte = "OR",
                     data_sample = data_sample_boot,
                     data_out_of_sample = data_out_of_sample_boot,
                     params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                 method = "RF",
                                 tune_RF = FALSE,
                                 xgboost_params = list(CV_XGB = TRUE,
                                                       nfolds = 5,
                                                       nrounds = 50),
                                 type_model = "continuous"))$tau

  # XGB --------------------------------------------------------------------
  XGB_boot[b, ] <- hte(type_hte = "OR",
                  data_sample = data_sample_boot,
                  data_out_of_sample = data_out_of_sample_boot,
                  params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                                   method = "XGB",
                                   tune_RF = FALSE,
                                   xgboost_params = list(CV_XGB = FALSE,
                                                         nfolds = 5,
                                                         nrounds = 50),
                                   type_model = "continuous"))$tau
  }
  # Variance --------------------------------------------------------------------------------------
  EBLUP_OR_var <- colMeans((EBLUP_boot - EBLUP_OR) ^ 2)
  MQ_OR_var <- colMeans((MQ_var - MQ_OR) ^ 2)
  RF_OR_var <- colMeans((RF_var - RF_OR) ^ 2)
  XGB_OR_var <- colMeans((XGB_var - XGB_OR) ^ 2)

  list(EBLUP_OR = EBLUP_OR,
       MQ_OR = MQ_OR,
       RF_OR = RF_OR,
       XGB_OR = XGB_OR,

       EBLUP_OR_var = EBLUP_OR_var,
       MQ_OR_var = MQ_OR_var,
       RF_OR_var = RF_OR_var,
       XGB_OR_var = XGB_OR_var)
}
stopImplicitCluster()
b = Sys.time()
results_all1 <- resEBLUP


  ###########################################################################
  # Store results in the list - standard for baobab.                         #
  ############################################################################
  Results = list(EBLUP_OR = EBLUP_OR,
                 MQ_OR = MQ_OR,
                 RF_OR = RF_OR,
                 XGB_OR = XGB_OR,

                 EBLUP_OR_var = EBLUP_OR_var,
                 MQ_OR_var = MQ_OR_var,
                 RF_OR_var = RF_OR_var,
                 XGB_OR_var = XGB_OR_var)

  outputName = paste("sim_OR_", a, ".RData",sep="")
  outputPath = file.path("/home/reluga/Comp", outputName)
  #outputPath=file.path("C:/Users/katar/Documents/Paper_3/sim_P_30u1e1",outputName)
  save("Results", file = outputPath)


#  registerDoParallel(14)  # use multicore, set to the number of our cores
#  x <- iris[which(iris[,5] != "setosa"), c(1,5)]
#  trials <- 100000
#  system.time({
#    r <- foreach(icount(trials), .combine=rbind) %dopar% {
#      ind <- sample(100, 100, replace=TRUE)
#      result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
#      coefficients(result1)
#    }
#  })
#  stopImplicitCluster()

#registerDoParallel(4)  # use multicore, set to the number of our cores
#foreach (i=1:3) %dopar% {
#  sqrt(i)
#}
#stopImplicitCluster()
#numCores <- detectCores()
#registerDoParallel(14)  # use multicore, set to the number of our cores
#x <- iris[which(iris[,5] != "setosa"), c(1,5)]
#trials <- 100000
#system.time({
#  r <- foreach(icount(trials), .combine=rbind) %dopar% {
#    ind <- sample(100, 100, replace=TRUE)
#    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
#    coefficients(result1)
#  }
#})
#stopImplicitCluster()
