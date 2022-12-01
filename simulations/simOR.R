# Model based simulations
setwd("./causalSAE")
devtools::load_all()

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

tau_true <- calculate_tau(list(populations), type_tau = "H")

a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# Simple checks of the code ------------------------------------------------------------------
#for (i in 1:NoSim) {
#a  = Sys.time()
#  print(i)
  set.seed(a * 2022)

  subpopulation <- sample_subpopulations(populations, frac_nc = 0.05, frac_nt = 0.05)
  data_sample <- data.frame(populations[subpopulation, ])
  data_out_of_sample <- populations[-subpopulation, ]


  # OR -------------------------------------------------------------------
  # EBLUP OR
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

  # MQ OR
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

  # RF OR

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

  # EBLUP XGB

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
                                                seed = 2 * a)

  EBLUP_var = matrix(0, nrow = n_boot, ncol = m)
  MQ_var = matrix(0, nrow = n_boot, ncol = m)
  RF_var = matrix(0, nrow = n_boot, ncol = m)
  XGB_var = matrix(0, nrow = n_boot, ncol = m)
#}


  for (i in 1:n_boot) {
    #print(i)
    # Modulus operation
    if(i %% 10==0) {
      # Print on the screen some message
      cat(paste0("Bootstrap iteration: ", i, "\n"))
    }
    # Bootstrap samples ---------------------------------------------------------

    index_sample <- bootstrap_indices[[i]]$ind_sample
    data_sample_boot <- data_sample[index_sample, ]
    row.names(data_sample_boot) <- 1 : dim(data_sample_boot)[1]

    # Bootstrap out of sample ----------------------------------------------------

    index_out_of_sample <- bootstrap_indices[[i]]$ind_population
    data_out_of_sample_boot <- data_out_of_sample[index_out_of_sample, ]
    row.names(data_out_of_sample_boot) <- 1 : dim(data_out_of_sample_boot)[1]

    # Create object of class hte --------------------------------------------------

    #EBLUP OR
    obj_hte_boot <- list(data_sample = data_sample_boot,
                         data_out_of_sample = data_out_of_sample_boot)

    class(obj_hte_boot) <- "OR"

    # EBLUP --------------------------------------
    EBLUP_var[i, ] <-  hte(type_hte = "OR",
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
    MQ_var[i, ] <-  hte(type_hte = "OR",
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
    RF_var[i, ] <- hte(type_hte = "OR",
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
    XGB_var[i, ] <- hte(type_hte = "OR",
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
#  b = Sys.time()
  EBLUP_OR_var <- colMeans((EBLUP_var - EBLUP_OR) ^ 2)
  MQ_OR_var <- colMeans((MQ_var - MQ_OR) ^ 2)
  RF_OR_var <- colMeans((RF_var - RF_OR) ^ 2)
  XGB_OR_var <- colMeans((XGB_var - XGB_OR) ^ 2)


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


