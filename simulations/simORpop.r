# Model based simulations
setwd("./causalSAE")
devtools::load_all()

b = Sys.time()
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

tau <- calculate_tau(list(populations), type_tau = "H")
c = Sys.time()

# ----------------------------------------
a = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
Results = list(tau = tau)

outputName = paste("sim_OR_", a, ".RData",sep="")
outputPath = file.path("/home/reluga/Comp", outputName)
#outputPath=file.path("C:/Users/katar/Documents/Paper_3/sim_P_30u1e1",outputName)
save("Results", file = outputPath)
