#' Bootstrap variance
#'
#' Calculate bootstrap variance
#'
#' @inheritParams hte
#' @param n_boot Number of bootstrap samples
#' @param estimated_tau Estimated value of parameter tau
#' @param seed Seed to repeat the simulations
#'
#' @return
#' \item{var_tau}{Bootstrap variance of parameter tau}
#'
#' @export
#'
#'
#' @examples
#'
#' m = 50
#' ni = rep(5, m)
#' Ni = rep(100, m)
#' N = sum(Ni)
#' n = sum(ni)
#'
#' X <- generate_X(
#'  n = N,
#'  p = 1,
#'  covariance_norm = NULL,
#'  cov_type = "unif",
#'  seed = 1
#' )
#'
#' X_outcome <- generate_X(
#'  n = N,
#'  p = 1,
#'  covariance_norm = NULL,
#'  cov_type = "lognorm",
#'  seed = 1
#' )
#'
#' populations <- generate_pop(X, X_outcome,
#' coeffs = get_default_coeffs(),
#' errors_outcome = get_default_errors_outcome(),
#' rand_eff_outcome = get_default_rand_eff_outcome(),
#' rand_eff_p_score = get_default_rand_eff_p_score(),
#' regression_type = "continuous",
#' Ni_size  = 100,
#' m = 50,
#' no_sim = 1,
#' seed = 1)
#'
#' samples <- generate_sample(populations, ni_size = 5,
#'                            sample_part = "sampled",
#'                            get_index = TRUE)
#'
#' data_sample <- data.frame(samples[[1]]$samp_data)
#' index_sample <- samples[[1]]$index_s
#' data_out_of_sample <- populations[-index_sample, ]
#'
#'
#' ht_OR <- hte(type_hte = "OR",
#'                       data_sample,
#'                       data_out_of_sample,
#'                       params_p_estimate = NULL,
#'                       params_OR = list(model_formula = y ~ X1 + Xo1 + A + (1 + A||group),
#'                                        method = "EBLUP",
#'                                        tune_RF = FALSE,
#'                                        xgboost_params = list(CV_XGB = TRUE,
#'                                                              nfolds = 5,
#'                                                             nrounds = 50),
#'                                        type_model = "gaussian"))
#'
#' estimated_tau <- ht_OR$tau
#'
##' boot_var <- bootstrap_variance(formula_y,
##'                                formula_p_score,
##'                                estimated_tau,
##'                                type_tau = "Hajek",
##'                                data_sample,
##'                                data_out_of_sample,
##'                                method_y = "EBLUP",
##'                                method_p_score = "EBLUP",
##'                                n_boot = 100,
##'                                seed = 1)
##'
#'
#'
#'
bootstrap_variance <- function(obj_hte = obj_hte,
                               params_p_score = params_p_score,
                               params_impute_y = params_impute_y,
                               params_OR = params_OR,
                               n_boot = 1000,
                               estimated_tau,
                               seed = 1) {

  # Retrieve sample sizes ----------------------------------------------------

  data_sample <- obj_hte$data_sample
  data_out_of_sample <- obj_hte$data_out_of_sample

  bootstrap_indices <- sample_bootstrap_indices(sample_sizes = as.data.frame(table(data_sample$group))$Freq,
                                                out_of_sample_sizes = as.data.frame(table(data_out_of_sample$group))$Freq,
                                                n_boot = n_boot,
                                                seed = seed)

  tau_boot = matrix(0, nrow = n_boot, ncol = length(estimated_tau))

  a = Sys.time()
  for (i in 1:n_boot) {
    print(i)
    # Bootstrap samples ---------------------------------------------------------

    index_sample <- bootstrap_indices[[i]]$ind_sample
    data_sample_boot <- data_sample[index_sample, ]
    row.names(data_sample_boot) <- 1 : dim(data_sample_boot)[1]

    # Bootstrap out of sample ----------------------------------------------------

    index_out_of_sample <- bootstrap_indices[[i]]$ind_population
    data_out_of_sample_boot <- data_out_of_sample[index_out_of_sample, ]
    row.names(data_out_of_sample_boot) <- 1 : dim(data_out_of_sample_boot)[1]

    # Create object of class hte --------------------------------------------------

    obj_hte_boot <- list(data_sample = data_sample_boot,
                         data_out_of_sample = data_out_of_sample_boot)

    class(obj_hte_boot) <- class(obj_hte)
    tau_boot[i, ] <- estimate_hte(obj_hte = obj_hte_boot,
                                  params_p_score = params_p_score,
                                  params_impute_y = params_impute_y,
                                  params_OR = params_OR)$tau


  }
  b = Sys.time()
  var_tau <- colMeans((tau_boot - estimated_tau) ^ 2)
  return(var_tau)

}

#' Sample bootstrap indices
#'
#' @param sample_sizes Data from the sample
#' @param out_of_sample_sizes Data from the population
#' @param n_boot Number of bootstrap samples
#' @param seed Seed to run simulations.
#'
#' @importFrom sampling strata
#'
#' @export
#'
#' @examples
#'
#' m = 50
#' ni = rep(5, m)
#' Ni = rep(100, m)
#' N = sum(Ni)
#' n = sum(ni)
#'
#' X <- generate_X(
#'  n = N,
#'  p = 1,
#'  covariance_norm = NULL,
#'  cov_type = "unif",
#'  seed = 1
#' )
#'
#'
#' X_outcome <- generate_X(
#'  n = N,
#'  p = 1,
#'  covariance_norm = NULL,
#'  cov_type = "lognorm",
#'  seed = 1
#' )
#'
#'
#' populations <- generate_pop(X, X_outcome,
#' coeffs = get_default_coeffs(),
#' errors_outcome = get_default_errors_outcome(),
#' rand_eff_outcome = get_default_rand_eff_outcome(),
#' rand_eff_p_score = get_default_rand_eff_p_score(),
#' regression_type = "continuous",
#' Ni_size  = 100,
#' m = 50,
#' no_sim = 1,
#' seed = 1)
#'
#' samples <- generate_sample(populations, ni_size = 10,
#'                            sample_part = "sampled",
#'                            get_index = TRUE)
#'
#' data_sample <- data.frame(samples[[1]]$samp_data)
#' index_sample <- samples[[1]]$index_s
#' data_out_of_sample <- populations[-index_sample, ]
#'
#' sample_sizes = as.data.frame(table(data_sample$group))$Freq
#' out_of_sample_sizes = as.data.frame(table(data_out_of_sample$group))$Freq
#'
#' bootstrap_indices <- sample_bootstrap_indices(sample_sizes,
#'                                                 out_of_sample_sizes,
#'                                                 n_boot = 500,
#'                                                 seed = 1)
#'

sample_bootstrap_indices <- function(sample_sizes,
                                     out_of_sample_sizes,
                                     n_boot = 50,
                                     seed = 1) {
  indices <- list()
  for (i in 1:n_boot) {
    seed_sim = seed * i
    indices[[i]] <- list(ind_population = sample_indices(sample_sizes = out_of_sample_sizes,
                                                         swr = TRUE,
                                                         seed = seed_sim),
                         ind_sample = sample_indices(sample_sizes = sample_sizes,
                                                     swr = TRUE,
                                                     seed = seed_sim * 2))
  }
  return(indices)

}

#' Sample indices
#'
#' @param sample_sizes Sample sizes
#' @param swr Simple random sampling with replacement? Default = \code{swr = TRUE}.
#' @param seed Seed to run simulations.
#'
#' @return Vector with indices.
#'
#'


sample_indices <- function(sample_sizes,
                           swr = TRUE,
                           seed = 10) {

  # Cumulative sum
  cum_sample_sizes <- cumsum(sample_sizes)
  list_indices <- list()
  list_indices[[1]] <- sample(1:cum_sample_sizes[1],
                              sample_sizes[1],
                              replace = swr)

  for (i in 2:length(sample_sizes)) {

    lower_bound <- cum_sample_sizes[i-1] + 1
    upper_bound <- cum_sample_sizes[i]

    list_indices[[i]] <- sample(lower_bound : upper_bound,
                                sample_sizes[i],
                                replace = swr)
  }

  indices <- sort(unlist(list_indices))
  return(indices)
}
