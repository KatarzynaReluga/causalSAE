#' Calculate causal effect
#'
#' This function calculates causal area effect tau
#'
#' @param populations List with one or several populations
#' @param type_tau Type of computation, \code{regular}: constant weights,
#' \code{adaptive}: adaptive weights
#'
#' @return
#' \item{tau_true}{List with taus accross small areas}
#'
#'
#'
#' @importFrom stats aggregate
#' @importFrom dplyr filter
#'
#' @export
#'
#' @examples
#'
#' #' m = 50
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
#'  start_seed = 1
#' )
#'
#'
#' X_outcome <- generate_X(
#'  n = N,
#'  p = 1,
#'  covariance_norm = NULL,
#'  cov_type = "lognorm",
#'  start_seed = 1
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
#' start_seed = 1)
#'
#'

calculate_tau <- function(populations, type_tau = c("regular", "adaptive")) {

  if (is.data.frame(populations)) {
    populations <- list(populations)
  }

  type_tau = match.arg(type_tau)

  population_treat <- function(one_pop, type_tau = type_tau) {

    one_pop_df <- data.frame(one_pop)
    group = one_pop_df$group
    group_name = unique(one_pop_df$group)
    list_group_name <- list()
    for (i in 1:length(group_name)) list_group_name[[i]] <- group_name[i]

    filter_by_group <- function(group_name) {
      data_group <- filter(one_pop_df, group == group_name)
      data_group
    }
    list_group <- lapply(list_group_name, filter_by_group)

#    list_group <- list()
#    m = length(unique(one_pop_df$group))
#    A_coef = unique(one_pop_df$coef_A_repeat)

#    for (i in unique(one_pop_df$group)) {
#      list_group[[i]] <- filter(one_pop_df, group == i)
#    }

    tau <- lapply(list_group, fct_treat_untreat, type_tau = type_tau)
    tau_subpopulations <- do.call(rbind.data.frame, tau)
    tau_subpopulations <-  cbind(tau_subpopulations, group_name)
    tau_subpopulations
  }

  tau_true <- lapply(populations, population_treat, type_tau = type_tau)

  return(tau_true)
}


#' Internal function to compute group means for the treated
#' and untreated
#' @param type_tau Type of computation, \code{regular}: constant weights,
#' \code{adaptive}: adaptive weights
#' @param subpopulation Subpopulation to compute causal effect
#'

fct_treat_untreat <- function(subpopulation, type_tau) {

  y  = subpopulation$y
  A = subpopulation$A
  p_score = subpopulation$p_score

  tau_treat_ind  <- (y * A) / p_score
  weights_treat <- A / p_score

  tau_untreat_ind <- (y * (1 - A)) / (1 - p_score)
  weights_untreat <- (1 - A) / (1 - p_score)

  if (type_tau == "regular") {
    tau_treat <- mean(tau_treat_ind)
    tau_untreat <- mean(tau_untreat_ind)
    tau <- tau_treat - tau_untreat
 } else {
    tau_treat <- sum(tau_treat_ind) / sum(weights_treat)
    tau_untreat <- sum(tau_untreat_ind) / sum(weights_untreat)
    tau <- tau_treat - tau_untreat
 }

  output <- data.frame(tau_treat = tau_treat,
                       tau_untreat = tau_untreat,
                       tau = tau)

  output
}


