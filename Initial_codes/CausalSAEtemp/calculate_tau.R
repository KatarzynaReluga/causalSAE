#' Calculate true causal effect
#' 
#' This function calculates true causal area effect tau
#' from generated populations
#'
#' @param populations List with generate populations
#'
#' @return 
#' \item{tau_true}{List with true taus accross small areas}
#'
#' @examples
#' 
#' populations <- generate_population(params = list(no_sim = 2), 
#'                                   regression_type = "continuous",
#'                                   start_seed = 1)
#' tau_true <- calculate_tau(populations)
#' 
#' @importFrom stats aggregate
#' 
#' @export
#'
#'

calculate_tau <- function(populations) {
  
  # Internal function to aggregate data from one population 
  aggregate_treat <- function(data_example) {
    temp <- aggregate(y ~ group, data = data_example, 
                      fct_treat_untreat, 
                      p_score = data_example$p_score,
                      treat_status = data_example$Tr, 
                      simplify = F)[[2]]
    do.call(rbind.data.frame, temp)
  }
  
  tau_true <- lapply(populations, aggregate_treat)

  return(tau_true)
}


#' Internal function to compute group means for the treated
#' and untreated
#'
#' @param y Outcome variable
#' @param p_score Propensity score
#' @param treat_status Treatment indicators
#'

fct_treat_untreat <- function(y, p_score, treat_status) {
  
  tau_treat <- mean((y * treat_status) / p_score)
  tau_untreat <- mean((y * (1 - treat_status)) / (1 - p_score))
  tau <- tau_treat - tau_untreat
  
  output <- data.frame(tau_treat = tau_treat, 
                       tau_untreat = tau_untreat, 
                       tau = tau)
  
  output
}

