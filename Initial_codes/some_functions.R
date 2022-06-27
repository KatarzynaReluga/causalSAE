populations <- generate_population(params = list(no_sim = 2), 
                                   regression_type = "continuous",
                                   start_seed = 1)


fct_treat(y = populations[[1]]$y, 
          p_score = populations[[1]]$p_score, 
          treat_status = populations[[1]]$Tr, 
          treated = T)

res = fct_treat_untreat(y = populations[[1]]$y, 
                        p_score = populations[[1]]$p_score, 
                        treat_status = populations[[1]]$Tr)



aggregate_treat2(data_example = populations[[1]])
lapply(populations, aggregate_treat)



#' Internal function to compute the effect for the treated
#'
#' @param y Outcome variable
#' @param p_score Propensity score
#' @param treat_status Treatment indicators
#' @param treated Tau component for treated or non treated, default: TRUE.  
#'

fct_treat <- function(y, p_score, treat_status, treated) {
  if (treated) {
    tau_component <- mean((y * treat_status) / p_score)
  } else {
    tau_component <- mean((y * (1-treat_status)) / (1 - p_score))
  }
  return(tau_component)
}


#' 
# #' @param data_example A data table with generated data
#' 

# aggregate_treat(data_example = populations[[1]])



