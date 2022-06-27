#'
#' Generate basic elements
#'
#' This is an internal function of gp_continuous
#'
#' @param start_seed Starting seed to reproduce 
#' @param Ni_size
#' @param m
#' @param additional_params
#' @param Xt
#' @param coef_p_score
#' @param hte
#' @param start_seed = 1
#' 
#' @return Data frame with common elements to generate
#' populations
#' 
#' @importFrom stats rbinom rnorm
#' 
#' @export
#'

generate_common_elements <-
  function(Ni_size, m, 
           additional_params, 
           Xt, coef_p_score,
           hte,
           start_seed = 1) {

    ## Set seed
    set.seed(start_seed)
    
    ## Number of units in areas and the total number of units
    Ni = rep(Ni_size, m)
    N  = sum(Ni)
    
    ## Group indicator
    group = rep(1:m, times = Ni)
    
    ## Generate random effects with outliers for outcome regression
    re = generate_re(
      m,
      var_re = additional_params$re$var_re,
      mean_re = 0,
      frac_out = additional_params$re$frac_out_re,
      var_re_out = additional_params$re$var_re_out,
      mean_re_out = additional_params$re$mean_re_out
    )
    
    re_repeat  = rep(re, times = Ni)
    
    ## Generate treatment status
    
    # Random effects for treatment status
    re_teat <- rnorm(m, 0, sqrt(additional_params$var_re_treat))
    re_teat_repeat <- rep(re_teat, times = Ni)
    
    # Treatment status
    exp_p_score = exp(coef_p_score$intercept + Xt$Xreg + re_teat_repeat)
    p_score = exp_p_score * (1 + exp_p_score)^(-1)
    A <- rbinom(N, 1, p_score)
  

    ## Generate heterogeneous treatment effect
    tau <- rnorm(m, hte$mean, 
                 sd = sqrt(hte$var))
    tau_repeat <- rep(tau, times = Ni) 
    

    ## Create a data frame with a population 
    common_elements <- list(group = group, A = A, 
                           tau_repeat = tau_repeat,
                           p_score = p_score, 
                           re_repeat = re_repeat)
    
    
    return(common_elements)
  }


