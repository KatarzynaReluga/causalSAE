#' Generate populations
#'
#' This function generates populations which can be used to
#' test the performance of our methods in simulations
#'
#' @inheritParams check_params_pop
#' @param start_seed Seed to reproduce simulations
#'
#' @return A data frame with generated population
#'
#' @examples
#' 
#' population <- generate_population(params = get_default_params(regression_type = "continuous"),
#'                                   regression_type = "continuous",
#'                                   start_seed = 1)
#'  
#' @export
#'
generate_population <-
  function(params = get_default_params(regression_type = "continuous"),
           regression_type = c("continuous",
                               "binary",
                               "poisson",
                               "nb"),
           start_seed = 1) {
    
    ## Check params and define a regression class
    params <- check_params_pop(params)
    regression_type <- match.arg(regression_type)
    class(regression_type) <- regression_type

    ## Generate no_sim populations  
    generate_pop_apply <- function(sim_seed) {
      gen_pop(regression_obj = regression_type,
              params = get_default_params(regression_type), 
              start_seed = sim_seed)
    }
    
    no_sim <- params$no_sim
    
    seed_sim <- list()
    for (j in 1:no_sim) seed_sim[[j]] <- start_seed * j
    
      populations <- lapply(seed_sim, generate_pop_apply)

    return(populations)
  }


#' Generate one population
#'
#' This is a generic function to generate outcome variables
#' 
#' @param regression_obj Object defining regression type
#' @param params List of parameters to generate a population
#' @param start_seed Start seed
#' @param ... Additional parameters 
#'
#' @return A data frame with generated population
#' 
#' @export
#'
gen_pop <- function(...) UseMethod("gen_pop")

#'
#' @describeIn gen_pop Generate population with normally distributed outcomes
#' @export
#' 
gen_pop.continuous <- function(regression_obj, 
                               params = get_default_params(regression_type = "continuous"),
                               start_seed = 1, ...) {
  
  ## One population
  basic_pop <- generate_basic_elements(params = get_default_params(regression_type = "continuous"),
                                       start_seed = 1)
  set.seed(start_seed)
  
  ## Generate errors with outliers for outcome regression
  
  e = generate_e(n = length(basic_pop$group),
                 var_e = params$var_e_y,
                 mean_e = 0,
                 frac_out = params$frac_e_out_y,
                 var_e_out = params$var_e_out_y,
                 mean_e_out = params$mean_e_out_y)
  
  
  ## Generate final outcome with different area specific effect
  y = params$offset_y + basic_pop$X_logn %*% params$coef_cov_lnorm + 
    basic_pop$X_unif %*% params$coef_cov_unif +  
    basic_pop$tau_repeat * basic_pop$Tr + basic_pop$re_repeat + e 
  

  ## Create a data frame with a population 
  
  basic_pop[["re_repeat"]] <- NULL
  basic_pop[["tau_repeat"]] <- NULL
  
  population_matrix <- data.frame(y, basic_pop)
  
}

#'
#' @describeIn gen_pop Generate population with Poisson distributed outcomes
#' @export
#'
gen_pop.binary <- function(regression_obj, 
                           params = get_default_params(regression_type = "binary"),
                           start_seed = 1,  ...) {
  
  ## One population
  basic_pop <- generate_basic_elements(params = get_default_params(regression_type = "binary"),
                                       start_seed = 1)
  set.seed(start_seed)
  exp_y = exp(params$offset_y + basic_pop$X_logn %*% params$coef_cov_lnorm + 
                basic_pop$X_unif %*% params$coef_cov_unif +  
                basic_pop$tau_repeat * basic_pop$Tr + basic_pop$re_repeat)
  
  p_y =  exp_y / (1 + exp_y)
  
  y  <- generate_binary(n = length(basic_pop$group),
                  p = p_y,
                  shift = 0.6,
                  frac_out = 0.02)
  
  ## Create a data frame with a population 
  
  basic_pop[["re_repeat"]] <- NULL
  basic_pop[["tau_repeat"]] <- NULL
  
  population_matrix <- data.frame(y, basic_pop)
  
}

#'
#' @describeIn gen_pop Generate population with binomial distributed outcomes
#' @export
#' 
gen_pop.poisson <- function(regression_obj, 
                           params = get_default_params(regression_type = "poisson"),
                           start_seed = 1,  ...) {
    
    ## One population
    basic_pop <- generate_basic_elements(params = get_default_params(regression_type = "poisson"),
                                         start_seed = 1)
    set.seed(start_seed)
    exp_y = exp(0.5 + basic_pop$X_logn %*% params$coef_cov_lnorm + 
                  basic_pop$X_unif %*% params$coef_cov_unif +  
                  basic_pop$tau_repeat * basic_pop$Tr + basic_pop$re_repeat)
    
    
    y = generate_poisson(n = length(basic_pop$group),
                         mu = exp_y, shift = params$shift_contam_y, 
                         frac_out = params$frac_contam_y) 
    
    
    ## Create a data frame with a population 
    
    basic_pop[["re_repeat"]] <- NULL
    basic_pop[["tau_repeat"]] <- NULL
    
    population_matrix <- data.frame(y, basic_pop)
    
  }


#'
#' Generate one population with continuous outcome
#'
#' This is an internal function of gp_continuous
#'
#' @inheritParams generate_population
#' @param start_seed Starting seed to reproduce 
#' 
#' @return A data frame with generated population
#' 
#' @importFrom stats rbinom rnorm
#'

generate_basic_elements <-
  function(params = get_default_params(regression_type = "continuous"),
           start_seed = 1) {
    set.seed(start_seed)
    
    ## Number of units in areas and the total number of units
    m  = params$m
    Ni = rep(params$Ni_size, m)
    N  = sum(Ni)
    
    ## Group indicator
    group = rep(1:m, times = Ni)
    
    ## Generate random effects with outliers for outcome regression
    re = generate_re(
      m,
      var_re = params$var_re_y,
      mean_re = 0,
      frac_out = params$frac_re_out_y,
      var_re_out = params$var_re_out_y,
      mean_re_out = params$mean_re_out_y
    )
    
    re_repeat  = rep(re, times = Ni)
    
    ## Generate random effects for treatment regression
    re_teat <- rnorm(m, 0, sqrt(params$var_re_treat))
    re_teat_repeat <- rep(re_teat, times = Ni)
    
    ## Generate treatment status
    X_unif  = generate_X_unif(p_cov = params$no_cov_unif, 
                              n_cov = N)
    
    if (params$no_cov_unif_treat == 0) {
      
      X_unif_treat  = generate_X_unif(p_cov = params$no_cov_unif_treat, 
                                      n_cov = N)
      
      exp_T <- exp(params$offset_treat + 
                     X_unif %*% params$coef_cov_unif +
                     re_teat_repeat)
    } else {
      exp_T <- exp(params$offset_treat + 
                     X_unif %*% params$coef_cov_unif + 
                     X_unif_treat %*% params$coef_cov_unif_treat + 
                     re_teat_repeat)
    }
    
    p_score =  exp_T / (1 + exp_T)
    Tr <- rbinom(N, 1, p_score)
    
    ## Generate heterogeneous treatment effect
    tau <- rnorm(m, params$mean_treat_effect, 
                 sd = sqrt(params$var_treat_effect))
    tau_repeat <- rep(tau, times = Ni) 
    
    
    ## Generate covariates for y 
    X_logn <- generate_X_logn(p_cov = params$no_cov_lnorm,
                              n_cov = N)
    
    ## Create a data frame with a population 
    
    if (params$no_cov_unif_treat == 0) {
      basic_population <- list(X_unif = X_unif, 
                               X_logn  = X_logn, 
                               group = group, 
                               Tr = Tr, 
                               tau_repeat = tau_repeat,
                               p_score = p_score, 
                               re_repeat = re_repeat)
    } else {
      basic_population <- list(X_unif = X_unif, 
                               X_unif_treat = X_unif_treat,
                               X_logn  = X_logn, 
                               group = group, 
                               Tr = Tr,
                               tau_repeat = tau_repeat,
                               p_score = p_score, 
                               re_repeat = re_repeat)
    }
    return(basic_population)
  }




#' Generate random effects with outliers
#'
#' This is an internal function to generate random effects from
#' a normal distribution with outliers. It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use this function separately.
#'
#' @param m Number of random effects
#' @param var_re Variance of random effects
#' @param mean_re Mean of random effects
#' @param frac_out Fraction of outlying random effects
#' @param var_re_out Variance of outlying random effects
#' @param mean_re_out Variance of outlying random effects
#'
#' @return \item{re}{Vector of random effects with outliers}
#'

generate_re <- function(m = 100,
                        var_re = 1,
                        mean_re = 0,
                        frac_out = 0,
                        var_re_out = 9,
                        mean_re_out = 20) {
  if (frac_out > 0) {
    no_out = as.integer(frac_out * m)
    re  = c(rnorm(m - no_out, mean_re, sqrt(var_re)),
            rnorm(no_out, mean_re_out, sqrt(var_re_out)))
  } else {
    re = rnorm(m, mean_re, sqrt(var_re))
  }
  
  return(re)
}


#' Generate errors with outliers
#'
#' This is an internal function to generate errors from
#' a normal distribution with outliers. It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use this function separately.
#'
#' @param n Number of errors (residuals)
#' @param var_e Variance of random effects
#' @param mean_e Mean of random effects
#' @param frac_out Fraction of outlying random effects
#' @param var_e_out Variance of outlying errors
#' @param mean_e_out Variance of outlying errors
#'
#' @return \item{re}{Vector of random effects with outliers}
#' 
#' @importFrom stats rbinom rnorm
#'

generate_e <- function(n = 10,
                       var_e = 1,
                       mean_e = 0,
                       frac_out = 0,
                       var_e_out = 9,
                       mean_e_out = 20) {
  
  
  select_out_e = rbinom(n, 1, frac_out)
  
  e <- (1 - select_out_e) * rnorm(n, mean_e, sqrt(var_e)) +
    select_out_e * rnorm(n, mean_e_out, sqrt(var_e_out))
  
  return(e)
}


#' Generate Poisson distributed rv with outliers
#'
#' This is an internal function to generate errors from
#' a normal distribution with outliers. It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use this function separately.
#'
#' @param n Number of errors (residuals)
#' @param mu Mean of Poisson rv
#' @param shift Shift parameter to generate outliers 
#' @param frac_out Fraction of outlying random effects
#'
#' @return \item{y_pois}{Vector of random effects with outliers}
#' 
#' @importFrom stats rbinom rpois
#'

generate_poisson <- function(n = 10,
                             mu = 1,
                             shift = 5,
                             frac_out = 0) {
  
  
  select_out = rbinom(n, 1, frac_out)
  
  y_pois <- (1 - select_out) * rpois(n, mu) +
    select_out * rpois(n, mu + shift)
  
  return(y_pois)
}


#' Generate negative binomial distributed rv with outliers
#'
#' This is an internal function to generate errors from
#' a normal distribution with outliers. It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use this function separately.
#'
#' @param n Number of errors (residuals)
#' @param mu Mean of Poisson rv
#' @param size Parameter
#' @param shift Shift parameter to generate outliers 
#' @param frac_out Fraction of outlying random effects
#'
#' @return \item{y_pois}{Vector of random effects with outliers}
#' 
#' @importFrom stats rbinom rnbinom
#'

generate_nb <- function(n = 10,
                        mu = 3,
                        size = 1,
                        shift = 5,
                        frac_out = 0) {
  
  
  select_out = rbinom(n, 1, frac_out)
  
  y_pois <- (1 - select_out) * rnbinom(n = n, mu = mu, size = size) +
    select_out * rnbinom(n = n, mu = mu + shift, size = size)
  
  return(y_pois)
}


#' Generate binomial distributed random variables with outliers
#'
#' This is an internal function to generate errors from
#' a normal distribution with outliers. It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use this function separately.
#'
#' @param n Number of errors (residuals)
#' @param p Vector of probabilities 
#' @param shift Vector of outlying probabilites 
#' @param frac_out Fraction of outlying random effects
#'
#' @return \item{y_pois}{Vector of random effects with outliers}
#' 
#' @importFrom stats rbinom rpois
#'

generate_binary <- function(n = 10,
                            p = 0.5,
                            shift = 0.6,
                            frac_out = 0) {
  
  
  select_out = rbinom(n, 1, frac_out)
  
  y_binary <- (1 - select_out) * rbinom(n, 1, p) +
    select_out * rbinom(n, 1, shift)
  #TODO: This shift needs to be fixed here
  
  return(y_binary)
}



#'
#' Generate uniformly distributed covariates
#'
#' This is an internal function to generate uniformly distributed covariates.
#' It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use it separately.
#'
#' @param p_cov Number of covariates to generate
#' @param n_cov Length of covariates
#'
#' @return \item{X_logn}{Matrix of log normally distributed covariates}
#' 
#' @importFrom stats runif
#'

generate_X_unif <- function(p_cov = 2, n_cov = 10) {
  X_unif <- matrix(runif(p_cov * n_cov, min = 0, max = 1),
                   nrow = n_cov,
                   ncol = p_cov)
  return(X_unif)
}

#'
#' Generate log normally distributed covariates
#'
#' This is an internal function to generate log normally distributed covariates.
#' It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use it separately.
#'
#' @param p_cov Number of covariates to generate
#' @param n_cov Length of covariates
#'
#' @return \item{X_logn}{Matrix of log normally distributed covariates}
#' 
#' @importFrom stats rlnorm 
#' 

generate_X_logn <- function(p_cov = 2, n_cov = 10) {
  X_logn <- matrix(rlnorm(p_cov * n_cov, log(4.5) - 0.5, 0.5),
                   nrow = n_cov,
                   ncol = p_cov)
  return(X_logn)
}

