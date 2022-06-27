#  > 0 
#> 
#> no_sim = 1,
#> ni = 5, 
#> Ni = 100, 
#> m = 50,
#> 
#> 
#> var_re_y
#> var_re_out_y
#> var_treat_effect
#> 
#> 
#> =>0
#> frac_re_out_y = 0.2,
#> no_cov_unif_y = 1,
#> coef_cov_unif_y = 2,
#> no_cov_lnorm_y = 1, 
#> 
#> coef_cov_lnorm_y = 2,
#> no_cov_unif_treat = 1,
#> coef_cov_unif_treat = 2,
#> no_cov_lnorm_treat = 1, 
#> coef_cov_lnorm_treat = 2,

pp <-get_default_params("continuous")


list(ni_size=5, #OK
     Ni_size=100, #OK
     NoSim=1, #OK
     sigma2u=3, # variance of random effect outcome model 
     sigma2e=6, #variance of the residuals outcome model
     sigma2B = 1, # variance of the treatment effects
     treat.mean=10,  # mean of treatment effects
     sigma2V = 0.25, # variance of random effects for treatment status W
     m=50, # number of areas
     out.u=0.2, #percentace of outliers in u
     out.mean.u= 9, # mean for outlier in u
     out.v.u= 20, #variance for outliers in u
     out.e=0.03,#percentace of outliers in e
     out.mean.e= 20,
     out.v.e= 150,
     seed.start=2# seed to simulation with
)



if (type == "ggboost") {
  
  list(cp = 0.01,
       eta = 0.1,
       niter = 500,
       seed = 1,
       normalize = TRUE,
       update_method = "descent_ascent")
  
} else if (type == "ggboost_cv") {
  
  params <- get_default_params("ggboost")
  params$eta <- c(0.1, 0.5)
  params$cp <- c(0.001, 0.005)
  
  params
  
  basic_params <- list(
    no_sim = 1,
    ni_size = 5,
    Ni_size = 100,
    m = 50,
    
    # Generate outcomes Y
    var_re_y = 3,
    frac_re_out_y = 0.2,
    var_re_out_y = 20,
    mean_re_out_y = 9,
    
    offset_y = 100,
    
    no_cov_unif = 1,
    coef_cov_unif = 1,
    
    no_cov_lnorm = 1,
    coef_cov_lnorm = 2,
    
    # Generate treatment T (1/0)
    var_re_treat = 0.25,
    offset_treat = -1,
    subset_cov_unif = c(1),
    
    no_cov_unif_treat = 0,
    coef_cov_unif_treat = 0,
    
    # Generate heterogeneous treatment effect
    mean_treat_effect = 10,
    var_treat_effect = 1, 
    frac_e_out_y = 0.2
  )
  
params = check_params_pop(params = basic_params, 
                 regression_type = "continuous")



gp_cont <-
  function(params = get_default_params(regression_type = "continuous"),
           start_seed = 1) {
    
    basic_pop <- gp_basic(params = get_default_params(regression_type = "continuous"),
                          start_seed = 1)
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
    
    ## Generate errors with outliers for outcome regression
    
    e = generate_e(n = N,
                   var_e = params$var_e_y,
                   mean_e = 0,
                   frac_out = params$frac_e_out_y,
                   var_e_out = params$var_e_out_y,
                   mean_e_out = params$mean_e_out_y)
    
    ## Generate random effects for treatment regression
    re_teat <- rnorm(m, 0, sqrt(params$var_re_treat))
    re_teat_repeat <- rep(re_teat, times = Ni)
    
    ## Generate treatment status
    X_unif  = generate_X_unif(p_cov = params$no_cov_unif, 
                              n_cov = N)
    
    if (params$no_cov_unif_treat == 0) {
      
      X_unif_treat  = generate_X_unif(p_cov = params$no_cov_unif_treat, 
                                      n_cov = N)
      
      exp_p <- exp(params$offset_treat + 
                     X_unif %*% params$coef_cov_unif +
                     re_teat_repeat)
    } else {
      exp_p <- exp(params$offset_treat + 
                     X_unif %*% params$coef_cov_unif + 
                     X_unif_treat %*% params$coef_cov_unif_treat + 
                     re_teat_repeat)
    }
    
    pi =  exp_p / (1 + exp_p)
    Tr <- rbinom(N, 1, pi)
    
    ## Generate heterogeneous treatment effect
    tau <- rnorm(m, params$mean_treat_effect, 
                 sd = sqrt(params$var_treat_effect))
    tau_repeat <- rep(tau, times = Ni) 
    
    
    ## Generate final outcome with different area specific effect
    X_logn <- generate_X_logn(p_cov = params$no_cov_lnorm,
                              n_cov = N)
    y = params$offset_y + X_logn %*% params$coef_cov_lnorm + 
      X_unif %*% params$coef_cov_unif + tau_repeat * Tr + 
      re_repeat + e 
    
    ## Create a data frame with a population 
    
    if (params$no_cov_unif_treat == 0) {
      population_matrix <- data.frame(y, X_unif, X_logn, group, Tr, pi)
    } else {
      population_matrix <- data.frame(y, X_unif, X_logn, 
                                      X_unif_treat, group, Tr, pi)   
    }
    
    return(population_matrix)
  }


params <- check_params_pop(params)

population <- switch(
  regression_type,
  continuous = gp_continuous(
    params = get_default_params(regression_type = "continuous"),
    start_seed = 1
  ),
  binary = gp_binary(
    params = get_default_params(regression_type = "binary"),
    start_seed = 1
  ),
  poisson = gp_poisson(
    params = get_default_params(regression_type = "poisson"),
    start_seed = 1
  ),
  nb = gp_nb(
    params = get_default_params(regression_type = "nb"),
    start_seed = 1
  )
)




generate_pop_cont <-
  function(params = get_default_params(regression_type = "continuous"),
           start_seed = 1) {
    
    no_sim <- params$no_sim
    seed_sim <- list()
    for (j in 1:no_sim) seed_sim[[j]] <- start_seed * j
    
    populations <- lapply(seed_sim, generate_one_pop_cont, 
                          params = get_default_params(regression_type = "continuous"))
    
  }

#'
#' Generate one population with continuous outcome
#'
#' This is an internal function of gp_continuous
#'
#' @inheritParams generate_population
#'
#' @return A data frame with generated population
#' 
#' @importFrom stats rbinom rnorm
#'

generate_pop_cont <-
  function(params = get_default_params(regression_type = "continuous"),
           start_seed = 1) {
    
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
    
    #y = params$offset_y + X_logn %*% params$coef_cov_lnorm + 
    #  X_unif %*% params$coef_cov_unif + tau_repeat * Tr + 
    #  re_repeat + e 
    
    
    ## Create a data frame with a population 
    
    basic_pop[["re_repeat"]] <- NULL
    basic_pop[["tau_repeat"]] <- NULL
    
    population_matrix <- data.frame(y, basic_pop)
    
    ## Sim populations
  }

#'
#' Generate one population with continuous outcome
#'
#' This is an internal function of gp_continuous
#'
#' @inheritParams generate_population
#'
#' @return A data frame with generated population
#' 
#' @importFrom stats rbinom rnorm
#'

generate_pop_binary <-
  function(params = get_default_params(regression_type = "binary"),
           start_seed = 1) {
    
    ## One population
    basic_pop <- generate_basic_elements(params = get_default_params(regression_type = "binary"),
                                         start_seed = 1)
    set.seed(start_seed)
    exp_y = exp(params$offset_y + basic_pop$X_logn %*% params$coef_cov_lnorm + 
                  basic_pop$X_unif %*% params$coef_cov_unif +  
                  basic_pop$tau_repeat * basic_pop$Tr + basic_pop$re_repeat)
    
    pi =  exp_y / (1 + exp_y)
    
    Tr <- rbinom(N, 1, pi)
    for (i in 1:n) y_sample[i] = rbinom(1, 1 , p_dj[i])
    
    
    
    ## Generate final outcome with different area specific effect
    y = params$offset_y + basic_pop$X_logn %*% params$coef_cov_lnorm + 
      basic_pop$X_unif %*% params$coef_cov_unif +  
      basic_pop$tau_repeat * basic_pop$Tr + basic_pop$re_repeat + e 
    
    #y = params$offset_y + X_logn %*% params$coef_cov_lnorm + 
    #  X_unif %*% params$coef_cov_unif + tau_repeat * Tr + 
    #  re_repeat + e 
    
    
    ## Create a data frame with a population 
    
    basic_pop[["re_repeat"]] <- NULL
    basic_pop[["tau_repeat"]] <- NULL
    
    population_matrix <- data.frame(y, basic_pop)
    
    ## Sim populations
  }

#'
#' Generate one population with continuous outcome
#'
#' This is an internal function of gp_continuous
#'
#' @inheritParams generate_population
#'
#' @return A data frame with generated population
#' 
#' @importFrom stats rbinom rnorm
#'

generate_pop_poisson <-
  function(params = get_default_params(regression_type = "poisson"),
           start_seed = 1) {
    
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
    
    
    pi =  exp_y / (1 + exp_y)
    
    Tr <- rbinom(N, 1, pi)
    for (i in 1:n) y_sample[i] = rbinom(1, 1 , p_dj[i])
    
    
    
    ## Generate final outcome with different area specific effect
    y = params$offset_y + basic_pop$X_logn %*% params$coef_cov_lnorm + 
      basic_pop$X_unif %*% params$coef_cov_unif +  
      basic_pop$tau_repeat * basic_pop$Tr + basic_pop$re_repeat + e 
    
    #y = params$offset_y + X_logn %*% params$coef_cov_lnorm + 
    #  X_unif %*% params$coef_cov_unif + tau_repeat * Tr + 
    #  re_repeat + e 
    
    
    ## Create a data frame with a population 
    
    basic_pop[["re_repeat"]] <- NULL
    basic_pop[["tau_repeat"]] <- NULL
    
    population_matrix <- data.frame(y, basic_pop)
    
    ## Sim populations
  }




#'
#' Generate one population with continuous outcome
#'
#' This is an internal function of gp_continuous
#'
#' @inheritParams generate_population
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
    
    pi =  exp_T / (1 + exp_T)
    Tr <- rbinom(N, 1, pi)
    
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
                               pi = pi, 
                               re_repeat = re_repeat)
    } else {
      basic_population <- list(X_unif = X_unif, 
                               X_unif_treat = X_unif_treat,
                               X_logn  = X_logn, 
                               group = group, 
                               Tr = Tr,
                               tau_repeat = tau_repeat,
                               pi = pi, 
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


#' Generate errors with outliers
#'
#' This is an internal function to generate errors from
#' a normal distribution with outliers. It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use this function separately.
#'
#' @param n Number of errors (residuals)
#' @param mu Variance of random effects
#' @param shift Mean of random effects
#' @param frac_out Fraction of outlying random effects
#' @param var_e_out Variance of outlying errors
#' @param mean_e_out Variance of outlying errors
#'
#' @return \item{re}{Vector of random effects with outliers}
#' 
#' @importFrom stats rbinom rpois
#'

generate_poisson <- function(n = 10,
                             mu = 1,
                             shift = 0,
                             frac_out = 0) {
  
  
  select_out_e = rbinom(n, 1, frac_out)
  
  y_pois <- (1 - select_out_e) * rpois(n, mu) +
    select_out_e * rpois(n, mu)
  
  return(y_pois)
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

##################################################
##################################################


#' Internal function to aggregate data 
#' 
#' @param data_example A data table with generated data
#' 

aggregate_treat <- function(data_example) {
  aggregate(y + pi + Tr ~ group, data = data_example, fct_treat, 
            treat_status = Tr, pi_treat = pi, treated = T)
}


#' Internal function to aggregate data 
#' 
#' @param data_example A data table with generated data
#' 
aggregate_untreat <- function(data_example) {
  aggregate(y + pi + Tr ~ group, data = ll, fct_treat, 
            treat_status = Tr, pi_treat = pi, treated = F)
  
}


aggregate_treat <- function(data_example) {
  aggregate(y + pi + Tr ~ group, data = data_example, fct_treat, 
            treat_status = Tr, pi_treat = pi, treated = T)
}

aggregate_untreat <- function(data_example) {
  #gab_cv_int(obj = obj,
  #           fold_label = fold_label,
  #           nfolds = nfolds,
  #           params = params_list, 
  #           progress = 0)
  aggregate(y  ~ group, data = ll, fct_treat, 
            treat_status = ll$Tr, pi_treat = ll$pi, treated = F)
  
}


population <- generate_population(params = list(no_sim = 2),
                                  regression_type = "continuous",
                                  start_seed = 1)



ll <- population[[1]]

subsetll <- ll[names(ll) %in% c("y", "group", "Tr", "pi")]

fct_treat <- function(y, pi_treat, treat_status, treated) {
  if (treated) {
    tau_treat <- mean((y * treat_status) / pi_treat)
  } else {
    tau_treat <- mean((y * (1-treat_status)) / (1 - pi_treat))
  }
  return(tau_treat)
}



treated$y-untreated$y

#do.call(aggregate(subsetll, by  = list(group), fct_treat, 
#                  treat_status = Tr, pi_treat = pi)
        
#        do.call("rbind", aggregate(. ~ id1 + id2, data = x, 
#                                   FUN = function(x) 
#                                     data.frame(m = mean(x), n = length(x))))
        
#        (temp3a <- aggregate(len ~ ., data = ToothGrowth,
#                             function(x) cbind(sum(x), mean(x))))
        
#        (temp2a <- aggregate(cbind(ncases, ncontrols) ~ alcgp + tobgp, data = esoph,
#                             function(x) cbind(sum(x), mean(x))))
        
#        (temp2a <- aggregate(. ~ group, data = subsetll,
#                             function(x) cbind(fct_treat(x), mean(x))))
        
        
#        [tau_treat = fct_treat(y = y, 
#                               pi_treat = pi, 
#                               treat_status = Tr),
#          tau_untreat = fct_treat(y = y, 
#                                  pi_treat = pi, 
#                                  treat_status = Tr,
#                                  treated = FALSE)]
        
#        distr.estimate <- subsetll[tau_treat = fct_treat(y = y, 
#                                     pi_treat = pi, 
#                                 treat_status = Tr),
#                                   tau_untreat = fct_treat(y = y, 
#                                                           pi_treat = pi, 
#                                                           treat_status = Tr,
#                                                           treated = FALSE)]
        
#        saddle_fun_fold_cp_eta <- cbind(t(rbindlist(list(saddle_fun_fold_cp_eta))), fold_cp_eta)
#        saddle_sum_fold_cp_eta <-  data.table(aggregate(. ~ cp + eta, data = saddle_fun_fold_cp_eta, FUN = "sum"))
        
        
