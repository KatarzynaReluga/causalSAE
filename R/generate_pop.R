#' Generate populations 
#' 
#'
#' @param X Covariates for generating outcome and propensity score
#' @param X_outcome Covariates for generating outcome
#' @param X_p_score Covariates for generating propensity score
#' @param coeffs List of coefficients for generating outcomes and propensity scores, tat is:
#' \itemize{
#' \item intercept_outcome
#' \item intercept_p_score
#' \item coef_outcome
#' \item coef_p_score
#' \item mean_A
#' \item var_A
#' } 
#' @param errors_outcome List of parameters to generate errors, that is:
#' \itemize{
#'  \item var_e
#'  \item mean_e
#'  \item frac_out
#'  \item var_e_out
#'  \item mean_e_out 
#'  \item disturbance_outcome
#' }
#' @param rand_eff_outcome List of parameters to generate random effects for outcome model, that is:
#' \itemize{
#'  \item var_re
#'  \item mean_re 
#'  \item frac_out 
#'  \item var_re_out 
#'  \item mean_re_out
#' }
#' @param rand_eff_p_score List of parameters to generate random effects for propensity score model, that is:
#' \itemize{
#'  \item var_re
#'  \item mean_re 
#'  \item frac_out 
#'  \item var_re_out 
#'  \item mean_re_out
#' }
#' @param regression_type Type of outcomes
#' @param Ni_size Vector of subpopulations sizes
#' @param m Number of subpopulations
#' @param no_sim Number of simulations
#' @param start_seed Seed to replicate simulations
#'
#' @return \item{populations}{List composed of data frames which contain covariates X, 
#' treatment status A, group indicator group, propensity score p_score, and outcome y}
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
#' no_sim = 100,
#' start_seed = 1)
#' 
#' 
#'


generate_pop <- function(X = matrix(),
                         X_outcome = NULL,
                         X_p_score = NULL,
 
                         coeffs = get_default_coeffs(),
                         
                         errors_outcome = get_default_errors_outcome(),
                         
                         rand_eff_outcome = get_default_rand_eff_outcome(),
                         rand_eff_p_score = get_default_rand_eff_p_score(),

                         regression_type = c("continuous",
                                             "binary",
                                             "poisson",
                                             "nb"),
                         Ni_size  = 100,
                         m = 50,
                         no_sim = 100,
                         start_seed = 1) {
  
  ## Check coefficients --------------------------------------------------------
  check_coeffs(coeffs)
  check_errors_outcome(errors_outcome)
  check_rand_eff(rand_eff_outcome)
  check_rand_eff(rand_eff_p_score)
  
  ## Retrieve number of units in areas and the total number of units ------------
  if (length(Ni_size) == 1) {
    Ni = rep(Ni_size, m)
  } else {
    Ni = Ni_size
  }
  
  N  = sum(Ni)
  
  ## Retrieve group indicator ---------------------------------------------------
  group = rep(1:m, times = Ni)
  
  ## Define a regression class --------------------------------------------------
  regression_type <- match.arg(regression_type)
  class(regression_type) <- regression_type
  
  
  ## Wrapper function ----------------------------------------------------------
  generate_pop_apply <- function(sim_seed) {
    
    ## Set seed ------------------------------------------------------------------
    set.seed(sim_seed)

    ## Generate random effects with outliers for outcome regression --------------
    random_effects <- generate_re(
      n = m,
      type_re = "random_effects",
      var_re = rand_eff_outcome$var_re,
      mean_re = rand_eff_outcome$mean_re,
      frac_out = rand_eff_outcome$frac_out,
      var_re_out = rand_eff_outcome$var_re_out,
      mean_re_out = rand_eff_outcome$mean_re_out, 
      start_seed = sim_seed
    )
    
    re_repeat  = rep(random_effects, times = Ni)
    
    ## Generate treatment status -------------------------------------------
    
    # Random effects for treatment status
    re_treat <- generate_re(
      n = m,
      type_re = "random_effects",
      var_re = rand_eff_p_score$var_re,
      mean_re = rand_eff_p_score$mean_re,
      frac_out = rand_eff_p_score$frac_out,
      var_re_out = rand_eff_p_score$var_re_out,
      mean_re_out = rand_eff_p_score$mean_re_out, 
      start_seed = sim_seed
    )
    
    re_treat_repeat <- rep(re_treat, times = Ni)
    
    # Treatment status
    if (is.null(X_p_score)) {
      if (length(coeffs$coef_p_score) != ncol(X)) {
        stop("Length of X coefficients for propensity score must be equal to the 
           number of propensity score covariates")
      }
      Xreg_p_score <- X %*% matrix(coeffs$coef_p_score)
    } else {
      if (length(coeffs$coef_p_score) != ncol(cbind(X, X_p_score))) {
        stop("Length of X coefficients for propensity score must be equal to the 
           number of propensity score covariates")
      }
      Xreg_p_score <- cbind(X, X_p_score) %*% matrix(coeffs$coef_p_score)
    }
    
    exp_p_score = exp(coeffs$intercept_p_score + Xreg_p_score + re_treat_repeat)
    p_score = exp_p_score * (1 + exp_p_score)^(-1)
    A <- rbinom(N, 1, p_score)
    
    ## Define Xreg_outcome -------------------------------------------------------
    
    if (is.null(X_outcome)) {
      if (length(coeffs$coef_outcome) != ncol(X)) {
        stop("Length of X coefficients for outcome must be equal to the 
           number of outcome covariates")
      }
      if (length(coeffs$coef_outcome) == 1) {
        Xreg_outcome <- X %*% matrix(coeffs$coef_outcome)
      } else {
        Xreg_outcome <- X %*% coeffs$coef_outcome
      } 
    } else {
      if (length(coeffs$coef_outcome) != ncol(cbind(X, X_outcome))) {
        stop("Length of X coefficients for outcome must be equal to the 
           number of outcome covariates")
      }
      
      Xreg_outcome <-  cbind(X, X_outcome) %*% coeffs$coef_outcome
       
    }
    
    ## Generate outcome -----------------------------------------------------
    gen_outcome <- gen_outcome(regression_obj = regression_type,
                               Xreg_outcome = Xreg_outcome,
                               errors_outcome = errors_outcome,
                               coeffs = coeffs,
                               A = A,
                               Ni = Ni, m = m, 
                               rand_eff_outcome = rand_eff_outcome,
                               start_seed = start_seed)
    
    ## Return data frame 
    Xcov <- data.frame(X)
    names(Xcov) <-  paste0("X", 1:ncol(Xcov))
    
    if (!is.null(X_p_score)){
      X_p_score <- data.frame(X_p_score)
      names(X_p_score) <-  paste0("Xp", 1:ncol(X_p_score))
      X_cov <- data.frame(Xcov, X_p_score)
    } 
    
    if (!is.null(X_outcome)){
      X_outcome <- data.frame(X_outcome)
      names(X_outcome) <-  paste0("Xo", 1:ncol(X_outcome))
      X_cov <- data.frame(Xcov, X_outcome)
    } 
    
    
    
    population <- data.frame(X_cov, A, group, p_score,
                             y = gen_outcome)
    return(population)
  }
  
  
  ## Generate no_sim populations 
  
  if (no_sim == 1) {
    populations <- lapply(start_seed, generate_pop_apply)
    populations <- data.frame(populations)
  } else {
    sim_seed <- list()
    for (j in 1:no_sim) sim_seed[[j]] <- start_seed * j
    
    populations <- lapply(sim_seed, generate_pop_apply)
    
  }
  
  return(populations)

}



#' Generate one population
#'
#' This is a generic function to generate outcome variables
#' 
#' @inheritParams generate_pop
#' @param regression_obj Object defining regression type
#' @param Xreg_outcome Outcome covariates times coefficients \code{X beta} (without intercept)
#' @param A Treatment status
#' @param Ni Number of units in areas 
#' @param ... Additional parameters
#'
#' @return A data frame with generated population
#'
#' @export
#'
gen_outcome <- function(...)
  UseMethod("gen_outcome")

#'
#' @describeIn gen_outcome Generate population with normally distributed outcomes
#' @export
#'
gen_outcome.continuous <- function(regression_obj,
                                   Xreg_outcome,
                                   errors_outcome,
                                   coeffs, 
                                   A,
                                   Ni, m, 
                                   rand_eff_outcome,
                                   start_seed = 1,
                                   ...) {
  set.seed(start_seed)
  
  ## Generate errors with outliers for outcome regression --------------------
  e = generate_re(
    n = length(A),
    type_re = "errors",
    var_re = errors_outcome$var_e,
    mean_re = errors_outcome$mean_e,
    frac_out = errors_outcome$frac_out,
    var_re_out = errors_outcome$var_e_out,
    mean_re_out = errors_outcome$mean_e_out, 
    start_seed = start_seed
  )
  
  ## Generate random effects with outliers for outcome regression --------------------
  
  re = generate_re(
    n = m,
    type_re = "random_effects",
    var_re = rand_eff_outcome$var_re,
    mean_re = rand_eff_outcome$mean_re,
    frac_out = rand_eff_outcome$frac_out,
    var_re_out = rand_eff_outcome$var_re_out,
    mean_re_out = rand_eff_outcome$mean_re_out, 
    start_seed = start_seed
  )
  
  re_repeat = rep(re, times = Ni)
  
  ## Generate coefficients for treatment heterogeneous treatment effect
  coef_A <- rnorm(m, coeffs$mean_A, 
                  sd = sqrt(coeffs$var_A))
  coef_A_repeat <- rep(coef_A, times = Ni) 

  ## Generate outcome with different area specific effect
  y = coeffs$intercept_outcome + Xreg_outcome + coef_A * A + re_repeat + e
  
  return(y)
  
}

#'
#' @describeIn gen_outcome Generate population with Poisson distributed outcomes
#' @export
#'
gen_outcome.binary <- function(regression_obj,
                               Xreg_outcome,
                               errors_outcome,
                               coeffs, 
                               A,
                               Ni, m, 
                               rand_eff_outcome,
                               start_seed = 1,
                               ...) {
  set.seed(start_seed)
  
  ## Generate random effects with outliers for outcome regression --------------------
  
  re = generate_re(
    n = length(A),
    type_re = "random_effects",
    var_re = rand_eff_outcome$var_re,
    mean_re = rand_eff_outcome$mean_re,
    frac_out = rand_eff_outcome$frac_out,
    var_re_out = rand_eff_outcome$var_re_out,
    mean_re_out = rand_eff_outcome$mean_re_out, 
    start_seed = start_seed
  )
  
  re_repeat = rep(re, times = Ni)
  
  ## Generate coefficients for treatment heterogeneous treatment effect
  coef_A <- rnorm(m, coeffs$mean_A, 
                  sd = sqrt(coeffs$var_A))
  coef_A_repeat <- rep(coef_A, times = Ni) 
 
  ## Generate outcomes
  exp_outcome = exp(
    coeffs$intercept_outcome + Xreg_outcome + coef_A * A + re_repeat
  )
  p_outcome = exp_outcome * (1 + exp_outcome) ^ (-1)
  y  <- generate_binary(
    n = length(Xreg_outcome),
    p = p_outcome,
    disturbance = errors_outcome$disturbance_outcome,
    frac_out = errors_outcome$frac_out,
    start_seed = start_seed
  )
  
  
  return(y)
  
}

#'
#' @describeIn gen_outcome Generate population with binomial distributed outcomes
#' @export
#'
gen_outcome.poisson <- function(regression_obj,
                                Xreg_outcome,
                                errors_outcome,
                                coeffs, 
                                A,
                                Ni, m, 
                                rand_eff_outcome,
                                start_seed = 1,
                                ...) {
  set.seed(start_seed)
  
  ## Generate random effects with outliers for outcome regression --------------------
  
  re = generate_re(
    n = length(A),
    type_re = "random_effects",
    var_re = rand_eff_outcome$var_re,
    mean_re = rand_eff_outcome$mean_re,
    frac_out = rand_eff_outcome$frac_out,
    var_re_out = rand_eff_outcome$var_re_out,
    mean_re_out = rand_eff_outcome$mean_re_out, 
    start_seed = start_seed
  )
  
  re_repeat = rep(re, times = Ni)
  
  ## Generate coefficients for treatment heterogeneous treatment effect
  coef_A <- rnorm(m, coeffs$mean_A, 
                  sd = sqrt(coeffs$var_A))
  coef_A_repeat <- rep(coef_A, times = Ni) 
  
  ## Generate outcomes
  exp_outcome = exp(
    coeffs$intercept_outcome + Xreg_outcome + coef_A * A + re_repeat
  )
  
  y <- generate_poisson(
    n = length(Xreg_outcome),
    mu = exp_outcome,
    disturbance = errors_outcome$disturbance_outcome,
    frac_out = errors_outcome$frac_out,
    start_seed = start_seed
  )
  
  return(y)
  
}
# coeffs
check_coeffs <- function(coeffs) {
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  stopifnot("coeffs must be a list" = identical(class(coeffs), "list"))

  for (name in c(names(coeffs))) {
    if (name %!in% c(
      "intercept_outcome", 
      "intercept_p_score",
      "coef_outcome",
      "coef_p_score",
      "mean_A",
      "var_A")) {
      stop(
        "Define one of the following coefficients: intercept_outcome, 
          intercept_p_score, coef_outcome, coef_p_score, mean_A, var_A."
      )
    }
  }

}

#errors_outcome
check_errors_outcome <- function(errors_outcome) {
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  stopifnot("coeffs must be a list" = identical(class(errors_outcome), "list"))
  
  for (name in c(names(errors_outcome))) {
    if (name %!in% c(
      "var_e",
      "mean_e",
      "frac_out",
      "var_e_out",
      "mean_e_out",
      "disturbance_outcome")) {
      stop(
        "Generate outcomes with disturbances using: var_e, 
          mean_e, frac_out, var_e_out, mean_e_out, disturbance_outcome"
      )
    }
  }

}

# check random effects 
check_rand_eff <- function(rand_eff) {
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  stopifnot("rand_eff must be a list" = identical(class(rand_eff), "list"))
  
  for (name in c(names(rand_eff))) {
    if (name %!in% c(
      "var_re",
      "mean_re", 
      "frac_out", 
      "var_re_out", 
      "mean_re_out")) {
      stop(
        "Define one of the following coefficients: intercept_outcome, 
          intercept_p_score, coef_outcome, coef_p_score, mean_A, var_A."
      )
    }
  }
  
}



#'
#' Get default coefficients
#'
#' @export
#'

get_default_coeffs <- function() {
  
  params_list <- list(
    intercept_outcome = 100,
    intercept_p_score = - 1,
    coef_outcome = c(1, 2),
    coef_p_score = 0.5, 
    mean_A = 10, 
    var_A = 1
    )
  return(params_list)
}

#'
#' Get default parameters to generate outcomes with possible disturbances
#'
#' @export
#'

get_default_errors_outcome <- function() {
  
  errors_outcome = list(
    var_e = 1,
    mean_e = 0,
    frac_out = 0.1,
    var_e_out = 10,
    mean_e_out = 0, 
    disturbance_outcome = 5
  )
  return(errors_outcome)
}


#'
#' Get default parameters to generate random effects with possible disturbances for outcome model
#'
#' @export
#'

get_default_rand_eff_outcome <- function() {
  
  rand_effects = list(
    var_re = 3,
    mean_re = 0,
    frac_out = 0,
    var_re_out = 0,
    mean_re_out = 0
  )
  return(rand_effects)
}

#'
#' Get default parameters to generate random effects with possible disturbances for propensity score model
#'
#' @export
#'

get_default_rand_eff_p_score <- function() {
  
  rand_effects = list(
    var_re = 0.25,
    mean_re = 0,
    frac_out = 0,
    var_re_out = 0,
    mean_re_out = 0
  )
  return(rand_effects)
}
