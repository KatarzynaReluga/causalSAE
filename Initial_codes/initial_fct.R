n = 100
measurement.error = 0.02
# Add measurement error to response variable with
n.error = floor(n * measurement.error)
n.clean = n - n.error
error = sd(data$Y) * error.strength
error.var = sample(c(rep(error, n.error), rep(0, n.clean)))
data$Y = data$Y +  replicate(n, sample(c(-1, 1), 1)) * error.var


' Generate one population
#'
#' This is a generic function to generate outcome variables
#' 
#' @param regression_obj Object defining regression type
#' @param ... Additional parameters 
#'
#' @return A data frame with generated population
#' 
#' @export
#'
select_par_reg <- function(...) UseMethod("select_par_reg")

#'
#' @describeIn gen_pop Generate population with normally distributed outcomes
#' @export
#' 
select_par_reg.continuous <- function(regression_obj, 
                                      additional_par, 
                                      ...){
  par_reg <- additional_par$par_continuous
  return(par_reg)
}

#'
#' @describeIn gen_pop Generate population with normally distributed outcomes
#' @export
#' 
select_par_reg.binary <- function(regression_obj, 
                                  additional_par, 
                                  ...){
  par_reg <- additional_par$par_binary
  return(par_reg)
}

#'
#' @describeIn gen_pop Generate population with normally distributed outcomes
#' @export
#' 
select_par_reg.nb <- function(regression_obj, 
                              additional_par, 
                              ...){
  par_reg <- additional_par$par_nb
  return(par_reg)
}

#'
#' @describeIn gen_pop Generate population with normally distributed outcomes
#' @export
#' 
select_par_reg.poisson <- function(regression_obj, 
                                   additional_par, 
                                   ...){
  par_reg <- additional_par$par_poisson
  return(par_reg)
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
#' @export
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
#' @export
#' 

generate_X_logn <- function(p_cov = 2, n_cov = 10) {
  X_logn <- matrix(rlnorm(p_cov * n_cov, log(4.5) - 0.5, 0.5),
                   nrow = n_cov,
                   ncol = p_cov)
  return(X_logn)
}

