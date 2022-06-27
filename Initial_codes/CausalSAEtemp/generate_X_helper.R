#'  
#' Default values to generate covariates
#'
get_default_cov <- function(type = c("outcome_p_score", 
                                     "outcome", 
                                     "p_score")) {
  
  type <- match.arg(type)
  
  if (type == "outcome_p_score") {
    default_cov = list(coef_unif_outcome = 1, 
                       coef_unif_p_score = 0.5)
    
  } else if (type == "outcome") {
    default_cov = list(intercept = 100, 
                       coef_lognorm = 2)
  } else {
    default_cov = list(intercept = -1)
  }
  
  return(default_cov)
}

reg_coef <- default_cov


####################
#' Generate normally distributed covariates  
#' 
#' @param n Number of units
#' @param coef_norm Number of units to generate
#' @param cov_norm Variance - covariance matrix to generate correlated covariates,
#'
#' @importFrom mgcv rmvn 
#' 
#' @examples 
#' X_norm = get_X_norm(n = 10, coef_norm = c(1, 1), cov_norm = NULL)
#' 
#' 


get_X_norm <- function(n = 10, 
                       coef_norm = c(1, 1), 
                       cov_norm = NULL) {
  
  if (is.null(cov_norm)) {
    X_n = rnorm(n * length(coef_norm), mean = 0, sd = 1)
    dim(X_n) = c(n, length(coef_norm))
  } else {
    
    if (length(coef_norm) != length(diag(cov_norm))) {
      stop("Number of coefficients does not match the dimansion
             of variance-covariance matrix.") }
    
    if (!isSymmetric(cov_norm)) {
      stop("Varince covaraince matrix needs to be symmetric") }
    
    X_n = rmvn(n, mu = rep(0, length(coef_norm)), V = cov_norm)
  }
  
  return(X_n)
}
get_X_norm()
################## 
#' Generate log normally distributed covariates  
#' 
#' @param n number of units
#' @param coef_lognorm Number of units to generate
#'
#' 
#' @examples 
#' X_ln = matrix(0.25, nrow = 2, ncol = 2)
#' diag(cov_norm) <- rep(1, 2)
#' reg_coef  = default_cov
#' 
#coef_outcome = list(intercept = 1, coef_norm = 1)
#reg_coef = coef_outcome


get_X_lognorm <- function(n = 10, 
                          coef_lognorm = c(1, 1)) {
  
  ####
  #      X_n = rlnorm(N * coef_lognorm, log(4.48169) - 0.5, 0.5)
  X_ln = rlnorm(N * length(coef_lognorm), log(4.5) - 0.5, 0.5)      
  dim(X_ln) = c(N, length(coef_lognorm))
  
  return(X_ln)
}
get_X_lognorm()
################## 
#' Generate uniformly distributed covariates  
#' 
#' @param n number of units
#' @param coef_unif Number of units to generate
#'
#' 
#' @examples 
#' X_ln = matrix(0.25, nrow = 2, ncol = 2)
#' diag(cov_norm) <- rep(1, 2)
#' reg_coef  = default_cov
#' 
#coef_outcome = list(intercept = 1, coef_norm = 1)
#reg_coef = coef_outcome


get_X_unif <- function(n = 10, 
                       coef_unif = c(1, 1)) {
  
  X_u = runif(N * length(coef_unif), min = 0, max = 1)      
  dim(X_u) = c(N, length(coef_unif))
  
  return(X_u)
}
get_X_unif()
