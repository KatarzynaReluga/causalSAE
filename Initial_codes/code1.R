#' Generate population
#'
#'
#'
#' @param X
#' @param cov_outcome List of parameters to generate cov to model the outcome
#' @param frac_out
#' @param cov_p_score List of parameters to generate cov to model the treatment
#' @param hte List of parameters to generate heteregouns area effects
#' @param regression_type 
#' @param additional_params List of additional parameters
#' @param cov_norm Variance - covariance matrix to generate correlated covariates,
#'  default: cov_norm = NULL
#'
#'
#' @export 
#' 
generate_pop <- function(coef_outcome_p_score = get_default_cov("outcome_p_score"),
                         coef_outcome = get_default_cov("outcome"),
                         coef_p_score = get_default_cov("p_score"),
                         hte = list(mean = 10, var = 1),
                         Ni_size  = 1000,
                         m = 50,
                         regression_type = c("continuous",
                                             "binary",
                                             "poisson",
                                             "nb"), 
                         cov_norm = NULL, 
                         additional_params = list(frac_out = 0.05,
                                                  re = list(var_re = 3,
                                                            frac_out_re = 0.25,
                                                            var_re_out = 10, 
                                                            mean_re_out = 15),
                                                  var_re_treat  = 0.25,
                                                  params_cont = list(var_e = 6,
                                                               var_e_out = 150,
                                                               mean_e_out = 20), 
                                                  disturbance = 5), 
                         start_seed = 1){
  set.seed(start_seed)
  
  check_coef(coef_outcome, coef_p_score, coef_outcome_p_score)
  
#  stopifnot("coef_outcome must be a list" = identical(class(coef_outcome),
#                                                      "list"))
#  stopifnot("coef_p_score must be a list" = identical(class(coef_p_score),
#                                                      "list"))
#  stopifnot("coef_outcome_p_score must be a list" = identical(class(coef_outcome_p_score),
#                                                      "list"))

#  for (name in c(names(coef_outcome), names(coef_p_score), names(coef_outcome_p_score))) {
#    if (name %!in% c("intercept","coef_norm",
#                     "coef_unif", "coef_lognorm") ) {
#      stop("Generate confounders from followoing types: intercept, coef_norm, coef_unif, coef_lognorm.")
#    }
#    if (name %in% c("coef_norm", "coef_unif", "coef_lognorm") && 
#        all(length(coef_outcome_p_score[[name]]) == 0)) {
#      stop("Generate at least one confounder for outcome.")
#    }
#  }
  
  ## Number of units in areas and the total number of units
  Ni = rep(Ni_size, m)
  N  = sum(Ni)
  
  ## Group indicator
  group = rep(1:m, times = Ni)
  
  ## Check cov_outcome
  Xo <- generate_X(N = N, reg_coef = default_cov2, 
                   type = "outcome_p_score", 
                   cov_norm = NULL)
 
  Xt <- generate_X(N = N, reg_coef = coef_p_score, 
                   type = "tretment", 
                   cov_norm = NULL)

  ## Define a regression class
  #  params <- check_params_pop(params)
  #  TODO: this needs to be checked 
  # additional_params = get_default_params()
  regression_type <- match.arg(regression_type)
  class(regression_type) <- regression_type
  
  ## Generate treatment status
  
  # Random effects for treatment status
#  re_teat <- rnorm(m, 0, sqrt(additional_params$var_re_treat))
#  re_teat_repeat <- rep(re_teat, times = Ni)
  
  # Treatment status
#  exp_p_score = exp(coef_p_score$intercept + Xt$Xreg + re_teat_repeat)
#  p_score = exp_p_score * (1 + exp_p_score)^(-1)
#  A <- rbinom(N, 1, p_score)
  
  ## Generate heterogeneous treatment effect
#  tau <- rnorm(m, params$mean_treat_effect, 
#               sd = sqrt(params$var_treat_effect))
#  tau_repeat <- rep(tau, times = Ni) 
  
  ## Generate common parameters 
  common_elements <- generate_common_elements(Ni_size, m, additional_params, 
                                            Xt, coef_p_score,
                                            hte,
                                            start_seed = 1) 
  
  
  

  
  ## Generate no_sim populations 

}
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

  
#' Not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

#'
#' Generate X 
#' 
#' @param N Number of units to generate
#' @param reg_coef List of parameters
#' @param type Type of covariate 
#' @param cov_norm Variance-covariance matrix for 
#'
#'
#' @importFrom mgcv rmvn 
#' 
#' @examples 
#' cov_norm = matrix(0.25, nrow = 2, ncol = 2)
#' diag(cov_norm) <- rep(1, 2)
#' reg_coef  = default_cov
#' 
#coef_outcome = list(intercept = 1, coef_norm = 1)
#reg_coef = coef_outcome

generate_X <- function(N, reg_coef,
                       cov_norm = NULL, 
                       start_seed = 1) {
  
  set.seed(start_seed)
  
  type <- match.arg(type)
  
if (type == "outcome_p_score") {
  if (length(reg_coef$coef_norm_outcome) != length(reg_coef$coef_norm_p_score)) {
    stop("For covariate type outcome_p_score, length of vectors coef_norm_outcome
      and coef_norm_p_score should be the same.") }
  
  if (length(reg_coef$coef_unif_outcome) != length(reg_coef$coef_unif_p_score)) {
    stop("For covariate type outcome_p_score, length of vectors coef_unif_outcome
      and coef_unif_p_score should be the same.") }
  
  if (length(reg_coef$coef_lognorm_outcome) != length(reg_coef$coef_lognorm_p_score)) {
    stop("For covariate type outcome_p_score, length of vectors coef_lognorm_outcome
      and coef_lognorm_p_score should be the same.") }
}
  
  
  Xcov <- NULL
  Xreg_outcome <- NULL
  Xreg_p_score <- NULL
  
  
  if (length(reg_coef$coef_norm_outcome) != 0) {
    ### Here should be a function 
###    
    if (is.null(cov_norm)) {
      X_n = rnorm(N * length(reg_coef$coef_norm_outcome), 
                  mean = 0, sd = 1)
      dim(X_n) = c(N, length(reg_coef$coef_norm_outcome))
      
    } else {
      
      if (length(reg_coef$coef_norm_outcome) != length(diag(cov_norm))) {
        stop("Number of coefficients does not match the dimansion
             of variance-covariance matrix.") }
      
      if (!isSymmetric(cov_norm)) {
        stop("Varince covaraince matrix needs to be symmetric") }
      
      X_n = rmvn(N, mu = rep(0, length(reg_coef$coef_norm_outcome)), 
                 V = cov_norm)
    }
###    
    
    Xcov <- X_n
    Xreg_outcome = X_n %*% reg_coef$coef_norm_outcome
    Xreg_p_score = X_n %*% reg_coef$coef_norm_p_score
  }  
  
  if (length(reg_coef$coef_lognorm)!= 0) {
    ####
    
    #      X_n = rlnorm(N * coef_lognorm, log(4.48169) - 0.5, 0.5)
    X_ln = rlnorm(N * length(reg_coef$coef_lognorm), log(4.5) - 0.5, 0.5)      
    dim(X_ln) = c(N, length(reg_coef$coef_lognorm))
    
    if (is.null(Xcov)) {
      Xcov <- X_ln 
    } else {
      Xcov <- cbind(Xcov, X_ln)
    }
  ####  
    if (is.null(Xreg)) {
      Xreg <- X_ln %*% reg_coef$coef_lognorm 
    } else {
      Xreg <- Xreg + X_ln %*% reg_coef$coef_lognorm 
    }
    
  } 
  
  if (length(reg_coef$coef_unif)!= 0) {
    ###
    X_u = runif(N * length(reg_coef$coef_unif), min = 0, max = 1)      
    dim(X_u) = c(N, length(reg_coef$coef_unif))
    
    if (is.null(Xcov)) {
      Xcov <- X_u 
    } else {
      Xcov <- cbind(Xcov, X_u)
    }
    ###
    if (is.null(Xreg)) {
      Xreg <- X_u %*% reg_coef$coef_unif 
    } else {
      Xreg <- Xreg + X_u %*% reg_coef$coef_unif  
    }
  } 
  
  Xcov <- data.frame(Xcov)
  
  if (type == "outcome") {
    names(Xcov) <-  paste0("Xo", 1:ncol(Xcov))  
  } else { 
    names(Xcov) <-  paste0("Xt", 1:ncol(Xcov))  
  }
  output <- list(Xcov = Xcov, 
                 Xreg = Xreg)
  
  
}
#'
#' Generate X 
#' 
#' @param N Number of units to generate
#' @param reg_coef List of parameters
#' @param type Type of covariate 
#' @param cov_norm Variance-covariance matrix for 
#'
#'
#' @importFrom mgcv rmvn 
#' 
#' @examples 
#' cov_norm = matrix(0.25, nrow = 2, ncol = 2)
#' diag(cov_norm) <- rep(1, 2)
#' reg_coef  = default_cov
#' 
#coef_outcome = list(intercept = 1, coef_norm = 1)
#reg_coef = coef_outcome

generate_X <- function(N, reg_coef, 
                       type = c("outcome", "p_score"), 
                       cov_norm = NULL, 
                       start_seed = 1) {
  set.seed(start_seed)
  
  type <- match.arg(type)
  
  Xcov <- NULL
  Xreg <- NULL
  

  if (length(reg_coef$coef_norm) != 0) {
    
    if (is.null(cov_norm)) {
      X_n = rnorm(N * length(reg_coef$coef_norm), 
                  mean = 0, sd = 1)
      dim(X_n) = c(N, length(reg_coef$coef_norm))
      
    } else {
      
      if (length(reg_coef$coef_norm) != length(diag(cov_norm))) {
        stop("Number of coefficients does not match the dimansion
             of variance-covariance matrix.") }
      
      if (!isSymmetric(cov_norm)) {
        stop("Varince covaraince matrix needs to be symmetric") }
      
      X_n = rmvn(N, mu = rep(0, length(reg_coef$coef_norm)), 
                 V = cov_norm)
    }
    Xcov <- X_n
    Xreg = X_n %*% reg_coef$coef_norm
  }  
    
  if (length(reg_coef$coef_lognorm)!= 0) {
    #      X_n = rlnorm(N * coef_lognorm, log(4.48169) - 0.5, 0.5)
    X_ln = rlnorm(N * length(reg_coef$coef_lognorm), log(4.5) - 0.5, 0.5)      
    dim(X_ln) = c(N, length(reg_coef$coef_lognorm))
    
    if (is.null(Xcov)) {
      Xcov <- X_ln 
    } else {
      Xcov <- cbind(Xcov, X_ln)
    }
    
    if (is.null(Xreg)) {
      Xreg <- X_ln %*% reg_coef$coef_lognorm 
    } else {
      Xreg <- Xreg + X_ln %*% reg_coef$coef_lognorm 
    }

  } 
  
  if (length(reg_coef$coef_unif)!= 0) {
    X_u = runif(N * length(reg_coef$coef_unif), min = 0, max = 1)      
    dim(X_u) = c(N, length(reg_coef$coef_unif))
    
    if (is.null(Xcov)) {
      Xcov <- X_u 
    } else {
      Xcov <- cbind(Xcov, X_u)
    }
    
    if (is.null(Xreg)) {
      Xreg <- X_u %*% reg_coef$coef_unif 
    } else {
      Xreg <- Xreg + X_u %*% reg_coef$coef_unif  
    }
  } 

  Xcov <- data.frame(Xcov)
  
  if (type == "outcome") {
    names(Xcov) <-  paste0("Xo", 1:ncol(Xcov))  
  } else { 
    names(Xcov) <-  paste0("Xt", 1:ncol(Xcov))  
  }
  output <- list(Xcov = Xcov, 
                 Xreg = Xreg)
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
                               coef_outcome,
                               common_elements,
                               additional_params, frac_out,
                               Xo, Xt, 
                               start_seed = 1, ...) {
  
  set.seed(start_seed)
  
  ## Generate errors with outliers for outcome regression
  # params_cont = default_params$params_cont
  # frac_out = default_params$frac_out
  # params_re = default_params$re
  e = generate_e(n = length(Xo$Xreg),
                 var_e = additional_params$params_cont$var_e,
                 mean_e = 0,
                 frac_out = frac_out,
                 var_e_out = additional_params$params_cont$var_e_out,
                 mean_e_out = additional_params$params_cont$mean_e_out)
  

  ## Generate outcome with different area specific effect
  y = coef_outcome$intercept + Xo$Xreg + 
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


#' Check coeffcients of covariates
#'
#'
#' @param coef_outcome
#' @param coef_p_score
#' @param coef_outcome_p_score
#' 
#'  


check_coef <- function(coef_outcome, coef_p_score, coef_outcome_p_score) {
  stopifnot("coef_outcome must be a list" = identical(class(coef_outcome),
                                                      "list"))
  stopifnot("coef_p_score must be a list" = identical(class(coef_p_score),
                                                      "list"))
  stopifnot("coef_outcome_p_score must be a list" = identical(class(coef_outcome_p_score),
                                                              "list"))
  
  for (name in c(names(coef_outcome), names(coef_p_score), names(coef_outcome_p_score))) {
    if (name %!in% c("intercept","coef_norm",
                     "coef_unif", "coef_lognorm") ) {
      stop("Generate confounders from followoing types: intercept, coef_norm, coef_unif, coef_lognorm.")
    }
    if (name %in% c("coef_norm", "coef_unif", "coef_lognorm") && 
        all(length(coef_outcome_p_score[[name]]) == 0)) {
      stop("Generate at least one confounder for outcome.")
    }
  }
}

