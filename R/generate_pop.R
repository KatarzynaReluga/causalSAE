#' Generate populations
#'
#' Function to generate populations using parametric model assumptions.
#'
#' @param X Covariates for generating outcome and propensity score.
#' @param X_outcome Covariates for generating outcome.
#' @param X_p_score Covariates for generating propensity score.
#' @param coeffs List of coefficients for generating outcomes and propensity scores, that is:
#' \itemize{
#' \item intercept_outcome - intercept for outcome regression model.
#' \item intercept_p_score - intercept for propensity score model.
#' \item coef_outcome - coefficients for outcome regression model.
#' \item coef_p_score - coefficients for propensity score model.
#' \item mean_A - mean of coefficients of treatments across different areas.
#' \item var_A - variance of coefficients of treatments across different areas.
#' }
#' @param errors_outcome List of parameters to generate errors, that is:
#' \itemize{
#'  \item var_e - variance of errors,
#'  \item mean_e - mean of errors,
#'  \item frac_out - fraction of outlying errors,
#'  \item var_e_out - variance of outlying errors,
#'  \item mean_e_out - mean of outlying errors,
#'  \item disturbance_outcome - shift parameter to generate outlying observations.
#' }
#' @param rand_eff_outcome List of parameters to generate random effects for outcome model, that is:
#' \itemize{
#'  \item var_re - variance of random effects,
#'  \item mean_re - mean of random effects,
#'  \item frac_out - fraction of outlying random effects,
#'  \item var_re_out - variance of outlying random effects,
#'  \item mean_re_out - mean of outlying random effects.
#' }
#' @param rand_eff_p_score List of parameters to generate random effects for propensity score model, that is:
#' \itemize{
#'  \item var_re - variance of random effects,
#'  \item mean_re - mean of random effects,
#'  \item frac_out - fraction of outlying random effects,
#'  \item var_re_out - variance of outlying random effects,
#'  \item mean_re_out - mean of outlying random effects.
#' }
#' @param fct_form Functional form of covariates.
#' @param regression_type Type of outcomes.
#' @param Ni_size Vector of subpopulations sizes.
#' @param m Number of subpopulations.
#' @param no_sim Number of simulations.
#' @param seed Seed to replicate simulations.
#'
#' @return
#' List composed of data frames which contain covariates X,
#' treatment status A, group indicator group, propensity score p_score, and outcome y
#'
#' @export
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
#'  seed = 1
#' )
#'
#'
#' X_outcome <- generate_X(
#'  n = N,
#'  p = 1,
#'  covariance_norm = NULL,
#'  cov_type = "lognorm",
#'  seed = 1
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
#' seed = 1)
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
                         fct_form = c("sin", "exp", "sin-exp", "lin"),
                         Ni_size  = 100,
                         m = 50,
                         no_sim = 100,
                         seed = 1) {

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

  # Define a functional form -------------------------------------
  fct_form <- match.arg(fct_form)

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
      seed = sim_seed
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
      seed = sim_seed
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
                               fct_form = fct_form,
                               rand_eff_outcome = rand_eff_outcome,
                               seed = seed)

    ## Return data frame
    X_cov <- data.frame(X)
    names(X_cov) <-  paste0("X", 1:ncol(X_cov))

    if (!is.null(X_p_score)){
      X_p_score <- data.frame(X_p_score)
      names(X_p_score) <-  paste0("Xp", 1:ncol(X_p_score))
      X_cov <- data.frame(X_cov, X_p_score)
    }

    if (!is.null(X_outcome)){
      X_outcome <- data.frame(X_outcome)
      names(X_outcome) <-  paste0("Xo", 1:ncol(X_outcome))
      X_cov <- data.frame(X_cov, X_outcome)
    }

    population <- data.frame(X_cov, A, group, p_score,
                             y = gen_outcome$y,
                             y1 = gen_outcome$y1,
                             y0 = gen_outcome$y0)
    return(population)
  }

  ## Generate no_sim populations

  if (no_sim == 1) {
    populations <- lapply(seed, generate_pop_apply)
    populations <- data.frame(populations)
  } else {
    sim_seed <- list()
    for (j in 1:no_sim) sim_seed[[j]] <- seed * j

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
                                   Ni, m, fct_form,
                                   rand_eff_outcome,
                                   seed = 1,
                                   ...) {
  set.seed(seed)

  ## Generate errors with outliers for outcome regression --------------------
  e = generate_re(
    n = length(A),
    type_re = "errors",
    var_re = errors_outcome$var_e,
    mean_re = errors_outcome$mean_e,
    frac_out = errors_outcome$frac_out,
    var_re_out = errors_outcome$var_e_out,
    mean_re_out = errors_outcome$mean_e_out,
    seed = seed
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
    seed = seed
  )

  re_repeat = rep(re, times = Ni)

  ## Generate coefficients for treatment heterogeneous treatment effect
  coef_A <- rnorm(m, coeffs$mean_A,
                  sd = sqrt(coeffs$var_A))
  coef_A_repeat <- rep(coef_A, times = Ni)

  ## Generate outcome with different area specific effect

  # define function class --------------------------------------------------
  fct_obj <- list(intercept = coeffs$intercept_outcome,
                  Xreg_outcome = Xreg_outcome)

  class(fct_obj) <- fct_form
  fx <- fct_form(fct_obj)

  y = fx + coef_A_repeat * A + re_repeat + e
  y1 = fx + coef_A_repeat * 1 + re_repeat + e
  y0 = fx + coef_A_repeat * 0 + re_repeat + e


#  y = coeffs$intercept_outcome + Xreg_outcome + coef_A_repeat * A + re_repeat + e
#  y1 = coeffs$intercept_outcome + Xreg_outcome + coef_A_repeat * 1 + re_repeat + e
#  y0 = coeffs$intercept_outcome + Xreg_outcome + coef_A_repeat * 0 + re_repeat + e

  output <- list(y = y,
                 y1 = y1,
                 y0 = y0)

  return(output)

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
                               seed = 1,
                               fct_form,
                               ...) {
  set.seed(seed)

  ## Generate random effects with outliers for outcome regression --------------------

  re = generate_re(
    n = length(A),
    type_re = "random_effects",
    var_re = rand_eff_outcome$var_re,
    mean_re = rand_eff_outcome$mean_re,
    frac_out = rand_eff_outcome$frac_out,
    var_re_out = rand_eff_outcome$var_re_out,
    mean_re_out = rand_eff_outcome$mean_re_out,
    seed = seed
  )

  re_repeat = rep(re, times = Ni)

  ## Generate coefficients for treatment heterogeneous treatment effect
  coef_A <- rnorm(m, coeffs$mean_A,
                  sd = sqrt(coeffs$var_A))
  coef_A_repeat <- rep(coef_A, times = Ni)

  ## Generate outcomes
  # define function class --------------------------------------------------
  fct_obj <- list(intercept = coeffs$intercept_outcome,
                  Xreg_outcome = Xreg_outcome)

  class(fct_obj) <- fct_form
  fx <- fct_form(fct_obj)
  #y --------------------------------
#  exp_outcome = exp(
#    coeffs$intercept_outcome + Xreg_outcome + coef_A_repeat * A + re_repeat
#  )
    exp_outcome = exp(
      fx + coef_A_repeat * A + re_repeat
    )
  p_outcome = exp_outcome * (1 + exp_outcome) ^ (-1)
  y  <- generate_binary(
    n = length(Xreg_outcome),
    p = p_outcome,
    disturbance = errors_outcome$disturbance_outcome,
    frac_out = errors_outcome$frac_out,
    seed = seed
  )
  #y1 ------------------------------------------
  exp_outcome1 = exp(
    fx + coef_A_repeat * 1 + re_repeat
  )

#  exp_outcome1 = exp(
#    coeffs$intercept_outcome + Xreg_outcome + coef_A_repeat * 1 + re_repeat
#  )
  p_outcome1 = exp_outcome1 * (1 + exp_outcome1) ^ (-1)
  y1  <- generate_binary(
    n = length(Xreg_outcome),
    p = p_outcome1,
    disturbance = errors_outcome$disturbance_outcome,
    frac_out = errors_outcome$frac_out,
    seed = seed
  )
  #y0 ------------------------------------------
  #  exp_outcome0 = exp(
  #    coeffs$intercept_outcome + Xreg_outcome + coef_A_repeat * 0 + re_repeat
  #  )
  exp_outcome0 = exp(
    fx + coef_A_repeat * 0 + re_repeat
  )

  p_outcome0 = exp_outcome0 * (1 + exp_outcome0) ^ (-1)
  y0  <- generate_binary(
    n = length(Xreg_outcome),
    p = p_outcome0,
    disturbance = errors_outcome$disturbance_outcome,
    frac_out = errors_outcome$frac_out,
    seed = seed
  )

  output <- list(y = y,
                 y1 = y1,
                 y0 = y0)

  return(output)
#  return(y)

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
                                seed = 1,
                                ...) {
  set.seed(seed)

  ## Generate random effects with outliers for outcome regression --------------------

  re = generate_re(
    n = length(A),
    type_re = "random_effects",
    var_re = rand_eff_outcome$var_re,
    mean_re = rand_eff_outcome$mean_re,
    frac_out = rand_eff_outcome$frac_out,
    var_re_out = rand_eff_outcome$var_re_out,
    mean_re_out = rand_eff_outcome$mean_re_out,
    seed = seed
  )

  re_repeat = rep(re, times = Ni)

  ## Generate coefficients for treatment heterogeneous treatment effect
  coef_A <- rnorm(m, coeffs$mean_A,
                  sd = sqrt(coeffs$var_A))
  coef_A_repeat <- rep(coef_A, times = Ni)

  ## Generate outcomes
  fct_obj <- list(intercept = coeffs$intercept_outcome,
                  Xreg_outcome = Xreg_outcome)

  class(fct_obj) <- fct_form
  fx <- fct_form(fct_obj)

  #y -------------------------
#  exp_outcome = exp(
#    coeffs$intercept_outcome + Xreg_outcome + coef_A_repeat * A + re_repeat
#  )

  exp_outcome = exp(
    fx + coef_A_repeat * A + re_repeat
  )


  y <- generate_poisson(
    n = length(Xreg_outcome),
    mu = exp_outcome,
    disturbance = errors_outcome$disturbance_outcome,
    frac_out = errors_outcome$frac_out,
    seed = seed
  )

  #y1 ---------------------
#  exp_outcome1 = exp(
#    coeffs$intercept_outcome + Xreg_outcome + coef_A_repeat * 1 + re_repeat
#  )
  exp_outcome1 = exp(
    fx + coef_A_repeat * 1 + re_repeat
  )

  y1 <- generate_poisson(
    n = length(Xreg_outcome),
    mu = exp_outcome1,
    disturbance = errors_outcome$disturbance_outcome,
    frac_out = errors_outcome$frac_out,
    seed = seed
  )

  #y0 ---------------------------
#  exp_outcome0 = exp(
#    coeffs$intercept_outcome + Xreg_outcome + coef_A_repeat * 0 + re_repeat
#  )
  exp_outcome0 = exp(
    fx + coef_A_repeat * 0 + re_repeat
  )

  y0 <- generate_poisson(
    n = length(Xreg_outcome),
    mu = exp_outcome0,
    disturbance = errors_outcome$disturbance_outcome,
    frac_out = errors_outcome$frac_out,
    seed = seed
  )

  output <- list(y = y,
                 y1 = y1,
                 y0 = y0)

  return(output)

}

#' Generate one population
#'
#' This is a generic function to generate outcome variables
#'
#' @param fct_obj Function form object
#' @param ... Additional parameters
#'
#' @return A data vector
#'
#' @export
#'

fct_form <- function(...)
  UseMethod("fct_form")

#'
#' @describeIn fct_form Sine functional form
#' @export
#'

fct_form.sin <- function(fct_obj, ...) {

  intercept <- fct_obj$intercept
  Xreg_outcome <- fct_obj$Xreg_outcome

  fx <- sin(intercept + Xreg_outcome)

  return(fx)
}

#'
#' @describeIn fct_form Sine functional form
#' @export
#'

fct_form.exp <- function(fct_obj, ...) {

  intercept <- fct_obj$intercept
  Xreg_outcome <- fct_obj$Xreg_outcome

  fx <- exp(intercept + Xreg_outcome)

  return(fx)
}

#'
#' @describeIn fct_form Sine functional form
#' @export
#'

fct_form.sin_exp <- function(fct_obj, ...) {

  intercept <- fct_obj$intercept
  Xreg_outcome <- fct_obj$Xreg_outcome

  fx <- sin(exp(intercept + Xreg_outcome))

  return(fx)
}

fct_form.lin <- function(fct_obj, ...) {

  intercept <- fct_obj$intercept
  Xreg_outcome <- fct_obj$Xreg_outcome

  fx <- intercept + Xreg_outcome

  return(fx)
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
