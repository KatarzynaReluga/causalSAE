#' Generate random elements
#'
#' This is an internal function to generate random effects or errors from
#' a normal distribution with outliers. It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use this function separately.
#'
#' @param n Number of random effects.
#' @param type_re Type of the random elements: random_effects, errors.
#' @param var_re Variance of random effects.
#' @param mean_re Mean of random effects.
#' @param frac_out Fraction of outlying random effects.
#' @param var_re_out Variance of outlying random effects.
#' @param mean_re_out Mean of outlying random effects.
#' @param seed Seed to replicate simulations.
#'
#' @importFrom stats rnorm rbinom rpois rnbinom
#'
#' @return Vector of random effects or errors with outliers.
#'
#' @export
#'
#' @examples
#'
#' random_effects <- generate_re(n = 100,
#' type_re = "random_effects",
#' var_re = 2,
#' mean_re = 0,
#' frac_out = 0.05,
#' var_re_out = 20,
#' mean_re_out = 0)
#'
#'
#' random_effects <- generate_re(n = 100,
#' type_re = "errors",
#' var_re = 2,
#' mean_re = 0,
#' frac_out = 0.05,
#' var_re_out = 20,
#' mean_re_out = 0)
#'
#'

generate_re <- function(n,
                        type_re = c("random_effects", "errors"),
                        var_re,
                        mean_re,
                        frac_out,
                        var_re_out,
                        mean_re_out,
                        seed = 1) {
  # Set seed
  set.seed(seed)

  #Set type of random elements
  type_re  = match.arg(type_re)

  if (type_re == "random_effects") {

    # Outlying random effects only at the end
    if (frac_out > 0) {
      no_out = as.integer(frac_out * n)
      re  = c(rnorm(n - no_out, mean_re, sqrt(var_re)),
              rnorm(no_out, mean_re_out, sqrt(var_re_out)))
    } else {
      re = rnorm(n, mean_re, sqrt(var_re))
    }

  } else {

    # Outlying errors anywhere

      select_out_re = rbinom(n, 1, frac_out)

      re <- (1 - select_out_re) * rnorm(n, mean_re, sqrt(var_re)) +
        select_out_re * rnorm(n, mean_re_out, sqrt(var_re_out))

  }
  return(re)
}

#'
#' Default parameters to generate random effects
#'
#' @export
#'

get_default_re <- function() {
  params_list <- list(
    var = 3,
    mean = 0,
    frac_out = 0,
    var_out = 3,
    mean_out = 0
  )
  params_list

}

#'
#' Default parameters to generate errors
#'
#' @export
#'

get_default_e <- function() {
  params_list <- list(
    var = 6,
    mean = 0,
    frac_out = 0.03,
    var_out = 150,
    mean_out = 20
  )
  params_list

}


#' Generate Poisson distributed rv with outliers
#'
#' This is an internal function to generate errors from
#' a Poisson distribution with outliers. It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use this function separately.
#'
#' @param n Number of errors (residuals)
#' @param mu Mean of Poisson rv
#' @param disturbance Shift parameter to generate outliers
#' @param frac_out Fraction of outlying random effects
#' @param seed Seed to replicate simulations
#'
#' @return Vector of random effects with outliers.
#'
#' @importFrom stats rbinom rpois
#'
#' @export
#'

generate_poisson <- function(n = 10,
                             mu = 1,
                             disturbance = 5,
                             frac_out = 0,
                             seed = 1) {
  set.seed(seed)


  select_out = rbinom(n, 1, frac_out)

  y_pois <- (1 - select_out) * rpois(n, mu) +
    select_out * rpois(n, mu + disturbance)

  return(y_pois)
}


#' Generate negative binomial distributed rv with outliers
#'
#' This is an internal function to generate errors from
#' a negative binomial distribution with outliers. It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use this function separately.
#'
#' @param n Number of errors (residuals)
#' @param mu Mean of Poisson rv
#' @param size Parameter
#' @param disturbance Shift parameter to generate outliers
#' @param frac_out Fraction of outlying random effects
#' @param seed Seed to replicate simulations
#'
#' @return Vector of random effects with outliers.
#'
#' @importFrom stats rbinom rnbinom
#'
#' @export
#'

generate_nb <- function(n = 10,
                        mu = 3,
                        size = 1,
                        disturbance = 5,
                        frac_out = 0,
                        seed = 1) {
  set.seed(seed)

  select_out = rbinom(n, 1, frac_out)

  y_nb <-
    (1 - select_out) * rnbinom(n = n, mu = mu, size = size) +
    select_out * rnbinom(n = n,
                         mu = mu + disturbance,
                         size = size)

  return(y_nb)
}


#' Generate binomial distributed random variables with outliers
#'
#' This is an internal function to generate errors from
#' a binomial distribution with outliers. It is assumed that the parameters
#' passed into the function are already checked, so we do not recommend to
#' use this function separately.
#'
#' @param n Number of errors (residuals)
#' @param p Vector of probabilities
#' @param disturbance Vector of outlying probabilities
#' @param frac_out Fraction of outlying random effects
#' @param seed Seed to replicate simulations
#'
#' @return Vector of random effects with outliers.
#'
#' @importFrom stats rbinom rpois
#'
#' @export
#'

generate_binary <- function(n = 10,
                            p = 0.5,
                            disturbance = 0.6,
                            frac_out = 0,
                            seed = 1) {
  set.seed(seed)

  select_out = rbinom(n, 1, frac_out)

  y_binary <- (1 - select_out) * rbinom(n, 1, p) +
    select_out * rbinom(n, 1, disturbance)

  return(y_binary)
}
