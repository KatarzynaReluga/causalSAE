#' Generate covariates
#'
#' Wrapper function to generate covariates
#'
#' @param n Number of observations.
#' @param p Number of coefficients.
#' @param cov_type Type of coefficients.
#' @param covariance_norm Variance-covariance matrix for
#' correlated normally distributed covariates.
#' @param seed Seed to run simulations.
#'
#' @importFrom stats rnorm rlnorm runif
#' @importFrom mgcv rmvn
#'
#' @return Matrix of generated covariates.
#'
#' @export
#'
#' @examples
#'
#' # Generate X from normal distribution
#' X_norm <- generate_X(
#' n = 10,
#' p = 2,
#' covariance_norm = NULL,
#' cov_type = c("norm"),
#'               seed = 1)
#'
#'
#' # Generate X from log-normal distribution
#' X_norm <- generate_X(
#' n = 10,
#' p = 2,
#' covariance_norm = NULL,
#' cov_type = c("lognorm"),
#'               seed = 1)
#'
#'
#' # Generate X from uniform distribution
#' X_norm <- generate_X(
#' n = 10,
#' p = 2,
#' covariance_norm = NULL,
#' cov_type = c("unif"),
#'               seed = 1)
#'
#'
#'

generate_X <- function(n,
                       p = 2,
                       covariance_norm = NULL,
                       cov_type = c("norm",
                                    "lognorm",
                                    "unif"),
                       seed = 1) {
  ## Type of coefficient -------------------------------------------------------------------
  cov_type <- match.arg(cov_type)
  class(cov_type) <- cov_type

  ## Generate X ----------------------------------------------------------------------------
  X <- get_X(cov_type,
             n,
             p,
             covariance_norm,
             seed)
  return(X)

}

#' Get covariate \code{X}
#'
#' Generic function to generate covariate \code{X} from normal, log-normal or uniform distribution
#'
#' @inheritParams generate_X
#' @param ... Additional parameters
#'
#' @return A data frame with generated population
#'
#' @export
#'
#'
get_X <- function(...)
  UseMethod("get_X")

#'
#' @describeIn get_X Generate covariate from normal distribution.
#' @export
#'

get_X.norm <- function(cov_type,
                       n,
                       p,
                       covariance_norm,
                       seed,
                       ...) {
  ## Set seed -----------------------------------------------------------------------------
  set.seed(seed)

  ## Generate X ---------------------------------------------------------------------------
  if (is.null(covariance_norm)) {
    X_n = rnorm(n * p, mean = 0, sd = 1)
    dim(X_n) = c(n, p)
  } else {
    if (p != length(diag(covariance_norm))) {
      stop(
        "Number of coefficients does not match the dimension
             of variance-covariance matrix."
      )
    }

    if (!isSymmetric(covariance_norm)) {
      stop("Variance-covariance matrix needs to be symmetric")
    }

    X_n = rmvn(n, mu = rep(0, p), V = covariance_norm)
  }

  return(X_n)
}

#'
#' @describeIn get_X Generate covariate from log-normal distribution.
#' @export
#'


get_X.lognorm <- function(cov_type,
                          n,
                          p,
                          covariance_norm,
                          seed,
                          ...) {
  ## Set seed -----------------------------------------------------------------------------
  set.seed(seed)

  ## Generate X ---------------------------------------------------------------------------
  X_ln = rlnorm(n * p, log(4.5) - 0.5, 0.5)
  dim(X_ln) = c(n, p)

  return(X_ln)
}

#'
#' @describeIn get_X Generate covariate from uniform distribution.
#' @export
#'


get_X.unif <- function(cov_type,
                       n,
                       p,
                       covariance_norm,
                       seed,
                       ...) {
  ## Generate X --------------------------------------------------------------------------
  set.seed(seed)

  ## Generate X -------------------------------------------------------------------------
  X_u = runif(n * p, min = 0, max = 1)
  dim(X_u) = c(n, p)

  return(X_u)
}
