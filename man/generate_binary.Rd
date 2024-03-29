% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_rv_out.R
\name{generate_binary}
\alias{generate_binary}
\title{Generate binomial distributed random variables with outliers}
\usage{
generate_binary(n = 10, p = 0.5, disturbance = 0.6, frac_out = 0, seed = 1)
}
\arguments{
\item{n}{Number of errors (residuals)}

\item{p}{Vector of probabilities}

\item{disturbance}{Vector of outlying probabilities}

\item{frac_out}{Fraction of outlying random effects}

\item{seed}{Seed to replicate simulations}
}
\value{
Vector of random effects with outliers.
}
\description{
This is an internal function to generate errors from
a binomial distribution with outliers. It is assumed that the parameters
passed into the function are already checked, so we do not recommend to
use this function separately.
}
