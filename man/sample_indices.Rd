% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residual_bootstrap.R
\name{sample_indices}
\alias{sample_indices}
\title{Sample indices}
\usage{
sample_indices(sample_sizes, swr = TRUE, seed = 10, rand_clust = FALSE)
}
\arguments{
\item{sample_sizes}{Sample sizes}

\item{swr}{Simple random sampling with replacement? Default = \code{swr = TRUE}.}

\item{seed}{Seed to run simulations.}

\item{rand_clust}{Select clusters randomly. Default: \code{rand_clust = FALSE}.}
}
\value{
Vector with indices.
}
\description{
Sample indices
}
