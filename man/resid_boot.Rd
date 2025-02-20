% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residual_bootstrap.R
\name{resid_boot}
\alias{resid_boot}
\alias{resid_boot.br1}
\alias{resid_boot.br2}
\alias{resid_boot.br3}
\alias{resid_boot.br4}
\alias{resid_boot.br5}
\alias{resid_boot.br6}
\alias{resid_boot.br7}
\title{Bootstrap residuals}
\usage{
resid_boot(...)

\method{resid_boot}{br1}(obj_boot, ...)

\method{resid_boot}{br2}(obj_boot, ...)

\method{resid_boot}{br3}(obj_boot, ...)

\method{resid_boot}{br4}(obj_boot, ...)

\method{resid_boot}{br5}(obj_boot, ...)

\method{resid_boot}{br6}(obj_boot, ...)

\method{resid_boot}{br7}(obj_boot, ...)
}
\arguments{
\item{...}{Additional parameters}

\item{obj_boot}{Object to estimate residual bootstrap}
}
\description{
Internal generic function to bootstrap residuals
}
\section{Methods (by class)}{
\itemize{
\item \code{resid_boot(br1)}: Two-level block residual bootstrap without re-scaling (sample from level one and level two residuals)

\item \code{resid_boot(br2)}: One-level block residual bootstrap by treatment group

\item \code{resid_boot(br3)}: Two-level block residual bootstrap by treatment group

\item \code{resid_boot(br4)}: Two-level block residual bootstrap by treatment group

\item \code{resid_boot(br5)}: Two-level block residual bootstrap by treatment group

\item \code{resid_boot(br6)}: Two-level block residual bootstrap by treatment group

\item \code{resid_boot(br7)}: Two-level block residual bootstrap without re-scaling (sample from level one and level two residuals)

}}
