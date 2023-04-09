
# causalSAE

<!-- badges: start -->
<!-- badges: end -->

## Overview

Package **causalSAE** implements estimators for causal small area
estimation(CSAE). CSAE-ORpost-cAIC confidence intervals (CI) for mixed
and fixed parameters under linear mixed models. Post-cAIC CI accounts
for the data-driven model selection using conditional Akaike information
criterion (cAIC).

## Installation

You can install the most recent version of **causalSAE** from
[GitHub](https://github.com/) together with the package **mquantreg**
which is used in the implementation of our methods.

``` r
# install.packages("remotes")

remotes::install_github("FelixSkarke/mquantreg")
remotes::install_github("KatarzynaReluga/causalSAE")
```

## Example

This is a basic example which shows you how to estimate heterogenous
average treatment effect at the level of subpopulation.

First we fix parameters and generate covariates.
