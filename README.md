
# causalSAE

<!-- badges: start -->
<!-- badges: end -->

## Overview

Package **causalSAE** implements methods for causal small area
estimation (CSAE): CSAE outcome regression estimator (CSAE-OR), CSAE
inverse probability weighting estimator (CSAE-IPW), CSAE normalized
inverse probability weighting estimator (CSAE-NIPW), CSAE augmented
inverse probability weighting estimator (CSAE-AIPW) and a direct
estimator.

Reference: Reluga, Kong, Ranjbar, Salvati, van der Laan (2025). *The
impact of job stability on monetary poverty in Italy: causal small area
estimation*.

Available at [arXiv.org](https://www.arxiv.org/abs/2502.12376).

In what follows, we use the abbreviation R-2025 to refer to the article
by Reluga, Kong, Ranjbar, Salvati, van der Laan (2025).

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

This is a basic example which shows you how to estimate heterogeneous
average treatment effects across subpopulations in the small area
setting.

First we fix parameters and generate covariates.

``` r
library(causalSAE)

m = 50
ni = rep(10, m)
Ni = rep(200, m)
N = sum(Ni)
n = sum(ni)

X <- generate_X(
  n = N,
  p = 1,
  covariance_norm = NULL,
  cov_type = "unif",
  seed = 1
)

X_outcome <- generate_X(
  n = N,
  p = 1,
  covariance_norm = NULL,
  cov_type = "lognorm",
  seed = 1
)
```

Second, we generate a population and a sample.

``` r
populations <- generate_pop(X, X_outcome,
                            coeffs = get_default_coeffs(),
                            errors_outcome = get_default_errors_outcome(),
                            rand_eff_outcome = get_default_rand_eff_outcome(),
                            rand_eff_p_score = get_default_rand_eff_p_score(),
                            regression_type = "continuous",
                            Ni_size  = 200,
                            m = 50,
                            no_sim = 1,
                            seed = 10)

samples <- generate_sample(populations, ni_size = 10,
                           sample_part = "sampled",
                           get_index = TRUE)

data_sample <- data.frame(samples[[1]]$samp_data)
index_sample <- samples[[1]]$index_s
data_out_of_sample <- populations[-index_sample, ]
```

Third, we obtain CSAE-OR estimate, CSAE-NIPW estimate and CSAE-AIPW
estimate.

``` r
hte_OR <- hte(type_hte = "OR",
              data_sample,
              data_out_of_sample,
              params_OR = list(model_formula = y ~ X1 + Xo1 + (1|group),
                               method = "EBLUP",
                               type_model = "gaussian"))

hte_NIPW <- hte(type_hte = "NIPW",
                data_sample,
                data_out_of_sample,
                params_impute_y = list(model_formula = y ~ X1 + Xo1 + (1 +                                               A||group),
                                       method = "EBLUP",
                                       type_model = "gaussian"),
                params_p_score =  list(model_formula = A ~ X1 + Xo1 + (1|group),
                                       method = "EBLUP"))

hte_AIPW <- hte(type_hte = "AIPW",
                data_sample,
                data_out_of_sample,
                params_impute_y = list(model_formula = y ~ X1 + Xo1 + (1 +                                               A||group),
                                       method = "EBLUP",
                                       type_model = "gaussian"),
                params_p_score =  list(model_formula = A ~ X1 + Xo1 + (1|group),
                                       method = "EBLUP"),
                params_OR = list(model_formula = y ~ X1 + Xo1 + (1 + A||group),
                                 method = "MQ",
                                 type_model = "continuous"))
```

## Simulation study from R-2025

All simulations in Section 6 of **R-2025** can be reproduced using
functions in the R package **causalSAE**. Detailed instructions and all
code are available in the directory
[reluga2025_results](https://github.com/KatarzynaReluga/causalSAE/tree/main/simulations/reluga2025_results).

Alongside the R files and Shell scripts, the directory includes:

1.  Additional algorithms to be added to the R package
    [SuperLearner](https://github.com/ecpolley/SuperLearner/tree/master),
    which are used to estimate nuisance parameters.

2.  The R package **BinaryMQ** and the R script **QRLM.r**, both
    authored by Nicola Salvati, which were used to implement M-quantile
    models for survey data (see references below).

## References

#### Super learner

Package
[SuperLearner](https://github.com/ecpolley/SuperLearner/tree/master)

Polley E.C., van der Laan M.J. (2010) **Super Learner in Prediction**.
U.C. Berkeley Division of Biostatistics Working Paper Series. Paper 226.
<http://biostats.bepress.com/ucbbiostat/paper266/>

van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2007) **Super
Learner**. Statistical Applications of Genetics and Molecular Biology,
6, article 25.
<http://www.degruyter.com/view/j/sagmb.2007.6.issue-1/sagmb.2007.6.1.1309/sagmb.2007.6.1.1309.xml>

van der Laan, M. J., & Rose, S. (2011). **Targeted learning: causal
inference for observational and experimental data**. Springer Science &
Business Media.

#### M-quantile for survey data

Chambers, R., Chandra, H., Salvati, N., & Tzavidis, N. (2014). **Outlier
robust small area estimation**. Journal of the Royal Statistical Society
Series B: Statistical Methodology, 76(1), 47-69.

Chambers, R. and Tzavidis, N. (2006) **M-quantile models for small area
estimation.** Biometrika, 93, 255–268.
