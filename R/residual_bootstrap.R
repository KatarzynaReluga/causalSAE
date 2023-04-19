#' Residual bootstrap
#'
#' One-level and two-level random block residual bootstrap for variance calculations
#'
#' @param y Vector of true outcomes
#' @param y_hat Vector of estimated outcomes
#' @param data_sample Sample data with group indicator
#' @param boot_seed Seed to repeat bootstrap sampling
#' @param n_boot Number of bootstrap samples
#' @param type_boot Type of bootstrap: \code{br1} - one-level block residual bootstrap,
#' \code{br2} - two-level block residual bootstrap.
#' @param method_scale - Method to scale residuals: \code{sd}: standard deviation,
#'  \code{Qn}: robust and efficient estimate of scale as in Rousseeuw and Croux (1993)
#'  \code{mad}: median absolute deviation. Default: \code{method_scale = sd}.
#' @param method_center Method to center residuals: \code{mean}: mean, \code{median}: median.
#' Default: \code{method_center = mean}.
#' @param ... Additional parameters
#'
#' @export
#'
#'
#'
#'


residual_bootstrap <- function(y,
                               y_hat,
                               data_sample,
                               n_boot  = 250,
                               type_boot = c("br1", "br2"),
                               method_scale  = c("sd", "Qn", "mad"),
                               method_center  = c("mean", "median"),
                               boot_seed = 10, ...) {
  # Type boot
  type_boot <- match.arg(type_boot)
  method_scale <- match.arg(method_scale)
  method_center <- match.arg(method_center)


  #y <- data_sample$y
  #y_hat <- y_hat$y_hat_sample
  obj_boot <- list(
    y = y,
    y_hat = y_hat,
    n_boot = n_boot,
    data_sample = data_sample,
    method_scale = method_scale,
    method_center = method_center,
    boot_seed = boot_seed
  )

  class(obj_boot) <- type_boot

  y_boot <- resid_boot(obj_boot)

  return(y_boot)
}
# One-level block residual bootstrap


#' Bootstrap residuals
#'
#' Internal generic function to bootstrap residuals
#'
#' @inheritParams residual_bootstrap
#' @param obj_boot Object to estimate residual bootstrap
#'
#' @importFrom robustbase Qn
#' @importFrom stats mad median sd
#'
#'


resid_boot <- function(...)
  UseMethod("resid_boot")

#'
#' @describeIn resid_boot One-level block residual bootstrap
#' @export
#'


resid_boot.br1 <- function(obj_boot, ...) {
  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  data_sample <- obj_boot$data_sample
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  method_scale <- obj_boot$method_scale
  method_center <- obj_boot$method_center

  resid <- unname(y - y_hat)
  sample_sizes = as.data.frame(table(data_sample$group))$Freq


  boot_ind <- bootstrap_indices(
    sample_sizes = sample_sizes,
    n_boot = n_boot,
    seed = boot_seed
  )


  # Center and scale in each cluster?
  group_names <- as.data.frame(table(data_sample$group))$Var1
  index_group <- which(sample_sizes > 2)
  each_cluster_cs <- length(index_group) == length(group_names)

  df <- data.frame(resid  = resid,
                   group = data_sample$group,
                   y_hat = y_hat)


  if (each_cluster_cs) {

    mean_group <-
      aggregate(df$resid, list(df$group), FUN = method_center)$x
    scale_group <-
      aggregate(df$resid, list(df$group), FUN = method_scale)$x

    resid_cs <-
      (resid - rep(mean_group, sample_sizes)) / (rep(scale_group, sample_sizes))

  } else {
    group_cs <- group_names[index_group]
    group_ncs <- group_names[-index_group]

    # Centered and scaled
    df_cs <- df[is.element(df$group, group_cs), ]
    sample_sizes_cs <-  as.data.frame(table(df_cs$group))$Freq
    mean_group_cs <-
      aggregate(df_cs$resid, list(df_cs$group), FUN = method_center)$x
    scale_group_cs <-
      aggregate(df_cs$resid, list(df_cs$group), FUN = method_scale)$x

    resid_cs <-
      (df_cs$resid - rep(mean_group_cs, sample_sizes_cs)) / (rep(scale_group_cs, sample_sizes_cs))

    df_cs$resid_boot <- resid_cs
 #   y_boot_cs <- df_cs$y_hat + resid_cs
 #    df_cs$y_boot <- y_boot_cs

    # No centering and scaling
    df_ncs <- df[is.element(df$group, group_ncs), ]
    df_ncs$resid_boot <-  df_ncs$resid
    # Merging together and sorting
    dfb <- rbind(df_cs, df_ncs)
    dfbo <- dfb[order(dfb$group), ]

    resid_cs <- dfbo$resid_boot
    #y_boot[[i]]  <- dfbo
  }

  # Center and scale
  y_boot <- list()
  for (i in 1:n_boot) {
    boot_sample_index <- boot_ind[[i]]
    y_boot[[i]]  <- df$y_hat + resid_cs[boot_sample_index]

  }
  return(y_boot)
}


#'
#' @describeIn resid_boot Two-level block residual bootstrap
#' @export
#'


resid_boot.br2 <- function(obj_boot, ...) {
  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  data_sample <- obj_boot$data_sample
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  method_scale <- obj_boot$method_scale
  method_mean <- obj_boot$method_center

  resid <- unname(y - y_hat)
  sample_sizes = as.data.frame(table(data_sample$group))$Freq

  df <- data.frame(resid  = resid,
                   group = data_sample$group,
                   y_hat = y_hat)
  resid2 <- aggregate(df$resid, list(df$group), FUN = mean)$x
  resid1 <- df$resid - rep(resid2, sample_sizes)


  if (method_mean == "mean") {
    mean_r1 <- mean(resid1)
    mean_r2 <- mean(resid2)
  } else {
    mean_r1 <- median(resid1)
    mean_r2 <- median(resid2)
  }

  if (method_scale == "sd") {
    scale_r1 <- sd(resid1)
    scale_r2 <- sd(resid2)
  } else  if (method_scale == "Qn") {
    scale_r1 <- Qn(resid1)
    scale_r2 <- Qn(resid2)
  } else {
    scale_r1 <- mad(resid1)
    scale_r2 <- mad(resid2)
  }
  resid1 <- (resid1 - mean_r1)/scale_r1
  resid2 <- (resid2 - mean_r2)/scale_r2


  boot_ind1 <- sample_bresiduals(
    sample_size = length(resid1),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )

  boot_ind2 <- sample_bresiduals(
    sample_size = length(unique(data_sample$group)),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )

  # Center and scale
  y_boot <- list()
  for (i in 1:n_boot) {
#    boot_sample_index <- boot_ind[[i]]

    y_boot[[i]]  <- unname(y_hat) + resid1[boot_ind1[[i]]] + rep(resid2[boot_ind2[[i]]], sample_sizes)

    }
  return(y_boot)
}





#' Bootstrap indices from the sample
#'
#' @param sample_sizes Data from the sample
#' @param n_boot Number of bootstrap samples
#' @param seed Seed to run simulations.
#'
#' @importFrom sampling strata
#'
#' @return List with following parameters:
#' \item{ind_population}{Vector with population indices}
#' \item{ind_sample}{Vector with sample indices}
#'
#' @export
#'


bootstrap_indices <- function(sample_sizes,
                              n_boot = 50,
                              seed = 1) {
  indices <- list()
  for (i in 1:n_boot) {
    seed_sim = seed * i
    indices[[i]] <- sample_indices(sample_sizes = sample_sizes,
                                   swr = TRUE,
                                   seed = seed_sim * 2)
  }
  return(indices)

}


#' Sample bootstrap residuals
#'
#' @param sample_size Sample sizes
#' @param swr Simple random sampling with replacement? Default = \code{swr = TRUE}.
#' @param seed Seed to run simulations.
#' @param n_boot Number of bootstrap samples
#'
#' @return List with indices.
#'
#' @export


sample_bresiduals <- function(sample_size,
                              n_boot,
                              swr = TRUE,
                              seed = 10) {

  # Set seed
  set.seed(seed)

  list_indices <- list()
  for (i in 1:n_boot) {
    list_indices[[i]] <- sample(1 : sample_size,
                               sample_size,
                               replace = swr)
  }
  return(list_indices)
}


#' Sample bootstrap indices
#'
#' @param sample_sizes Data from the sample
#' @param out_of_sample_sizes Data from the population
#' @param n_boot Number of bootstrap samples
#' @param seed Seed to run simulations.
#' @param type_boot Type of bootstrap sampling
#'
#' @importFrom sampling strata
#'
#' @return List with following parameters:
#' \item{ind_population}{Vector with population indices}
#' \item{ind_sample}{Vector with sample indices}
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
#' no_sim = 1,
#' seed = 1)
#'
#' samples <- generate_sample(populations, ni_size = 10,
#'                            sample_part = "sampled",
#'                            get_index = TRUE)
#'
#' data_sample <- data.frame(samples[[1]]$samp_data)
#' index_sample <- samples[[1]]$index_s
#' data_out_of_sample <- populations[-index_sample, ]
#'
#' sample_sizes = as.data.frame(table(data_sample$group))$Freq
#' out_of_sample_sizes = as.data.frame(table(data_out_of_sample$group))$Freq
#'
#' bootstrap_indices <- sample_bootstrap_indices(sample_sizes,
#'                                                 out_of_sample_sizes,
#'                                                 type_boot = "sample",
#'                                                 n_boot = 500,
#'                                                 seed = 1)


sample_bootstrap_indices <- function(sample_sizes,
                                     out_of_sample_sizes,
                                     type_boot,
                                     n_boot = 50,
                                     seed = 1) {


  obj_boot <- list(sample_sizes = sample_sizes,
                   out_of_sample_sizes = out_of_sample_sizes)
  class(obj_boot) <- type_boot

  boot_indices <- boot_indices_internal(obj_boot,
                                        n_boot = n_boot,
                                        seed = seed)

  return(boot_indices)

}

#'
#' Bootstrap indices
#'
#' Internal generic function to obtain bootstrap indices
#'
#' @inheritParams sample_bootstrap_indices
#' @param obj_boot Object to estimate hte
#' @param ... Additional parameters
#'
#' @return List with following parameters:
#' \item{ind_population}{Vector with population indices}
#' \item{ind_sample}{Vector with sample indices}
#'
#'


boot_indices_internal <- function(...)
  UseMethod("boot_indices_internal")

#'
#' @describeIn boot_indices_internal Obtain bootstrap indices bootstrapping withing the sample
#' @export
#'


boot_indices_internal.sample <- function(obj_boot,
                                         n_boot = 50,
                                         seed = 1,
                                         ...) {

  sample_sizes = obj_boot$sample_sizes
  out_of_sample_sizes = obj_boot$out_of_sample_sizes


  indices <- list()
  for (i in 1:n_boot) {
    seed_sim = seed * i
    indices[[i]] <- list(ind_population = 1 : sum(out_of_sample_sizes),
                         ind_sample = sample_indices(sample_sizes = sample_sizes,
                                                     swr = TRUE,
                                                     seed = seed_sim * 2))
  }
  return(indices)

}


#'
#' @describeIn boot_indices_internal Obtain bootstrap indices bootstrapping withing the sample
#' @export
#'


boot_indices_internal.both <- function(obj_boot,
                                       n_boot = 50,
                                       seed = 1,
                                       ...) {

  sample_sizes = obj_boot$sample_sizes
  out_of_sample_sizes = obj_boot$out_of_sample_sizes


  indices <- list()
  for (i in 1:n_boot) {
    seed_sim = seed * i
    indices[[i]] <- list(ind_population = sample_indices(sample_sizes = out_of_sample_sizes,
                                                         swr = TRUE,
                                                         seed = seed_sim),
                         ind_sample = sample_indices(sample_sizes = sample_sizes,
                                                     swr = TRUE,
                                                     seed = seed_sim * 2))
  }
  return(indices)

}


#'
#' Sample indices
#'
#' @param sample_sizes Sample sizes
#' @param swr Simple random sampling with replacement? Default = \code{swr = TRUE}.
#' @param seed Seed to run simulations.
#'
#' @return Vector with indices.
#'
#' @export


sample_indices <- function(sample_sizes,
                           swr = TRUE,
                           seed = 10) {

  #
  set.seed(seed)

  # Cumulative sum
  cum_sample_sizes <- cumsum(sample_sizes)
  list_indices <- list()
  list_indices[[1]] <- sample(1:cum_sample_sizes[1],
                              sample_sizes[1],
                              replace = swr)

  for (i in 2:length(sample_sizes)) {

    lower_bound <- cum_sample_sizes[i-1] + 1
    upper_bound <- cum_sample_sizes[i]

    list_indices[[i]] <- sample(lower_bound : upper_bound,
                                sample_sizes[i],
                                replace = swr)
  }

  indices <- sort(unlist(list_indices))
  return(indices)
}

