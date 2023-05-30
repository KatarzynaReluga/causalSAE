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
#' @param rand_clust Random clusters
#' @param ... Additional parameters
#'
#' @export
#'
#'
#'
#'


residual_bootstrap <- function(y,
                               y_hat,
#                               y_hat_c = NULL,
                               data_sample,
                               n_boot  = 250,
                               type_boot = c("br1", "br2", "br3","br4", "br5"),
#                               method_scale  = c("sd", "Qn", "mad"),
#                               method_center  = c("mean", "median"),
                               boot_seed = 10,
                               rand_clust = FALSE, ...) {
  # Type boot
  type_boot <- match.arg(type_boot)
#  method_scale <- match.arg(method_scale)
#  method_center <- match.arg(method_center)
  sample_sizes = as.data.frame(table(data_sample$group))$Freq

  obj_boot <- list(
    y = y,
    y_hat = y_hat,
#    y_hat_c = y_hat_c,
    n_boot = n_boot,
    data_sample = data_sample,
#    method_scale = method_scale,
#    method_center = method_center,
    boot_seed = boot_seed,
    rand_clust = rand_clust,
    sample_sizes = sample_sizes
#    y_hat_c = y_hat_c
  )

  class(obj_boot) <- type_boot

  resid_res <- resid_boot(obj_boot)

  boot_ind1 <- resid_res$boot_ind1
  boot_ind2 <- resid_res$boot_ind2

  resid1 <- resid_res$resid1
  resid2 <- resid_res$resid2

  y_boot <- list()

  for (i in 1:n_boot) {
 #   if (is.null(y_hat_c)) {
      y_boot[[i]]  <- unname(y_hat) + resid1[boot_ind1[[i]]] + rep(resid2[boot_ind2[[i]]], sample_sizes)
#    } else {
#      resid1_c <- resid_res$resid1_c
#      resid2_c <- resid_res$resid2_c

#      y_hat <- unname(y_hat) + resid1[boot_ind1[[i]]] + rep(resid2[boot_ind2[[i]]], sample_sizes)
#      y_hat_c <- unname(y_hat_c) + resid1_c[boot_ind1[[i]]] + rep(resid2_c[boot_ind2[[i]]], sample_sizes)

#      y_boot[[i]] <- list(y_hat = y_hat,
#                          y_hat_c = y_hat_c)
#    }

  }

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
#' @describeIn resid_boot Two-level block residual bootstrap
#' @export
#'

#with centering and scaling, but not rescaling, sampling using carpenter
resid_boot.br1 <- function(obj_boot, ...) {
  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  data_sample <- obj_boot$data_sample
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  sample_sizes <- obj_boot$sample_sizes
#  y_hat_c <- obj_boot$y_hat_c

  # Residuals and sample sizes
  resid <- unname(y - y_hat)
#  sample_sizes = as.data.frame(table(data_sample$group))$Freq

  df <- data.frame(resid  = resid,
                   group = data_sample$group,
                   y_hat = y_hat)

  # Cluster residuals
  resid20 <- aggregate(df$resid, list(df$group), FUN = mean)$x
  # Individual residuals
  resid10 <- df$resid - rep(resid20, sample_sizes)


  mean_r1 <- mean(resid10)
  mean_r2 <- mean(resid20)

  scale_r1 <- sd(resid10)
  scale_r2 <- sd(resid20)

  resid1 <- (resid10 - mean_r1)/scale_r1
  resid2 <- (resid20 - mean_r2)/scale_r2

 # Sample from the pool of level-one residuals
  boot_ind1 <- sample_bresiduals(
    sample_size = length(resid1),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )
  # Sample from the pool of level-two residuals
  boot_ind2 <- sample_bresiduals(
    sample_size = length(resid2),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )

  output <- list(resid1 = resid1,
                 resid2 = resid2,

                 boot_ind1 = boot_ind1,
                 boot_ind2 = boot_ind2)



  return(output)

}



#'
#' @describeIn resid_boot Two-level block residual bootstrap
#' @export
#'
#with centering and scaling, and rescaling using Q
#two cases because of rand_clust
resid_boot.br2 <- function(obj_boot, ...) {


  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  data_sample <- obj_boot$data_sample
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
#  method_scale <- obj_boot$method_scale
#  method_mean <- obj_boot$method_center
  rand_clust <- obj_boot$rand_clust
  sample_sizes <- obj_boot$sample_sizes
#  y_hat_c <- obj_boot$y_hat_c

  resid <- unname(y - y_hat)
#  sample_sizes = as.data.frame(table(data_sample$group))$Freq
  group_names <- as.data.frame(table(data_sample$group))$Var1

  df <- data.frame(resid  = resid,
                   group = data_sample$group,
                   y_hat = y_hat)
  resid20 <- aggregate(df$resid, list(df$group), FUN = mean)$x
  resid10 <- df$resid - rep(resid20, sample_sizes)

  mean_r1 <- mean(resid10)
  mean_r2 <- mean(resid20)

  sd_r1 <- sd(resid10)
  sd_r2 <- sd(resid20)

  Q_r1 <- Qn(resid10)
  Q_r2 <- Qn(resid20)

  resid1 <- Q_r1/sd_r1 * (resid10 - mean_r1)
  resid2 <- Q_r2/sd_r2 * (resid20 - mean_r2)


  boot_ind2 <- sample_bresiduals(
    sample_size = length(resid2),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )


  boot_ind1 <- bootstrap_indices(
    sample_sizes = sample_sizes,
    n_boot = n_boot,
    seed = boot_seed,
    rand_clust = rand_clust
  )


  output <- list(resid1 = resid1,
                 resid2 = resid2,

                 boot_ind1 = boot_ind1,
                 boot_ind2 = boot_ind2)


#  if (!is.null(y_hat_c)) {
    # Residuals and sample sizes
#    resid_c <- unname(y - y_hat_c)
#    df_c <- data.frame(resid  = resid_c,
#                       group = data_sample$group,
#                       y_hat = y_hat_c)

    # Cluster residuals
#    resid20_c <- aggregate(df_c$resid, list(df_c$group), FUN = mean)$x
    # Individual residuals
#    resid10_c <- df$resid_c - rep(resid20_c, sample_sizes)

#    mean_r1_c <- mean(resid10_c)
#    mean_r2_c <- mean(resid20_c)

#    scale_r1_c <- sd(resid10_c)
#    scale_r2_c <- sd(resid20_c)

#    Q_r1_c <- Qn(resid10_c)
#    Q_r2_c <- Qn(resid20_c)

#    resid1 <- Q_r1_c/scale_r1_c * (resid10_c - mean_r1_c)
#    resid2 <- Q_r1_c/scale_r2_c * (resid20_c - mean_r2_c)

#    output$resid1_c <- resid1_c
#    output$resid2_c <- resid2_c
#  }

  return(output)

}


#'
#' @describeIn resid_boot Two-level block residual bootstrap
#' @export
#'

#with centering and scaling, and rescaling using Q, but sampling collectively from
# the pool of residulas at the level 1 and 2 (aka Cerpenter method)
resid_boot.br3 <- function(obj_boot, ...) {

  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  data_sample <- obj_boot$data_sample
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  sample_sizes <- obj_boot$sample_sizes
#  method_scale <- obj_boot$method_scale
#  method_mean <- obj_boot$method_center
#  rand_clust <- obj_boot$rand_clust
#  y_hat_c <- obj_boot$y_hat_c

  resid <- unname(y - y_hat)
#  sample_sizes = as.data.frame(table(data_sample$group))$Freq
  group_names <- as.data.frame(table(data_sample$group))$Var1

  df <- data.frame(resid  = resid,
                   group = data_sample$group,
                   y_hat = y_hat)
  resid20 <- aggregate(df$resid, list(df$group), FUN = mean)$x
  resid10 <- df$resid - rep(resid20, sample_sizes)


  #  if (method_mean == "mean") {
  mean_r1 <- mean(resid10)
  mean_r2 <- mean(resid20)

  sd_r1 <- sd(resid10)
  sd_r2 <- sd(resid20)

  Q_r1 <- Qn(resid10)
  Q_r2 <- Qn(resid20)

  resid1 <- Q_r1/sd_r1 * (resid10 - mean_r1)
  resid2 <- Q_r2/sd_r2 * (resid20 - mean_r2)


  boot_ind2 <- sample_bresiduals(
    sample_size = length(resid2),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )


  boot_ind1 <- sample_bresiduals(
    sample_size = length(resid1),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )


  output <- list(resid1 = resid1,
                 resid2 = resid2,

                 boot_ind1 = boot_ind1,
                 boot_ind2 = boot_ind2)


#  if (!is.null(y_hat_c)) {
    # Residuals and sample sizes
#    resid_c <- unname(y - y_hat_c)
#    df_c <- data.frame(resid  = resid_c,
#                       group = data_sample$group,
#                       y_hat = y_hat_c)

    # Cluster residuals
#    resid20_c <- aggregate(df_c$resid, list(df_c$group), FUN = mean)$x
    # Individual residuals
#    resid10_c <- df$resid_c - rep(resid20_c, sample_sizes)

#    mean_r1_c <- mean(resid10_c)
#    mean_r2_c <- mean(resid20_c)

#    scale_r1_c <- sd(resid10_c)
#    scale_r2_c <- sd(resid20_c)

#    Q_r1_c <- Qn(resid10_c)
#    Q_r2_c <- Qn(resid20_c)

#    resid1 <- Q_r1_c/scale_r1_c * (resid10_c - mean_r1_c)
#    resid2 <- Q_r1_c/scale_r2_c * (resid20_c - mean_r2_c)

#    output$resid1_c <- resid1_c
#    output$resid2_c <- resid2_c
#  }

  return(output)

}



#'
#' @describeIn resid_boot Two-level block residual bootstrap
#' @export
#'

# without centering and scaling, sampling aka Carpenter from the second level,
# and sampling cases from the the level-one residuals
# two cases because of rand_clust
resid_boot.br4 <- function(obj_boot, ...) {

  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  data_sample <- obj_boot$data_sample
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  rand_clust <- obj_boot$rand_clust
  sample_sizes <- obj_boot$sample_sizes
#  y_hat_c <- obj_boot$y_hat_c

  resid <- unname(y - y_hat)
#  sample_sizes = as.data.frame(table(data_sample$group))$Freq

  df <- data.frame(resid  = resid,
                   group = data_sample$group,
                   y_hat = y_hat)
  resid2 <- aggregate(df$resid, list(df$group), FUN = mean)$x
  resid1 <- df$resid - rep(resid2, sample_sizes)

  boot_ind2 <- sample_bresiduals(
    sample_size = length(resid2),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )


  boot_ind1 <- bootstrap_indices(
    sample_sizes = sample_sizes,
    n_boot = n_boot,
    seed = boot_seed,
    rand_clust = rand_clust
  )

  output <- list(resid1 = resid1,
                 resid2 = resid2,

                 boot_ind1 = boot_ind1,
                 boot_ind2 = boot_ind2)


#  if (!is.null(y_hat_c)) {
    # Residuals and sample sizes
#    resid_c <- unname(y - y_hat_c)
#    df_c <- data.frame(resid  = resid_c,
#                       group = data_sample$group,
#                       y_hat = y_hat_c)

    # Cluster residuals
#    resid2_c <- aggregate(df_c$resid, list(df$group), FUN = mean)$x
    # Individual residuals
#    resid1_c <- df$resid_c - rep(resid2_c, sample_sizes)

#    output$resid1_c <- resid1_c
#    output$resid2_c <- resid2_c
#  }
  return(output)
 }



#'
#' @describeIn resid_boot Two-level block residual bootstrap
#' @export
#'

# without centering and scaling, sampling aka Carpenter level-one and level-two

resid_boot.br5 <- function(obj_boot, ...) {

  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  data_sample <- obj_boot$data_sample
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  sample_sizes <- obj_boot$sample_sizes
#  rand_clust <- obj_boot$rand_clust
#  y_hat_c <- obj_boot$y_hat_c

  resid <- unname(y - y_hat)
#  sample_sizes = as.data.frame(table(data_sample$group))$Freq

  df <- data.frame(resid  = resid,
                   group = data_sample$group,
                   y_hat = y_hat)
  resid2 <- aggregate(df$resid, list(df$group), FUN = mean)$x
  resid1 <- df$resid - rep(resid2, sample_sizes)

  boot_ind2 <- sample_bresiduals(
    sample_size = length(resid2),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )


  boot_ind1 <- sample_bresiduals(
    sample_size = length(resid1),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )

  output <- list(resid1 = resid1,
                 resid2 = resid2,

                 boot_ind1 = boot_ind1,
                 boot_ind2 = boot_ind2)


#  if (!is.null(y_hat_c)) {
    # Residuals and sample sizes
#    resid_c <- unname(y - y_hat_c)
#    df_c <- data.frame(resid  = resid_c,
#                       group = data_sample$group,
#                       y_hat = y_hat_c)

    # Cluster residuals
#    resid2_c <- aggregate(df_c$resid, list(df$group), FUN = mean)$x
    # Individual residuals
#    resid1_c <- df$resid_c - rep(resid2_c, sample_sizes)

#    output$resid1_c <- resid1_c
#    output$resid2_c <- resid2_c
#  }
  return(output)
}




#' Bootstrap indices from the sample
#'
#' @param sample_sizes Data from the sample
#' @param n_boot Number of bootstrap samples
#' @param seed Seed to run simulations.
#' @param rand_clust Select clusters randomly. Default: \code{rand_clust = FALSE}.
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
                              seed = 1,
                              rand_clust) {
  indices <- list()
  for (i in 1:n_boot) {
    seed_sim = seed * i
    indices[[i]] <- sample_indices(sample_sizes = sample_sizes,
                                   swr = TRUE,
                                   seed = seed_sim * 2,
                                   rand_clust = rand_clust)
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



#'
#' Sample indices
#'
#' @param sample_sizes Sample sizes
#' @param swr Simple random sampling with replacement? Default = \code{swr = TRUE}.
#' @param seed Seed to run simulations.
#' @param rand_clust Select clusters randomly. Default: \code{rand_clust = FALSE}.
#'
#' @return Vector with indices.
#'
#' @export

# add random bootstrap clustering
sample_indices <- function(sample_sizes,
                           swr = TRUE,
                           seed = 10,
                           rand_clust = FALSE) {

  #
  set.seed(seed)

  # Cumulative sum
  cum_sample_sizes <- cumsum(sample_sizes)

  # List of sample sizes
  list_lower_upper <- list()
  list_lower_upper[[1]] <- 1:cum_sample_sizes[1]

  for (i in 2:length(sample_sizes)) {

    lower_bound <- cum_sample_sizes[i-1] + 1
    upper_bound <- cum_sample_sizes[i]
    list_lower_upper[[i]] <- lower_bound : upper_bound
  }


  if (rand_clust) {
    clusters <- sample(1:length(sample_sizes), replace = TRUE)
  } else {
    clusters <- c(1:length(sample_sizes))
  }
    # List of indices
    list_indices <- list()
    for (i in 1:length(sample_sizes)) {

      list_indices[[i]] <- sample(list_lower_upper[[clusters[i]]],
                                  sample_sizes[i],
                                  replace = swr)
    }



#  list_indices <- list()
#  list_indices[[1]] <- sample(1:cum_sample_sizes[1],
#                              sample_sizes[1],
#                              replace = swr)

#  for (i in 2:length(sample_sizes)) {

#    lower_bound <- cum_sample_sizes[i-1] + 1
#    upper_bound <- cum_sample_sizes[i]

#    list_indices[[i]] <- sample(lower_bound : upper_bound,
#                                sample_sizes[i],
#                                replace = swr)
#  }

    indices <- unlist(list_indices)
#  indices <- do.call(rbind, list_indices)
#  output <- list(indices = indices,
#                 clusters = clusters)
  return(indices)
}

