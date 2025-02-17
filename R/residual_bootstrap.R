#' Residual bootstrap
#'
#' One-level and two-level random block residual bootstrap for variance calculations
#'
# #' @param y Vector of true outcomes
#' @param y_hat Vector of estimated outcomes
# #' @param A treatment status 
#' @param data_sample Sample data 
#' @param boot_seed Seed to repeat bootstrap sampling
#' @param n_boot Number of bootstrap samples
#' @param type_boot Type of bootstrap: \code{br1} - one-level block residual bootstrap,
#' \code{br2} - two-level block residual bootstrap.
#' @param cent_resid Center residuals? Default: TRUE
#' @param c_Win Constant to Winsorize residuals at c_rob times robust scale estimate
#' @param ... Additional parameters
#'
#' @export
#'
#'
#'
#'


residual_bootstrap <- function(
    #y,
                             y_hat,
#                               A = NULL, 
                               data_sample,
                               n_boot  = 250,
                               type_boot = c("br1", "br2", "br3", "br4", "br5", "br6", "br7"),
#                               method_scale  = c("sd", "Qn", "mad"),
                               boot_seed = 10,
                               cent_resid = TRUE, 
                               c_Win = 3, ...) {
  # Type boot
  type_boot <- match.arg(type_boot)
  group_names = as.data.frame(table(data_sample$group))$Var1
  A = data_sample$A
  group = data_sample$group
#  method_scale <- match.arg(method_scale)
#  method_center <- match.arg(method_center)
#  sample_sizes = as.data.frame(table(data_sample$group))$Freq

  obj_boot <- list(
   y = data_sample$y,
   y_hat = y_hat,
   A = A,
   n_boot = n_boot,
   group = group,
#    data_sample = data_sample,
#    method_scale = method_scale,
#    method_center = method_center,
    boot_seed = boot_seed,
    sample_sizes = as.data.frame(table(data_sample$group))$Freq,
    cent_resid = cent_resid, 
    c_Win = c_Win
  )

  class(obj_boot) <- type_boot

  resid_res <- resid_boot(obj_boot)
  return(resid_res)
  
  # Create bootstrap versions of y
  
  # # Load indices from all bootstrap samples
  # boot_ind10 <- resid_res$boot_ind10
  # boot_ind11 <- resid_res$boot_ind11
  # 
  # boot_ind20 <- resid_res$boot_ind21
  # boot_ind21 <- resid_res$boot_ind20
  # 
  # # Load residuals
  # resid10 <- resid_res$resid10
  # resid11 <- resid_res$resid11
  # 
  # resid20 <- resid_res$resid20
  # resid21 <- resid_res$resid21
  # 
  # # Create y_boot by adding residuals to each y_hat
  # 
  # y_boot <- list()
  # if (is.null(resid20)) {
  #   for (i in 1:n_boot) {
  #     resid1 <- numeric(length(y_hat))
  #     resid1[A == 0] <- resid10[boot_ind10[[i]]]
  #     resid1[A == 1] <- resid11[boot_ind11[[i]]]
  #     
  #     y_boot[[i]]  <- unname(y_hat) + resid1
  #   }
  #   
  # } else { 
  #   
  #   for (i in 1:n_boot) {
  #     resid1 <- numeric(length(y_hat))
  #     resid1[A == 0] <- resid10[boot_ind10[[i]]]
  #     resid1[A == 1] <- resid11[boot_ind11[[i]]]
  #     
  #     resid2 <- numeric(length(y_hat))
  #     for (j in 1:length(group_names)) {
  #       which_group0 <- which(A == 0 & group == group_names[j])
  #       resid2[which_group0] <- resid20[boot_ind20[[i]][j]]
  #       which_group1 <- which(A == 1 & group == group_names[j])
  #       resid2[which_group1] <- resid21[boot_ind21[[i]][j]]
  #     }
  #     
  #     
  #     y_boot[[i]]  <- unname(y_hat) + resid1 + resid2
  #   }
  #   
  # }
  # 
  # list_results  = list(y_boot = y_boot, 
  #                      
  #                      resid10 = resid10,
  #                      resid11 = resid11,
  #                      
  #                      resid20 = resid20,
  #                      resid21 = resid21,
  #                      
  #                      boot_ind10 = boot_ind10,
  #                      boot_ind11 = boot_ind11,
  #                      
  #                      boot_ind20 = boot_ind20, 
  #                      boot_ind21 = boot_ind21 
  #   
  # )
  # return(list_results)
}


#' Bootstrap residuals
#'
#' Internal generic function to bootstrap residuals
#'
#' @inheritParams residual_bootstrap
#' @param obj_boot Object to estimate residual bootstrap
#'
#' @importFrom robustbase Qn
#' @importFrom stats mad median sd optim var
#'
#'


resid_boot <- function(...)
  UseMethod("resid_boot")



#'
#' @describeIn resid_boot Two-level block residual bootstrap without re-scaling (sample from level one and level two residuals)
#' @export
#'

resid_boot.br1 <- function(obj_boot, ...) {
  
  
  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  # data_sample <- obj_boot$data_sample
  group = obj_boot$group
  group_names <- as.data.frame(table(group))$group
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  sample_sizes <- obj_boot$sample_sizes
  cent_resid <- obj_boot$cent_resid
  c_Win <- obj_boot$c_Win
  
  resid <- unname(y - y_hat)
  
  df <- data.frame(resid  = resid,
                   group = group,
                   y_hat = y_hat)
  resid20 <- aggregate(df$resid, list(df$group), FUN = mean)$x
  resid10 <- df$resid - rep(resid20, sample_sizes)
  
  
  if (cent_resid) {
    mean_r1 <- mean(resid10)
    mean_r2 <- mean(resid20)
    
    resid1 <- resid10 - mean_r1
    resid2 <- resid20 - mean_r2
  } else {
    resid1 <- resid10
    resid2 <- resid20
  }
  
  
  
  # Sample indices level-1 residuals
  boot_ind1 <- sample_bresiduals(
    sample_size = length(resid1),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )
  
  # Sample indices level-2 residuals
  boot_ind2 <- sample_bresiduals(
    sample_size = length(resid2),
    n_boot = n_boot,
    seed = boot_seed * 2,
    swr = TRUE
  )
  
  
  # Load indices from all bootstrap samples
#  boot_ind1 <- resid_res$boot_ind1
#  boot_ind2 <- resid_res$boot_ind2
  
  # Load residuals
#  resid1 <- resid_res$resid1
#  resid2 <- resid_res$resid2
  
  # Create y_boot by adding residuals to each y_hat
  
  y_boot <- list()
  for (i in 1:n_boot) {
    y_boot[[i]]  <- unname(y_hat) + resid1[boot_ind1[[i]]] + rep(resid2[boot_ind2[[i]]], sample_sizes)
  }
  
  # Save indices and residulas
  
  list_results  = list(y_boot = y_boot, 
                       
                       resid1 = resid1,
                       resid2 = resid2,
                 
                       boot_ind1 = boot_ind1,
                       boot_ind2 = boot_ind2)
  
  
  
  return(list_results)
  
}





#'
#' @describeIn resid_boot One-level block residual bootstrap by treatment group 
#' @export
#'

resid_boot.br2 <- function(obj_boot, ...) {
  
  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  A <- obj_boot$A
  group = obj_boot$group
  group_names <- as.data.frame(table(group))$group
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  sample_sizes <- obj_boot$sample_sizes
  cent_resid <- obj_boot$cent_resid
  
  # Prepare residuals for bootstrapping
  ## Extract total estimation errors 
  resid <- unname(y - y_hat)
  
  ##  Random effects 
  A0_index <- which(A == 0)
  A1_index <- which(A == 1)
  
  ### Control, level 2 
  # df0 <- data.frame(resid  = resid[A0_index],
  #                   group = group[A0_index],
  #                   y_hat = y_hat[A0_index])
  # resid20_prelim <- aggregate(df0$resid, list(df0$group), FUN = mean)$x
  # 
  # ### Treated, level 2 
  # df1 <- data.frame(resid  = resid[A1_index],
  #                   group = group[A1_index],
  #                   y_hat = y_hat[A1_index])
  # resid21_prelim <- aggregate(df1$resid, list(df1$group), FUN = mean)$x
  # 
  # resid2_full <- numeric(length(A))
  # 
  # for (i in 1:length(group_names)) {
  #   which_group0 <- which(A == 0 & group == group_names[i])
  #   resid2_full[which_group0] <- resid20_prelim[i]
  #   which_group1 <- which(A == 1 & group == group_names[i])
  #   resid2_full[which_group1] <- resid21_prelim[i]
  # }
  
  ## Total error minus area random effects 
#  resid1_full <- resid - resid2_full
  
  ### Unit-level random effect
  resid10_prelim <- resid[A == 0] 
  resid11_prelim <- resid[A == 1]
  
  if (cent_resid) {
    #Level 1
    mean_r10 <- mean(resid10_prelim)
    mean_r11 <- mean(resid11_prelim)
    
    resid10 <- resid10_prelim - mean_r10
    resid11 <- resid11_prelim - mean_r11
    
    # Level 2
    # mean_r20 <- mean(resid20_prelim)
    # mean_r21 <- mean(resid21_prelim)
    # 
    # 
    # resid20 <- resid20_prelim - mean_r20
    # resid21 <- resid21_prelim - mean_r21
  } else {
    resid10 <- resid10_prelim
    resid11 <- resid11_prelim
    
    # resid20 <- resid20_prelim
    # resid21 <- resid21_prelim
  }
  
  
  
  # Sample indices level-1 residuals, control
  boot_ind10 <- sample_bresiduals(
    sample_size = length(resid10),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )
  
  # Sample indices level-1 residuals, treated
  boot_ind11 <- sample_bresiduals(
    sample_size = length(resid11),
    n_boot = n_boot,
    seed = boot_seed * 2,
    swr = TRUE
  )
  
  # Sample indices level-2 residuals, control
  # boot_ind20 <- sample_bresiduals(
  #   sample_size = length(resid20),
  #   n_boot = n_boot,
  #   seed = boot_seed * 3,
  #   swr = TRUE
  # )
  # 
  # # Sample indices level-2 residuals, treated
  # boot_ind21 <- sample_bresiduals(
  #   sample_size = length(resid21),
  #   n_boot = n_boot,
  #   seed = boot_seed * 4,
  #   swr = TRUE
  # )
  
  # Create bootstrap versions of y
  
  # Load indices from all bootstrap samples
  # boot_ind10 <- resid_res$boot_ind10
  # boot_ind11 <- resid_res$boot_ind11
  # 
  # boot_ind20 <- resid_res$boot_ind21
  # boot_ind21 <- resid_res$boot_ind20
  
  # Load residuals
  # resid10 <- resid_res$resid10
  # resid11 <- resid_res$resid11
  # 
  # resid20 <- resid_res$resid20
  # resid21 <- resid_res$resid21
  
  # Create y_boot by adding residuals to each y_hat
  
  y_boot <- list()
 # if (is.null(resid20)) {
    for (i in 1:n_boot) {
      resid1 <- numeric(length(y_hat))
      resid1[A == 0] <- resid10[boot_ind10[[i]]]
      resid1[A == 1] <- resid11[boot_ind11[[i]]]
      
      y_boot[[i]]  <- unname(y_hat) + resid1

    }
  
  list_results  = list(y_boot = y_boot, 
                       
                       resid10 = resid10,
                       resid11 = resid11,
                       
                       # resid20 = resid20,
                       # resid21 = resid21,
                       
                       boot_ind10 = boot_ind10,
                       boot_ind11 = boot_ind11)
                       
                       # boot_ind20 = boot_ind20, 
                       # boot_ind21 = boot_ind21 
                       # 

  return(list_results)
  
  
  # Save indices and residulas
  # output <- list(resid10 = resid10,
  #                resid11 = resid11,
  #                
  #                resid21 = NULL,
  #                resid20 = NULL,
  #                
  #                boot_ind10 = boot_ind10,
  #                boot_ind11 = boot_ind11,
  #                
  #                boot_ind20 = NULL,
  #                boot_ind21 = NULL)
  
  
  
 # return(list_results)
  
  
}




#'
#' @describeIn resid_boot Two-level block residual bootstrap by treatment group 
#' @export
#'

resid_boot.br3 <- function(obj_boot, ...) {
  
  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  A <- obj_boot$A
  group = obj_boot$group
  group_names <- as.data.frame(table(group))$group
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  sample_sizes <- obj_boot$sample_sizes
  cent_resid <- obj_boot$cent_resid
  
  # Prepare residuals for bootstrapping
  ## Extract total estimation errors 
  resid <- unname(y - y_hat)
  
  ##  Random effects 
  A0_index <- which(A == 0)
  A1_index <- which(A == 1)
  
  ### Control, level 2 
  df0 <- data.frame(resid  = resid[A0_index],
                   group = group[A0_index],
                   y_hat = y_hat[A0_index])
  resid20_prelim <- aggregate(df0$resid, list(df0$group), FUN = mean)$x
  ### Treated, level 2 
  df1 <- data.frame(resid  = resid[A1_index],
                    group = group[A1_index],
                    y_hat = y_hat[A1_index])
  resid21_prelim <- aggregate(df1$resid, list(df1$group), FUN = mean)$x
  
  resid2_full <- numeric(length(A))
  
  for (i in 1:length(group_names)) {
    which_group0 <- which(A == 0 & group == group_names[i])
    resid2_full[which_group0] <- resid20_prelim[i]
    which_group1 <- which(A == 1 & group == group_names[i])
    resid2_full[which_group1] <- resid21_prelim[i]
  }

  resid20_prelim <- Qn(resid20_prelim)/sd(resid20_prelim) * resid20_prelim
  resid21_prelim <- Qn(resid21_prelim)/sd(resid21_prelim) * resid21_prelim
  
  ## Total error minus area random effects 
  resid1_full <- resid - resid2_full
  
  ### Unit-level random effect
#  resid10_prelim <- resid1_full[A == 0] 
#  resid11_prelim <- resid1_full[A == 1] 

  resid10_prelim <- Qn(resid1_full[A == 0])/sd(resid1_full[A == 0]) * resid1_full[A == 0]
  resid11_prelim <- Qn(resid1_full[A == 1])/sd(resid1_full[A == 1]) * resid1_full[A == 1]
  
    
  if (cent_resid) {
    #Level 1
    mean_r10 <- mean(resid10_prelim)
    mean_r11 <- mean(resid11_prelim)
    
    resid10 <- resid10_prelim - mean_r10
    resid11 <- resid11_prelim - mean_r11
    
    # Level 2
    mean_r20 <- mean(resid20_prelim)
    mean_r21 <- mean(resid21_prelim)
    
    
    resid20 <- resid20_prelim - mean_r20
    resid21 <- resid21_prelim - mean_r21
  } else {
    resid10 <- resid10_prelim
    resid11 <- resid11_prelim
    
    resid20 <- resid20_prelim
    resid21 <- resid21_prelim
  }
  
  
  
  # Sample indices level-1 residuals, control
  boot_ind10 <- sample_bresiduals(
    sample_size = length(resid10),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )
  
  # Sample indices level-1 residuals, treated
  boot_ind11 <- sample_bresiduals(
    sample_size = length(resid11),
    n_boot = n_boot,
    seed = boot_seed * 2,
    swr = TRUE
  )
  
  # Sample indices level-2 residuals, control
  boot_ind20 <- sample_bresiduals(
    sample_size = length(resid20),
    n_boot = n_boot,
    seed = boot_seed * 3,
    swr = TRUE
  )
  
  # Sample indices level-2 residuals, treated
  boot_ind21 <- sample_bresiduals(
    sample_size = length(resid21),
    n_boot = n_boot,
    seed = boot_seed * 4,
    swr = TRUE
  )
  
  y_boot <- list()
  
  for (i in 1:n_boot) {
    resid1 <- numeric(length(y_hat))
    resid1[A == 0] <- resid10[boot_ind10[[i]]]
    resid1[A == 1] <- resid11[boot_ind11[[i]]]
    
    resid2 <- numeric(length(y_hat))
    for (j in 1:length(group_names)) {
      which_group0 <- which(A == 0 & group == group_names[j])
      resid2[which_group0] <- resid20[boot_ind20[[i]][j]]
      which_group1 <- which(A == 1 & group == group_names[j])
      resid2[which_group1] <- resid21[boot_ind21[[i]][j]]
    }
    
    
    y_boot[[i]]  <- unname(y_hat) + resid1 + resid2
  }
  
  # Save indices and residulas
  list_results <- list(y_boot = y_boot, 
                 
                       resid10 = resid10,
                   resid11 = resid11,
                 
                   resid21 = resid21,
                   resid20 = resid20,
                 
                   boot_ind10 = boot_ind10,
                   boot_ind11 = boot_ind11,
                   boot_ind20 = boot_ind20,
                   boot_ind21 = boot_ind21)
  
  
  
  return(list_results)
  
}



#'
#' @describeIn resid_boot Two-level block residual bootstrap by treatment group 
#' @export
#'

resid_boot.br4 <- function(obj_boot, ...) {
  
  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  A <- obj_boot$A
  group = obj_boot$group
  group_names <- as.data.frame(table(group))$group
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  sample_sizes <- obj_boot$sample_sizes
  cent_resid <- obj_boot$cent_resid
  
  # Prepare residuals for bootstrapping
  ## Extract total estimation errors 
  resid <- unname(y - y_hat)
  df <- data.frame(resid  = resid,
                   group = group,
                   y_hat = y_hat)
  resid20 <- aggregate(df$resid, list(df$group), FUN = mean)$x
#  resid10 <- df$resid - rep(resid20, sample_sizes)
  
  
  ### Control, level 2 
  # df0 <- data.frame(resid  = resid[A0_index],
  #                   group = group[A0_index],
  #                   y_hat = y_hat[A0_index])
  # resid20_prelim <- aggregate(df0$resid, list(df0$group), FUN = mean)$x
  # 
  # ### Treated, level 2 
  # df1 <- data.frame(resid  = resid[A1_index],
  #                   group = group[A1_index],
  #                   y_hat = y_hat[A1_index])
  # resid21_prelim <- aggregate(df1$resid, list(df1$group), FUN = mean)$x
  
  resid2_full <- numeric(length(A))
  
  for (i in 1:length(group_names)) {
    which_group <- which(group == group_names[i])
    resid2_full[which_group] <- resid20[i]
    # which_group1 <- which(A == 1 & group == group_names[i])
    # resid2_full[which_group1] <- resid21_prelim[i]
  }
  
  ## Total error minus area random effects 
  resid1_full <- resid - resid2_full
  
  
  ##  Random effects 
  A0_index <- which(A == 0)
  A1_index <- which(A == 1)
  
  
  ### Unit-level random effect
  resid10_prelim <- resid1_full[A == 0] 
  resid11_prelim <- resid1_full[A == 1] 
  
  if (cent_resid) {
    #Level 1
    mean_r10 <- mean(resid10_prelim)
    mean_r11 <- mean(resid11_prelim)
    
    resid10 <- resid10_prelim - mean_r10
    resid11 <- resid11_prelim - mean_r11
    
    # Level 2
    mean_r2 <- mean(resid20)
#    mean_r21 <- mean(resid21_prelim)
    
    
    resid2 <- resid20 - mean_r2
#    resid21 <- resid21_prelim - mean_r21
  } else {
    resid10 <- resid10_prelim
    resid11 <- resid11_prelim
    
    resid2 <- resid20
#    resid21 <- resid21_prelim
  }
  
  
  
  # Sample indices level-1 residuals, control
  boot_ind10 <- sample_bresiduals(
    sample_size = length(resid10),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )
  
  # Sample indices level-1 residuals, treated
  boot_ind11 <- sample_bresiduals(
    sample_size = length(resid11),
    n_boot = n_boot,
    seed = boot_seed * 2,
    swr = TRUE
  )
  
  # Sample indices level-2 residuals, control
  boot_ind2 <- sample_bresiduals(
    sample_size = length(resid2),
    n_boot = n_boot,
    seed = boot_seed * 3,
    swr = TRUE
  )
  
  # Sample indices level-2 residuals, treated
  # boot_ind21 <- sample_bresiduals(
  #   sample_size = length(resid21),
  #   n_boot = n_boot,
  #   seed = boot_seed * 4,
  #   swr = TRUE
  # )
  # 
  y_boot <- list()
  
  for (i in 1:n_boot) {
    resid1 <- numeric(length(y_hat))
    resid1[A == 0] <- resid10[boot_ind10[[i]]]
    resid1[A == 1] <- resid11[boot_ind11[[i]]]
    
    resid2b <- numeric(length(y_hat))
    for (j in 1:length(group_names)) {
      which_group <- which(group == group_names[j])
      resid2b[which_group] <- resid2[boot_ind2[[i]][j]]
 #     which_group1 <- which(A == 1 & group == group_names[j])
#      resid2[which_group1] <- resid21[boot_ind21[[i]][j]]
    }
    
    
    y_boot[[i]]  <- unname(y_hat) + resid1 + resid2b
  }
  
  # Save indices and residulas
  list_results <- list(y_boot = y_boot, 
                       
                       resid10 = resid10,
                       resid11 = resid11,
                       
                       resid2 = resid2,
#                       resid20 = resid20,
                       
                       boot_ind10 = boot_ind10,
                       boot_ind11 = boot_ind11,
                       boot_ind2 = boot_ind2)
#                       boot_ind21 = boot_ind21)
  
  
  
  return(list_results)
  
}



#'
#' @describeIn resid_boot Two-level block residual bootstrap by treatment group 
#' @export
#'

resid_boot.br5 <- function(obj_boot, ...) {
  
  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  A <- obj_boot$A
  group = obj_boot$group
  group_names <- as.data.frame(table(group))$group
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  sample_sizes <- obj_boot$sample_sizes
  cent_resid <- obj_boot$cent_resid
  
  # Prepare residuals for bootstrapping
  ## Extract total estimation errors 
  resid <- unname(y - y_hat)
  
  ##  Random effects 
  A0_index <- which(A == 0)
  A1_index <- which(A == 1)
  
  ### Control, level 2 
  df0 <- data.frame(resid  = resid[A0_index],
                    group = group[A0_index],
                    y_hat = y_hat[A0_index])
  resid20_prelim <- aggregate(df0$resid, list(df0$group), FUN = mean)$x
  
  ### Treated, level 2 
  df1 <- data.frame(resid  = resid[A1_index],
                    group = group[A1_index],
                    y_hat = y_hat[A1_index])
  resid21_prelim <- aggregate(df1$resid, list(df1$group), FUN = mean)$x
  
  resid2_full <- numeric(length(A))
  
  for (i in 1:length(group_names)) {
    which_group0 <- which(A == 0 & group == group_names[i])
    resid2_full[which_group0] <- resid20_prelim[i]
    which_group1 <- which(A == 1 & group == group_names[i])
    resid2_full[which_group1] <- resid21_prelim[i]
  }
  
  ## Total error minus area random effects 
  resid1_full <- resid - resid2_full
  
  ### Unit-level random effect
  resid1_prelim <- resid1_full
#  resid11_prelim <- resid1_full[A == 1] 
  
  if (cent_resid) {
    #Level 1
    mean_r1 <- mean(resid1_prelim)
#    mean_r11 <- mean(resid11_prelim)
    
    resid1 <- resid1_prelim - mean_r1
#    resid11 <- resid11_prelim - mean_r11
    
    # Level 2
    mean_r20 <- mean(resid20_prelim)
    mean_r21 <- mean(resid21_prelim)
    
    
    resid20 <- resid20_prelim - mean_r20
    resid21 <- resid21_prelim - mean_r21
  } else {
    resid1 <- resid1_prelim
#    resid11 <- resid11_prelim
    
    resid20 <- resid20_prelim
    resid21 <- resid21_prelim
  }
  
  
  
  # Sample indices level-1 residuals, control
  boot_ind1 <- sample_bresiduals(
    sample_size = length(resid1),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )
  
  # Sample indices level-1 residuals, treated
#  boot_ind11 <- sample_bresiduals(
#    sample_size = length(resid11),
#    n_boot = n_boot,
#    seed = boot_seed * 2,
#    swr = TRUE
#  )
  
  # Sample indices level-2 residuals, control
  boot_ind20 <- sample_bresiduals(
    sample_size = length(resid20),
    n_boot = n_boot,
    seed = boot_seed * 3,
    swr = TRUE
  )
  
  # Sample indices level-2 residuals, treated
  boot_ind21 <- sample_bresiduals(
    sample_size = length(resid21),
    n_boot = n_boot,
    seed = boot_seed * 4,
    swr = TRUE
  )
  
  y_boot <- list()
  
  for (i in 1:n_boot) {
    resid1b <- numeric(length(y_hat))
    resid1b <- resid1[boot_ind1[[i]]]
#    resid1[A == 1] <- resid11[boot_ind11[[i]]]
    
    resid2 <- numeric(length(y_hat))
    for (j in 1:length(group_names)) {
      which_group0 <- which(A == 0 & group == group_names[j])
      resid2[which_group0] <- resid20[boot_ind20[[i]][j]]
      which_group1 <- which(A == 1 & group == group_names[j])
      resid2[which_group1] <- resid21[boot_ind21[[i]][j]]
    }
    
    
    y_boot[[i]]  <- unname(y_hat) + resid1b + resid2
  }
  
  # Save indices and residulas
  list_results <- list(y_boot = y_boot, 
                       
                       resid1 = resid1,
#                       resid11 = resid11,
                       
                       resid21 = resid21,
                       resid20 = resid20,
                       
                       boot_ind1 = boot_ind1,
#                       boot_ind11 = boot_ind11,
                       boot_ind20 = boot_ind20,
                       boot_ind21 = boot_ind21)
  
  
  
  return(list_results)
  
  
}


#'
#' @describeIn resid_boot Two-level block residual bootstrap by treatment group 
#' @export
#'
#'
#Different sampling sample from normal

resid_boot.br6 <- function(obj_boot, ...) {
  
  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  A <- obj_boot$A
  group = obj_boot$group
  group_names <- as.data.frame(table(group))$group
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  sample_sizes <- obj_boot$sample_sizes
  cent_resid <- obj_boot$cent_resid
  
  # Prepare residuals for bootstrapping
  ## Extract total estimation errors 
  resid <- unname(y - y_hat)
  
  ##  Random effects 
  A0_index <- which(A == 0)
  A1_index <- which(A == 1)
  
  ### Control, level 2 
  df0 <- data.frame(resid  = resid[A0_index],
                    group = group[A0_index],
                    y_hat = y_hat[A0_index])
  resid20_prelim <- aggregate(df0$resid, list(df0$group), FUN = mean)$x
  var20 <- var(resid20_prelim)
  ### Treated, level 2 
  df1 <- data.frame(resid  = resid[A1_index],
                    group = group[A1_index],
                    y_hat = y_hat[A1_index])
  resid21_prelim <- aggregate(df1$resid, list(df1$group), FUN = mean)$x
  var21 <- var(resid21_prelim)
  
  resid2_full <- numeric(length(A))
  
  for (i in 1:length(group_names)) {
    which_group0 <- which(A == 0 & group == group_names[i])
    resid2_full[which_group0] <- resid20_prelim[i]
    which_group1 <- which(A == 1 & group == group_names[i])
    resid2_full[which_group1] <- resid21_prelim[i]
  }
  
  ## Total error minus area random effects 
  resid1_full <- resid - resid2_full
  
  ### Unit-level random effect
  resid10_prelim <- resid1_full[A == 0] 
  resid11_prelim <- resid1_full[A == 1] 
  
  if (cent_resid) {
    #Level 1
    mean_r10 <- mean(resid10_prelim)
    mean_r11 <- mean(resid11_prelim)
    
    resid10 <- resid10_prelim - mean_r10
    resid11 <- resid11_prelim - mean_r11
    
    # Level 2
#    mean_r20 <- mean(resid20_prelim)
#    mean_r21 <- mean(resid21_prelim)
    
    
#    resid20 <- resid20_prelim - mean_r20
#    resid21 <- resid21_prelim - mean_r21
  } else {
    resid10 <- resid10_prelim
    resid11 <- resid11_prelim
    
#   resid20 <- resid20_prelim
#    resid21 <- resid21_prelim
  }
  
  
  
  # Sample indices level-1 residuals, control
  boot_ind10 <- sample_bresiduals(
    sample_size = length(resid10),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )
  
  # Sample indices level-1 residuals, treated
  boot_ind11 <- sample_bresiduals(
    sample_size = length(resid11),
    n_boot = n_boot,
    seed = boot_seed * 2,
    swr = TRUE
  )
  
  # # Sample indices level-2 residuals, control
  # boot_ind20 <- sample_bresiduals(
  #   sample_size = length(resid20),
  #   n_boot = n_boot,
  #   seed = boot_seed * 3,
  #   swr = TRUE
  # )
  # 
  # # Sample indices level-2 residuals, treated
  # boot_ind21 <- sample_bresiduals(
  #   sample_size = length(resid21),
  #   n_boot = n_boot,
  #   seed = boot_seed * 4,
  #   swr = TRUE
  # )
  # 
  y_boot <- list()
  boot_ind21 <- list()
  boot_ind20 <- list() 
  
  for (i in 1:n_boot) {
    resid1 <- numeric(length(y_hat))
    resid1[A == 0] <- resid10[boot_ind10[[i]]]
    resid1[A == 1] <- resid11[boot_ind11[[i]]]
    
    resid2 <- numeric(length(y_hat))
    for (j in 1:length(group_names)) {
      resid20 <- rnorm(0, 41, sqrt(var20))
      which_group0 <- which(A == 0 & group == group_names[j])
      resid2[which_group0] <- resid20[[j]]
      resid21 <- rnorm(0, 41, sqrt(var21))
      which_group1 <- which(A == 1 & group == group_names[j])
      resid2[which_group1] <- resid21[[j]]
    }
    boot_ind21[[i]] <- resid21
    boot_ind20[[i]] <- resid20
    
    y_boot[[i]]  <- unname(y_hat) + resid1 + resid2
  }
  
  # Save indices and residulas
  list_results <- list(y_boot = y_boot, 
                       
                       resid10 = resid10,
                       resid11 = resid11,
                       
                       resid21 = resid21,
                       resid20 = resid20,
                       
                       boot_ind10 = boot_ind10,
                       boot_ind11 = boot_ind11,
                      boot_ind20 = boot_ind20,
                       boot_ind21 = boot_ind21)
  
  
  
  return(list_results)
  
}



#'
#' @describeIn resid_boot Two-level block residual bootstrap without re-scaling (sample from level one and level two residuals)
#' @export
#'

resid_boot.br7 <- function(obj_boot, ...) {
  
  
  y <- obj_boot$y
  y_hat <- obj_boot$y_hat
  # data_sample <- obj_boot$data_sample
  A <- obj_boot$A
  group = obj_boot$group
  group_names <- as.data.frame(table(group))$group
  n_boot <- obj_boot$n_boot
  boot_seed <- obj_boot$boot_seed
  sample_sizes <- obj_boot$sample_sizes
  cent_resid <- obj_boot$cent_resid
  c_Win <- obj_boot$c_Win
  
  resid <- unname(y - y_hat)
  
  df <- data.frame(resid  = resid,
                   group = group,
                   y_hat = y_hat)
  
  
  ##  Indicators for treated and control
  A0_index <- which(A == 0)
  A1_index <- which(A == 1)
  
  # Treatment random effects
  resid_t0 <- numeric(length(group_names))
  
  # Mean of resid by group
  resid20 <- numeric(length(group_names))
  resid_gm <- numeric(length(group_names))
  
  for (i in 1:length(group_names)) {
    
    # Group indicator
    which_group <- which(group == group_names[i])
    
    # Treatment and resid by group
    A_g <- A[which_group]
    resid_g <- resid[which_group]
    
    # Means by group
    A_gm <- mean(A_g)
    resid_gm[i] <- mean(resid_g)

    resid_t0[i] <- sum((A_g-A_gm) * (resid_g-resid_gm[i]))/sum((A_g-A_gm)^2)
    
#    gamma_hat[i] <- sum((A[which_group]-mean(A[which_group]))*(resid[which_group]-mean(resid[which_group])))/sum(A[which_group]-mean(A[which_group])^2)
 #   gamma_hat[i] <- sum((samp$w[samp$area==i]-mean(samp$w[samp$area==i]))*(r_hat[samp$area==i]-mean(r_hat[samp$area==i])))/sum(samp$w[samp$area==i]-mean(samp$w[samp$area==i])^2)
    
    if( resid_t0[i] == "NaN"){
      resid_t0[i] = 0
     }
    resid20[i] <-  resid_gm[i] - resid_t0[i] * A_gm
    
  }
  #  e_hat <- r_hat - samp$w*rep(gamma_hat,times=ni) - rep(u_hat,times=ni)
  resid10 <- df$resid - rep(resid20, sample_sizes) - A * rep(resid_t0, sample_sizes)  
 
  
  # resid10 <- Qn(resid10)/sd(resid10) * resid10
  # resid20 <- Qn(resid20)/sd(resid20) * resid20
  # resid_t <- Qn(resid_t0)/sd(resid_t0) * resid_t0
  
  
  if (cent_resid) {
    mean_r1 <- mean(resid10)
    mean_r2 <- mean(resid20)
    mean_t <- mean(resid_t0)
    
    resid1 <- resid10 - mean_r1
    resid2 <- resid20 - mean_r2
    resid_t0 <- resid_t0 - mean_t
  } else {
    resid1 <- resid10
    resid2 <- resid20
    resid_t <- resid_t0
  }
  
  
  
  # Sample indices level-1 residuals
  boot_ind1 <- sample_bresiduals(
    sample_size = length(resid1),
    n_boot = n_boot,
    seed = boot_seed,
    swr = TRUE
  )
  
  # Sample indices level-2 residuals
  boot_ind2 <- sample_bresiduals(
    sample_size = length(resid2),
    n_boot = n_boot,
    seed = boot_seed * 2,
    swr = TRUE
  )
  
  # Sample indices level-2 residuals
  boot_ind_t <- sample_bresiduals(
    sample_size = length(resid_t),
    n_boot = n_boot,
    seed = boot_seed * 2,
    swr = TRUE
  )
  
  
  # Load indices from all bootstrap samples
  #  boot_ind1 <- resid_res$boot_ind1
  #  boot_ind2 <- resid_res$boot_ind2
  
  # Load residuals
  #  resid1 <- resid_res$resid1
  #  resid2 <- resid_res$resid2
  
  # Create y_boot by adding residuals to each y_hat
  
  y_boot <- list()
  for (i in 1:n_boot) {
    y_boot[[i]]  <- unname(y_hat) + resid1[boot_ind1[[i]]] + rep(resid2[boot_ind2[[i]]], sample_sizes) + A * rep(resid_t[boot_ind_t[[i]]], sample_sizes)
  }
  
  # Save indices and residulas
  
  list_results  = list(y_boot = y_boot, 
                       
                       resid1 = resid1,
                       resid2 = resid2,
                       resid_t  = resid_t, 
                       
                       boot_ind1 = boot_ind1,
                       boot_ind2 = boot_ind2, 
                       boot_ind_t = boot_ind_t)
  
  
  
  return(list_results)
  
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

  list_indices <- list()
  for (i in 1:n_boot) {
    # Set seed
    set.seed(seed * i)
    
    #Sample indices
    list_indices[[i]] <- sample(1 : sample_size,
                               sample_size,
                               replace = swr)
  }
  return(list_indices)
}

# Most likeley redundant

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

#################################################################
################################################################
### Another bootstrap type

# #'
# #' @describeIn resid_boot One-level block residual bootstrap by treatment group 
# #' @export
# #'
# 
# resid_boot.br1 <- function(obj_boot, ...) {
#   
#   y <- obj_boot$y
#   y_hat <- obj_boot$y_hat
#   A <- obj_boot$A
#   group = obj_boot$group
#   group_names <- as.data.frame(table(group))$group
#   n_boot <- obj_boot$n_boot
#   boot_seed <- obj_boot$boot_seed
#   sample_sizes <- obj_boot$sample_sizes
#   cent_resid <- obj_boot$cent_resid
#   
#   # Prepare residuals for bootstrapping
#   ## Extract total estimation errors 
#   resid <- unname(y - y_hat)
#   
#   ##  Random effects 
#   A0_index <- which(A == 0)
#   A1_index <- which(A == 1)
#   
#   ### Control, level 2 
#   # df0 <- data.frame(resid  = resid[A0_index],
#   #                   group = group[A0_index],
#   #                   y_hat = y_hat[A0_index])
#   # resid20_prelim <- aggregate(df0$resid, list(df0$group), FUN = mean)$x
#   # 
#   # ### Treated, level 2 
#   # df1 <- data.frame(resid  = resid[A1_index],
#   #                   group = group[A1_index],
#   #                   y_hat = y_hat[A1_index])
#   # resid21_prelim <- aggregate(df1$resid, list(df1$group), FUN = mean)$x
#   # 
#   # resid2_full <- numeric(length(A))
#   # 
#   # for (i in 1:length(group_names)) {
#   #   which_group0 <- which(A == 0 & group == group_names[i])
#   #   resid2_full[which_group0] <- resid20_prelim[i]
#   #   which_group1 <- which(A == 1 & group == group_names[i])
#   #   resid2_full[which_group1] <- resid21_prelim[i]
#   # }
#   
#   ## Total error minus area random effects 
#   #  resid1_full <- resid - resid2_full
#   
#   ### Unit-level random effect
#   resid10_prelim <- resid[A == 0] 
#   resid11_prelim <- resid[A == 1] 
#   
#   llik = function(x, par){
#    m = par[1]
#     s = par[2]
#     n = length(x)
#     # log of the normal likelihood
#     # -n/2 * log(2*pi*s^2) + (-1/(2*s^2)) * sum((x-m)^2)
#    ll = -(n/2)*(log(2*pi*s^2)) + (-1/(2*s^2)) * sum((x-m)^2)
#     # return the negative to maximize rather than minimize
#     return(-ll)
#     
#   }
#  
#   MLr10 = optim(par = c(mean(resid10), sd(resid10)), llik, x = resid10)
#   MLr11 = optim(par = c(mean(resid11), sd(resid11)), llik, x = resid11)
#  
#  
#   if (cent_resid) {
#    #Level 1
#     mean_r10 <- mean(resid10_prelim)
#     mean_r11 <- mean(resid11_prelim)
#     
#     resid10 <- resid10_prelim - mean_r10
#     resid11 <- resid11_prelim - mean_r11
#     
#     # Level 2
#     # mean_r20 <- mean(resid20_prelim)
#     # mean_r21 <- mean(resid21_prelim)
#     # 
#     # 
#     # resid20 <- resid20_prelim - mean_r20
#     # resid21 <- resid21_prelim - mean_r21
#  } else {
#     resid10 <- resid10_prelim
#     resid11 <- resid11_prelim
#     
#     # resid20 <- resid20_prelim
#     # resid21 <- resid21_prelim
#   }
#   
#   
#   
#   # Sample indices level-1 residuals, control
#   boot_ind10 <- sample_bresiduals(
#     sample_size = length(resid10),
#     n_boot = n_boot,
#     seed = boot_seed,
#     swr = TRUE
#   )
#   
#   # Sample indices level-1 residuals, treated
#   boot_ind11 <- sample_bresiduals(
#     sample_size = length(resid11),
#     n_boot = n_boot,
#     seed = boot_seed * 2,
#     swr = TRUE
#   )
#   
#   # Sample indices level-2 residuals, control
#   # boot_ind20 <- sample_bresiduals(
#   #   sample_size = length(resid20),
#   #   n_boot = n_boot,
#   #   seed = boot_seed * 3,
#   #   swr = TRUE
#   # )
#   # 
#   # # Sample indices level-2 residuals, treated
#   # boot_ind21 <- sample_bresiduals(
#   #   sample_size = length(resid21),
#   #   n_boot = n_boot,
#   #   seed = boot_seed * 4,
#   #   swr = TRUE
#   # )
#   
#   # Save indices and residulas
#   output <- list(resid10 = resid10,
#                  resid11 = resid11,
#                  
#                  resid21 = NULL,
#                  resid20 = NULL,
#                  
#                  boot_ind10 = boot_ind10,
#                  boot_ind11 = boot_ind11,
#                  
#                  boot_ind20 = NULL,
#                  boot_ind21 = NULL)
#   
#   
#   
#   return(output)
#   
# }

