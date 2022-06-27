tau_treat_w_ind <- numeric(length(tau_treat_ind))
tau_treat_w_ind[which(weights_treat != 0)] <- 
  tau_treat_ind[which(weights_treat != 0)]/weights_treat[which(weights_treat != 0)]


tau_untreat_w_ind <- numeric(length(tau_untreat_ind))
tau_untreat_w_ind[which(weights_untreat != 0)] <- 
  tau_untreat_ind[which(weights_untreat != 0)]/weights_untreat[which(weights_untreat != 0)]


tau_ind <- tau_treat_ind - tau_untreat_ind

tau_ind_weight <- tau_treat_w_ind - tau_untreat_w_ind

tau_var <- (length(tau_treat_ind) - 1) ^ (-1) *  sum((tau_ind - tau)^2)
tau_var <- (length(tau_treat_ind) - 1) ^ (-1) *  sum((tau_ind_weight - tau)^2)
