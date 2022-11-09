#' Sample subpopulations
#'
#' @param one_pop Population to sample from
#' @param frac_nc Fraction of controls in each subpopulation
#' @param frac_nt Fraction of treated in each subpopulation
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
#' X_outcome <- generate_X(
#'  n = N,
#'  p = 1,
#'  covariance_norm = NULL,
#'  cov_type = "lognorm",
#'  seed = 1
#' )
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
#' subpopulation <- sample_subpopulations(populations,
#' frac_nc = 0.1, frac_nt = 0.1)
#'

sample_subpopulations <- function(one_pop, frac_nc, frac_nt) {

  one_pop_df <- data.frame(one_pop)
  one_pop_df$index <- 1:dim(one_pop_df)[1]
  group <- one_pop_df$group
  group_name <- unique(one_pop_df$group)


  list_group_name <- list()
  for (i in 1:length(group_name)) list_group_name[[i]] <- group_name[i]

  filter_by_group <- function(group_name) {
    data_group <- filter(one_pop_df, group == group_name )
    data_group
  }
  list_group <- lapply(list_group_name, filter_by_group)

  Nc = as.numeric(table(one_pop$group[one_pop$A == 0]))
  Nt = as.numeric(table(one_pop$group[one_pop$A == 1]))

  nic <- round(frac_nc * Nc)
  nic[which(nic == 0)] <- 1

  nit <- round(frac_nt * Nt)
  nit[which(nit == 0)] <- 1

  sample_group <- list()

  for (i in 1:length(list_group)) {
    group_data <- list_group[[i]]

    index_control <- group_data$index[group_data$A == 0]

    if (length(index_control) == 1) {
      sample_control <- index_control
    } else {
      sample_control <- sample(index_control, nic[i])
    }

    index_treat <- group_data$index[group_data$A == 1]

    if (length(index_treat) == 1) {
      sample_treat <- index_treat
    } else {
      sample_treat <- sample(index_treat, nit[i])
    }

    sample_group[[i]] <- sort(c(sample_control, sample_treat))

  }
  sample_subpopulations <- unlist(sample_group, use.names = FALSE)
  sample_subpopulations
}
