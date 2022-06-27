#' Sample from the population
#' 
#' @param populations List with populations
#' @param ni_size Sample size 
#' @param sample_part Select sampled, non_sampled or both parts 
#' @param get_index Get the index of sampled observations
#' or population, default: sample_part = sampled  
#'
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
#'  start_seed = 1
#' )
#'
#'
#' X_outcome <- generate_X(
#'  n = N,
#'  p = 1,
#'  covariance_norm = NULL,
#'  cov_type = "lognorm",
#'  start_seed = 1
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
#' no_sim = 10,
#' start_seed = 1)
#' 
#' samples <- generate_sample(populations,  ni_size = 5,
#' sample_part = "sampled")
#' 


generate_sample <- function(populations, ni_size = 5,
                            sample_part = "sampled", 
                            get_index = F){
  
  if (is.data.frame(populations)) {
    populations <- list(populations)
  }
  
  if (ni_size > length(populations[[1]]$group)) {
    stop("ni_size cannot be bigger than the size of the population.")
  }
  
  m = length(unique(populations[[1]]$group))
  ni = rep(ni_size, m)
  
  samples <- lapply(populations, generate_one_sample, 
                    size_samples = ni, sample_part, get_index)
  return(samples)
}

#'
#' Generate one sample from population 
#'
#' @param population Data frame with a generated population
#' @param size_samples Sample sizes in each subpopulation/area
#' @param sample_part Select sampled, non_sampled or both parts 
#' or population, default: sample_part = sampled  
#' @param get_index Get the index of sampled observations
#' 
#' @importFrom sampling strata
#'

generate_one_sample <- function(population, 
                                size_samples, 
                                sample_part = "sampled", 
                                get_index = F) {
  
  s <- strata(population, stratanames = "group", 
              size = size_samples, 
              method = "srswor") 

  samp_data <- population[s[, 2] ,]
  non_samp_data <- population[ - s[, 2], ]
  
  if (sample_part == "sampled") {
    if (get_index) {
      output <- list(samp_data = samp_data, 
                     index_s = s[, 2])
    } else {
      output <- samp_data 
    }

  } else {
    if (get_index) {
      output <- list(samp_data = samp_data,
                     non_samp_data = non_samp_data,
                     index_s = s[, 2])
    } else {
      output <- list(samp_data = samp_data, 
                     non_samp_data = non_samp_data)  
    }
  }
  output
}


