#' Sample from the population
#' 
#' @param populations List with populations
#' @param ni_size Sample size 
#' @param sample_part Select sampled, non_sampled or both parts 
#' or population, default: sample_part = sampled  
#'
#' @examples 
#' 
#' gen_pops <- generate_pop(coef_outcome_p_score = get_default_cov("outcome_p_score"),
#'                          coef_outcome = get_default_cov("outcome"),
#'                          coef_p_score = get_default_cov("p_score"),
#'                          hte = list(mean = 10, var = 1),
#'                          Ni_size  = 1000,
#'                          m = 50,
#'                          no_sim = 10,
#'                          regression_type = "continuous", 
#'                          additional_params = list(frac_out = 0.03,
#'                          re = list(var_re = 3,
#'                                    frac_out_re = 0.2,
#'                                    var_re_out = 9, 
#'                                    mean_re_out = 20),
#'                                    var_re_treat  = 0.25, 
#'                                    params_cont = list(var_e = 6,
#'                                                       var_e_out = 150,
#'                                                       mean_e_out = 20), 
#'                                    disturbance = 5), 
#'                                    start_seed = 1)
#' samples <- generate_sample(gen_pops)
#' 
#' @export
#' 


generate_sample <- function(populations, ni_size = 5,
                            sample_part = "sampled"){
  
  if (ni_size > length(populations[[1]]$group)) {
    stop("ni_size cannot be bigger than the size of the population.")
  }
  
  m = length(unique(populations[[1]]$group))
  ni = rep(ni_size, m)
  
  samples <- lapply(populations, generate_one_sample, 
                    size_samples = ni, sample_part)
  return(samples)
}

#'
#' Generate one sample from population 
#'
#' @param population Data frame with a generated population
#' @param size_samples Sample sizes in each subpopulation/area
#' @param sample_part Select sampled, non_sampled or both parts 
#' or population, default: sample_part = sampled  
#' 
#' @importFrom sampling strata
#'

generate_one_sample <- function(population, 
                                size_samples, 
                                sample_part = "sampled") {
  
  s <- strata(population, stratanames = "group", 
              size = size_samples, 
              method = "srswor") 

  samp_data <- population[s[, 2] ,]
  non_samp_data <- population[ - s[, 2], ]
  
  if (sample_part == "sampled") {
    output <- samp_data
  } else if (sample_part == "non_sampled") {
    output <- samp_data
  } else {
    output <- list(samp_data = samp_data, 
                   non_samp_data = non_samp_data)  
  }
  output
}


