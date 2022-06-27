#' Check or initialise the default parameters to generate populations
#'
#' This function checks or initialises parameter used to generate populations
#' to carry out simulation studies
#'
#' @param params List of parameters
#' @param regression_type Type of regression model we aim to generate:
#' \itemize{
#' \item \code{"continuous"} Continuous regression
#' \item \code{"binary"} Binary regression
#' \item \code{"poisson"} Poisson regression
#' \item \code{"nb"} Negative binomial regression
#' }  
#' 
#' @examples
#'
#' # Check the default parameters for different types of 
#' # regression problems. If some parameters are not set, they are set to defult
#' # values.
#' 
#' check_continuous <- check_params_pop(params = list(no_sim = 2,
#'                                                    ni_size = 5,
#'                                                    Ni_size = 100,
#'                                                    m = 50, 
#'                                                    frac_re_out_y = 0.2,
#'                                                    no_cov_unif = 1,
#'                                                    coef_cov_unif = 1), 
#'                                      regression_type = "continuous")
#'  
#' @export
#' 
#' 

check_params_pop <- function(params = NULL,
                             regression_type = c("continuous",
                                                 "binary",
                                                 "poisson",
                                                 "nb")) {
  stopifnot("params must be a list" = identical(class(params), "list"))
  
  regression_type <- match.arg(regression_type)
  
  default_params <- get_default_params(regression_type)
  
  for (name in names(default_params)) {
    if (is.null(params[[name]])) {
      params[[name]] <- default_params[[name]]
    } else {
      if (mode(params[[name]]) != mode(default_params[[name]])) {
        stop(sprintf(
          "%s is not of the expected type %s.",
          name,
          mode(default_params[name])
        ))
      }
    }
    if (name %in% c(
#      "no_sim",
      "ni_size",
      "Ni_size",
      "m",
      "var_re_y",
      "var_re_out_y",
      "var_treat_effect",
      "subset_cov_unif"
    ) && all(params[[name]] <= 0)) {
      stop(paste(name, "should be positive."))
    }
    if (name %in% c("no_sim", "ni_size", "Ni_size", "m") &&
        identical(params[[name]], as.integer(params[[name]]))) {
      stop(paste(name, "should be integer"))
    }
    if (name %in% c("frac_re_out_y",
                    "no_cov_unif",
                    "no_cov_lnorm",
                    "no_cov_unif_treat") && all(params[[name]] < 0)) {
      stop(paste(name, "should not be negative"))
    }
    if (name %in% c("frac_re_out_y") && all(params[[name]] >= 1)) {
      stop(paste(name, "should be a positve fraction"))
    }
    
    if (regression_type == "continuous") {
      if (name %in% c("var_e_y", "var_e_out_y") && all(params[[name]] <= 0)) {
        stop(paste(name, "should be positive."))
      }
      if (name == "frac_e_out_y" & (all(params[[name]] < 0)|all(params[[name]] > 1))) {
        stop(paste(name, "should be a positive fraction."))
      }
    } else {
      if (name == "frac_contam_y" & (all(params[[name]] < 0)|all(params[[name]] > 1))) {
        stop(paste(name, "should be a positive fraction."))
      }
    }
    
  }
#  if (params$ni_size > params$Ni_size) {
#    stop("ni_size should be bigger than Ni_size.")
#  }
  
  if (max(params$subset_cov_unif) > params$no_cov_unif) {
    stop("Confounders of treatment should be a subset of
         covariates.")
  }
  
  stopifnot(
    "Covariates and offset cannot be 0! No regression model for y" = (
      identical(params$coef_cov_unif_y,
                params$coef_cov_lnorm_y, 0) &&
        !identical(params$offset_y, 0)
    )
  )
  return(params)
}


#' Default parameters for different regression types
#' @inheritParams check_params_pop
#' 
#' @return A list the default parameters
#' 
#' @export
#' 
#' @examples 
#' # To obtain a list with defualt parameters, run
#' 
#' list1 <- get_default_params(regression_type = c("continuous"))
#' list2 <- get_default_params(regression_type = c("binary"))
#' list3 <- get_default_params(regression_type = c("poisson"))
#' list4 <- get_default_params(regression_type = c("nb"))  
#' 
#' 
get_default_params <- function(regression_type = c("continuous",
                                                   "binary",
                                                   "poisson",
                                                   "nb")) {
  regression_type <- match.arg(regression_type)
  
  
  basic_params <- list(
    no_sim = 1,
#    ni_size = 5,
    Ni_size = 100,
    m = 50,
    
    # Generate covariates for Y
    frac_re_out_y = 0.2,
    no_cov_unif = 1,
    coef_cov_unif = 1,
    
    no_cov_lnorm = 1,
    coef_cov_lnorm = 2,
    
    # Generate treatment T (1/0)
    var_re_treat = 0.25,
    offset_treat = -1,
    subset_cov_unif = c(1),
    
    no_cov_unif_treat = 0,
    coef_cov_unif_treat = 0,
    
    # Generate heterogeneous treatment effect
    mean_treat_effect = 10,
    var_treat_effect = 1
  )
  
  
  if (regression_type == "continuous") {
    params = basic_params
    
    #Generate y
    params$var_re_y = 3
    params$var_re_out_y = 20
    params$mean_re_out_y = 9
    params$offset_y = 100
    
    # Generate errors
    params$var_e_y = 6
    params$frac_e_out_y = 0.03
    params$var_e_out_y = 150
    params$mean_e_out_y = 20
    
  } else {
    #Generate y
    params = basic_params
    
    params$var_re_y = 0.25
    params$var_re_out_y = 1.5
    params$mean_re_out_y = 0.8
    params$offset_y = 0.5
    
    params$frac_contam_y = 0.03
    params$shift_contam_y = 5
    
  }
  
  return(params)
}
