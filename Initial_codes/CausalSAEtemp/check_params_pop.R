#' Default parameters for different regression types
#' 
#' @return List of parameters, that is
#' \item{frac_out}{Fraction of outliers, default: frac_out = 0}
#' \item{re}{List of parameters to generate random effects (re) with outliers, 
#' default: re = list(var_re = 3, frac_out_re = 0.25, var_re_out = 10, mean_re_out = 15)}
#' \item{params_cont}{List of parameters to generate errors for continuous outcome, default:
#' params_cont = list(var_e = 6, var_e_out = 150, mean_e_out = 20)}
#' \item{disturbance}{Disturbance parameter to generate non continuous outcome}
#' 
#' @export
#' 
#' @examples 
#' # To obtain a list with defualt parameters, run
#' 
#' list1 <- get_default_params()
#'
get_default_params <- function() {
  default_params <- list(frac_out = 0.05,
                         re = list(var_re = 3,
                                   frac_out_re = 0.25,
                                   var_re_out = 10, 
                                   mean_re_out = 15),
                         var_re_treat  = 0.25,
                         params_cont = list(var_e = 6,
                                            var_e_out = 150,
                                            mean_e_out = 20),
                         disturbance = 5)
  return(default_params)
}

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
#' check_continuous <- check_params_pop(params = get_default_params(), 
#'                                      regression_type = "continuous")
#'  
#' @export
#' 
#' 

check_params_pop <- function(params = get_default_params()) {
  
  stopifnot("params must be a list" = identical(class(params), "list"))
  
  default_params <- get_default_params()
  
  for (par_name in names(default_params)) {
    if (is.null(params[[par_name]])) {
      params[[par_name]] <- default_params[[par_name]]
    } else {
      if (mode(params[[par_name]]) != mode(default_params[[par_name]])) {
        stop(sprintf("%s is not of the expected type %s.", par_name, mode(default_params[name])))
      }
    }
  }
  
  #  default_params <- list(cov_norm = NULL, 
  #                         frac_out = 0.05,
  #                         re = list(var_re = 3,
  #                                   frac_out_re = -0.25,
  #                                   var_re_out = 10, 
  #                                   mean_re_out = 15),
  #                         params_cont = list(var_e = 6,
  #                                            var_e_out = 150,
  #                                            mean_e_out = 20),
  #                         disturbance = 5)
  
  params  = default_params
  
  ## Check if fraction of outliers is valid
  for (par_name in c(params$frac_out, params$re$frac_out_re))  {
    if (any(par_name < 0) | any (par_name >= 1)) {
      stop("Choosen faction of outliers is not valid. Choose a fraction  0 <= fraction < 1.")
    }
    
  }
  
  #  if (any(frac_out < 0) | any (frac_out >= 1)) {
  #    stop("Choosen faction of outliers is not valid. One should select 0 <= frac_out < 1.")
  #  }
  #  if (any(re$frac_out_re < 0) | any (re$frac_out_re >= 1)) {
  #    stop("Choosen faction of outliers is not valid. One should select 0 <= frac_out_re < 1.")
  #  }
  
  ## Check if variance parameter is valid
  for (par_name in c(params$params_cont$var_e,
                     params$params_cont$var_e_out,
                     params$re$var_re, params$re$var_re_out, 
                     params$var_re_treat))  {
    if (any(par_name < 0)) {
      stop("Choosen faction of outliers is not valid. Choose a fraction  0 <= fraction < 1.")
    }
  }
  
  
  return(params)
}
