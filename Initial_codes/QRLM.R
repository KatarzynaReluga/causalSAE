## General comment: this can be a small package per se
## we have a main function and at least 4-5 functions
## I think we should remove dots . and replace them with 
## _ (dots are used to represent methods of generic functions, 
## it might be confusing to leave to use them)
##
#' What does function do?
#'
#' More details about the function
#'
#'
#' @param x Design matrix including the intercept.
#' Matrix frame of covariates?
#' 
#' @param y The response variable.
#' @param case.weights An optional vector of 
#' weights to be used in the fitting process.
#' Default: \code{case.weights = rep(1, nrow(x))}
#' Possible sample wights, you can use any wights 
#' @param var.weights An optional vector of weights of 
#' the variance to be used in the fitting process.
#' Default: \code{ var.weights = rep(1, nrow(x))}
#' Weights to compute the variance you can have different wights
#' for var and for fitting
#' @param init ? Default: \code{init = "ls"}
#' Linear regression for the starting point, starting beta for comupitng the residuslas
#' @param psi ? Default: \code{psi = psi.huber}
#' MASS  --> I need to take it from there
#' 
#' @param scale.est ? Avaialble choices:
#' different choice of the estimation of scale, MAD or different robust methods
#' @param k2 
#' In the psi function -- you need k 1.345
#' @param k = 1.345, psi.huber check in MASS
#'  Default: \code{k2 = 1.345}
#'  is a parameter in the function to compute the residuals
#'  they can be the same or diff
#' @param method ? Available choices: M estimation or MM estimation 
#' @param maxit An integer giving the number of maximum iterations
#' Default: \code{maxit = 20}
#' @param acc ? Positive convergence 
#' tolerance  Default: \code{acc = 1e-04}
#' @param test.vec Test vector in residuals 
#' Method to test the convergence from one iter to the second 
#' Default: \code{test.vec = "resid"}
#' @param q A vector of M-quantiles. Default: \code{q = 0.5}
#'
#' @return List with elements:
#' \item{fitted.values}{}
#' \item{residuals}{}
#' \item{q.values}{}
#' \item{q.weights}{}
#' \item{coef}{}
#' \item{qscale}{}
#'
#' @details
#' Function QRLM implements (what method) published
#' RLM - robust linear model
#' Q - quantile
#' in (what journal) Chambers and Tzavidis (2006) Small area estiamtion ..., Biometrika
#' 
#' @importFrom MASS lqs psi.huber
#'

QRLM <- function (x,
                  y,
                  case.weights = rep(1, nrow(x)),
                  var.weights = rep(1, nrow(x)),
                  ...,
                  w = rep(1, nrow(x)),
                  init = "ls",
                  psi = psi.huber,
                  scale.est = c("MAD", "Huber", "proposal 2"),
                  k2 = 1.345,
                  method = c("M", "MM"),
                  maxit = 20,
                  acc = 1e-04,
                  test.vec = "resid",
                  q = 0.5)
{
  ###############################################
  ## This is a different function, it should be
  ## outside this file, helper utitlites
  ###################################
  #' What does the function do?
  #'
  #' @param old ? old coef 
  #' @param new ? new coef
  irls.delta <-
    function(old, new)
      sqrt(sum((old - new) ^ 2) / max(1e-20, sum(old ^ 2)))
  # consequtive values of the coeff.
  ##############################################
  ## A different function, should be outside
  ## this file
  ################################################
  #'What does it do?
  #'
  #'@param 
  #'
  #'
  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix(
      r * w, 1, length(r)
    ) %*% x) / sqrt(matrix(w, 1, length(
      r
    )) %*% (x ^ 2)))) / sqrt(sum(w * r ^ 2))
  }
  # if you use residual, the other function, if you use other methods
  # the fucntion to evaluate the convergence (apart from residuals)
  method <- match.arg(method)
  #########################################
  #What does nmx stand for? Names of x?
  ########################################
  
  nmx <- deparse(substitute(x))
  
  ##########################################
  # Checks of data? 
  # This should be a separate function
  # check_x <- function(x) {
  #
  #}
  # or 
  # check_paramters <- function(x, test.vec, ...)
  ########################################
  if (is.null(dim(x))) {
    x <- as.matrix(x)
    colnames(x) <- nmx
  } else {
    x <- as.matrix(x)
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste("X", seq(ncol(x)), sep = "")
  }
  
  if (qr(x)$rank < ncol(x)){
    stop("x is singular: singular fits are not implemented in rlm")
  }
  
  #######################################
  # Check test.vec
  # Why cannot we write:
  # test.vec = c("resid", "coef", "w", "NULL")
  # test.vec <- match.arg(test.vec)
  #########################################
  if (!(any(test.vec == c("resid", "coef", "w", "NULL")) ||
        is.null(test.vec))) {
    stop("invalid testvec")
  }
  
  ############################################
  ## Checks of var.weights 
  ## we can use   stopifnot
  #stopifnot("Length of var.weights must equal number of observations" = 
  #            length(var.weights) == nrow(x))
  if (length(var.weights) != nrow(x))
    stop("Length of var.weights must equal number of observations")
  if (any(var.weights < 0))
    stop("Negative var.weights value")
  
  #########################################
  ## The same as above
  if (length(case.weights) != nrow(x))
    stop("Length of case.weights must equal number of observations")
  
  #### Part of the checks
  w <- (w * case.weights) / var.weights
  
  ### Checks once again
  if (method == "M") {
    scale.est <- match.arg(scale.est)
    ### Where is psi.huber defined? 
    if (!is.function(psi))
      psi <- get(psi, mode = "function")
    arguments <- list(...)
    if (length(arguments)) {
      pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
      if (any(pm == 0))
        warning(paste("some of ... do not match"))
      pm <- names(arguments)[pm > 0]
      formals(psi)[pm] <- unlist(arguments[pm])
    }
    if (is.character(init)) {
      if (init == "ls")
        #Linear weighted model, not sure what
        # it does
        # why method here not a generic?
        temp <- lm.wfit(x, y, w, method = "qr")
      else if (init == "lts")
        # lqs --> resistant regression ?
        # is it robust regression?
        # why method here not a generic?
        temp <- lqs.default(x, y, intercept = FALSE, nsamp = 200)
      else
        stop("init method is unknown")
      coef <- temp$coef
      resid <- temp$resid
    }
    else {
      if (is.list(init))
        coef <- init$coef
      else
        coef <- init
      resid <- y - x %*% coef
    }
  }
  else if (method == "MM") {
    scale.est <- "MM"
    # Once again method not generic?
    temp <-
      lqs.default(x,
                  y,
                  intercept = FALSE,
                  method = "S",
                  k0 = 1.548)
    #k0 fixed 
    #MM psi b-square
    coef <- temp$coef
    resid <- temp$resid
    psi <- psi.bisquare
    if (length(arguments <- list(...)))
      if (match("c", names(arguments), nomatch = FALSE)) {
        c0 <- arguments$c
        if (c0 > 1.548) {
          psi$c <- c0
        }
        else
          warning("c must be at least 1.548 and has been ignored")
      }
    scale <- temp$scale
  }
  else
    stop("method is unknown")
  
  # to stop the alg.
  done <- FALSE
  conv <- NULL #convergence?
  n1 <- nrow(x) - ncol(x)
  if (scale.est != "MM")
    scale <- mad(resid / sqrt(var.weights), 0)
  #############################################
  ##########################################
  # why is it commented
  #TO be removed
  #theta <- 2 * pnorm(k2) - 1
  #gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
  # expect to compute, in MAD median of the residuals /0.sth
  # this sth obtained by using gamma
  ######################################
  ### This function returns
  ###############################
  #fitted.values = qfit,
  #residuals = qres,
  #q.values = q,
  #q.weights = qwt,
  #coef = qest,
  #qscale = qscale
  #estimated?
  qest <- matrix(0, nrow = ncol(x), ncol = length(q))
  #q.weights
  qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
  # fitted.values
  qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
  # residuals Why cannot we call them residuals?
  qres <- matrix(0, nrow = nrow(x), ncol = length(q))
  ####################
  # No idea?
  # This has to be explained to me
  ##########################
  qscale <- NULL
  for (i in 1:length(q)) {
    for (iiter in 1:maxit) {
      if (!is.null(test.vec))
        testpv <- get(test.vec)
      if (scale.est != "MM") {
        if (scale.est == "MAD")
          scale <- median(abs(resid / sqrt(var.weights))) / 0.6745
        else {
          gamma <-
            4 * k2 ^ 2 * (1 - pnorm(k2)) * ((1 - q[i]) ^ 2 + q[i] ^ 2) - 4 * k2 * dnorm(k2) *
            ((1 - q[i]) ^ 2 + q[i] ^ 2) + 4 * (1 - q[i]) ^ 2 * (pnorm(0) - (1 - pnorm(k2))) + 4 *
            q[i] ^ 2 * (pnorm(k2) - pnorm(0))
          scale <-
            sqrt(sum(pmin(
              resid ^ 2 / var.weights, (k2 * scale) ^ 2
            )) / (n1 * gamma))
        }
        if (scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid / (scale * sqrt(var.weights))) * case.weights
      ww <- 2 * (1 - q[i]) * w
      ww[resid > 0] <- 2 * q[i] * w[resid > 0]
      w <- ww
      temp <- lm.wfit(x, y, w, method = "qr")
      coef <- temp$coef
      resid <- temp$residuals
      if (!is.null(test.vec))
        convi <- irls.delta(testpv, get(test.vec))
      else
        convi <- irls.rrxwr(x, wmod, resid)
      conv <- c(conv, convi)
      done <- (convi <= acc)
      if (done)
        break
    }
    if (!done)
      warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
    qest[, i] <- coef
    qscale[i] <- scale
    qwt[, i] <- w
    qfit[, i] <- temp$fitted.values
    qres[, i] <- resid
  }
  list(
    fitted.values = qfit,
    residuals = qres,
    q.values = q,
    q.weights = qwt,
    coef = qest,
    qscale = qscale
  )
}

#########################################
### This should be outside the main function
### I don't know what this function do?
#### It computes order of quantiles?
######################################
# COMPUTING OF THE QUANTILE-ORDERS
"zerovalinter" <- function(y, x)
{
  if (min(y) > 0) {
    xmin <- x[y == min(y)]
    if (length(xmin) > 0)
      xmin <- xmin[length(xmin)]
    xzero <- xmin
  } 
  
  
  
  else {
    if (max(y) < 0) {
      xmin <- x[y == max(y)]
      if (length(xmin) > 0)
        xmin <- xmin[1]
      xzero <- xmin
    }
    else {
      y1 <- min(y[y > 0])
      if (length(y1) > 0)
        y1 <- y1[length(y1)]
      y2 <- max(y[y < 0])
      if (length(y2) > 0)
        y2 <- y2[1]
      x1 <- x[y == y1]
      if (length(x1) > 0)
        x1 <- x1[length(x1)]
      x2 <- x[y == y2]
      if (length(x2) > 0)
        x2 <- x2[1]
      xzero <- (x2 * y1 - x1 * y2) / (y1 - y2)
      xmin <- x1
      if (abs(y2) < y1)
        xmin <- x2
    }
  }
  resu <-  xzero
  resu
}

### This should be outside 

# Function for Finding the Quantile Orders by Linear Interpolation
# Assumes that "zerovalinter" function has been already loaded
# For each y  we want to get the quantile coef
# we have to find M-quant order (Chambers, Salvati paper)
#' @param y observed values
#' @param fitted_y noexpec fitted values fro all the quatiles 
#' @param Q quant estiamted in the M-quantile function (28-29 quantlies)
"gridfitinter" <- function(y, expectile, Q)
  # computing of the expectile-order of each observation of y by interpolation
{
  nq <- length(Q)
  diff <- y %*% t(as.matrix(rep(1, nq))) - expectile
  vectordest <- apply(diff, 1, zerovalinter, Q)
  #print(vectordest)
  #qord<-list(ord=c(vectordest))
  #qord
}
