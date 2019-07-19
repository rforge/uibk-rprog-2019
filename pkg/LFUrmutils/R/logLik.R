## ## ## ##
## logLik methods ----
## ## ## ##

logLik.fGARCH <- function(object, ...){
  structure(-object@fit$llh, df = length(object@fit$par), nobs = length(residuals(object)), class = "logLik")
}

logLik.UnivVola <- function(object, ...){
  structure(-object$llh, df = 1, nobs = length(object$Returns), class = "logLik")
}

logLik.MultiEWMA <- function(object, ...){
  structure(-object$llh, df = 1, nobs = dim(object$Returns)[1], class = "logLik")
}


