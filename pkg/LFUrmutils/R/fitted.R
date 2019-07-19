## ## ## ##
## Fitted values ----
## ## ## ## 

## Univariate volatility models ----
fitted.UnivVola <- function(object, ...){
  fit <- object$Variance
  return(fit)
}

# Multivariate EWMA model
fitted.MultiEWMA <- function(object, ...){
  fit <- object$Variances
  return(fit)
}


