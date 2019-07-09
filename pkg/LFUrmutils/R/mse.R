## ## ## ##
## Mean squared error ----
## ## ## ## 

## Generic function ----
mse <- function(object, ...){
  UseMethod("mse")
}

## Default function ----
# mse.default <- function(object) {
#   cat("MSE is currently only defined for objects of the following classes:")
#   methods(mse)
# }


## Univariate volatility models ----
mse.UnivVola <- function(object, ...){
  merg <- merge(object$Variance, object$Returns)
  merg <- na.trim(merg)
  mse <- mean((merg[, 2]^2 - merg[, 1])^2)
  return(mse)
}


## fGARCH model ----
mse.fGARCH <- function(object, ...){
  # Extract dates, variance and returns
  timestamps <- as.Date(attr(object@data, "names"))
  sigSQ <- zoo(x = object@h.t, order.by = timestamps)
  ret <- zoo(x = object@data, order.by = timestamps)
  
  # Merge conditional volatility with squared return series
  merg <- merge(sigSQ, ret)
  
  # Delete trailing missing values
  merg <- na.trim(merg)
  
  # Compute and return mean squared error
  mse <- mean((merg[, 2]^2 - merg[, 1])^2)
  return(mse)
}


# Multivariate EWMA model
mse.MultiEWMA <- function(object, ...){
    n <- dim(object$Variances)[1]
    c <- sqrt(dim(object$Variances)[2])
    
    # Extract diagonal entries of variance-covariance matrix
    volat <- matrix(NA, nrow = n, ncol = c)
    for (k in 1:c) {
      volat[, k] <- object$Variances[, grep( paste0(k,k), colnames(object$Variances) )]
    }
    volat <- zoo(volat, index(object$Variances))
    
    # Find dimension of merged objects
    merg <- merge(volat[, 1], object$Returns[, 1])
    merg <- na.trim(merg)
    m <- dim(merg)[1]
    
    # Initialize object to store squared errors
    errorsSQ <- matrix(NA, nrow = m, ncol = c)
    
    # Compute (mean) squared error
    for(i in 1:c){
      merg <- merge(volat[, i], object$Returns[, i])
      merg <- na.trim(merg)
      errorsSQ[, i] <- (merg[, 1] - merg[, 2]^2)^2
    }
    mse <- mean(errorsSQ)

  return(mse)
}


