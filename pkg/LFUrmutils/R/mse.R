## ## ## ##
## Mean squared error ----
## ## ## ## 

# Generic function
MSE <- function(object, ...){
  UseMethod("MSE")
}

# Default function
MSE.default <- function(object) {
  cat("MSE is only defined for objects of the following classes:")
  methods(MSE)
}

# Univariate moving average model
MSE.UniMA <- function(cvar, y = NULL, squared = FALSE, volatility = FALSE){
  # cvar: sigma_t^2 estimates
  # y: (daily) log return series
  if (squared == FALSE) {
    y <- y^2
  }
  
  if (volatility == TRUE){
    cvar <- cvar^2
  }
  
  merg <- merge(cvar, y)
  merg <- na.trim(merg)
  
  mse <- mean((merg[, 2] - merg[, 1])^2)
  return(mse)
}

# Univariate EWMA model
MSE.UniEWMA <- function(object, y = NULL, squared = FALSE){
  x <- 10
}

# Multivariate EWMA model
# Univariate EWMA model
MSE.MultEWMA <- function(object, y = NULL, squared = FALSE){
  x <- 10
}

# fGARCH
MSE.fGARCH <- function(model, ret, squared = FALSE){
  # model: object of class fGARCH
  # ret: return series
  # squared: set to TRUE, if returns are already squared
  
  # Compute squared returns, if necessary
  if (squared == FALSE) {
    ret <- ret^2
  }
  
  # fitGarch computes conditional volatility by default
  # Therefore, we have to compute the conditional variance manually
  sigma_t2 <- zoo(model@sigma.t^2, time(ret))
  
  # Merge conditional volatility with squared return series
  merg <- merge(sigma_t2, ret)
  
  # Delete trailing missing values
  merg <- na.trim(merg)
  
  # Compute and return mean squared error
  mse <- mean((merg[, 2] - merg[, 1])^2)
  return(mse)
  
}

