## ## ## ##
## Univariate volatility models ----
## ## ## ##

# UnivMA <- function(returns, width = 30, center = FALSE){
#   # Center returns for comparability reasons
#   if (center == TRUE){
#     returns <- scale(returns, center = TRUE, scale = FALSE)
#   }
#   
#   # Compute variance
#   sigma_t2 <- rollmeanr(returns^2, k = width, fill = NA)
#   
#   # Add empty row
#   NBD <- NextBusinessDay(returns)
#   newrow <- zoo(NA, NBD)
#   sigma_t2 <- suppressWarnings(rbind(sigma_t2, newrow))
#   colnames(sigma_t2) <- "Variance"
#   
#   # Lead series by one observation
#   sigma_t2 <- lag(sigma_t2, -1)
#   
#   # Create and return output object
#   object <- list(Variance = sigma_t2, Returns = returns, width = width, center = center)
#   attr(object, "class") <- c("UnivMA", "UnivVola")
#   return(object)
# }

UnivVola <- function(returns, width = 30, lambda = 0.94, type = c("RiskMetrics", "WeightedAverage", "MovingAverage"), center = FALSE){
  # Check validity of model variant
  "%!in%" <- Negate("%in%")
  if (type %!in% c("RiskMetrics", "WeightedAverage", "MovingAverage")){
    stop("type must be either \"RiskMetrics\", \"WeightedAverage\" or \"MovingAverage\".")
  }
  
  # Center returns for comparability reasons
  if (center == TRUE){
    returns <- scale(returns, center = TRUE, scale = FALSE)
  }
  
  if (type %in% c("RiskMetrics", "WeightedAverage")){
    if(lambda <= 0 | lambda >= 1){
      stop("lambda must be between 0 and 1.")
    }
    # Compute 1-lambda for efficiency reasons
    nlambda <- 1 - lambda
  }
  
  if (type == "RiskMetrics") {
    # Initialize output object
    obs <- length(returns) + 1
    sigma_t2 <- rep(NA, obs)
    sigma_t2[1] <- var(returns)
    
    # Compute conditional variances
    for(i in 2:obs){
      sigma_t2[i] <- nlambda * returns[i-1]^2 + lambda * sigma_t2[i-1]
    }
    
    # Add empty row
    NBD <- NextBusinessDay(returns)
    newrow <- zoo(NA, NBD)
    newindex <- merge(returns, newrow)
    sigma_t2 <- zoo(sigma_t2, index(newindex))
    
    # Set width to NA, whatever the original value is
    width <- NA
  }
    
  if (type == "WeightedAverage") {
    # Compute correction factor and weights
    CF <- nlambda/(lambda*(1-lambda^width))
    s <- c(width:1)
    lambda <- lambda^s
    
    # Compute conditional variances
    sigma_t2 <- rollapplyr(returns^2, width = width, FUN = function(x) CF * sum(lambda * x), fill = NA)
    
    # Add empty row
    NBD <- NextBusinessDay(returns)
    newrow <- zoo(NA, NBD)
    sigma_t2 <- suppressWarnings(rbind(sigma_t2, newrow))
    
    # Lead series by one observation
    sigma_t2 <- lag(sigma_t2, -1, na.pad = TRUE)
  }
  
  if (type =="MovingAverage"){
    # Compute variance
    sigma_t2 <- rollmeanr(returns^2, k = width, fill = NA)
    
    # Add empty row
    NBD <- NextBusinessDay(returns)
    newrow <- zoo(NA, NBD)
    sigma_t2 <- suppressWarnings(rbind(sigma_t2, newrow))
    colnames(sigma_t2) <- "Variance"
    
    # Lead series by one observation
    sigma_t2 <- lag(sigma_t2, -1, na.pad = TRUE)
    
    # Set lambda to NA, whatever the original input is
    lambda <- NA
  }
  
  # Create and return output object
  object <- list(Variance = sigma_t2, Returns = returns, variant = match.arg(type), width = width, lambda = lambda, centered = center)
  attr(object, "class") <- c("UnivVola")
  
  return(object)
}


## ## ## ##
## Multivariate volatility models ----
## ## ## ##

MultiEWMA <- function(returns, lambda = 0.94, center = FALSE){
  # Extract necessary information from the return series
  x <- coredata(returns)
  n <- dim(x)[1] + 1
  d <- dim(x)[2]
  
  # Center returns for comparability reasons
  if (center == TRUE){
    x <- scale(x, center = TRUE, scale = FALSE)
  }
  
  # Create names
  names <- matrix(NA, d, d)
  for (i in 1:d){
    for (j in 1:d) {
      names[j,i] <- paste0("Sigma", j,i)
    }
  }
  names <- c(t(names))
  
  # Initialize output
  SIGMA_t <- matrix(nrow = n, ncol = d^2)
  colnames(SIGMA_t) <- names
  sigma <- cov(x)
  SIGMA_t[1, ] <- as.vector(sigma)
  nlambda <- 1-lambda
  
  # Core of the function: Compute EWMA
  for(i in 2:n){
    sigma <- lambda * sigma + nlambda * tcrossprod(x[(i-1), ])
    SIGMA_t[i, ] <- as.vector(sigma)
  }
  SIGMA_t <- zoo(SIGMA_t, index(returns))
  
  # Create new index
  NBD <- NextBusinessDay(returns)
  newrow <- zoo(NA, NBD)
  newindex <- rbind(returns[, 1], newrow)
  SIGMA_t <- zoo(SIGMA_t, index(newindex))
  
  object <- list(Variances = SIGMA_t, Returns = returns, lambda = lambda, centered = center)

  attr(object, "class") <- "MultiEWMA"
  return(object)
}
