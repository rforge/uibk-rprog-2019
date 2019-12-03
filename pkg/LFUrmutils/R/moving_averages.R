## ## ## ##
## Utility functions to estimate lambda in EWMA models ----
## ## ## ##

## COPYRIGHT INFORMATION
## The original code for the two following functions is due to David Ruppert and David S. Matteson.

nllik.ewma <- function(lambda0, innov){
  clambda0 <- 1-lambda0
  Sigma.hat <- var(innov)
  invSigma.hat <- chol2inv(chol(Sigma.hat))
  detSigma.hat <- det(Sigma.hat)
  
  llik <- -0.5 * log(detSigma.hat) - 0.5 * crossprod(innov[1,],invSigma.hat) %*% innov[1,]
  llik <- llik - 0.5 * log(detSigma.hat) - 0.5 * crossprod(innov[2,],invSigma.hat) %*% innov[2,]
  
  n <- dim(innov)[1]
  
  for(i in 3:n){
    atm1 <- innov[(i-1),]
    at <- innov[i,]
    denom <- 1 - lambda0^(i-1)
    Sigma.hat <- (clambda0/denom) * tcrossprod(atm1) + (lambda0 * (1-lambda0^(i-2)) / denom) * Sigma.hat
    invSigma.hat <- chol2inv(chol(Sigma.hat))
    detSigma.hat <- det(Sigma.hat)
    llik <- llik - 0.5 * (log(detSigma.hat) + crossprod(at,invSigma.hat) %*% at)
  }
  
  nllik <- -llik
  
  return(nllik)
}

est.ewma <- function(lambda, innov){
  out <- optim(lambda, nllik.ewma, lower = 0.001, upper = 0.999, innov = innov, method = "L-BFGS-B", hessian = TRUE)
  llh <- out$value
  lambda.hat <- out$par
  lambda.hat.se <- 1 / sqrt(out$hessian)
  output <- list(lambda.hat = lambda.hat, lambda.hat.se = lambda.hat.se, llh = llh)
}


## ## ## ##
## Univariate volatility models ----
## ## ## ##

UnivVola <- function(returns, width = 30, lambda = 0.94, type = c("RiskMetrics", "WeightedAverage", "MovingAverage"), center = FALSE, exchange = "UnitedStates/NYSE"){
  # Check validity of model variant
  type <- match.arg(type, c("RiskMetrics", "WeightedAverage", "MovingAverage"))
  
  # Set lambda.se to NA
  lambda.se <- NA
  llh <- NA
  
  # Center returns for comparability reasons
  if (center){
    returns <- scale(returns, center = TRUE, scale = FALSE)
  }
  
  # Estimate lambda if desired
  if (type %in% c("RiskMetrics", "WeightedAverage")){
    
    if(lambda == 0 || lambda >= 1){
      stop("lambda must be between 0 and 1. If a negative value is supplied, lambda will be estimated from the data.")
    }
    
    if(lambda < 0){
      lambda.est <- est.ewma(0.95, as.matrix(returns))
      lambda <- lambda.est$lambda.hat
      lambda.se <- lambda.est$lambda.hat.se
      llh <- lambda.est$llh
      
      tval <- lambda/lambda.se
      matcoef <- cbind(lambda, lambda.se, tval, 2 * (1 - pnorm(abs(tval))))
      dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
      cat("\nCoefficient(s):\n")
      printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
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
    NBD <- NextBusinessDay(returns, exchange = exchange)
    newrow <- zoo(NA, NBD)
    newindex <- suppressWarnings(rbind(returns, newrow))
    sigma_t2 <- zoo(sigma_t2, index(newindex))
    names(sigma_t2) <- "Variance"
    
    # Set width to NA, whatever the original value is
    width <- NA
  }
  
  if (type == "WeightedAverage") {
    # Compute correction factor and weights
    CF <- nlambda/(lambda*(1-lambda^width))
    s <- c(width:1)
    lambda <- lambda^s
    
    # Compute conditional variances
    sigma_t2 <- rollapplyr(returns^2, width = width, FUN = function(x) CF * sum(lambda * x), fill = NA, by.column = FALSE)
    
    # Add empty row
    NBD <- NextBusinessDay(returns, exchange = exchange)
    newrow <- zoo(NA, NBD)
    newindex <- suppressWarnings(rbind(returns, newrow))
    sigma_t2 <- zoo(sigma_t2, index(newindex))
    names(sigma_t2) <- "Variance"
    
    # Lead series by one observation
    sigma_t2 <- lag(sigma_t2, -1, na.pad = TRUE)
  }
  
  if (type =="MovingAverage"){
    # Compute variance
    sigma_t2 <- rollapplyr(returns^2, width = width, FUN = mean, fill = NA, by.column = FALSE)
    
    # Add empty row
    NBD <- NextBusinessDay(returns, exchange = exchange)
    newrow <- zoo(NA, NBD)
    newindex <- suppressWarnings(rbind(returns, newrow))
    sigma_t2 <- zoo(sigma_t2, index(newindex))
    names(sigma_t2) <- "Variance"
    
    # Lead series by one observation
    sigma_t2 <- lag(sigma_t2, -1, na.pad = TRUE)
    
    # Set lambda to NA, whatever the original input is
    lambda <- NA
  }
  
  # Create and return output object
  object <- list(Variance = sigma_t2, Returns = returns, variant = match.arg(type), width = width, 
                 lambda = lambda, lambda.se = lambda.se, llh = llh, centered = center)
  attr(object, "class") <- c("UnivVola")
  
  return(object)
}


## ## ## ##
## Multivariate volatility models ----
## ## ## ##

MultiEWMA <- function(returns, lambda = 0.94, center = FALSE, exchange = "UnitedStates/NYSE"){
  # Extract necessary information from the return series
  x <- coredata(returns)
  n <- dim(x)[1] + 1
  d <- dim(x)[2]
  lambda.se <- NA
  llh <- NA
  
  # Center returns for comparability reasons
  if (center){
    x <- scale(x, center = TRUE, scale = FALSE)
  }
  
  # Check lambda
  if(lambda == 0 || lambda >= 1){
    stop("lambda must be between 0 and 1. If a negative value is supplied, lambda will be estimated from the data.")
  }
  
  if(lambda < 0){
    lambda.est <- est.ewma(0.95, as.matrix(returns))
    lambda <- lambda.est$lambda.hat
    lambda.se <- lambda.est$lambda.hat.se
    llh <- lambda.est$llh
    
    tval <- lambda/lambda.se
    matcoef <- cbind(lambda, lambda.se, tval, 2 * (1 - pnorm(abs(tval))))
    dimnames(matcoef) = list(names(tval), c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
    cat("\nCoefficient(s):\n")
    printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
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
  NBD <- NextBusinessDay(returns, exchange = exchange)
  newrow <- zoo(NA, NBD)
  newindex <- rbind(returns[, 1], newrow)
  SIGMA_t <- zoo(SIGMA_t, index(newindex))
  x <- zoo(x, index(returns))
  
  object <- list(Variances = SIGMA_t, Returns = x, lambda = lambda, lambda.se = lambda.se, llh = llh, centered = center)

  attr(object, "class") <- "MultiEWMA"
  return(object)
}
