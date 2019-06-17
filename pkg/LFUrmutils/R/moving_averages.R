## ## ## ##
## Univariate volatility models ----
## ## ## ##

MovingAverage <- function(returns, width = 30, center = FALSE){
  # y: (daily) log returns

  # Center returns for comparability reasons
  if (center == TRUE){
    returns <- scale(returns, center = TRUE, scale = FALSE)
  }
  
  # Compute variance
  sigma_t2 <- rollmeanr(returns^2, k = width, fill = NA)
  
  # Add empty row
  NBD <- NextBusinessDay(y)
  newrow <- zoo(NA, NBD)
  sigma_t2 <- rbind(sigma_t2, newrow)
  
  # Lead series by one observation
  sigma_t2 <- lag(sigma_t2, -1)
  
  # Create and return output object
  object <- list(Variance = sigma_t2, Returns = returns)
  attr(object, "class") <- "UnivMA"
  return(object)
}


EWMA <- function(returns, width = 30, lambda = 0.94, type = "RiskMetrics", center = FALSE){
  "%!in%" <- Negate("%in%")
  
  if (type %!in% c("RiskMetrics", "WeightedAverage")){
    stop("Error: Type must be either RiskMetrics or WeightedAverage.")
  }
  
  # Center returns for comparability reasons
  if (center == TRUE){
    returns <- scale(returns, center = TRUE, scale = FALSE)
  }
  
  # Compute 1-lambda for efficiency reasons
  nlambda <- 1 - lambda
  
  
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
    newindex <- merge(y, newrow)
    sigma_t2 <- zoo(sigma_t2, index(newindex))
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
    sigma_t2 <- rbind(sigma_t2, newrow)
    
    # Lead series by one observation
    sigma_t2 <- lag(sigma_t2, -1)
    
  }
  
  # Create and return output object
  object <- list(Variance = sigma_t2, Returns = returns, lambda = lambda, centered = center)
  attr(object, "class") <- "UnivEWMA"
  
  if (type == "RiskMetrics"){
    attr(object, "variant") <- "RiskMetrics"
  }
  
  if (type == "WeightedAverage"){
    attr(object, "variant") <- "WeightedAverage"
  }
  
  return(object)
}



## ## ## ##
## Multivariate volatility models ----
## ## ## ##

EWMAmult <- function(returns, lambda = 0.94, center = FALSE){
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
  nlambda <- 1-lambda
  SIGMA_t <- matrix(nrow = n, ncol = d^2)
  colnames(SIGMA_t) <- names
  sigma <- cov(x)
  SIGMA_t[1, ] <- as.vector(sigma)
  
  # Core of the function: Compute EWMA
  for(i in 2:n){
    sigma <- lambda * sigma + nlambda * tcrossprod(x[i-1, ])
    SIGMA_t[i, ] <- as.vector(sigma)
  }
  SIGMA_t <- zoo(SIGMA_t, index(returns))
  
  # Create new index
  NBD <- NextBusinessDay(y)
  newrow <- zoo(NA, NBD)
  newindex <- rbind(returns[, 1], newrow)
  SIGMA_t <- zoo(SIGMA_t, index(newindex))
  
  # SIGMA_t <- lag(SIGMA_t, -1)
  object <- list(Variances = SIGMA_t, Returns = returns, lambda = lambda, centered = center)

  attr(object, "class") <- "MultiEWMA"
  return(object)
}




##################################################################

library(tseries)
library(zoo)

p1 <- get.hist.quote("msft", "2000-01-01", "2009-12-31",
                     quote = "Adjusted", quiet = TRUE)
p2 <- get.hist.quote("ibm", "2000-01-01", "2009-12-31",
                     quote = "Adjusted", quiet = TRUE)

# Combine both series
p <- merge.zoo(p1,p2)


# Compute log returns
y <- diff(log(p)) # * 100
colnames(y) <- c("msft", "ibm")

# Compute de-meaned returns
# Since y is a matrix, we have to use the scale function to compute column-wise means
y <- scale(y, center = TRUE, scale = FALSE)

ma <- MovingAverage(y$msft)

head(ma$Variance, 10)
tail(ma$Variance, 10)

ewma <- EWMA(y$msft, type = "RiskMetrics")

plot(ewma$Variance)

head(ewma$Variance, 10)
tail(ewma$Variance, 10)

head(ewma$Returns, 10)
tail(ewma$Returns, 10)

length(ewma$Returns)
length(ewma$Variance)

dim(ewma$Returns)
dim(ewma$Variances)

View(ewma[["Variances"]])

ewma <- EWMAmult(y)

ewma <- MTS::EWMAvol(y, lambda = 0.94)

tail(ewma$Sigma.t)
dim(ewma$Sigma.t)

rbind.zoo

weekdays(end(y), abbreviate = TRUE)


RQuantLib::isHoliday("UnitedStates/NYSE", end(y)+2)

NextBusinessDay <- function(x, calendar = "UnitedStates/NYSE"){
  i <- FALSE
  NBD <- end(x)
  while(i == FALSE){
    i <- RQuantLib::isBusinessDay(calendar, NBD + 1)
    NBD <- NBD + 1
  }
  return(NBD)
}

NextBusinessDay(y)


res_ma <- residuals.UnivMA(ma, standardize = FALSE)

head(res_ma)
tail(res_ma)

plot(res_ma)

res_ewma <- residuals(ewma, standardize = FALSE)

head(res_ewma)
tail(res_ewma)

plot(res_ewma)

##################################################################
