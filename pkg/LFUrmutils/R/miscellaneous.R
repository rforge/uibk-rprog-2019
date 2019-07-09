## ## ## ##
## Determine next business day ----
## ## ## ##

NextBusinessDay <- function(x){
  i <- FALSE
  NBD <- end(x) + 1
  while(i == FALSE){
    i <- RQuantLib::isHoliday("UnitedStates/NYSE", NBD)
    NBD <- NBD + 1
  }
  return(NBD)
}


## ## ## ##
## Matrix manipulation ----
## ## ## ##

# Matrix square-root
matrix.sqrt <- function(A){
  sva <- svd(A)
  if (min(sva$d)>=0){
    Asqrt <- t(tcrossprod(sva$v) * sqrt(sva$d))
  } else {
    stop("Matrix square root is not defined")
  }
  return(Asqrt)
}

# Inverse of the matrix square-root
matrix.sqrt.inv <- function(A){
  sva <- svd(A)
  if (min(sva$d)>=0)
    Asqrt <- t(tcrossprod(sva$v) / sqrt(sva$d))
  else
    stop("Matrix square root is not defined")
  return(Asqrt)
}

## ## ## ##
## Likelihood ratio tests ----
## ## ## ##

lrtest <- function(object, ...){
  UseMethod("lrest")
}

# fGARCH
lrtest.fGARCH <- function(object1, object2){
  x_llh <- as.numeric(abs(object1@fit$llh))
  y_llh <- as.numeric(abs(object2@fit$llh))
  
  x_par <- length(object1@fit$par)
  y_par <- length(object2@fit$par)
  
  par <- matrix(nrow = 1, ncol = 3)
  colnames(par) <- c("LR statistic", "Restrictions", "Pr(>Chisq)")
  
  par[1, 1] <- 2 * abs(x_llh - y_llh)
  par[1, 2] <- abs(x_par - y_par)
  par[1, 3] <- round(pchisq(par[1], par[2], lower.tail = FALSE), 3)
  par[1, 1] <- round(par[1, 1], 3)
  
  return(par)
}