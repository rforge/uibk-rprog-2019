## ## ## ##
## Determine next business day
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
