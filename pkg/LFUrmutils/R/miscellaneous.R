## ## ## ##
## Determine next business day ----
## ## ## ##

NextBusinessDay <- function(x, exchange = "UnitedStates/NYSE"){
  i <- FALSE
  NBD <- end(x) + 1
  while(!i){
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
## fGARCH methods ----
## ## ## ##

"$.fGARCH" <- function(object, x) {
  wi <- pmatch(x, colnames(object))
  if (is.na(wi)) 
    NULL
  else object[, wi]
}
