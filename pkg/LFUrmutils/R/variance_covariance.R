## ## ## ##
## Variance-covariance matrices ----
## ## ## ##

## Generic function ----

varcov <- function(object, ...){
  UseMethod("varcov")
}


## Multivariate EWMA model ----

# Variance-Covariance matrix (Sigma)
varcov.MultiEWMA <- function(object, offdiagonal = TRUE, duplicates = TRUE){
  
  # Extract variance-covariance matrix
  varcov <- object$Variances

  # Delete off-diagonal elements if not necessary
  if(offdiagonal == FALSE){
    n <- dim(varcov)[1]
    c <- sqrt(dim(varcov)[2])
    diago <- matrix(NA, nrow = n, ncol = c)
    for (k in 1:c){
      # Keep only diagonal elements
      diago[, k] <- varcov[, grep(paste0(k,k), colnames(varcov))]
    }
    varcov <- diago
    # Update names
    names <- matrix(NA, nrow = c, ncol = 1)
    for(k in 1:c){
      names[k, ] <- paste0("Sigma", k,k)
    }
    colnames(varcov) <- names
  }
  
  # Delete duplicate columns
  if(offdiagonal == TRUE){
    if(duplicates == FALSE){
      corematrix <- matrix(coredata(varcov), nrow = dim(varcov)[1], ncol = dim(varcov)[2])
      duplicate <- duplicated(round(t(corematrix), 10))
      varcov <- varcov[, !duplicate]
    }
  }
  return(varcov)
}



## OGARCH model ----

    #### TO BE DONE/UPDATED

# Variance-Covariance matrix (Sigma)
OGvarcov <- function(OGobject, offdiagonal = TRUE, duplicates = TRUE){
  VarCov <- OGobject@H
  n <- length(VarCov)
  c <- sqrt(length(VarCov[[1]]))
  
  names <- matrix(NA, c, c)
  for (i in 1:c){
    for (j in 1:c) {
      names[j,i] <- paste0("Sigma", j,i)
    }
  }
  names <- c(t(names))
  
  og.varcov <- matrix(NA, n, c^2)
  colnames(og.varcov) <- names
  
  for(i in 1:n){
    og.varcov[i, ] <- as.vector(VarCov[[i]])
  }
  
  og.varcov <- zoo(og.varcov, 
                   as.Date( attr(OGobject@ica$S, "dimnames")[[1]] ))
  
  if(offdiagonal == FALSE){
    OG.diag <- matrix(NA, nrow = n, ncol = c)
    for (k in 1:c){
      OG.diag[, k] <- og.varcov[, grep(paste0(k,k), colnames(og.varcov))]
    }
    og.varcov <- OG.diag
    names <- matrix(NA, nrow = c, ncol = 1)
    for(k in 1:c){
      names[k, ] <- paste0("Sigma", k,k)
    }
    colnames(og.varcov) <- names
  }
  
  if(offdiagonal == TRUE){
    if(duplicates == FALSE){
      corematrix <- matrix(coredata(og.varcov), nrow = dim(og.varcov)[1], ncol = dim(og.varcov)[2])
      duplicate <- duplicated(round(t(corematrix), 10))
      og.varcov <- og.varcov[, !duplicate]
    }
  }
  
  return(og.varcov)
}



## Dynamic conditional volatility model ----

    #### TO BE DONE/UPDATED

# Variance-Covariance matrix (Sigma)
DCCvarcov <- function(DCCobject, offdiagonal = TRUE, duplicates = TRUE){
  
  varVec <- DCCobject$h
  corMat <- DCCobject$DCC
  
  c <- dim(varVec)[2]
  n <- dim(varVec)[1]
  
  names <- matrix(NA, c, c)
  for (i in 1:c){
    for (j in 1:c) {
      names[j,i] <- paste0("Sigma", j,i)
    }
  }
  names <- c(t(names))
  
  DCC.varcov <- matrix(NA, n, c^2)
  for (i in 1:n){
    sdMat <- diag(sqrt(varVec[i, ]))
    DCC.varcov[i, ] <- c(sdMat %*% matrix(corMat[i, ], c, c, byrow = TRUE) %*% t(sdMat))
  }
  DCC.varcov <- zoo(DCC.varcov, as.Date(attr(DCCobject$std.resid, which = "dimnames")[[1]]))
  colnames(DCC.varcov) <- names
  
  if(offdiagonal == FALSE){
    DCC.diag <- matrix(NA, nrow = n, ncol = c)
    for (k in 1:c){
      DCC.diag[, k] <- DCC.varcov[, grep(paste0(k,k), colnames(DCC.varcov))]
    }
    DCC.varcov <- DCC.diag
    names <- matrix(NA, nrow = c, ncol = 1)
    for(k in 1:c){
      names[k, ] <- paste0("Sigma", k,k)
    }
    colnames(DCC.varcov) <- names
  }
  if(offdiagonal == TRUE){
    if(duplicates == FALSE){
      corematrix <- matrix(coredata(DCC.varcov), nrow = dim(DCC.varcov)[1], ncol = dim(DCC.varcov)[2])
      duplicate <- duplicated(round(t(corematrix), 10))
      DCC.varcov <- DCC.varcov[, !duplicate]
    }
  }
  return(DCC.varcov)
}

