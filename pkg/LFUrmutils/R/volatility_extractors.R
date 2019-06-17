## ## ## ##
## Volatilities ----
## ## ## ##

## Generic function ----

vola <- function(object, ...){
  UseMethod("vola")
}

## Univariate MA model ----

    ### TO BE DONE


## Univariate EWMA model ----

    ### TO BE DONE


## Multivariate EWMA model ----

# Volatilities, i.e. Sigma^(1/2)
vola.MultiEWMA <- function(object, offdiagonal = TRUE, duplicates = TRUE){
  
    sig <- object$Variances
    # Sig <- coredata(object$Variances)
    # TimeStamps <- index(object$Variances)
  
  # Initialize output
  # n <- dim(Sig)[1]
  c <- sqrt(dim(sig)[2])
  # names <- matrix(NA, c, c)
  # for (i in 1:c){
  #   for (j in 1:c) {
  #     names[j,i] <- paste0("sigma", j,i)
  #   }
  # }
  # names <- c(t(names))
  # volat <- matrix(NA, nrow = n, ncol = c^2)
  # colnames(volat) <- names
  
  # Core of the function
  for (i in 1:c){
    # Adjust diagonal elements
    sig[, grep(paste0(i,i), colnames(sig))] <- sqrt(sig[, grep(paste0(i,i), colnames(sig))])
  }
  
  # Delete off-diagonal elements if not necessary
  if(offdiagonal == FALSE){
    n <- dim(sig)[1]
    c <- sqrt(dim(sig)[2])
    diago <- matrix(NA, nrow = n, ncol = c)
    for (i in 1:c){
      # Keep only diagonal elements
      diago[, i] <- sig[, grep(paste0(i,i), colnames(sig))]
    }
    sig <- diago
  }
  
  # Delete duplicate columns (e.g. keep only sigma_12, and delete sigma_21)
  if(offdiagonal == TRUE){
    if(duplicates == FALSE){
      corematrix <- matrix(coredata(sig), nrow = dim(sig)[1], ncol = dim(sig)[2])
      duplicate <- duplicated(round(t(corematrix), 10))
      sig <- sig[, !duplicate]
    }
  }
  colnames(sig) <- sub("Sigma", "Volatility", colnames(sig))
  
  return(sig)
}



## OGARCH model ----

    ## TO BE DONE

# Volatility, i.e. Sigma^(1/2)
OGvola <- function(OGobject, offdiagonal = TRUE, duplicates = TRUE){
  
  OG.sigma <- OGvarcov(OGobject)
  
  n <- dim(OG.sigma)[1]
  c <- sqrt(dim(OG.sigma)[2])
  TimeStamps <- as.Date( attr(OGobject@ica$S, "dimnames")[[1]] )
  names <- matrix(NA, c, c)
  for (i in 1:c){
    for (j in 1:c) {
      names[j,i] <- paste0("sigma", j,i)
    }
  }
  names <- c(t(names))
  OG.vola <- matrix(NA, nrow = n, ncol = c^2)
  colnames(OG.vola) <- names
  for (i in 1:n){
    # OG.vola[i, ] <- c(matrix.sqrt(matrix(OG.sigma[i, ], 
    #                               nrow = c, ncol = c, byrow = TRUE)))
    OG.vola[i, ] <- c(sqrt(matrix(OG.sigma[i, ], 
                                  nrow = c, ncol = c, byrow = TRUE)))
  }
  OG.vola <- zoo(OG.vola, TimeStamps)
  
  if(offdiagonal == FALSE){
    OG.diag <- matrix(NA, nrow = n, ncol = c)
    for (k in 1:c){
      OG.diag[, k] <- OG.vola[, grep(paste0(k,k), colnames(OG.vola))]
    }
    OG.vola <- OG.diag
    names <- matrix(NA, nrow = c, ncol = 1)
    for(k in 1:c){
      names[k, ] <- paste0("sigma", k,k)
    }
    colnames(OG.vola) <- names
  }
  
  if(offdiagonal == TRUE){
    if(duplicates == FALSE){
      corematrix <- matrix(coredata(OG.vola), nrow = dim(OG.vola)[1], ncol = dim(OG.vola)[2])
      duplicate <- duplicated(round(t(corematrix), 10))
      OG.vola <- OG.vola[, !duplicate]
    }
  }
  return(OG.vola)
}


## Dynamic conditional volatility model ----

    ## TO BE DONE

# Volatility, i.e. Sigma^(1/2)
DCCvola <- function(DCCobject, offdiagonal = TRUE, duplicates = TRUE){
  DCC.sigma <- DCCvarcov(DCCobject)
  n <- dim(DCC.sigma)[1]
  c <- sqrt(dim(DCC.sigma)[2])
  TimeStamps <- index(DCC.sigma)
  names <- matrix(NA, c, c)
  for (i in 1:c){
    for (j in 1:c) {
      names[j,i] <- paste0("sigma", j,i)
    }
  }
  names <- c(t(names))
  DCC.vola <- matrix(NA, nrow = n, ncol = c^2)
  colnames(DCC.vola) <- names
  for (i in 1:n){
    # DCC.vola[i, ] <- c(matrix.sqrt(matrix(DCC.sigma[i, ], 
    #                                nrow = c, ncol = c, byrow = TRUE)))
    DCC.vola[i, ] <- c(sqrt(matrix(DCC.sigma[i, ], 
                                   nrow = c, ncol = c, byrow = TRUE)))
  }
  DCC.vola <- zoo(DCC.vola, TimeStamps)
  
  if(offdiagonal == FALSE){
    DCC.diag <- matrix(NA, nrow = n, ncol = c)
    for (k in 1:c){
      DCC.diag[, k] <- DCC.vola[, grep(paste0(k,k), colnames(DCC.vola))]
    }
    DCC.vola <- DCC.diag
    names <- matrix(NA, nrow = c, ncol = 1)
    for(k in 1:c){
      names[k, ] <- paste0("sigma", k,k)
    }
    colnames(DCC.vola) <- names
  }
  if(offdiagonal == TRUE){
    if(duplicates == FALSE){
      corematrix <- matrix(coredata(DCC.vola), nrow = dim(DCC.vola)[1], ncol = dim(DCC.vola)[2])
      duplicate <- duplicated(round(t(corematrix), 10))
      DCC.vola <- DCC.vola[, !duplicate]
    }
  }
  return(DCC.vola)
}

