## ## ## ##
## Volatilities ----
## ## ## ##

## Generic function ----
vola <- function(object, ...){
  UseMethod("vola")
}

## Default function ----
# volatility.default <- function(object) {
#   cat("volatility is currently only defined for objects of the following classes:")
#   methods(volatility)
# }

## Univariate volatility models ----
vola.UnivVola <- function(object){
  sig <- sqrt(object$Variance)
  return(sig)
}

## fGARCH model ----
vola.fGARCH <- function(object){
  sig <- zoo(x = object@sigma.t, 
             order.by = as.Date(attr(object@data, "names")))
  return(sig)
}


## Multivariate EWMA model ----
vola.MultiEWMA <- function(object, offdiagonal = TRUE, duplicates = TRUE){
  
  sig <- object$Variances
  colnames(sig) <- sub("Sigma", "Volatility", colnames(sig))

  d <- as.integer(sqrt(dim(sig)[2]))
  
  # Core of the function
  for (i in 1:d){
    # Adjust diagonal elements
    sig[, grep(paste0(i,i), colnames(sig))] <- sqrt(sig[, grep(paste0(i,i), colnames(sig))])
  }
  
  if(offdiagonal == FALSE){
    # Delete off-diagonal elements if not necessary
      # Create sequence of 11, 22, etc.
      diago <- seq.int(1,d) + 10L*seq.int(1,d)
    sig <- sig[, grep(paste(diago, collapse="|"), colnames(sig))]
  } else {
    # Otherwise keep off-diagonal elements
    for (i in 1:d){
      for (j in seq.int(i,d)){
        if (i != j){
          # Compute volatility of off-diagonal elements (only upper elements)
          sig[, grep(paste0(i, j), colnames(sig))] <- sig[, grep(paste0(i,j), colnames(sig))]/sqrt(sig[, grep(paste0(i,i), colnames(sig))] * sig[, grep(paste0(j,j), colnames(sig))])
        }
    }
  }
    if(duplicates == TRUE){
      for (i in 1:d){
        for (j in seq.int(i,d)){
          if (i != j){
            # Copy upper elements to lower elements
            sig[, grep(paste0(j, i), colnames(sig))] <- sig[, grep(paste0(i, j), colnames(sig))]  
          }
        }
      }
    } else {
      for (i in 1:d){
        for (j in seq.int(i,d)){
          if (i != j){
            # Otherwise delete lower elements
            sig <- sig[, -grep(paste0(j, i), colnames(sig))]
          }
        }
      }
    }
  }
  return(sig)
}


# ## OGARCH model ----
# 
#     ## TO BE DONE
# 
# # Volatility, i.e. Sigma^(1/2)
# volatility.OG <- function(object, offdiagonal = TRUE, duplicates = TRUE){
#   
#   sigma <- varcov(object)
#   
#   n <- dim(sigma)[1]
#   c <- sqrt(dim(sigma)[2])
#   TimeStamps <- as.Date( attr(object@ica$S, "dimnames")[[1]] )
#   names <- matrix(NA, c, c)
#   for (i in 1:c){
#     for (j in 1:c) {
#       names[j,i] <- paste0("sigma", j,i)
#     }
#   }
#   names <- c(t(names))
#   volat <- matrix(NA, nrow = n, ncol = c^2)
#   colnames(volat) <- names
#   for (i in 1:n){
#     # volat[i, ] <- c(matrix.sqrt(matrix(OG.sigma[i, ], 
#     #                               nrow = c, ncol = c, byrow = TRUE)))
#     volat[i, ] <- c(sqrt(matrix(sigma[i, ], 
#                                   nrow = c, ncol = c, byrow = TRUE)))
#   }
#   volat <- zoo(volat, TimeStamps)
#   
#   if(offdiagonal == FALSE){
#     diago <- matrix(NA, nrow = n, ncol = c)
#     for (k in 1:c){
#       diago[, k] <- volat[, grep(paste0(k,k), colnames(volat))]
#     }
#     volat <- diago
#     names <- matrix(NA, nrow = c, ncol = 1)
#     for(k in 1:c){
#       names[k, ] <- paste0("sigma", k,k)
#     }
#     colnames(volat) <- names
#   }
#   
#   if(offdiagonal == TRUE){
#     if(duplicates == FALSE){
#       corematrix <- matrix(coredata(volat), nrow = dim(volat)[1], ncol = dim(volat)[2])
#       duplicate <- duplicated(round(t(corematrix), 10))
#       volat <- volat[, !duplicate]
#     }
#   }
#   return(volat)
# }
# 
# 
# ## Dynamic conditional volatility model ----
# 
#     ## TO BE DONE
# 
# # Volatility, i.e. Sigma^(1/2)
# volatility.DCC <- function(object, offdiagonal = TRUE, duplicates = TRUE){
#   sigma <-  varcov(object)
#   n <- dim(sigma)[1]
#   c <- sqrt(dim(sigma)[2])
#   TimeStamps <- index(sigma)
#   names <- matrix(NA, c, c)
#   for (i in 1:c){
#     for (j in 1:c) {
#       names[j,i] <- paste0("sigma", j,i)
#     }
#   }
#   names <- c(t(names))
#   volat <- matrix(NA, nrow = n, ncol = c^2)
#   colnames(volat) <- names
#   for (i in 1:n){
#     # volat[i, ] <- c(matrix.sqrt(matrix(DCC.sigma[i, ], 
#     #                                nrow = c, ncol = c, byrow = TRUE)))
#     volat[i, ] <- c(sqrt(matrix(DCC.sigma[i, ], 
#                                    nrow = c, ncol = c, byrow = TRUE)))
#   }
#   volat <- zoo(volat, TimeStamps)
#   
#   if(offdiagonal == FALSE){
#     diago <- matrix(NA, nrow = n, ncol = c)
#     for (k in 1:c){
#       diago[, k] <- volat[, grep(paste0(k,k), colnames(volat))]
#     }
#     volat <- diago
#     names <- matrix(NA, nrow = c, ncol = 1)
#     for(k in 1:c){
#       names[k, ] <- paste0("sigma", k,k)
#     }
#     colnames(volat) <- names
#   }
#   if(offdiagonal == TRUE){
#     if(duplicates == FALSE){
#       corematrix <- matrix(coredata(volat), nrow = dim(volat)[1], ncol = dim(volat)[2])
#       duplicate <- duplicated(round(t(corematrix), 10))
#       volat <- volat[, !duplicate]
#     }
#   }
#   return(volat)
# }

