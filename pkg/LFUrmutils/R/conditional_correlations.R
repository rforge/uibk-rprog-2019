## ## ## ##
## Conditional correlations ----
## ## ## ##

## Generic function ----
ccor <- function(object, ...){
  UseMethod("ccor")
}

## Default function ----
# ccor.default <- function(object) {
#   cat("ccor is only defined for objects of the following classes:")
#   methods(ccor)
# }


## Multivariate EWMA model ----
ccor.MultiEWMA <- function(object, diagonal = TRUE, duplicates = TRUE){
  
  # Initialize output object
  correlations <- object$Variances
  
  n <- dim(correlations)[1]
  c <- sqrt(dim(correlations)[2])
  
  for (t in 1:n){
    correlations[t, ] <- c(cov2cor(matrix(correlations[t, ], nrow = c, ncol = c, byrow = TRUE)))
  }
  
  if(diagonal == FALSE){
    for (k in 1:c){
      correlations <- correlations[, -grep(paste0(k,k), colnames(correlations))]
    }
  }
  
  if(duplicates == FALSE){
    # Extract and sort column names, and keep only unique ones
    nams <- colnames(correlations)
    nams <- unname(sapply(nams, FUN = function(x) paste(sort(unlist(strsplit(x, ""))), collapse = "")))
    correlations <- correlations[, !duplicated(nams), drop = FALSE]
  }
  
  # Update column names
  colnames(correlations) <- sub("Sigma", "rho", colnames(correlations))
  return(correlations)
}



## OGARCH model ----

    ## TO BE DONE

# OGccor <- function(OGobject, diagonal = TRUE, duplicates = TRUE){
#   OG.sigma <- OGvarcov(OGobject)
#   Sig.t <- coredata(OG.sigma)
#   TimeStamps <- as.Date( attr(OGobject@ica$S, "dimnames")[[1]] )
#   n <- dim(Sig.t)[1]
#   c <- sqrt(dim(Sig.t)[2])
#   names <- matrix(NA, c, c)
#   for (i in 1:c){
#     for (j in 1:c) {
#       names[j,i] <- paste0("rho", j,i)
#     }
#   }
#   names <- c(t(names))
#   OG.cor <- matrix(NA, nrow = n, ncol = c^2)
#   colnames(OG.cor) <- names
#   for (t in 1:n){
#     OG.cor[t, ] <- c(cov2cor(matrix(Sig.t[t, ], 
#                                     nrow = c, ncol = c, byrow = TRUE)))
#   }
#   OG.cor <- zoo(OG.cor, TimeStamps)
#   if(diagonal == FALSE){
#     for (k in 1:c){
#       OG.cor <- OG.cor[, -grep(paste0(k,k), colnames(OG.cor))]
#     }
#   }
#   if(duplicates == FALSE){
#     corematrix <- matrix(coredata(OG.cor), nrow = dim(OG.cor)[1], ncol = dim(OG.cor)[2])
#     duplicate <- duplicated(round(t(corematrix), 10))
#     OG.cor <- OG.cor[, !duplicate]
#   }
#   return(OG.cor)
# }


## Dynamic conditional volatility model ----

    ## TO BE DONE

# DCCccor <- function(DCCobject, diagonal = TRUE, duplicates = TRUE){
#   DCC.ccor <- zoo( DCCobject$DCC, 
#                    as.Date( attr(DCCobject$std.resid, "dimnames")[[1]] ) )
#   n <- dim(DCC.ccor)[1]
#   c <- sqrt(dim(DCC.ccor)[2])
#   names <- matrix(NA, c, c)
#   for (i in 1:c){
#     for (j in 1:c) {
#       names[j,i] <- paste0("rho", j,i)
#     }
#   }
#   names <- c(t(names))
#   colnames(DCC.ccor) <- names
#   
#   if(diagonal == FALSE){
#     for (k in 1:c){
#       DCC.ccor <- DCC.ccor[, -grep(paste0(k,k), colnames(DCC.ccor))]
#     }
#   }
#   
#   if(duplicates == FALSE){
#     corematrix <- matrix(coredata(DCC.ccor), nrow = dim(DCC.ccor)[1], ncol = dim(DCC.ccor)[2])
#     duplicate <- duplicated(round(t(corematrix), 10))
#     DCC.ccor <- DCC.ccor[, !duplicate]
#   }
#   return(DCC.ccor)
# }
