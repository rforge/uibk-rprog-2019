## ## ## ##
## Residuals ----
## ## ## ##

## Generic function ----

# residuals <- function(object, ...){
#   UseMethod("residuals")
# }



## Univariate moving average ----
residuals.UnivMA <- function(object, standardize = TRUE) {
  if (standardize == TRUE) {
    merged <- merge(object$Returns, sqrt(object$Variance))
    resids <- merged[, 1]/merged[, 2]
  } else {
    resids <- object$Returns
  }
  return(resids)
}


## Univariate EWMA model ----

residuals.UniEWMA <- function(object, standardize = TRUE) {
  if (standardize == TRUE) {
    merged <- merge(object$Returns, sqrt(object$Variance))
    resids <- merged[, 1]/merged[, 2]
  } else {
    resids <- object$Returns
  }
  return(resids)
}

## Multivariate EWMA model ----

residuals.MultiEWMA <- function(object, standardize = TRUE) {
  if (standardize == TRUE) {
    n <- dim(object$Variances)[1]
    c <- sqrt(dim(object$Variances)[2])
    
    # Compute relevant standard deviations
    StdDev <- matrix(NA, nrow = n, ncol = c)
    for (k in 1:c) {
      StdDev[, k] <- object$Variances[, grep( paste0(k,k), colnames(object$Variances) )]
      StdDev[, k] <- sqrt(StdDev[, k])
    }
    StdDev <- zoo(StdDev, index(object$Variances))
    
    # Compute residuals
    resids <- matrix(NA, nrow = n, ncol = c)
    for (k in 1:c) {
      merged <- merge(object$Returns[, k], StdDev[, k])
      resids[, k] <- merged[, 1]/merged[, 2]
    }
    resids <- zoo(resids, index(object$Variances))
    colnames(resids) <- colnames(object$Returns)
  } else {
    message("Note that non-standardized residuals simply correspond to the orginial returns.")
    resids <- object$Returns
  }
  return(resids)
}



## OGARCH model ----

#### TO BE DONE

# Residuals
OGresiduals <- function(object){
  resid <- zoo( residuals(object), 
                  as.Date( attr(object@ica$S, "dimnames")[[1]] ) )
  colnames(resid) <- attr(object@ica$X, "dimnames")[[2]]
  return(resid)
}



## Dynamic conditional volatility model ----

#### TO BE DONE

# Residuals 
DCCresiduals <- function(DCCobject){
  dcc.res <- zoo(DCCobject$std.resid, 
                 as.Date( attr(DCCobject$std.resid, "dimnames")[[1]] ) )
  # dcc.res <- zoo( residuals(DCCobject), 
  #                 as.Date( attr(DCCobject$std.resid, "dimnames")[[1]] ) )
  colnames(dcc.res) <- attr(DCCobject$std.resid, "dimnames")[[2]]
  return(dcc.res)
}

