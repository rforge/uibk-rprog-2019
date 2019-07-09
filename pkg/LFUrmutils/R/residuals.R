## ## ## ##
## Residuals ----
## ## ## ##

## Univariate volatility models ----
residuals.UnivVola <- function(object, standardize = TRUE, na.action = "na.pass", ...){
  if (standardize == TRUE) {
    merged <- merge(object$Returns, sqrt(object$Variance))
    if(na.action == "na.trim"){
      merged <- na.trim(merged)
    }
    res <- merged[, 1]/merged[, 2]
  } else {
    message("Note that non-standardized residuals simply correspond to the orginial returns.")
    res <- object$Returns
  }
  return(res)
}

## Multivariate EWMA model ----
residuals.MultiEWMA <- function(object, standardize = TRUE, na.action = "na.pass", ...){
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
    if (na.action == "na.trim"){
      resids <- na.trim(resids)
    }
    colnames(resids) <- colnames(object$Returns)
  } else {
    message("Note that non-standardized residuals simply correspond to the orginial returns.")
    resids <- object$Returns
  }
  return(resids)
}



## OGARCH model ----

#### TO BE DONE

# # Residuals
# OGresiduals <- function(object){
#   resid <- zoo( residuals(object), 
#                   as.Date( attr(object@ica$S, "dimnames")[[1]] ) )
#   colnames(resid) <- attr(object@ica$X, "dimnames")[[2]]
#   return(resid)
# }
# 
# 
# 
# ## Dynamic conditional volatility model ----
# 
# #### TO BE DONE
# 
# # Residuals 
# DCCresiduals <- function(DCCobject){
#   dcc.res <- zoo(DCCobject$std.resid, 
#                  as.Date( attr(DCCobject$std.resid, "dimnames")[[1]] ) )
#   # dcc.res <- zoo( residuals(DCCobject), 
#   #                 as.Date( attr(DCCobject$std.resid, "dimnames")[[1]] ) )
#   colnames(dcc.res) <- attr(DCCobject$std.resid, "dimnames")[[2]]
#   return(dcc.res)
# }

