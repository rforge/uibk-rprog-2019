## ## ## ##
## Residuals ----
## ## ## ##

## Univariate volatility models ----
residuals.UnivVola <- function(object, standardize = TRUE, na.action = "na.pass", ...){
  if (standardize) {
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

# Multivariate EWMA model ----
# residuals.MultiEWMA <- function(object, standardize = TRUE, na.action = "na.pass", ...){
#   if (standardize) {
#     n <- dim(object$Variances)[1]
#     cols <- sqrt(dim(object$Variances)[2])
# 
#     # Merge returns with square root of diagonal elements of variance-covariance matrix
#     diagonals <- seq.int(from = 1, to = cols^2, by = cols+1)
#     merged <- merge(object$Returns, sqrt(object$Variances[, diagonals]))
#     
#     # Determine corresponding columns
#     retcols <- seq.int(from = 1, to = cols, by = 1)
#     sdcols <- seq.int(from = cols+1, to = dim(merged)[2], by = 1)
#     
#     # Compute residuals
#     resids <- merged[, retcols] / merged[, sdcols]
# 
#     # Delete NA at the end if desired
#     if (na.action == "na.trim"){
#       resids <- na.trim(resids)
#     }
#     
#     # Format output
#     colnames(resids) <- colnames(object$Returns)
#     
#   } else { # end if(standardize == TRUE)
#     message("Note that non-standardized residuals simply correspond to the orginial returns.")
#     resids <- object$Returns
#   }
#   return(resids)
# }

## Ruppert and Matteson (2015, p. 431 -- 443)
residuals.MultiEWMA <- function(object, standardize = TRUE, na.action = "na.pass", ...){
  if (standardize) {
    n <- dim(object$Variances)[1]
    cols <- sqrt(dim(object$Variances)[2])
    
    # Merge returns and variances and determine corresponding columns
    merged <- merge(object$Returns, object$Variances)
    retcols <- seq.int(from = 1, to = cols)
    varcols <- seq.int(from = cols+1, to = dim(merged)[2])
    
    # Function to compute residuals (using a singular value decomposition)
    residualcomputer <- function(x, retcols = retcols, varcols = varcols, nrow = nrow, ncol = ncol){
      periodresidual <- x[retcols] %*% matrix.sqrt.inv(matrix(x[varcols], byrow = TRUE, nrow = cols, ncol = cols)) 
      return(periodresidual)
    }
    
    # Compute residuals (using a singular value decomposition)
    resids <- rollapplyr(merged, width = 1, FUN = residualcomputer, by.column = FALSE, 
                         # Additional arguments
                         retcols = retcols, varcols = varcols, nrow = cols, ncol = cols)
    
    if (na.action == "na.trim"){
      resids <- na.trim(resids)
    }
    
    # Format output
    colnames(resids) <- colnames(object$Returns)
    
  } else { # end if(standardize == TRUE)
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

