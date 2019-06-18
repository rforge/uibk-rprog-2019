## ## ## ##
## Mean squared error ----
## ## ## ## 

# Generic function
MSE <- function(object, ...){
  UseMethod("MSE")
}

# Default function
MSE.default <- function(object) {
  cat("MSE is only defined for objects of the following classes:")
  methods(MSE)
}

# Univariate moving average model
MSE.UnivMA <- function(cvar, y = NULL, squared = FALSE, volatility = FALSE){
  # cvar: sigma_t^2 estimates
  # y: (daily) log return series
  if (squared == FALSE) {
    y <- y^2
  }
  
  if (volatility == TRUE){
    cvar <- cvar^2
  }
  
  merg <- merge(cvar, y)
  merg <- na.trim(merg)
  
  mse <- mean((merg[, 2] - merg[, 1])^2)
  return(mse)
}

# Univariate EWMA model
MSE.UnivEWMA <- function(object, y = NULL, squared = FALSE){
  x <- 10
}


# fGARCH
MSE.fGARCH <- function(model, ret, squared = FALSE){
  # model: object of class fGARCH
  # ret: return series
  # squared: set to TRUE, if returns are already squared
  
  # Compute squared returns, if necessary
  if (squared == FALSE) {
    ret <- ret^2
  }
  
  # fitGarch computes conditional volatility by default
  # Therefore, we have to compute the conditional variance manually
  sigma_t2 <- zoo(model@sigma.t^2, time(ret))
  
  # Merge conditional volatility with squared return series
  merg <- merge(sigma_t2, ret)
  
  # Delete trailing missing values
  merg <- na.trim(merg)
  
  # Compute and return mean squared error
  mse <- mean((merg[, 2] - merg[, 1])^2)
  return(mse)
  
}


## ## ## ##
## Likelihood ratio test ----
## ## ## ##

my_LRtest <- function(object1, object2){
  # x, y: two objects of class fGARCH
  
  # Extract log likelihoods
  x_llh <- as.numeric(abs(object1@fit$llh))
  y_llh <- as.numeric(abs(object2@fit$llh))
  
  # Determine length of parameter vector
  x_par <- length(object1@fit$par)
  y_par <- length(object2@fit$par)
  
  # Set up a matrix to store the results
  par <- matrix(nrow = 1, ncol = 3)
  colnames(par) <- c("LR statistic", "Restrictions", "Pr(>Chisq)")
  
  # Compute and save LR statistic
  par[1, 1] <- 2 * abs(x_llh - y_llh)
  
  # Determine and save number of restrictions
  par[1, 2] <- abs(x_par - y_par)
  
  # Compute p-value
  par[1, 3] <- round(pchisq(par[1], par[2], lower.tail = FALSE), 3)
  
  # Additionally round LR statistic
  par[1, 1] <- round(par[1, 1], 2)
  
  # Return results
  return(par)
  
}


## ## ## ##
## Residuals ----
## ## ## ##

# MA and EWMA residuals
residuals.ma <- function(x) {
  x <- 10
}

residuals.ewma <- function(x) {
  x <- 10
}


## ewmaRM_res <- SPret/sqrt(ewmaRM)


# Annualized volatilities
# Define function
AnnualizedVola.fGARCH <- function(x, h = 252, rounded = 3){

  # Get constant (omega)
  numerator <- sum(coef(x)[grepl("omega", names(coef(x)))])
  
  # Get sum of alphas and betas
  denominator <- sum(coef(x)[grepl(paste(c("alpha", "beta"), 
                                         collapse="|"), names(coef(x)))])
  
  # Set up a matrix to store the results
  V <- matrix(NA, nrow = 3, ncol = 1)
  rownames(V) <- c("Unconditional variance", 
                   "Unconditional volatility", 
                   "Annualized volatility")
  
  # Compute the unconditional variance
  V[1,1] <- numerator/(1-denominator)
  
  # Compute the unconditional volatility
  V[2,1] <- sqrt(V[1,1])
  
  # Compute the annualized volatility
  V[3,1] <- V[2,1]*sqrt(h)
  
  if (rounded != FALSE) {
    return(V)
  } else {
    return(round(V, rounded))
  }
  
} # end of function



## Multivariate EWMA model ----

# Variance-covariance matrix
EWMAmult <- function(ReturnSeries, lambda = 0.94, center = FALSE){
  # Extract necessary information from the return series
  x <- coredata(ReturnSeries)
  n <- dim(x)[1] + 1
  c <- dim(x)[2]^2
  
  # If non-centered return series is supplied to the function
  if (center == TRUE){
    x <- scale(x, center = TRUE, scale = FALSE)
  }
  
  # Create names
  names <- matrix(NA, sqrt(c), sqrt(c))
  for (i in 1:sqrt(c)){
    for (j in 1:sqrt(c)) {
      names[j,i] <- paste0("Sigma", j,i)
    }
  }
  names <- c(t(names))
  
  # Initialize output
  SIGMA.t <- matrix(nrow = n, ncol = c)
  colnames(SIGMA.t) <- names
  sigma <- cov(x)
  SIGMA.t[1, ] <- as.vector(sigma)
  
  # Core of the function: Compute EWMA
  for(i in 2:n){
    sigma <- lambda*sigma + (1-lambda)*x[i-1, ] %*% t(x[i-1, ])
    SIGMA.t[i, ] <- as.vector(sigma)
  }
  
  # Adjust output, and add return series
  SIGMA.t <- zoo(SIGMA.t, index(y))
  EWMA <- list(SIGMA.t = SIGMA.t, Returns = ReturnSeries)
  
  return(EWMA)
}


# Residuals
EWMAresiduals <- function(EWMAobject){
  
  # Get variance-covariance matrix
  if(names(EWMAobject)[1] == "Sigma.t"){
    # If EWMA is estimated with EWMAvol() from pkg MTS
    y <- coredata(EWMAobject$return)
    Sig.t <- coredata(EWMAobject$Sigma.t)
    TimeStamps <- index(EWMAobject$return)
    names <- attr(EWMAobject$return, "dimnames")[[2]]
  } else {
    # If EWMA is estimated with the user-defined function
    y <- coredata(EWMAobject$Returns)
    Sig.t <- coredata(EWMAobject$SIGMA.t)
    TimeStamps <- index(EWMAobject$SIGMA.t)
    names <- attr(EWMAobject$Returns, "dimnames")[[2]]
  }
  
  # Initialize loop and output object
  n <- dim(Sig.t)[1]
  c <- sqrt(dim(Sig.t)[2])
  stdResid <- matrix(NA, n, c)
  colnames(stdResid) <- names
  
  # Core of the function: Compute residuals
  for(i in 1:n){
    stdResid[i, ] <- y[i, ] %*% matrix.sqrt.inv(matrix(Sig.t[i, ], 
                                                       nrow = c, ncol = c, byrow = TRUE))
  }
  # Return estimated residuals
  stdResid <- zoo(stdResid, TimeStamps)
  return(stdResid)
} 



# Variance-Covariance matrix (Sigma)
EWMAvarcov <- function(EWMAobject, offdiagonal = TRUE, duplicates = TRUE){
  if(names(EWMAobject)[1] == "Sigma.t"){
    # If EWMAvol() from MTS is used
    EWMA.varcov <- EWMAobject$Sigma.t
    n <- dim(EWMA.varcov)[1]
    c <- sqrt(dim(EWMA.varcov)[2])
    TimeStamps <- index(EWMAobject$return)
    EWMA.varcov <- zoo(EWMA.varcov, TimeStamps)
    names <- matrix(NA, c, c)
    for (i in 1:c){
      for (j in 1:c) {
        names[j,i] <- paste0("Sigma", j,i)
      }
    }
    names <- c(t(names))
    colnames(EWMA.varcov) <- names
  } else {
    # If user-defined EWMA function is used
    EWMA.varcov <- EWMAobject$SIGMA.t
    # n <- dim(EWMA.varcov)[1]
    # c <- sqrt(dim(EWMA.varcov)[2])
    # names <- matrix(NA, c, c)
    # for (i in 1:c){
    #   for (j in 1:c) {
    #     names[j,i] <- paste0("Sigma", j,i)
    #   }
    # }
    # names <- c(t(names))
    # colnames(EWMA.varcov) <- names
  }
  
  # Delete off-diagonal elements if not necessary
  if(offdiagonal == FALSE){
    n <- dim(EWMA.varcov)[1]
    c <- sqrt(dim(EWMA.varcov)[2])
    EWMA.diag <- matrix(NA, nrow = n, ncol = c)
    for (k in 1:c){
      EWMA.diag[, k] <- EWMA.varcov[, grep(paste0(k,k), colnames(EWMA.varcov))]
    }
    EWMA.varcov <- EWMA.diag
    names <- matrix(NA, nrow = c, ncol = 1)
    for(k in 1:c){
      names[k, ] <- paste0("Sigma", k,k)
    }
    colnames(EWMA.varcov) <- names
  }
  
  # Delete duplicate columns (e.g. keep only sigma_12, and delete sigma_21)
  if(offdiagonal == TRUE){
    if(duplicates == FALSE){
      corematrix <- matrix(coredata(EWMA.varcov), nrow = dim(EWMA.varcov)[1], ncol = dim(EWMA.varcov)[2])
      duplicate <- duplicated(round(t(corematrix), 10))
      EWMA.varcov <- EWMA.varcov[, !duplicate]
    }
  }
  return(EWMA.varcov)
}



# Volatilities, i.e. Sigma^(1/2)
EWMAvola <- function(EWMAobject, offdiagonal = TRUE, duplicates = TRUE){
  if(names(EWMAobject)[1] == "Sigma.t"){
    # If EWMAvol() from MTS is used
    Sig.t <- coredata(EWMAobject$Sigma.t)
    TimeStamps <- index(EWMAobject$Sigma.t)
  } else {
    # If user-defined function is used
    Sig.t <- coredata(EWMAobject$SIGMA.t)
    TimeStamps <- index(EWMAobject$SIGMA.t)
  }
  
  # Initialize output
  n <- dim(Sig.t)[1]
  c <- sqrt(dim(Sig.t)[2])
  names <- matrix(NA, c, c)
  for (i in 1:c){
    for (j in 1:c) {
      names[j,i] <- paste0("sigma", j,i)
    }
  }
  names <- c(t(names))
  EWMA.vola <- matrix(NA, nrow = n, ncol = c^2)
  colnames(EWMA.vola) <- names
  
  # Core of the function
  for (i in 1:n){
    # EWMA.vola[i, ] <- c(matrix.sqrt(matrix(Sig.t[i, ], 
    #                                 nrow = c, ncol = c, byrow = TRUE)))
    EWMA.vola[i, ] <- c(sqrt(matrix(Sig.t[i, ], 
                                    nrow = c, ncol = c, byrow = TRUE)))
  }
  EWMA.vola <- zoo(EWMA.vola, TimeStamps)
  
  # Delete off-diagonal elements if not necessary
  if(offdiagonal == FALSE){
    EWMA.diag <- matrix(NA, nrow = n, ncol = c)
    for (k in 1:c){
      EWMA.diag[, k] <- EWMA.vola[, grep(paste0(k,k), colnames(EWMA.vola))]
    }
    EWMA.vola <- EWMA.diag
    names <- matrix(NA, nrow = c, ncol = 1)
    for(k in 1:c){
      names[k, ] <- paste0("sigma", k,k)
    }
    colnames(EWMA.vola) <- names
    EWMA.vola <- zoo(EWMA.vola, TimeStamps)
  }
  
  # Delete duplicate columns (e.g. keep only sigma_12, and delete sigma_21)
  if(offdiagonal == TRUE){
    if(duplicates == FALSE){
      corematrix <- matrix(coredata(EWMA.vola), nrow = dim(EWMA.vola)[1], ncol = dim(EWMA.vola)[2])
      duplicate <- duplicated(round(t(corematrix), 10))
      EWMA.vola <- EWMA.vola[, !duplicate]
    }
  }
  return(EWMA.vola)
}



# Conditional correlations
EWMAccor <- function(EWMAobject, diagonal = TRUE, duplicates = TRUE){
  if(names(EWMAobject)[1] == "Sigma.t"){
    Sig.t <- coredata(EWMAobject$Sigma.t)
    TimeStamps <- index(EWMAobject$Sigma.t)
  } else {
    Sig.t <- coredata(EWMAobject$SIGMA.t)
    TimeStamps <- index(EWMAobject$SIGMA.t)
  }
  n <- dim(Sig.t)[1]
  c <- sqrt(dim(Sig.t)[2])
  
  names <- matrix(NA, c, c)
  for (i in 1:c){
    for (j in 1:c) {
      names[j,i] <- paste0("rho", j,i)
    }
  }
  names <- c(t(names))
  EWMA.cor <- matrix(NA, nrow = n, ncol = c^2)
  colnames(EWMA.cor) <- names
  for (t in 1:n){
    EWMA.cor[t, ] <- c(cov2cor(matrix(Sig.t[t, ], 
                                      nrow = c, ncol = c, byrow = TRUE)))
  }
  EWMA.cor <- zoo(EWMA.cor, TimeStamps)
  if(diagonal == FALSE){
    for (k in 1:c){
      EWMA.cor <- EWMA.cor[, -grep(paste0(k,k), colnames(EWMA.cor))]
    }
  }
  if(duplicates == FALSE){
    corematrix <- matrix(coredata(EWMA.cor), nrow = dim(EWMA.cor)[1], ncol = dim(EWMA.cor)[2])
    duplicate <- duplicated(round(t(corematrix), 10))
    EWMA.cor <- EWMA.cor[, !duplicate]
  }
  return(EWMA.cor)
}



## OGARCH model ----

# Residuals    
OGresiduals <- function(OGobject){
  dcc.res <- zoo( residuals(OGobject), 
                  as.Date( attr(OGobject@ica$S, "dimnames")[[1]] ) )
  colnames(dcc.res) <- attr(OGobject@ica$X, "dimnames")[[2]]
  return(dcc.res)
}

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


# Conditional correlations    
OGccor <- function(OGobject, diagonal = TRUE, duplicates = TRUE){
  OG.sigma <- OGvarcov(OGobject)
  Sig.t <- coredata(OG.sigma)
  TimeStamps <- as.Date( attr(OGobject@ica$S, "dimnames")[[1]] )
  n <- dim(Sig.t)[1]
  c <- sqrt(dim(Sig.t)[2])
  names <- matrix(NA, c, c)
  for (i in 1:c){
    for (j in 1:c) {
      names[j,i] <- paste0("rho", j,i)
    }
  }
  names <- c(t(names))
  OG.cor <- matrix(NA, nrow = n, ncol = c^2)
  colnames(OG.cor) <- names
  for (t in 1:n){
    OG.cor[t, ] <- c(cov2cor(matrix(Sig.t[t, ], 
                                    nrow = c, ncol = c, byrow = TRUE)))
  }
  OG.cor <- zoo(OG.cor, TimeStamps)
  if(diagonal == FALSE){
    for (k in 1:c){
      OG.cor <- OG.cor[, -grep(paste0(k,k), colnames(OG.cor))]
    }
  }
  if(duplicates == FALSE){
    corematrix <- matrix(coredata(OG.cor), nrow = dim(OG.cor)[1], ncol = dim(OG.cor)[2])
    duplicate <- duplicated(round(t(corematrix), 10))
    OG.cor <- OG.cor[, !duplicate]
  }
  return(OG.cor)
}


## Dynamic conditional volatility model ----

# Residuals 
DCCresiduals <- function(DCCobject){
  dcc.res <- zoo(DCCobject$std.resid, 
                 as.Date( attr(DCCobject$std.resid, "dimnames")[[1]] ) )
  # dcc.res <- zoo( residuals(DCCobject), 
  #                 as.Date( attr(DCCobject$std.resid, "dimnames")[[1]] ) )
  colnames(dcc.res) <- attr(DCCobject$std.resid, "dimnames")[[2]]
  return(dcc.res)
}

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


# Conditional correlations  
DCCccor <- function(DCCobject, diagonal = TRUE, duplicates = TRUE){
  DCC.ccor <- zoo( DCCobject$DCC, 
                   as.Date( attr(DCCobject$std.resid, "dimnames")[[1]] ) )
  n <- dim(DCC.ccor)[1]
  c <- sqrt(dim(DCC.ccor)[2])
  names <- matrix(NA, c, c)
  for (i in 1:c){
    for (j in 1:c) {
      names[j,i] <- paste0("rho", j,i)
    }
  }
  names <- c(t(names))
  colnames(DCC.ccor) <- names
  
  if(diagonal == FALSE){
    for (k in 1:c){
      DCC.ccor <- DCC.ccor[, -grep(paste0(k,k), colnames(DCC.ccor))]
    }
  }
  
  if(duplicates == FALSE){
    corematrix <- matrix(coredata(DCC.ccor), nrow = dim(DCC.ccor)[1], ncol = dim(DCC.ccor)[2])
    duplicate <- duplicated(round(t(corematrix), 10))
    DCC.ccor <- DCC.ccor[, !duplicate]
  }
  return(DCC.ccor)
}
