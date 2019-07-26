"bcVAR" <- function(data = data, p = 1, type = c("const", "none"), ...){
  
  ## dimensions
  T = nrow(data)
  K = ncol(data)
  
  ## concise matrix notation of the reduced-form VAR model
  y.orig <- as.matrix(data)
  y <- t(as.matrix(data))
  
  if (any(is.na(y))) 
    stop("\nNAs in y.\n")

  Y <- y[, p:T]
  
  if (p < 2) {
    Y <- Y
  } else {
    for (i in 1:(p-1)){
      Y <- rbind(Y, y[, (p-i):(T-i)])
    }
  }
  
  if (type == "none"){
    X <- Y[, 1:(T-p)]
  } else {
    X <- rbind(matrix(1, 1, T-p), Y[, 1:(T-p)])
  }
  Y <- Y[, 2:(T-p+1)]
  
  ## for datamat (Pfaff)
  p1 = p + 1
  datamat <- y[, p1:T]
  for (i in 1:(p1-1)){
    datamat <- rbind(datamat, y[, (p1-i):(T-i)])
  } 
  if (type == "none"){
    datamat <- datamat
  } else {
    datamat <- rbind(datamat, matrix(1, 1, T-p))
    rownames(datamat)[K*p1+1] <- "const"
  }
  datamat <- as.data.frame(t(datamat))
  
  ## OLS (reduced-form)
  # Ahat <- Y%*%t(X) %*% solve(X%*%t(X))
  Ahat <- lm.fit(t(X), t(Y))$coefficients
  Ahat <- t(Ahat)

  ## split up Ahat (constant, coeff)
  if (type == "none"){
    A.companion.hat <- Ahat[, 1:dim(Ahat)[2]]
  } else {
    const <- Ahat[1:K, 1]
    A.companion.hat <- Ahat[, 2:dim(Ahat)[2]]
  }

  ## estiamted reduced-form residuals
  Uhat <- Y - Ahat%*%X
  residuals <- t(Uhat[1:K, ])
  
  ## cov-var matrix reduced-form residuals
  SIGMA <- Uhat%*%t(Uhat)*(T-p-K*p-1)^(-1)
  SIGMA <- SIGMA[1:K, 1:K]
  
  ## Bias-corrected LS ##
  tA <- t(A.companion.hat)
  SIGMAU <- Uhat%*%t(Uhat)*(T-p-K*p-1)^(-1)
  I <- diag(K*p)
  
  vecSIGMAY <- solve(diag((K*p)^2) - kronecker(A.companion.hat, A.companion.hat)) %*% vec(SIGMAU)
  dim(vecSIGMAY) <- c(K*p, K*p)          
  GAMMAY <- vecSIGMAY
  
  e <- eigen(A.companion.hat)
  
  list <- list()
  for (i  in 1:(K*p)){
    list[[i]] <- e$values[i]*solve(I-e$values[i]*tA)
  }
  sumeigen <- Re(as.matrix(Reduce('+', list))) 
  
  BA <- SIGMAU%*%(solve(I-tA) + tA%*%solve(I-tA%*%tA) + sumeigen)%*%solve(GAMMAY) 
  
  correction <- BA/(T-p)
  bcLS <- A.companion.hat + correction
  
  if (type == "none"){
    coeff <- bcLS
  } else {
    coeff <- cbind(bcLS, const)
  }
  
  ## rename the lags
  temp1 <- NULL
  for (i in 1:p) {
    temp <- paste(colnames(data), ".l", i, sep = "")
    temp1 <- c(temp1, temp)
  }
  names <- c(temp1, "const")
  names(datamat)[(K+1):(K + K*p + 1)] <- names

  ## save the results like the vars package (Pfaff)
  varresult <- list()
  for (i in 1:K){
    varresult[[i]] <- list()
    class(varresult[[i]]) <- "lm"
    varresult[[i]]$coefficients <- coeff[i, ]
    names(varresult[[i]]$coefficients) <- names
    varresult[[i]]$residuals <- residuals[, i]
  }
  names(varresult) <- c(colnames(data))
  
  call <- match.call()
  
  ## return everything
  result <- list(
    varresult = varresult,
    datamat = datamat, 
    y = y.orig,
    type = type, 
    p = p, 
    K = K,
    obs = T - p,
    totobs = T,
    restrictions = NULL,
    call = call
  )
  
  class(result) <- "varest"
  return(result)
  
}





