"bcVAR" <- function(data = data, p = 1, constant = TRUE, ...){
  
  ## dimensions
  T = nrow(data)
  K = ncol(data)
  
  ## concise matrix notation of the reduced-form VAR model
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
  
  if (constant == FALSE){
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
  if (constant == FALSE){
    datamat <- datamat
  } else {
    datamat <- rbind(datamat, matrix(1, 1, T-p))
    rownames(datamat)[K*p1+1] <- "const"
  }
  
  ## OLS (reduced-form)
  Ahat <- lm.fit(t(X), t(Y))$coefficients
  Ahat <- t(Ahat)

  ## split up Ahat (constant, coeff)
  if (constant == FALSE){
    A.companion.hat <- Ahat[, 1:dim(Ahat)[2]]
  } else {
    V <- Ahat[1:K, 1]
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
  
  if (constant == FALSE){
    coeff <- bcLS
  } else {
    coeff <- cbind(bcLS, V)
  }
  
  ## save the results like the vars package (Pfaff)
  varresult <- list()
  for (i in 1:K){
    varresult[[i]] <- list()
    class(varresult[[i]]) <- "lm"
    varresult[[i]]$coefficients <- coeff[i, ]
    varresult[[i]]$residuals <- residuals[, i]
  }
  names(varresult) <- c(names(data))
  
  call <- match.call()
  
  ## return everything
  result <- list(
    varresult = varresult,
    datamat = t(datamat), 
    y = data,
    p = p, 
    K = K,
    obs = T,
    toobs = T - p,
    restrictions = NULL,
    call = call
  )
  
  class(result) <- "varest"
  return(result)
  
}





