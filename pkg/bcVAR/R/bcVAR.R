"bcVAR" <- function(data = data, p = 1, constant = TRUE, ...){
  
  ## dimensions
  T = nrow(data)
  K = ncol(data)
  
  ## concise matrix notation of the reduced-form VAR model
  y <- t(as.matrix(data))
  
  if (any(is.na(y))) 
    stop("\nNAs in y.\n")

  Y <- y[, p:T]
  
  if (p < 2){
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
  
  ## cov-var matrix reduced-form residuals
  SIGMA <- Uhat%*%t(Uhat)*(T-p-K*p-1)^(-1)
  SIGMA <- SIGMA[1:K, 1:K]
  
  ## Bias-corrected LS ##
  tA <- t(A.companion.hat)
  SIGMAU <- Uhat%*%t(Uhat)*(T-p-K*p-1)^(-1)
  I <- diag(K*p)
  
  vecSIGMAY <- solve(diag((K*p)^2) - kronecker(A.companion.hat, A.companion.hat)) %*% matrixcalc::vec(SIGMAU)
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
  
  ## return everything
  result <- list(
    "V" = if (constant == FALSE){
      "V" = c("no constant")
    } else {
      "V" = V
    },
    "A" = unname(A.companion.hat),
    "bcLS" = unname(bcLS),
    "SIGMA" = SIGMA,
    "Y" = Y,
    "Uhat" = Uhat,
    "sumstat" = c(
      "T" = T,
      "K" = K,
      "p" = p
    )
  )
  class(result) <- "varest"
  return(result)
}




