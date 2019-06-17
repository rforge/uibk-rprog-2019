imom <- function(data, use_names=TRUE, moneyness){
  mn <- moneyness/100 - 1
  df <- as.matrix(data)
  ret <- list()
  r2 <- r <- c()
  nn <- nrow(df)
  mom <- matrix(NA, nn, 3)
  colnames(mom) <- c("sigma", "skew", "kurt")
  x <- cbind(1, mn, mn^2)
  for(i in 1:nn){
    y <- df[i,]
    mdl <- lm.fit(x, y)
    f <- mdl$fitted.values
    m <- mean(y)
    r2[i] <- sum((f - m)^2) / sum((y -m)^2)
    mom[i,] <- mdl$coefficients
    # hier kann evt runtime verbessert werden:
    r <- c(r, mdl$residuals)
  }
  ret$r.squared <- r2
  ret$mom <- mom
  ret$residuals <- r
  ret$method <- "mom_gramchalier"
  ret$last_mdl <- mdl$fitted.values
  class(ret) <- "ivol"
  return(ret)
}

ihurst <- function(data, maturities){
  tau <- maturities
  df <- log(as.matrix(data))
  ret <- list()
  h <- sigma <- r2 <- r <- c()
  x <- cbind(1, log(tau))
  df <- as.matrix(df)
  for(i in 1:nrow(df)){
    y <- df[i,]
    mdl <- lm.fit(x, y)
    f <- mdl$fitted.values
    m <- mean(y)
    r2[i] <- sum((f - m)^2) / sum((y -m)^2)
    h[i] <- mdl$coefficients[2] + 0.5
    sigma[i] <- exp(mdl$coefficients[1])
    # hier kann evt runtime verbessert werden:
    r <- c(r, mdl$residuals)
  }
  ret$r.squared <- r2
  ret$h <- h
  ret$sigma_f <- sigma
  ret$residuals <- r
  ret$method <- "h_lm"
  ret$last_mdl <- mdl$fitted.values
  class(ret) <- "ivol"
  return(ret)
}

summary.ivol <- function(object){
  z <- object
  r2 <- z$r.squared
  r <- z$residuals
  met <- z$method
  z$r.squared <- z$residuals <- z$method <- NULL
  qq <- quantile(r2, c(0, 0.25, 0.5, 0.75, 1))
  rr <- quantile(r, c(0, 0.25, 0.5, 0.75, 1))
  names(qq) <- names(rr) <- c("Min", "1Q", "Median", "3Q", "Max")
  m <- sapply(z, mean)
  s <- sapply(z, sd)
  ss <- round(cbind(m, s), 3)
  colnames(ss) <- c("avg.", "sd.")
  
  if(met == "h_lm"){
    cat("\nCall:\nimplied Hurst time series: lm approach\n\nCoefficients: \n")
  } else {
    cat("\nCall:\nimplied Moments time series: quadratic approach\n\nCoefficients: \n")
  }
  print(ss)
  cat("\nResiduals:\n")
  print(round(rr, 3))
  cat("\nR squared: \n")
  round(qq, 3)
}



plot.ivol <- function(object){
  z <- object
  if(z$method == "h_lm"){
    par(mfrow=c(2, 1))
    plot(z$h, type='l', ylab = "implied Hurst", xlab = "time", ylim=c(0,1))
    abline(0.5, 0)
    plot(z$sigma_f, type='l', ylab = "fractal Volatility", xlab = "time", ylim=c(0, max(z$sigma_f)*1.05))
    par(mfrow=c(1, 1))
  }
}


predict.ivol <- function(object, newdata, horizon){
  z <- object
  if(missing(horizon) || is.null(horizon)){
    hz <- names(z$last_mdl)
    hz <- as.numeric(gsub("T", "", hz))
  } else {
    hz <- horizon
  }
  
  if(missing(newdata) || is.null(newdata)){
    y <- exp(z$last_mdl[1]) 
    h <- z$h[length(z$h)]
  } else {
    y <-  newdata[1]
    h <- newdata[2]
  }
  pr <- y * hz^(h-0.5)
  names(pr) <- hz
  pr
}

