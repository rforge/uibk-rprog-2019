imom <- function(data, moneyness=NULL){
  if(is.null(moneyness)) moneyness <- colnames(data)
  if(is.character(moneyness)){
    moneyness <- gsub("[^0-9.]", "",  moneyness)
    moneyness <- as.numeric(moneyness)
  } 
  if(!is.numeric(moneyness)) stop("unknown specification of 'moneyness'")
  mn <- moneyness
  if(min(mn) > 10)  mn <- mn/100 
  mn <- mn - 1
  
  t <- NULL
  if(inherits(data, c("xts", "zoo"))) t <- index(data)
  
  df <- as.matrix(data)
  ret <- list()
  x <- cbind(1, mn, mn^2)
  mdl <- lm.fit(x, t(df))
  co <- t(mdl$coefficients)
  colnames(co) <- c("iVol", "iSkew", "iKurt")
  f <- t(mdl$fitted.values)
  m <- colMeans(df)
  ss_tot <- rowSums((df-m)^2)
  ss_res <- rowSums((f - df)^2)
  r2 <- 1 - ss_res/ss_tot
  co <- cbind(co, r2)

  ret$coefficients <- co
  ret$residuals <- t(mdl$residuals)
  ret$fitted.values <- mdl$fitted.values
  if(!is.null(t) & is.zoo(data)) ret <- lapply(ret, function(x) zoo(x, t))
  if(!is.null(t) & is.xts(data)) ret <- lapply(ret, function(x) xts(x, t))
  ret$method <- "gramchalier"
  ret$call <- match.call()
  class(ret) <- c("imom", "ivol")
  return(ret)
}

ihurst <- function(data, maturities=NULL){
  if(is.null(maturities)) maturities <- colnames(data)
  if(is.character(maturities)){ 
    maturities <- gsub("[^0-9.]", "",  maturities)
    maturities <- as.numeric(maturities)
  } 
  if(!is.numeric(maturities)) stop("unknown specification of 'maturities'")
  tau <- maturities
  
  t <- NULL
  if(inherits(data, c("xts", "zoo"))) t <- index(data)
  
  df <- log(as.matrix(data))
  ret <- list()
  x <- cbind(1, log(tau))
  mdl <- lm.fit(x, t(df))
  co <-t(mdl$coefficients)
  colnames(co) <- c("fVola", "iHurst")
  co[,1] <- exp(co[,1])
  co[,2] <- 0.5 + co[,2]
  f <- t(mdl$fitted.values)
  res <- t(mdl$residuals)
  m <- colMeans(df)
  ss_tot <- rowSums((df-m)^2)
  ss_res <- rowSums((f - df)^2)
  r2 <- 1 - ss_res/ss_tot
  co <- cbind(co, r2)
  
  ret$coefficients <- co
  ret$residuals <- res
  ret$fitted.values <- f
  if(!is.null(t) && is.zoo(data)) ret <- lapply(ret, function(x) zoo(x, t))
  if(!is.null(t) && is.xts(data)) ret <- lapply(ret, function(x) xts(x, t))
  ret$method <- "hu_oskendal"
  ret$call <- match.call()
  class(ret) <-  c("ihurst", "ivol")
  return(ret)
}

summary.ivol <- function(object, ...){
  qq <- quantile(object$coefficients[,"r2"], c(0, 0.25, 0.5, 0.75, 1))
  rr <- quantile(object$residuals, c(0, 0.25, 0.5, 0.75, 1))
  names(qq) <- names(rr) <- c("Min", "1Q", "Median", "3Q", "Max")
  m <- colMeans(object$coefficients)
  s <- apply(object$coefficients, 2, sd)
  ss <- round(cbind(m, s), 3)
  colnames(ss) <- c("avg.", "sd.")
  
  if(class(object)[1] == "ihurst") k <- "Hurst" else k <- "Moments"
  cat(paste("\nCall:\nimplied", k, "time series | Approach:", object$method, "\n\nCoefficients: \n"))
  print(ss)
  cat("\nResiduals:\n")
  print(round(rr, 3))
  cat("\nR squared: \n")
  round(qq, 3)
}

plot.ivol <- function(x, time_format = "%Y-%m-%d", ...){
  df <- coef(x)
  if(inherits(df, c("xts", "zoo"))){
    t <- index(df)
  } else {
    t <- try(as.Date(rownames(df), time_format), silent = TRUE)
    if(inherits(t, "try-error")) t <- 1:nrow(df)
  }
  
  n <- if(inherits(x, "ihurst")) 3 else 4
  ask <- prod(par("mfcol")) < n && dev.interactive()
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  
  dict <- c("implied Hurst", "fractal Volatility", "Goodness of Fit", "implied Vola", "implied Skew", "implied Kurt")
  names(dict) <- c("iHurst", "fVola", "r2", "iVol", "iSkew", "iKurt")
  nm <- colnames(df)
  
  for(i in nm){
    if(i %in% c("iHurst", "r2")) ylim <- c(0,1) else ylim <- NULL
    plot(t, df[,i], type = 'l', xlab="time", ylab = i, ylim = ylim, main = dict[i])
    if(i == "iHurst") abline(0.5, 0)
  }
  cat(paste("\ngenerated", n, "plots\n"))
}

print.ivol <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)
  cat("\nMethod:\n")
  print(x$method)
  cat("\nLast Coefficients:\n")
  print(tail(coef(x)))
  invisible(x)
}
