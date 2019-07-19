## ## ## ##
## Generic plot function for multivariate EWMA models ----
## ## ## ##

plot.MultiEWMA <- function (x, which = "ask", ...){
  .interactiveMultiEWMAPlot(
    x, 
    choices = c(
      "Time Series", 
      "Conditional SD", 
      "Series with 1.96 Conditional SD Superimposed", 
      "ACF of Observations", 
      "ACF of Squared Observations", 
      "Cross Correlation Plot between (Squared) Returns",
      "Standardized Residuals", 
      "ACF of Standardized Residuals", 
      "ACF of Squared Standardized Residuals", 
      "Cross Correlation Plot between (Squared) Residuals",
      "Conditional Variance-Covariance Matrix",
      "Diagonal Elements of the Conditional Variance-Covariance Matrix", 
      "Conditional Correlations", 
      "Snapshot of the Model"
    ), 
    plotFUN = paste(".plot_MultiEWMA", 1:14, sep = "_"), 
    which = which, ...)
# Return value
invisible(x)
}


## ## ## ##
## Interactive plot function for multivariate EWMA models ----
## ## ## ##
  
.interactiveMultiEWMAPlot <- function(object, choices, plotFUN, which, ...){
  # Some checks
  if (length(choices) != length(plotFUN))
    stop("Arguments choices and plotFUN must be of same length.")
  if (length(which) > length(choices))
    stop("Arguments which has incorrect length.")
  if (length(which) > length(plotFUN))
    stop("Arguments which has incorrect length.")
  
  # Plot:
  if (is.numeric(which)) {
    Which = rep(FALSE, times = length(choices))
    Which[which] = TRUE
  }
  
  if (which[1] == "all") {
    Which = rep(TRUE, times = length(choices))
  }
  
  if (which[1] == "ask") {
    .multMultiEWMAPlot(object, choices, plotFUN, ...)
  } else {
    for ( i in 1:length(choices) ) {
      FUN = match.fun(plotFUN[i])
      if (Which[i]) FUN(object)
    }
  }
  
  # Return Value:
  invisible(object)
}


## ## ## ##
## Function to select plots for multivariate EWMA models ----
## ## ## ##  

.multMultiEWMAPlot <- function (object, choices, ...){
    pick = 1
    while (pick > 0) {
      pick = menu (
        choices = paste(" ", choices),
        title = "\nMake a plot selection (or 0 to exit):")
      # up to 19 plot functions ...
      switch(pick,
             .plot_MultiEWMA_1(object),
             .plot_MultiEWMA_2(object),
             .plot_MultiEWMA_3(object),
             .plot_MultiEWMA_4(object),
             .plot_MultiEWMA_5(object),
             .plot_MultiEWMA_6(object),
             .plot_MultiEWMA_7(object),
             .plot_MultiEWMA_8(object), 
             .plot_MultiEWMA_9(object), 
             .plot_MultiEWMA_10(object),
             .plot_MultiEWMA_11(object),
             .plot_MultiEWMA_12(object), 
             .plot_MultiEWMA_13(object), 
             .plot_MultiEWMA_14(object))
    }
}


## ## ## ##
## Functions for specific plots of multivariate EWMA models ----
## ## ## ##

# Return series
.plot_MultiEWMA_1 <- function(object, ...){
  rets <- object$Returns
  
  my.panel <- function(x, ...){
    lines(x, ...)
    abline(h = 0, col = "grey", lty = 3)
    grid()
  }
  
  plot(rets, col = "steelblue", main = "Time Series", 
       panel = my.panel)
}

# Conditional volatilities
.plot_MultiEWMA_2 <- function(object, ...){
  condvola <- vola(object)
  
  my.panel <- function(x, ...){
    lines(x, ...)
    grid()
  }
  
  plot(condvola, col = "steelblue", main = "Conditional SD", 
       panel = my.panel)
}

# Return with 1.96 Conditional SD superimposed
.plot_MultiEWMA_3 <- function(object, ...){
  ci <- qnorm(0.025)
  rets <- object$Returns
  condvola <- vola(object)
  
  my.panel <- function(x, ...){
    lines(x, ...)
    panel.number <- parent.frame()$panel.number
    lines(+ ci * condvola[, panel.number], col = "grey")
    lines(- ci * condvola[, panel.number], col = "grey")
    abline(h = 0, col = "grey", lty = 3)
    grid()
  }
  
  plot(rets, panel = my.panel, type = "l", col = "steelblue", 
       main = "Series with 1.96 Conditional SD Superimposed")
}


# ACF of Observations
.plot_MultiEWMA_4 <- function(object, ...){
  rets <- object$Returns
  n <- length(rets)
  cols <- dim(rets)[2]
  lag.max <- as.integer(10*log10(n))
  par(mfrow = rev(n2mfrow(cols)))
  for (i in 1:cols){
  acf(rets[, i], lag.max = lag.max, xlab = "Lags", col = "steelblue",
      main = "ACF of Observations", plot = TRUE, 
      na.action = na.pass)
  }
  par(mfrow = c(1, 1))
}


# ACF of Squared Observations
.plot_MultiEWMA_5 <- function(object, ...){
  rets <- object$Returns
  rets2 <- rets^2
  n <- length(rets2)
  cols <- dim(rets)[2]
  lag.max <- as.integer(10*log10(n))
  par(mfrow = rev(n2mfrow(cols)))
  for (i in 1:cols){
  acf(rets2[, i], lag.max = lag.max, xlab = "Lags", col = "steelblue",
      main = "ACF of Squared Observations", plot = TRUE, 
      na.action = na.pass)
  }
  par(mfrow = c(1, 1))
}


# Cross correlation between squared returns and returns
.plot_MultiEWMA_6 <- function(object, ...){
  rets <- object$Returns
  rets2 <- rets^2
  n <- length(rets)
  cols <- dim(rets)[2]
  lag.max <- as.integer(10*log10(n))
  par(mfrow = rev(n2mfrow(cols)))
  for(i in 1:cols){
  ccf(rets2[, i], rets[, i], lag.max = lag.max, xlab = "Lags",
      main = "Cross Correlation", plot = TRUE, col = "steelblue", 
      na.action = na.pass)
  }
  par(mfrow = c(1, 1))
}


# Standardized Residuals
.plot_MultiEWMA_7 <- function(object, ...){
  resids <- residuals(object)
  
  my.panel <- function(x, ...){
    lines(x, ...)
    abline(h = 0, col = "grey", lty = 3)
    grid()
  }
  
  plot(resids, type = "l", main = "Standardized Residuals", col = "steelblue", panel = my.panel, ...)
  #abline(h = 0, lty = 3)
}


# ACF of Standardized Residuals
.plot_MultiEWMA_8 <- function(object, ...){
  resids <- residuals(object)
  n <- length(resids)
  cols <- dim(resids)[2]
  lag.max <- as.integer(10*log10(n))
  par(mfrow = rev(n2mfrow(cols)))
  for(i in 1:cols){
  acf(resids[, 1], lag.max = lag.max, xlab = "Lags", col = "steelblue",
      main = "ACF of Standardized Residuals", plot = TRUE, 
      na.action = na.pass)
  }
}


# ACF of Squared Standardized Residuals
.plot_MultiEWMA_9 <- function(object, ...){
  resids <- residuals(object)
  resids2 <- resids^2
  n <- length(resids2)
  cols <- dim(resids2)[2]
  lag.max <- as.integer(10*log10(n))
  par(mfrow = rev(n2mfrow(cols)))
  for(i in 1:cols){
  acf(resids2[, i], lag.max = lag.max, xlab = "Lags", col = "steelblue",
      main = "ACF of Squared Standardized Residuals", plot = TRUE, 
      na.action = na.pass)
  }
}

# Cross correlation between squared returns and returns
.plot_MultiEWMA_10 <- function(object, ...){
  resids <- residuals(object)
  resids2 = resids^2
  n = length(resids)
  cols <- dim(resids2)[2]
  lag.max <- as.integer(10*log10(n))
  par(mfrow = rev(n2mfrow(cols)))
  for(i in 1:cols){
  ccf(resids2[, i], resids[, i], lag.max = lag.max, xlab = "Lags",
      main = "Cross Correlation", plot = TRUE, col = "steelblue", 
      na.action = na.pass)
  }
}


# Conditional variance-covariance matrix
.plot_MultiEWMA_11 <- function(object, ...){
  variances <- fitted(object)
  
  my.panel <- function(x, ...){
    lines(x, ...)
    grid()
  }
  
  plot(variances, col = "steelblue", main = "Conditional variance-covariance matrix", 
       panel = my.panel)
}

# Diagonal elements of the conditional variance-covariance matrix
.plot_MultiEWMA_12 <- function(object, ...){
  diagonal <- varcov(object, offdiagonal = FALSE)
  
  my.panel <- function(x, ...){
    lines(x, ...)
    grid()
  }
  
  plot(diagonal, col = "steelblue", main = "Diagonal elements of the\n conditional variance-covariance matrix", 
       panel = my.panel)
}

# Conditional correlations
.plot_MultiEWMA_13 <- function(object, ...){
  condcor <- ccor(object, diagonal = FALSE, duplicates = FALSE)
  
  my.panel <- function(x, ...){
    lines(x, ...)
    grid()
  }
  
  plot(condcor, col = "steelblue", main = "Conditional correlations", 
       panel = my.panel)
}

# Snapshot
.plot_MultiEWMA_14 <- function(object, ...){
  volat <- vola(object, offdiagonal = FALSE)
  varcovs <- fitted(object)
  condcor <- ccor(object, diagonal = FALSE, duplicates = FALSE)
  
  par(mfrow = c(2,2))
  plot(volat[, 1], ylab = "", main = expression(hat(sigma)["1,t"]))
  plot(varcovs[, 2], ylab = "", main = expression(hat(sigma)["12,t"]))
  plot(condcor[, 1], ylab = "", main = expression(hat(rho)["12,t"]))
  plot(volat[, 2], ylab = "", main = expression(hat(sigma)["2,t"]))
  par(mfrow = c(1,1))
}