## ## ## ##
## Generic plot function for univariate volatility models ----
## ## ## ##

## COPYRIGHT INFORMATION
## Parts of the following code were originally developed by Diethelm Wuertz.

plot.UnivVola <- function (x, which = "ask", ...){
  .interactiveUnivVolaPlot(
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
      "Cross Correlation Plot between (Squared) Residuals"
    ), 
    plotFUN = paste(".plot_UnivVola", 1:10, sep = "_"), 
    which = which, ...)
# Return value
invisible(x)
}


## ## ## ##
## Interactive plot function for univariate volatility models ----
## ## ## ##
  
.interactiveUnivVolaPlot <- function(object, choices, plotFUN, which, ...){
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
    .multUnivVolaPlot(object, choices, plotFUN, ...)
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
## Function to select plots for univariate volatility models ----
## ## ## ##  

.multUnivVolaPlot <-
  function (object, choices, ...)
  {
    
    pick = 1
    while (pick > 0) {
      pick = menu (
        choices = paste(" ", choices),
        title = "\nMake a plot selection (or 0 to exit):")
      # up to 19 plot functions ...
      switch(pick,
             .plot_UnivVola_1(object),
             .plot_UnivVola_2(object),
             .plot_UnivVola_3(object),
             .plot_UnivVola_4(object),
             .plot_UnivVola_5(object),
             .plot_UnivVola_6(object),
             .plot_UnivVola_7(object),
             .plot_UnivVola_8(object), 
             .plot_UnivVola_9(object), 
             .plot_UnivVola_10(object))
    }
  }


## ## ## ##
## Functions for specific plots of univariate volatility models ----
## ## ## ##

# Return series
.plot_UnivVola_1 <- function(object, ...){
  rets <- object$Returns
  plot(rets, col = "steelblue",
       main = "Time Series")
  abline(h = 0, col = "grey", lty = 3)
  grid()
}

# Conditional volatility
.plot_UnivVola_2 <- function(object, ...){
  condvola <- vola(object)
  plot(condvola, col = "steelblue", main = "Conditional SD")
  abline(h = 0, col = "grey", lty = 3)
  grid()
}

# Return with 1.96 Conditional SD superimposed
.plot_UnivVola_3 <- function(object, ...){
  rets <- object$Returns
  condvola <- vola(object)
  ci <- qnorm(0.025)
  plot(rets, type = "l", col = "steelblue", ylab = "x",
       main = "Series with 1.96 Conditional SD Superimposed")
  lines(mean(rets) + ci * condvola, col = "grey")
  lines(mean(rets) - ci * condvola, col = "grey")
  abline(h = 0, col = "grey", lty = 3)
  grid()
}

# ACF of Observations
.plot_UnivVola_4 <- function(object, ...){
  rets <- as.vector(object$Returns)
  n <- length(rets)
  lag.max <- as.integer(10*log10(n))
  acf(rets, lag.max = lag.max, xlab = "Lags", col = "steelblue",
      main = "ACF of Observations", plot = TRUE, 
      na.action = na.pass)
}


# ACF of Squared Observations
.plot_UnivVola_5 <- function(object, ...){
  rets <- as.vector(object$Returns)
  rets2 <- rets^2
  n <- length(rets)
  lag.max <- as.integer(10*log10(n))
  acf(rets2, lag.max = lag.max, xlab = "Lags", col = "steelblue",
      main = "ACF of Squared Observations", plot = TRUE, 
      na.action = na.pass)
}


# Cross correlation between squared returns and returns
.plot_UnivVola_6 <- function(object, ...){
  rets <- as.vector(object$Returns)
  rets2 <- rets^2
  n <- length(rets)
  lag.max <- as.integer(10*log10(n))
  ccf(rets2, rets, lag.max = lag.max, xlab = "Lags",
      main = "Cross Correlation", plot = TRUE, col = "steelblue", 
      na.action = na.pass)
}


# Standardized Residuals
.plot_UnivVola_7 <- function(object, ...){
  resids <- residuals(object)
  plot(resids, type = "l", main = "Standardized Residuals", col = "steelblue", ...)
  abline(h = 0, lty = 3)
  grid()
}


# ACF of Standardized Residuals
.plot_UnivVola_8 <- function(object, ...){
  resids <- residuals(object)
  n <- length(resids)
  lag.max <- as.integer(10*log10(n))
  acf(resids, lag.max = lag.max, xlab = "Lags", col = "steelblue",
      main = "ACF of Standardized Residuals", plot = TRUE, 
      na.action = na.pass)
}


# ACF of Squared Standardized Residuals
.plot_UnivVola_9 <- function(object, ...){
  resids <- residuals(object)
  resids2 <- resids^2
  n <- length(resids2)
  lag.max = as.integer(10*log10(n))
  acf(resids2, lag.max = lag.max, xlab = "Lags", col = "steelblue",
      main = "ACF of Squared Standardized Residuals", plot = TRUE, 
      na.action = na.pass)
}

# Cross correlation between squared returns and returns
.plot_UnivVola_10 <- function(object, ...){
  resids <- residuals(object)
  resids2 = resids^2
  n = length(resids)
  lag.max = as.integer(10*log10(n))
  ccf(resids2, resids, lag.max = lag.max, xlab = "Lags",
      main = "Cross Correlation", plot = TRUE, col = "steelblue", 
      na.action = na.pass)
}