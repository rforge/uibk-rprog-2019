\name{plot.UnivVola}
\alias{plot.UnivVola}
\title{
Plot univariate volatility models
}
\description{
Plot methods for univariate volatility models.
}
\usage{
\method{plot}{UnivVola}(x, which = "ask", ...)
}
\arguments{
  \item{x}{An object of class "\code{UnivVola}"}
  \item{which}{A character string denoting which plot should be displayed.}
  \item{...}{Optional Arguments to be passted.}
}
\details{
The generic function \code{plot} allows to display 10 graphs. These are the

Time Series Plot
Conditional Standard Deviation (Volatility) Plot
Series Plot with 1.96 Conditional SD Superimposed
Autocorrelation Function Plot of Observations
Autocorrelation Function Plot of Squared Observations
Cross Correlation Plot between (Squared) Observations
Standardized Residuals Plot
ACF Plot of Standardized Residuals
ACF Plot of Squared Standardized Residuals
Cross Correlation Plot between (Standardized Residuals)
}
\seealso{
\code{\link{UnivVola}}
}
\author{
Bernhard Eder
}
\examples{
# Load returns
data("CRSPday_zoo")
y <- CRSPday_zoo

# Compute univariate volatility model
ma <- UnivVola(y[, 1], type = "MovingAverage", center = TRUE)

# Plot volatility
plot(ma, which = 2)

\dontrun{
# Interactive plot
plot(ma)
}
}