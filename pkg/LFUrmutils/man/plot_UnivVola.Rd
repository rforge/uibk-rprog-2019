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
The generic function \code{plot} allows to display 10 graphs. These are the \cr

Time Series Plot, \cr
Conditional Standard Deviation (Volatility) Plot, \cr
Series Plot with 1.96 Conditional SD Superimposed, \cr
Autocorrelation Function Plot of Observations, \cr
Autocorrelation Function Plot of Squared Observations, \cr
Cross Correlation Plot between (Squared) Observations, \cr
Standardized Residuals Plot, \cr
ACF Plot of Standardized Residuals, \cr
ACF Plot of Squared Standardized Residuals, \cr
Cross Correlation Plot between (Standardized Residuals), \cr
QQ-Plot of Standardized Residuals, \cr
Time Series Plot of the Simple Errors, \cr
Time Series Plot of the Squared Errors \cr
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
\keyword{plot}