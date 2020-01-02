\name{plot.MultiEWMA}
\alias{plot.MultiEWMA}
\title{
Plot multivariate EWMA models
}
\description{
Plot methods for multiviarate EWMA models models.
}
\usage{
\method{plot}{MultiEWMA}(x, which = "ask", ...)
}
\arguments{
  \item{x}{An object of class "\code{MultiEWMA}"}
  \item{which}{A character string denoting which plot should be displayed.}
  \item{...}{Optional Arguments to be passed.}
}
\details{
The generic function \code{plot} allows to display 14 graphs. These are the \cr

Time Series, \cr
Conditional SD, \cr
Series with 1.96 Conditional SD Superimposed, \cr
ACF of Observations, \cr
ACF of Squared Observations, \cr
Cross Correlation Plot between (Squared) Returns, \cr
Standardized Residuals, \cr
ACF of Standardized Residuals, \cr
ACF of Squared Standardized Residuals, \cr
Cross Correlation Plot between (Squared) Residuals, \cr
Conditional Variance-Covariance Matrix, \cr
Diagonal Elements of the Conditional Variance-Covariance Matrix,\cr
Conditional Correlations, \cr
Snapshot of the Model\cr
}

\note{
The cross correlation plots return only the corresponding plots for the first and second firm.

The "Snapshot of the Model" plots the conditional volatilities of the first two return series, as well as their conditional covariance, and their conditional correlation. See Ruppert and Matteson (2015, p. 428 -- 443) for further details.
}

\seealso{
\code{\link{MultiEWMA}}
}
\author{
Bernhard Eder
}
\references{
Ruppert and Matteson (2015). Statistics and Data Analysis for Financial Engineering, 2nd. Springer. New York.
}
\examples{
# Load returns
data("CRSPday_zoo")
y <- CRSPday_zoo

# Compute multivariate volatility model
EWMA <- MultiEWMA(y)

# Plot conditional standard deviations
plot(EWMA, which = 2)

\dontrun{
# Interactive plot
plot(EWMA)
}
}
\keyword{hplot}