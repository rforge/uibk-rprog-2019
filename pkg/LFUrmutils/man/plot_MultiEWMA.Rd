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
  \item{...}{Optional Arguments to be passted.}
}
\details{
The generic function \code{plot} allows to display 14 graphs. These are the

Time Series, 
Conditional SD, 
Series with 1.96 Conditional SD Superimposed, 
ACF of Observations, 
ACF of Squared Observations, 
Cross Correlation Plot between (Squared) Returns,
Standardized Residuals, 
ACF of Standardized Residuals, 
ACF of Squared Standardized Residuals, 
Cross Correlation Plot between (Squared) Residuals,
Conditional Variance-Covariance Matrix,
Diagonal Elements of the Conditional Variance-Covariance Matrix, 
Conditional Correlations, 
Snapshot of the Model

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

\dontrun{
# Interactive plot
plot(EWMA)
}
}