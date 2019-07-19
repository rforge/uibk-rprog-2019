\name{fitted.UnivVola}
\alias{fitted.UnivVola}
\alias{fitted.MultiEWMA}
\title{
Extract Model Fitted Values
}
\description{
These functions are based on the generic \code{fitted} method for objects of classes \code{UnivVola} and \code{MultiEWMA}, and extract the conditional variances.
}
\usage{
\method{fitted}{UnivVola}(object, ...)
\method{fitted}{MultiEWMA}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{UnivVola} or \code{MultiEWMA}.}
  \item{...}{Arguments passed to methods (unused at the moment).}
}
\value{
For univariate volatility models, the function extracts the conditional variances. For the multivariate EWMA model, the function extracts the conditional variance-covariance matrices. In this case, each row in the output corresponds to the complete conditional variance-vovariance matrix of the corresponding day. 
}
\note{
For the multivariate EWMA model, we recommend that end users use the function \code{\link{varcov}}, which also allows to extract the conditional variance-covariance matrix, but provides additional arguments to control the output.
}
\seealso{
\code{\link{UnivVola}}, \code{\link{MultiEWMA}}, \code{\link{varcov}}.
}
\references{
Danielsson (2011). Financial Risk Forecasting. Wiley. Chichester.
Jorion (2007). Value at Risk, 3rd. McGraw-Hill. New York.
Ruppert and Matteson (2015). Statistics and Data Analysis for Financial Engineering, 2nd. Springer. New York.
}
\author{
Bernhard Eder
}
\examples{
# Load returns
data("CRSPday_zoo")
y <- CRSPday_zoo

# Compute univariate volatility model
RM <- UnivVola(y[, 1], type = "RiskMetrics", lambda = -1, center = TRUE)

# Extract fitted values
cond_vola <- fitted(RM)
}
