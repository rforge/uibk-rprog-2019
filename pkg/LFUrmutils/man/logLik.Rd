\name{logLik}
\alias{logLik}
\alias{logLik.UnivVola}
\alias{logLik.MultiEWMA}
\alias{logLik.fGARCH}
\title{
Extract Log-Likelihood
}
\description{
These functions are based on the generic \code{logLik} method for objects of classes \code{UnivVola}, \code{MultiEWMA}, and \code{fGARCH}.
}
\usage{
\method{logLik}{UnivVola}(object, ...)
\method{logLik}{MultiEWMA}(object, ...)
\method{logLik}{fGARCH}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{UnivVola}, \code{MultiEWMA} or \code{fGARCH}.}
  \item{...}{Arguments passed to methods (unused at the moment).}
}
\value{
Returns an object of class \code{logLik}.
}
\seealso{
\code{\link{UnivVola}}, \code{\link{MultiEWMA}}, \code{\link[fGarch:garchFit]{garchFit}}, \code{\link[fGarch:fGARCH-class]{fGARCH-class}}
}
\references{
Jorion (2007). Value at Risk, 3rd. McGraw-Hill. New York.
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

# Extract log likelihood
logLik(RM)

# Compute AIC
AIC(RM)
}
\concept{log likelihood}
\concept{AIC}