\name{residuals.UnivVola}
\alias{residuals.UnivVola}
\alias{residuals.MultiEWMA}
\title{
Residuals of volatility models
}
\description{
Function to compute the residuals of univariate or multivariate volatility models. Currently only implemented for objects of class \code{UnivVola} or \code{MultiEWMA}.
}
\usage{
\method{residuals}{UnivVola}(object, standardize = TRUE, na.action = "na.pass", ...)
\method{residuals}{MultiEWMA}(object, standardize = TRUE, na.action = "na.pass", ...)
}
\arguments{
  \item{object}{Object of class \code{UnivVola} or \code{MultiEWMA}.}
  \item{standardize}{Should standardized residuals be computed (\code{TRUE}, the default) or not (\code{FALSE}).}
  \item{na.action}{\code{NA} handling. See \link[=stats]{na.action} for further details.}
  \item{...}{Arguments passed to methods.}
}
\details{
By default standardized residuals are computed. For comparability reasons with the function \code{residuals()} from package \code{fGarch} non-standardized residuals can be computed as well. However, these correspond simply to the return series used to compute the conditional volatilities.
}
\value{
A \code{zoo} object containing the residuals for each return series.
}
\references{
Jorion (2007). Value at Risk, 3rd. McGraw-Hill. New York.
}
\author{
Bernhard Eder
}
\seealso{
\code{\link{UnivVola}}, \code{\link{MultiEWMA}}.
}
\examples{
# Load financial prices and compute returns
library("tseries")
library("zoo")
p1 <- get.hist.quote("msft", "2000-01-01", "2009-12-31",
                     quote = "Adjusted", quiet = TRUE)
p2 <- get.hist.quote("ibm", "2000-01-01", "2009-12-31",
                     quote = "Adjusted", quiet = TRUE)
y1 <- diff(log(p1)) * 100
y2 <- diff(log(p2)) * 100
y <- merge(y1, y2)

# Compute multivariate EWMA model (RiskMetrics)
EWMA <- MultiEWMA(y, center = TRUE)

# Compute and plot standardized residuals
EWMAres <- residuals(EWMA)
plot(EWMAres)
}

