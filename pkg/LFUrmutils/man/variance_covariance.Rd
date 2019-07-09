\name{varcov}
\alias{varcov}
\alias{varcov.MultiEWMA}
\title{
Variance-covariance matrix
}
\description{
Function to extract the variance-covariance matrix of multivariate volatility models. Currently only implemented for models of class \code{MultiEWMA}.
}
\usage{
varcov(object, ...)
\method{varcov}{MultiEWMA}(object, offdiagonal = TRUE, duplicates = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Output of a multivariate volatility model.}
  \item{offdiagonal}{Whether diagonal elements should be kept or not. \code{TRUE} (the default) or \code{FALSE}.}
  \item{duplicates}{Whether duplicate off-diagonal elements should be kept or not. \code{TRUE} (the default) or \code{FALSE}.}
  \item{...}{Arguments passed to methods.}
}
\details{
Returns the variance-covariance matrix of a multivariate volatility model. 
}
\value{
A multivariate \code{zoo} object. Each row of the output objects represents the correlation matrix (or elements thereof) of the corresponding day.
}
\references{
Jorion (2007). Value at Risk, 3rd. McGraw-Hill. New York.
}
\author{
Bernhard Eder
}
\seealso{
\code{\link{MultiEWMA}}, \code{\link{ccor}}
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

# Compute multivariate EWMA model
ewma <- MultiEWMA(y, center = TRUE)

# Extract variance covariance matrix
EWMAvc <- varcov(ewma)
plot(EWMAvc)
}

