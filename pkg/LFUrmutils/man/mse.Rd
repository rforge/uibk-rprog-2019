\name{mse}
\alias{mse}
\alias{mse.UnivVola}
\alias{mse.MultiEWMA}
\alias{mse.fGARCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Squared Error
}
\description{
Computes the mean squared error for various models.
}
\usage{
mse(object, ...)
}
\arguments{
  \item{object}{An object of class \code{UnivVola}, \code{MultiEWMA} or \code{fGARCH} (package \code{fGarch}).}
  \item{...}{Arguments passed to methods (unused at the moment).}
}
\value{
MSE of the corresponding model.
}
\references{
Jorion (2007). Value at Risk, 3rd. McGraw-Hill. New York.
}
\author{
Bernhard Eder
}
\examples{
# Load financial prices and compute returns
library("tseries")
library("zoo")
price <- get.hist.quote("db", "2010-01-01", "2019-07-08",
                     quote = "Adjusted", quiet = TRUE)
return <- diff(log(price))

# Compute univariate volatility model
ma <- UnivVola(return, type = "MovingAverage", center = TRUE)

# Compute MSE
mse(ma)
}
