\name{imom}
\alias{imom}
\title{Implied Moments}
\description{Estimates ex-ante skewness and kurtosis from option implied volatilities over degrees of moneyness. Uses Gram-Chalier expansion.}
\usage{imom(data, use_names=TRUE, moneyness=NA)}
\arguments{
  \item{data}{Time series of implied volatilities with columns of different moneyness. Either in matrix, \code{data.frame} or \code{xts} format.}
  \item{use_names}{Logical. If set to \code{TRUE} then moneyness is imported from column names of the data and values passed to \code{moneyness} are ignored.}
  \item{moneyness}{Degrees of moneyness (in percent) of the data columns.}
}
\details{\code{imom} Estimates implied base-volatility, skewness and kurtosis. Can be used as risk-neutral moments.}
\value{Implied Volatility, Skewness and Kurtosis as time series of class \code{ivol}}
\seealso{
  \code{\link{ihurst}}
}
\examples{
## estimate implied Moments:
# library(ivol)
# data(otm)
# x <- imom(otm)
# plot(x)
# x$mom
}
\keyword{options, volatility}

