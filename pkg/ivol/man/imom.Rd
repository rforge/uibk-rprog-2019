\name{imom}
\alias{imom}
\title{Implied Moments}
\description{Estimates ex-ante skewness and kurtosis from option implied volatilities over degrees of moneyness. Uses Gram-Chalier (19..) expansion.}
\usage{imom(data, moneyness)}
\arguments{
  \item{data}{Time series of implied volatilities with columns of different moneyness. Either in matrix, \code{data.frame} or \code{xts} format.}
  \item{moneyness}{Moneyness (in percent) of the data columns.}
}
\details{\code{imom} Estimates implied base-volatility, skewness and kurtosis. Can be used as risk-neutral moments.}
\value{Implied Volatility, Skewness and Kurtosis as time series of class \code{ivol}}
\seealso{
  \code{\link{ihurst}}
}
\examples{
## estimate implied Hurst exponent
im <- ihurst(data = implied_volatilities, moneyness)
summary(im)

}
\keyword{options, volatility}

