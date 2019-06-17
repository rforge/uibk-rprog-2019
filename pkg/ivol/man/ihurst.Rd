\name{ihurst}
\alias{ihurst}
\title{Implied Hurst Exponent}
\description{Estimates the ex-ante Hurst exponent from option implied volatilities via a decomposition over the term structure. Builds on fractal option pricing of Hu and Oksendal (2003)}
\usage{ihurst(data, maturities)}
\arguments{
  \item{data}{Time series of implied volatilities with columns of different time-to-maturity. Either in matrix, \code{data.frame} or \code{xts} format.}
  \item{maturities}{Maturities (in days) of the data columns.}
}
\details{\code{ihurst} Estimates the implied Hurst exponent. Can be used as risk-neutral momentum or sentiment analysis.}
\value{Implied Hurst exponent as time series of class \code{ivol}}
\seealso{
  \code{\link{imom}}
}
\examples{
## estimate implied Hurst exponent FIXME
## ih <- ihurst(data = implied_volatilities, maturities)
## summary(ih)
}
\keyword{options, volatility}

