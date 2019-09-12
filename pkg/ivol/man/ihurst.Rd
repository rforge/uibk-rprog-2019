\name{ihurst}
\alias{ihurst}
\title{Implied Hurst Exponent}

\description{Estimates the ex-ante Hurst exponent and fractal volatility from option implied volatilities via a decomposition over the term structure. Builds on fractal option pricing of Hu and Oksendal (2003).}

\usage{ihurst(data, maturities=NULL)}

\arguments{
  \item{data}{Time series of implied volatilities with columns of different time-to-maturity. \code{data} can be of any class, as long as it is transformable to a matrix when calling \code{\link{as.matrix()}}.}
  \item{maturities}{Maturities (in days) of the data columns. If not specified, then maturities are tried to be imported from \code{data}'s column names.}
}

\details{\code{ihurst} Estimates the implied Hurst exponent, which can be used as risk-neutral momentum or for sentiment analysis. Methods \code{summary}, \code{plot} and \code{print} are defined for objects of class \code{ivol}.
}

\value{\code{ihurst} returns an \code{list} object of \code{\link{class}} "\code{c(ihurst, ivol)}". List elements are as follows:
\item{coefficients}{a \code{data.frame} of estimated Hurst exponents and fractal volatility.}
\item{residuals}{the residuals, that is response minus fitted values.}
\item{method}{underlying estimation method. Currently only "\code{h_lm}" available.}
\item{r.squared}{vector of goodness of fit.}
\item{call}{the matched call.}
}

\seealso{
  \code{\link{imom}}
}
\examples{
# estimate implied Hurst exponent

library(ivol)
data(atm)
x <- ihurst(atm)
plot(x)
x$h
}
\keyword{options, volatility}

