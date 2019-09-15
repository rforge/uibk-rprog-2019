\name{imom}
\alias{imom}
\title{Implied Moments}

\description{Estimates ex-ante skewness and kurtosis from option implied volatilities over degrees of moneyness. Uses Gram-Chalier expansion.}

\usage{imom(data, moneyness=NULL)}

\arguments{
  \item{data}{Time series of implied volatilities with columns of different moneyness. \code{data} can be of any class, as long as it is transformable to a matrix when calling \code{as.matrix()}.}
  \item{moneyness}{Degrees of moneyness (in percent) of the data columns. If not specified, moneyness is tried to be imported from \code{data}'s column names.}
}

\details{\code{imom} Estimates implied base-volatility, skewness and kurtosis. Can be used as risk-neutral moments. Methods \code{summary}, \code{plot} and \code{print} are defined for objects of class \code{ivol}.}

\value{
\code{imom} returns a \code{list} object of \code{\link{class}} "\code{c(imom, ivol)}". List elements are as follows:
\item{coefficients}{A matrix of estimated implied volatility, skewness, kurtosis plus goodness of fit (r2).}
\item{residuals}{The residuals, that is response minus fitted values.}
\item{method}{Underlying estimation method. Currently only "\code{gramchalier}" available.}
\item{call}{The matched call.}
\item{fitted.values}{Regression's fitted values}
}

\seealso{
  \code{\link{ihurst}}
}

\examples{
# estimate implied Moments:

library(ivol)
data(otm)
x <- imom(otm)
summary(x)
plot(x)
coef(x)
}
\keyword{options, volatility}

