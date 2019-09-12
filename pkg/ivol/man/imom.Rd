\name{imom}
\alias{imom}
\title{Implied Moments}

\description{Estimates ex-ante skewness and kurtosis from option implied volatilities over degrees of moneyness. Uses Gram-Chalier expansion.}

\usage{imom(data, moneyness=NULL)}

\arguments{
  \item{data}{Time series of implied volatilities with columns of different moneyness. \code{data} can be of any class, as long as it is transformable to a matrix when calling \code{\link{as.matrix()}}.}
  \item{moneyness}{Degrees of moneyness (in percent) of the data columns. If not specified, moneyness is tried to be imported from \code{data}'s column names.}
}

\details{\code{imom} Estimates implied base-volatility, skewness and kurtosis. Can be used as risk-neutral moments. Methods \code{summary}, \code{plot} and \code{print} are defined for objects of class \code{ivol}.}

\value{
\code{imom} returns an \code{list} object of \code{\link{class}} "\code{c(imom, ivol)}". List elements are as follows:
\item{coefficients}{a \code{data.frame} of estimated implied volatility, skewness and kurtosis.}
\item{residuals}{the residuals, that is response minus fitted values.}
\item{method}{underlying estimation method. Currently only "\code{mom_gramchalier}" available.}
\item{r.squared}{vector of goodness of fit.}
\item{call}{the matched call.}}

\seealso{
  \code{\link{ihurst}}
}

\examples{
# estimate implied Moments:

library(ivol)
data(otm)
x <- imom(otm)
plot(x)
x$mom
}
\keyword{options, volatility}

