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
\method{mse}{UnivVola}(object, ...)
\method{mse}{MultiEWMA}(object, ...)
\method{mse}{fGARCH}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{UnivVola}, \code{MultiEWMA} or \code{fGARCH}.}
  \item{...}{Arguments passed to methods (unused at the moment).}
}
\value{
MSE of the corresponding model.
}
\seealso{
\code{\link{UnivVola}}, \code{\link{MultiEWMA}}, \code{\link[fGarch:garchFit]{garchFit}}
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
ma <- UnivVola(y[, 1], type = "MovingAverage", center = TRUE)

# Compute MSE
mse(ma)
}
