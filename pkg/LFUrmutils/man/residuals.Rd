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
By default standardized residuals are computed. For comparability reasons with the function \code{residuals()} from package \code{fGarch}, non-standardized residuals can be computed as well. Depending on the data used to compute the volatility models, the non-standardized residuals will either correspond to the return series or the residuals of the conditional mean model.
}
\value{
A \code{zoo} object containing the residuals for each return series.
}
\references{
Danielsson (2011). Financial Risk Forecasting. Wiley. Chichester.
Jorion (2007). Value at Risk, 3rd. McGraw-Hill. New York.
Ruppert and Matteson (2015). Statistics and Data Analysis for Financial Engineering, 2nd. Springer. New York.
}
\author{
Bernhard Eder
}
\seealso{
\code{\link{UnivVola}}, \code{\link{MultiEWMA}}.
}
\examples{
# Load returns
data("CRSPday_zoo")
y <- CRSPday_zoo

# Compute multivariate EWMA model (RiskMetrics)
EWMA <- MultiEWMA(y, center = TRUE)

# Compute and plot standardized residuals
EWMAres <- residuals(EWMA)
plot(EWMAres)
}
\concept{standardized residuals}
