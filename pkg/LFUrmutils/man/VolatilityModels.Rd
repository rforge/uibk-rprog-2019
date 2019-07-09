\name{VolatilityModels}
\alias{UnivVola}
\alias{MultiEWMA}

\title{Volatility models}
\description{
Computes univariate (moving average, weighted average, Risk Metrics) and multivariate (Risk Metrics) volatility models.
}
\usage{
UnivVola(returns, width = 30, lambda = 0.94, 
    type = c("RiskMetrics", "WeightedAverage", "MovingAverage"), 
    center = FALSE)
MultiEWMA(returns, lambda = 0.94, center = FALSE)
}
\arguments{
  \item{returns}{A univariate (multivariate) return series in the case of univariate (multivariate) volatility models.}
  \item{width}{How many past observations should be use to calculate the predicted volatility for the next day? Default is 30.}
  \item{lambda}{Weight of previous day's volatility. Ignored for models of type "MovingAverage".}
  \item{type}{Specifies the type of volatility model.}
  \item{center}{Should returns be centered before computing the volatility? Included for comparability reasons in \code{EWMAvol} in package \code{MTS}.}
}
\details{
\code{UnivVola} computes univariate volatility models. Currently simple moving average models as well as both types of exponentially weighted moving average models are implemented. \code{MultiEWMA} computes the multivariate exponentially weighted moving average model in RiskMetrics form.
}
\value{
\item{Variance }{\code{zoo} object of computed variances. In the multivariate model each row of the output objects represents the variance-covariance matrix of the corresponding day.}
\item{Returns }{\code{zoo} object of used return series.}
\item{variant }{Character string indicating the type of volatility model used. \code{NA} in the case of multivariate volatility models.}
\item{width }{Numeric value indicating the number of past observations that have been used to calculate the predicted volatility. \code{NA} in the case of \code{RiskMetrics} models.}
\item{lambda }{Value of the lambda parameter used. Only available for univariate and multivariate models of type "RiskMetrics".}
\item{centered }{Whether the return series has been centered within the function or not.}
}
\references{
Jorion (2007). Value at Risk, 3rd. McGraw-Hill. New York.
}
\author{
Bernhard Eder
}
\seealso{
\code{\link{residuals}}, \code{\link{mse}}, \code{\link{varcov}}, \code{\link{vola}}.
}
\examples{
# Load financial prices and compute returns
library("tseries")
library("zoo")
price <- get.hist.quote("db", "2010-01-01", "2019-07-08",
                     quote = "Adjusted", quiet = TRUE)
return <- diff(log(price))

# Compute univariate volatility models
ma <- UnivVola(return, type = "MovingAverage", center = TRUE)
wa <- UnivVola(return, type = "WeightedAverage", center = TRUE)

# Access variances and compare returns
head(ma$Variances)
head(ma$Returns)
head(return)
}
