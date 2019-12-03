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
    center = FALSE, exchange = "UnitedStates/NYSE")
MultiEWMA(returns, lambda = 0.94, center = FALSE, exchange = "UnitedStates/NYSE")
}
\arguments{
  \item{returns}{A univariate (multivariate) return series in the case of univariate (multivariate) volatility models.}
  \item{width}{How many past observations should be use to calculate the predicted volatility for the next day? Default is 30.}
  \item{lambda}{Weight of previous day's volatility. Ignored for models of type "MovingAverage". The value must be between 0 and 1. If a negative value is supplied to the function, lambda is estimated from the data.}
  \item{type}{Specifies the type of volatility model.}
  \item{center}{Should returns be centered before computing the volatility? Included for comparability reasons in \code{EWMAvol} in package \code{MTS}.}
  \item{exchange}{Stock exchange at which the asset is traded. Default is NYSE. See the details of \code{Calendars} in package \code{RQuantLib} for further options.}
}
\details{
\code{UnivVola} computes univariate volatility models. Currently simple moving average models as well as both types of exponentially weighted moving average models are implemented. \code{MultiEWMA} computes the multivariate exponentially weighted moving average model in RiskMetrics form. The default value of \code{lambda} is 0.94, as suggested for daily data by RiskMetrics (J.P.Morgan/Reuters, 1996). However, when a negative value for lambda is supplied, an optimal lambda will be estimated from the data, using a conditional maximum likelihood approach (see Ruppert and Matteson (2015) for further details).
}
\note{
Note that the models as implemented here implicitly assume that the time series of returns have a conditional expectation of zero, following the usual risk management conventions. Where this is not the case in practical applications, you may estimate a model for the conditional expected value before using the functions herein and use the residuals of that model instead of the original returns as inputs.
}
\value{
\item{Variance }{\code{zoo} object of computed variances. In the multivariate model each row of the output objects represents the variance-covariance matrix of the corresponding day.}
\item{Returns }{\code{zoo} object of used return series.}
\item{variant }{Character string indicating the type of volatility model used. \code{NA} in the case of multivariate volatility models.}
\item{width }{Numeric value indicating the number of past observations that have been used to calculate the predicted volatility. \code{NA} in the case of \code{RiskMetrics} models.}
\item{lambda }{Value of the lambda parameter used. Only available for univariate and multivariate models of type "RiskMetrics" or "WeightedAverage".}
\item{lambda.se }{Standard error of lambda, if estimated.}
\item{llh }{Value of the log-likelihood, if lambda is estimated.}
\item{centered }{Whether the return series has been centered within the function or not.}
}
\references{
Danielsson (2011). Financial Risk Forecasting. Wiley. Chichester. 
J.P.Morgan/Reuters (1996). RiskMetrics - Technical Document. 4th. New York. 
Jorion (2007). Value at Risk, 3rd. McGraw-Hill. New York. 
Ruppert and Matteson (2015). Statistics and Data Analysis for Financial Engineering, 2ed. Springer. New York. 
}
\author{
Bernhard Eder
}
\seealso{
\code{\link{residuals}}, \code{\link{mse}}, \code{\link{varcov}}, \code{\link{vola}}, \code{\link{NextBusinessDay}}.
}
\examples{
# Load returns
data("CRSPday_zoo")
y <- CRSPday_zoo

# Compute univariate volatility models
ma <- UnivVola(y[, 1], type = "MovingAverage", center = TRUE)
wa <- UnivVola(y[, 1], type = "WeightedAverage", center = TRUE)

# Access variances and compare returns
head(ma$Variances)
head(ma$Returns)
head(y[, 1])
}
\keyword{models}
\keyword{ts}
\concept{moving average}
\concept{exponential weighted moving average}
\concept{EWMA}
\concept{RiskMetrics}
\concept{univariate volatility models}
\concept{multivariate volatility models}