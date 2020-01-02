\name{vola}
\alias{vola}
\alias{vola.UnivVola}
\alias{vola.MultiEWMA}
\alias{vola.fGARCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conditional volatilities
}
\description{
A function to compute the conditional volatilities of univariate and multivariate volatility models.
}
\usage{
\method{vola}{UnivVola}(object, ...)
\method{vola}{fGARCH}(object, ...)
\method{vola}{MultiEWMA}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{UnivVola}, \code{fGARCH} or \code{MultiEWMA}.}
  \item{...}{Arguments passed to methods.}
}
\details{
Returns conditional volatilities of univariate and multivariate volatility models. In the latter case, the values in each row correspond to the square-root of the diagonal elements of the corresponding conditional variance-covariance matrix.
}
\value{
A \code{zoo} object containing the volatilities.
}
\references{
Danielsson (2011). Financial Risk Forecasting. Wiley. Chichester.\cr

Jorion (2007). Value at Risk, 3rd. McGraw-Hill. New York.\cr

Ruppert and Matteson (2015). Statistics and Data Analysis for Financial Engineering, 2nd. Springer. New York.
}
\author{
Bernhard Eder
}
\seealso{
\code{\link{UnivVola}}, \code{\link{MultiEWMA}}, \code{\link{varcov}}, \code{\link{ccor}}.
}
\examples{
library("zoo")

# Load returns
data("CRSPday_zoo")
y <- CRSPday_zoo

# Compute multivariate EWMA model (RiskMetrics)
EWMA <- MultiEWMA(y, center = TRUE)

# Compute variance-covariance matrix based on volatilities and correlation matrix
  # Compute relevant volatilities
  EWMAvola <- vola(EWMA)
  head(EWMAvola)
  
  # Compute conditional correlations
  EWMAcor <- ccor(EWMA, duplicates = TRUE, diagonal = TRUE)
  
  # Initialize output object
  TT <- dim(EWMAcor)[1]
  c <- sqrt(dim(EWMAcor)[2])
  VarCov <- matrix(NA, dim(EWMAcor)[1], dim(EWMAcor)[2])

  # Compute variance-covariance matrix
  for(i in 1:TT){
    VarCov[i, ] <- c(diag(as.numeric(EWMAvola[i, ])) 
                      \%*\% matrix(as.numeric(EWMAcor[i, ]), c, c,  byrow = TRUE) 
                      \%*\% diag(as.numeric(EWMAvola[i, ]))
                      )
  }
  
# Format output
VarCov <- zoo(VarCov, order.by = index(EWMAvola))

# Compare variance-covariance matrices
head(VarCov)
head(EWMA$Variances)
}
\concept{conditional volatility}
\concept{standard deviation}