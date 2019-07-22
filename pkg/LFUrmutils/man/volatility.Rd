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
\method{vola}{MultiEWMA}(object, offdiagonal = FALSE, duplicates = FALSE, ...)
}
\arguments{
  \item{object}{An object of class \code{UnivVola}, \code{fGARCH} or \code{MultiEWMA}.}
  \item{offdiagonal}{Whether diagonal elements should be kept or not. \code{FALSE} (the default) or \code{TRUE}.}
  \item{duplicates}{Whether duplicate off-diagonal elements should be kept or not. \code{TRUE} (the default) or \code{FALSE}.}
  \item{...}{Arguments passed to methods.}
}
\details{
Depending on the application, diagonal elements of multivariate volatility models may be printed. However, they are usually not required in practice, and may not have a meaningful interpretation. However, if they are printed, duplicated off-diagonal elements may be deleted, due to the symmetry of matrix.
}
\value{
A \code{zoo} object containing the volatilities. In the multivariate model each row of the output objects represents the volatility matrix (or elements thereof) of the corresponding day.
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
\code{\link{UnivVola}}, \code{\link{MultiEWMA}}, \code{\link{varcov}}, \code{\link{ccor}}.
}
\examples{
library("zoo")

# Load returns
data("CRSPday_zoo")
y <- CRSPday_zoo

# Compute multivariate EWMA model (RiskMetrics)
EWMA <- MultiEWMA(y, center = TRUE)

# Compute variance-covariance matrix
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