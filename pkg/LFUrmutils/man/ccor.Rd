\name{ccor}
\alias{ccor}
\alias{ccor.MultiEWMA}

\title{Conditional correlations}
\description{
Computes conditional correlations based on multivariate volatility models. Currently only implemented for objects of class \code{MultiEWMA}.
}
\usage{
ccor(object, ...)
\method{ccor}{MultiEWMA}(object, diagonal = TRUE, duplicates = TRUE, ...)
}
\arguments{
  \item{object}{Object of class \code{MultiEWMA}.}
  \item{diagonal}{Whether diagonal elements should be kept or not. \code{TRUE} (the default) or \code{FALSE}.} 
  \item{duplicates}{Whether duplicate off-diagonal elements should be kept or not. \code{TRUE} (the default) or \code{FALSE}.}
  \item{...}{Arguments passed to methods.}
}
\details{
Depending on the application, diagonal elements may be deleted, since they are always equal to one. Likewise, duplicated off-diagonal elements may be deleted due to the symmetry of the variance-covariance matrix. However, keep in mind that both diagonal elements as well as duplicated off-diagonal elements may be necessary to compute the variance-covariance matrix. 
}
\note{
Note that the function is currently only implemented for objects of class \code{MultiEWMA}. Implementation for other multivariate volatility models may be added in future releases.
}
\value{
\item{rho }{A \code{zoo} object containing the conditional correlations between the volatilities of the corresponding assets. Each row of the output object represents the correlation matrix (or elements thereof) of the corresponding day.}
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
\code{\link{MultiEWMA}}, \code{\link{varcov}}, \code{\link{vola}}, \code{\link{ccor}}.
}
\examples{
library("zoo")

# Load returns
data("CRSPday_zoo")
y <- CRSPday_zoo

# Compute multivariate EWMA model (RiskMetrics)
EWMA <- MultiEWMA(y, center = TRUE)

# Compute conditional correlations
EWMAcor <- ccor(EWMA, duplicates = TRUE, diagonal = TRUE)
head(EWMAcor)

# Compute conditional variance-covariance matrix
  # Compute relevant volatilities
  EWMAvola <- vola(EWMA, duplicates = TRUE, offdiagonal = FALSE)
  
  # Initialize output object
  TT <- dim(EWMAcor)[1]
  c <- sqrt(dim(EWMAcor)[2])
  VarCov <- matrix(NA, dim(EWMAcor)[1], dim(EWMAcor)[2])

# Compute conditional variance-covariance matrix
  for(i in 1:TT){
    VarCov[i, ] <- c(diag(as.numeric(EWMAvola[i, ])) 
                      \%*\% matrix(as.numeric(EWMAcor[i, ]), c, c,  byrow = TRUE) 
                      \%*\% diag(as.numeric(EWMAvola[i, ]))
                      )
  }
  
# Format output
VarCov <- zoo(VarCov, order.by = index(EWMAvola))

# Compare conditional variance-covariance matrices
head(VarCov)
head(EWMA$Variances)
}
