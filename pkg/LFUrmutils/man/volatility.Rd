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
\method{vola}{MultiEWMA}(object, offdiagonal = TRUE, duplicates = TRUE, ...)
}
\arguments{
  \item{object}{An object of class \code{UnivVola}, \code{fGARCH} or \code{MultiEWMA}.}
  \item{offdiagonal}{Whether diagonal elements should be kept or not. \code{TRUE} (the default) or \code{FALSE}.}
  \item{duplicates}{Whether duplicate off-diagonal elements should be kept or not. \code{TRUE} (the default) or \code{FALSE}.}
  \item{...}{Arguments passed to methods.}
}
\details{
Depending on the application, off-diagonal elements of multivariate volatility models may be deleted. Likewise, duplicated off-diagonal elements may be deleted as well, due to the symmetry of the volatility matrix.
}
\value{
A \code{zoo} object containing the volatilities. In the multivariate model each row of the output objects represents the correlation matrix (or elements thereof) of the corresponding day.
}
\references{
Jorion (2007). Value at Risk, 3rd. McGraw-Hill. New York.
}
\author{
Bernhard Eder
}
\seealso{
\code{\link{UnivVola}}, \code{\link{MultiEWMA}}, \code{\link{varcov}}, \code{\link{ccor}}.
}
\examples{
# Load financial prices and compute returns
library("tseries")
library("zoo")
p1 <- get.hist.quote("msft", "2000-01-01", "2009-12-31",
                     quote = "Adjusted", quiet = TRUE)
p2 <- get.hist.quote("ibm", "2000-01-01", "2009-12-31",
                     quote = "Adjusted", quiet = TRUE)
y1 <- diff(log(p1)) * 100
y2 <- diff(log(p2)) * 100
y <- merge(y1, y2)

# Compute multivariate EWMA model (RiskMetrics)
EWMA <- MultiEWMA(y, center = TRUE)

# Compute variance-covariance matrix
  # Compute relevant volatilities
  EWMAvola <- vola(EWMA, duplicates = TRUE, offdiagonal = FALSE)
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
