\name{CRSPday_zoo}
\alias{CRSPday_zoo}
\docType{data}
\title{
Daily Returns from the CRSP Database
}
\description{
Daily observations from 1969-01-03 to 1998-12-31. \cr
Number of Observations: 2528. \cr
Observation: returns. \cr
Country: United States. \cr
}
\usage{data("CRSPday_zoo")}
\format{
  A multivariate \code{zoo} object containing:
  \describe{
    \item{\code{ge}}{the return for General Electric, Permno 12060}
    \item{\code{ibm}}{the return for IBM, Permno 12490}
    \item{\code{mobil}}{the return for Mobil Corporation, Permno 15966}
  }
}
\source{
Dataset based on \code{\link[Ecdat:CRSPday]{CRSPday}} in package \pkg{Ecdat}. Original source: Center for Research in Security Prices, Graduate School of Business, University of Chicago, 725 South Wells - Suite 800, Chicago, Illinois 60607, \url{http://www.crsp.com}.
}
\references{
Davidson, R. and James G. MacKinnon (2004) Econometric Theory and Methods, New York, Oxford University Press, \url{http://www.econ.queensu.ca/ETM/}, chapter 7, 9 and 15.
}
\keyword{datasets}
\concept{financial returns}