\name{bcVAR}
\alias{bcVAR}


\title{Bias-corrected LS VAR}

\encoding{latin1}

\description{
  Estimation of the Bias-corrected LS of a VAR(p) model.
}

\usage{
bcVAR(data = data, p = 1, constant = TRUE, ...)
}

\arguments{
  \item{data}{Data item containing the endogenous variables.}
  \item{p}{Integer for the lag order (default is p = 1).}
  \item{constant}{If constant = TRUE, a intercept vector is included in the model.}
  \item{\dots}{Not used.}
}

\details{
  Estimates a VAR by OLS and applies the bias correction of Pope, A. L. (1990). The model is of the following form: 
  
  \deqn{
    \bold{y}_t = [\bold{v}, A_1, \ldots, A_p] + \bold{u}_t,
  }

  where \eqn{\bold{y}_t} is a \eqn{K \times 1} vector of endogenous
  variables and \eqn{\bold{u}_t} is assumed to be iid white noise. The coefficient matrices \eqn{A_1, \ldots, A_p} are   of dimension \eqn{K \times K}. The intercept vector \eqn{\bold{v}} is of dimension \eqn{K \times 1}. The LS estimator   is
  
   \deqn{
    \hat{A} = [\hat{\bold{v}}, \hat{A}_1, \ldots, \hat{A}_p] = YZ'(ZZ')^-1,
  }
  
  where...
  }

\value{
  A list with class attribute \sQuote{\code{varest}} holding the
  following elements:\cr
  
  \item{varresult}{list of \sQuote{\code{lm}} objects with the bias-corrected LS of the VAR(p) process.}
  \item{datamat}{The data matrix of the endogenous and explanatory variables.}
  \item{y}{The data matrix of the endogenous variables.}
  \item{p}{An integer specifying the lag order.}
  \item{K}{An integer specifying the dimension of the VAR.}
  \item{obs}{An integer specifying the number of used observations.}
  \item{totobs}{An integer specifying the total number of observations.}
  \item{restrictions}{Always \code{NULL}. Zero restrictions of the VAR(p) are not implemented.}
  \item{call}{The \code{call} to \command{bcVAR()}.}
  
}

\author{Simon Röck}

\references{
  Pope, A. L. (1990),  \emph{Biases of estimators in multivariate non-gausssian autoregressions}, Journal of Time        Series   Analysis.

  Kilian, L., & Lütkepohl, H. (2017), \emph{Structural Vector Autoregressive Analysis}, Cambridge University Press,      Cambridge.
  
}


\examples{
## load data of package
data("USmacro", package = "bcVAR")

## detrend data (substract mean)
dataDT <- apply(USmacro, 2, function(y) y - mean(y))

## example for a bias-corrected LS VAR(4) model (see Chapter 2 of Kilian, L., & Lütkepohl, H. (2017))
bcVAR(dataDT, p = 4, constant = TRUE)

}

\keyword{regression}
\concept{VAR}
\concept{Vector autoregressive model}
\concept{Bias-corrected LS VAR}

