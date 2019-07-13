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
  \item{data}{Data item containing the endogenous variables}
  \item{p}{Integer for the lag order (default is p = 1).}
  \item{constant}{If constant = TRUE, a intercept is included in the model.}
  \item{\dots}{Not used.}
}

\details{
  Estimates a VAR by OLS. The model is of the following
  form:
}

\value{
  A list with class attribute \sQuote{\code{varest}} holding the
  following elements:\cr

  \item{V}{Constants of the endogenous variables.}
  \item{A}{Estiamted companion matrix of the VAR(p) process.}
  \item{bcLS}{Bias-corrected LS of the companion matrix of the VAR(p) process.}
  \item{SIGMA}{The covariance matrix of the innovations.}
  \item{Y}{The model matrix.}
  \item{Uhat}{The estimated residuals.}
  \item{sumstat}{Returns the values T, K and p.}
}

\author{Simon Röck}

\references{
  Pope, A. L. (1990),  \emph{Biases of estimators in multivariate non-gausssian autoregressions}, Journal of Time        Series   Analysis.

  Kilian, L., & Lütkepohl, H. (2017), \emph{Structural Vector Autoregressive Analysis}, Cambridge University Press,      Cambridge.
  
}



\examples{
data("USmacro", package = "bcVAR")

## detrend data (substract mean)
dataDT <- apply(USmacro, 2, function(y) y - mean(y))

## Example for a bias-corrected LS VAR(4) model (see Chapter 2 of Kilian, L., & Lütkepohl, H. (2017))
bcVAR(dataDT, p = 4, constant = TRUE)

}

\keyword{regression}
\concept{VAR}
\concept{Vector autoregressive model}
\concept{Bias-corrected LS VAR}

