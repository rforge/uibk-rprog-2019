\name{bcVAR}
\alias{bcVAR}

\title{Bias-corrected LS VAR}

\encoding{latin1}

\description{
  Estimation of the Bias-corrected LS of a VAR(p) model.
}

\usage{
bcVAR(data = data, p = 1, type = c("const", "none"), ...)
}

\arguments{
  \item{data}{Data item containing the endogenous variables.}
  \item{p}{Integer for the lag order (default is p = 1).}
  \item{type}{Type of deterministic regressors to include.}
  \item{\dots}{Currently not used.}
}

\details{
  Estimates a VAR by OLS and applies the bias correction proposed by Pope, A. L. (1990). The model is of the following form: 
  
  \deqn{
    Y_t = \bold{v} + \bold{A}Y_{t-1} + \bold{u}_t,
  }
where \eqn{Y_t} is a \eqn{Kp \times 1} vector of endogenous
  variables and \eqn{\bold{u}_t} with the same dimension is assumed to be iid white noise. The comapanion matrix \eqn{\bold{A}} is of dimension \eqn{Kp \times Kp}. The intercept vector \eqn{\bold{v}} is of dimension \eqn{Kp \times 1}. By setting the type argument to \code{const}, the intercept vector \eqn{\bold{v}} is included. The bias of the LS estimator \eqn{\hat{\bold{A}}} for \eqn{\bold{A}} is 
  
    \deqn{
    -B_{\bold{A}}/T + O(T^{(-3/2)},
  }
  
where
\deqn{
  B_{\bold{A}} = \Sigma_U\left[(I_{Kp}-\bold{A}')^{-1}+\bold{A}'(I_{Kp}-\bold{A}'^2)^{-1}+\sum_{\lambda} \lambda(I_{Kp}-\lambda \bold{A}')^{-1}\right] \Gamma_Y(0)^{-1},
}

\eqn{\Gamma_Y(0)^{-1} = E(Y_tY_t'), \Sigma_U = E(U_tU_t')} and the sum over the eigenvalues  \eqn{\lambda} of  \eqn{\bold{A}} weighted by their multiplicities. For more details regarding the bias-corrected LS see for example chapter 2 of Kilian, L., & Lütkepohl, H. (2017). 

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

## example for a bias-corrected LS VAR(4) model (see Chapter 2 of Kilian, L., & Luetkepohl, H. (2017))
bcVAR(dataDT, p = 4, type = "const")

}

\keyword{regression}
\concept{VAR}
\concept{Vector autoregressive model}
\concept{Bias-corrected LS VAR}

