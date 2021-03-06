\name{bcVAR}
\alias{bcVAR}

\title{Bias-Corrected Least Squares Estimator for VAR Models}

\encoding{latin1}

\description{
  Estimation of the bias-corrected least squares of a VAR(p) model.
}

\usage{
bcVAR(data = data, p = 1, type = c("const", "none"), ...)
}

\arguments{
  \item{data}{Data item containing the endogenous variables with the dimension \eqn{T \times K}, where T is the length of the time series and K is the number of endogenous variables.}
  \item{p}{Integer for the lag order (default is p = 1).}
  \item{type}{Type of deterministic regressors to include.}
  \item{\dots}{Currently not used.}
}

\details{
  Estimates a VAR by OLS and applies the bias correction proposed by Pope (1990). The model is of the following form: 
  
  \deqn{
    Y_t = \bold{v} + \bold{A}Y_{t-1} + \bold{u}_t,
  }
where \eqn{Y_t} is a \eqn{Kp \times 1} vector of endogenous
  variables and \eqn{\bold{u}_t} with the same dimension is assumed to be iid white noise. The companion matrix \eqn{\bold{A}} is of dimension \eqn{Kp \times Kp}. The intercept vector \eqn{\bold{v}} is of dimension \eqn{Kp \times 1}. By setting the type argument to \code{const}, the intercept vector \eqn{\bold{v}} is included. The bias of the least squares (LS) estimator \eqn{\hat{\bold{A}}} for \eqn{\bold{A}} is 
  
    \deqn{
    -B_{\bold{A}}/T + O(T^{-3/2}),
  }
  
where
\deqn{
  B_{\bold{A}} = \Sigma_U\left[(I_{Kp}-\bold{A}')^{-1}+\bold{A}'(I_{Kp}-\bold{A}'^2)^{-1}+\sum_{\lambda} \lambda(I_{Kp}-\lambda \bold{A}')^{-1}\right] \Gamma_Y(0)^{-1},
}

\eqn{\Gamma_Y(0)^{-1} = E(Y_tY_t'), \Sigma_U = E(U_tU_t')} and the sum over the eigenvalues  \eqn{\lambda} of  \eqn{\bold{A}} weighted by their multiplicities. Adding \eqn{B_{\hat{\bold{A}}}/T} to \eqn{\hat{\bold{A}}} yields the bias-corrected LS estimator \eqn{\hat{\bold{A}}^{BC}}. For more details regarding the bias-corrected LS see for example chapter 2 of Kilian, L., & L�tkepohl, H. (2017). The resulting object of \command{bcVAR()} has the same class attribute as the object of \command{VAR()} from the \sQuote{vars} package of Pfaff, B. (2008).
  }
  
  

\value{
  A list with class attribute \sQuote{\code{varest}} (class attribute of the \sQuote{vars} package) holding the
  following elements:\cr
  
  \item{varresult}{List of pseudo \sQuote{\code{lm}} objects with the bias-corrected LS of the VAR(p) process. Pseudo \sQuote{\code{lm}} objects, because the objects are not generated by \command{lm()}.}
  \item{datamat}{The data matrix of the endogenous and explanatory variables.}
  \item{y}{The data matrix of the endogenous variables.}
  \item{p}{An integer specifying the lag order.}
  \item{K}{An integer specifying the dimension of the VAR.}
  \item{obs}{An integer specifying the number of used observations.}
  \item{totobs}{An integer specifying the total number of observations.}
  \item{restrictions}{Always \code{NULL}. Zero restrictions of the VAR(p) are not implemented.}
  \item{call}{The \code{call} to \command{bcVAR()}.} 
}

\author{Simon R�ck}

\references{
  Pope, A. L. (1990). \dQuote{Biases of Estimators in Multivariate Non-Gausssian Autoregressions},
  \emph{Journal of Time Series Analysis}, \bold{11}(3), 249--258. \doi{10.1111/j.1467-9892.1990.tb00056.x}

  Kilian, L., & L�tkepohl, H. (2017). \emph{Structural Vector Autoregressive Analysis}, Cambridge University Press, Cambridge.
  
  Pfaff, B. (2008). \dQuote{VAR, SVAR and SVEC Models: Implementation within R Package vars},
  \emph{Journal of Statistical Software}, \bold{27}(4), 1--32. \doi{10.18637/jss.v027.i04}
}


\seealso{\code{\link[vars]{VAR}}} 


\examples{
## load data of package
data("USmacro", package = "bcVAR")

## detrend data (substract mean)
dataDT <- apply(USmacro, 2, function(y) y - mean(y))

## bias-corrected LS VAR(4) model (see Chapter 2 of Kilian, L., & Luetkepohl, H. (2017))
bcVAR(dataDT, p = 4, type = "const")
}

\keyword{regression}
\concept{VAR}
\concept{Vector autoregressive model}
\concept{Bias-corrected LS VAR}

