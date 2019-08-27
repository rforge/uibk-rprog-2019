\name{ppml}
\alias{ppml}

\title{Poisson Pseudo Maximum Likelihood Estimator}

\encoding{latin1}

\description{

Can be used to estimate gravity models using \code{\link[stats]{glm}} with the quasipoisson distribution, a log-link and heteroskedasticity robust standard errors.

}

\usage{
ppml(formula = formula, data = data, ...)
}

\arguments{
  \item{formula}{An expression of the form \eqn{y \tilde model}}
  \item{data}{Data containing the variables of interest}
  \item{\dots}{Currently not used.}
}

\details{

The maximum likelihood estimator for the Poisson distribution allows for consistent estimates provided that the conditional mean is correctly specified. Consistency requires the score vector to have expected value zero:
\deqn{
\sum_i^n E[(y_i - \exp(\bold{x}_i' \bold{\beta})) \bold{x}_i] = 0
}

This is only the case if
\deqn{
E[y_i|\bold{x}_i] = \exp(\bold{x}_i' \bold{\beta})
}
Therefore the PMLE can still be used even if the data generating process of \eqn{y_i} is not the Poisson (Poisson Pseudo Maximum Likelihood) (Cameron and Trivedi, 2013).

The assumption of equidispersion does not hold for the PPMLE as the nuisance parameter \eqn{\phi} from the GLM framework is not restricted to be 1 anymore:

\deqn{
Var[y_i|x_i] = \phi \frac{\partial^2 b(\theta_i)}{\partial \theta_i^2} = \phi \cdot e^{\theta_i} = \phi \mu_i
}

and
 \deqn{
  \frac{Var[y|x]}{E[y|x]} = \phi
}

This means that the PPMLE is optimal if the conditional variance is proportional to the conditional mean (Santos Silva, 2015).

As stated by Arvis and Shepherd (2013), the Pseudo Maximum Likelihood applicable to continuous variables can be derived from the log-likelihood function and the large number Stirling approximation for \eqn{\ln(y_i!)}:  \eqn{\ln{y_i!} \approx y_i \ln{y_i} - y_i}. The pseudo-log-likelihood function \eqn{\tilde{\ell}} of the Poisson PML is:

\deqn{
\tilde{\ell}(\bold{\beta}|y_i) \approx \sum_i^n\biggl[y_i \bold{x}_i' \bold{\beta} - \exp(\bold{x}_i' \bold{\beta})-y_i \ln{y_i}+y_i\biggr]
}

For the next step it is assumed that \eqn{y_{ij}} depicts the value of trade flows and \eqn{\hat{y}_{ij} = \exp(\bold{x}_{ij}' \bold{\beta})} is the modelled value regressed on the actual value \eqn{y_{ij}}. Maximizing the pseudo-log-likelihood function with respect to \eqn{\bold{\beta}}:

\deqn{
\frac{\partial \tilde{\ell}(\bold{\beta}|y_i)}{\partial \bold{\beta}} = \sum_i \sum_j \biggl[y_{ij} \bold{x}_{ij}' - \exp(\bold{x}_{ij}' \bold{\beta})\bold{x}_{ij}' \biggr] = 0
}

\deqn{
\sum_i \sum_j y_{ij} \bold{x}_{ij}' = \sum_i \sum_j \exp(\bold{x}_{ij}' \bold{\beta})\bold{x}_{ij}'
}

\deqn{
\sum_i \sum_j y_{ij} \bold{x}_{ij}' = \sum_i \sum_j \hat{y}_{ij} \bold{x}_{ij}'
}
The predicted total trade flows between i and j when a dummy is true equals the actual total trade flows for the same subset (Arvis and Shepherd, 2013). Using PPML with importer- and exporter fixed effects therefore yields the desirable property of preserving  the sums of the trade flows.

Another advantage of using \code{ppml} compared to the estimation of log-linearized gravity models via OLS is its ability to handle zero trade flows. As \eqn{\ln{0}} is not defined, the respective observations will be dropped. This loss of information can be avoided by estimating the multiplicative model via \code{ppml}.

}
  
  
\value{

 A list with class attribute \code{\link[stats]{glm}} and \code{ppml} holding the
  following elements:\cr
  
  \item{coefficients}{}
  \item{residuals}{}
  \item{fitted.values}{}
  \item{effects}{}
  \item{R}{}
  \item{rank}{}
  \item{qr}{}
  \item{family}{}
  \item{linear.predictors}{} 
  \item{deviance}{} 
  \item{aic}{} 
  \item{null.deviance}{} 
  \item{iter}{} 
  \item{weights}{} 
  \item{prior.weights}{} 
  \item{df.residual}{} 
  \item{df.null}{} 
  \item{y}{} 
  \item{converged}{} 
  \item{boundary}{} 
  \item{model}{} 
  \item{call}{} 
  \item{formula}{} 
  \item{terms}{}
  \item{data}{} 
  \item{offset}{} 
  \item{control}{} 
  \item{method}{} 
  \item{contrasts}{} 
  \item{xlevels}{} 
  
  
  See \code{\link[stats]{glm}} for more information for the elements. 
}

\seealso{
\code{\link{glm}}
}

\author{Lukas Kleinheinz}

\references{
  Arvis, J.-F. and Shepherd, B. (2013). \dQuote{The poisson quasi-maximum likelihood estimator: a solution to the ``adding up'' problem in gravity models}, \emph{Applied Economics Letters}, \bold{20}(6), 515--519.
  
  Cameron, A. C. and Trivedi, P.K. (2013). \dQuote{Regression Analysis of Count Data}, Econometric Society Monographs, 2 edn, Cambridge University Press

  Santos Silva, J. (2015). \dQuote{The log of gravity page}, \url{http://personal.lse.ac.uk/tenreyro/lgw.html}, Accessed: 2019-08-20


}

\examples{
## load data of package
data("TradeFlows", package = "ppml")

## gravity model with robust standard errors 
model <- ppml(formula = s ~ log_dist + contig + comlang_off + 
                colony + smctry + i_j + log_stri_i_j
              + rta + i + j, data = TradeFlows)
              
## summary of the model
summary(model)

}

\keyword{regression}



