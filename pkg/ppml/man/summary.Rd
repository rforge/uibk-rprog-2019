\name{summary}
\alias{summary}
\alias{summary.ppml}


\title{
Summary Table with Heteroskedasticity Robust Standard Errors
}

\description{
This method is based on the function \code{coeftest} from the \sQuote{lmtest} package and the generic function \code{\link[stats]{summary.glm}} for the object of class \code{ppml}.
}

\usage{
\S3method{summary}{ppml}(object, ...)
}

\arguments{
  \item{object}{An object of class \code{ppml}.}
  \item{...}{Currently not used.}
}

\value{
 A list with class attribute \code{\link[stats]{summary.glm}} holding the
  following elements:\cr
  
  \item{call}{}
  \item{terms}{}
  \item{family}{}
  \item{deviance}{}
  \item{aic}{}
  \item{contrasts}{}
  \item{df.residual}{}
  \item{null.deviance}{}
  \item{df.null}{} 
  \item{iter}{} 
  \item{deviance.resid}{} 
  \item{coefficients}{} 
  \item{aliased}{} 
  \item{dispersion}{} 
  \item{df}{} 
  \item{cov.unscaled}{} 
  \item{cov.scaled}{} 
  \item{attr}{} 

  See \code{\link[stats]{glm}} for more information for the elements. 

}

\author{
Lukas Kleinheinz
}

\examples{
## load data of package
data("TradeFlows", package = "ppml")

## gravity model with robust standard errors 
model <- ppml(formula = s ~ ., data = TradeFlows)

## summary
summary(model)


}

\concept{bread}

