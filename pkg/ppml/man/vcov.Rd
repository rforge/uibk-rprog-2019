\name{vcov}
\alias{vcov}
\alias{vcov.ppml}


\title{
Making Sandwiches with Bread and Meat
}

\description{
This method is based on the function \code{sandwich} from the \sQuote{sandwich} package for the object of class \code{ppml} to estimate the sandwich covariance matrix.
}

\usage{
\S3method{vcov}{ppml}(object, ...)
}

\arguments{
  \item{object}{An object of class \code{ppml}.}
  \item{...}{Currently not used.}
}

\value{
Returns the sandwich covariance matrix.
}


\references{
  Zeileis A (2006), Object-Oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://www.jstatsoft.org/v16/i09/}.
}

\author{
Lukas Kleinheinz
}

\examples{
## load data of package
data("TradeFlows", package = "ppml")

## gravity model with robust standard errors 
model <- ppml(formula = s ~ ., data = TradeFlows)

## sandwich matrix
vcov(model)


}

\concept{bread}

