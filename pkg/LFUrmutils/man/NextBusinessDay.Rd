\name{NextBusinessDay}
\alias{NextBusinessDay}
\title{
Next Business Day
}
\description{
Compute the next business day on the New York Stock Exchange (NYSE) that follows at the end of a time series. The function is mainly used internally by functions within this package, but may be useful for other purposes as well.
}
\usage{
NextBusinessDay(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{zoo}.}
}
\value{
An object of class \code{date}.
}
\author{
Bernhard Eder
}
\note{
Objects of other classes than \code{zoo} may also be supplied, but are not tested. The function may work for any class for which function \code{end()} is defined. Call \code{methods("end")} to see for which classes this function is defined.
}
\examples{
# Load financial prices and compute returns
library("tseries")
library("zoo")
prices <- get.hist.quote("msft", "2009-01-01", "2009-12-31",
                     quote = "Adjusted", quiet = TRUE)

# End of the sample
end(prices)

# Next business day
NextBusinessDay(prices)
}
