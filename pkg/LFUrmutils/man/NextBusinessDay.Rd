\name{NextBusinessDay}
\alias{NextBusinessDay}
\title{
Next Business Day
}
\description{
Compute the next business day on the New York Stock Exchange (NYSE) that follows at the end of a time series. The function is mainly used internally by functions within this package, but may be useful for other purposes as well.
}
\usage{
NextBusinessDay(x, exchange = "UnitedStates/NYSE")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{zoo}.}
  \item{exchange}{A string identifying one of the supported QuantLib calendars. Default is "UnitedStates/NYSE". See Details in \link[RQuantLib:Calendars]{Calendars} for further options.}
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
# Load returns
data("CRSPday_zoo")
y <- CRSPday_zoo

# End of the sample
end(y)

# Next business day
NextBusinessDay(y)
}
\concept{calendar}
\concept{holidays}