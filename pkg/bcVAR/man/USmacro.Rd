\name{USmacro}

\alias{USmacro}

\title{US: Macroeconomic Time Series}

\encoding{latin1}

\description{
US macroeconomic data for log of U.S. real GNP, the corresponding GNP deflator in logs and the federal funds rate, averaged by quarter. The sample range is from the 4thQ 1954 until 4thQ 2007. The dataset is also used in Kilian, L., and L�tkepohl, H. (2017). 
}

\usage{data("USmacro", package = "bcVAR")}

\format{
A data frame containing three time series variables with 213 observations from the 4thQ 1954 until 4thQ 2007.
\describe{
  \item{drgdp}{Log of U.S. real GNP.}
  \item{irate}{GNP deflator in logs.}
  \item{infl}{Federal funds rate, averaged by quarter.}
}
}

\source{
Data originally provided by FRED (\url{https://fred.stlouisfed.org/}).
}


\examples{
data("USmacro", package = "bcVAR")

plot.ts(USmacro)
}

\keyword{datasets}
