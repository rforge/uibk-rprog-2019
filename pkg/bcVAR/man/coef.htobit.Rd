\name{coef.htobit}

\alias{coef.htobit}
\alias{vcov.htobit}
\alias{terms.htobit}
\alias{model.matrix.htobit}

\title{Extractor Methods for htobit Objects}

\description{
  Methods for extracting information from fitted \code{htobit} objects.
}

\usage{
\S3method{coef}{htobit}(object, model = c("full", "location", "scale"), \dots)
\S3method{vcov}{htobit}(object, model = c("full", "location", "scale"), \dots)

\S3method{terms}{htobit}(x, model = c("location", "scale", "full"), \dots)
\S3method{model.matrix}{htobit}(object, model = c("location", "scale"), \dots)
}
\arguments{
  \item{object, x}{an object of class \code{"htobit"}.}
  \item{model}{character indicating (sub)model for which information
    should be extracted.}
  \item{\dots}{currently not used.}
}

\details{
  In addition to the methods above, a set of standard extractor functions for 
  \code{"htobit"} objects is available, see \code{\link{htobit}} for an overview.
}

\seealso{\code{\link{htobit}}}

\keyword{regression}
