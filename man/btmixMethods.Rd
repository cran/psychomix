\name{btmix-methods}
\docType{methods}
\alias{worth,btmix-method}
% \alias{labels<-}
% \alias{mscale}
% \alias{mscale<-}


\title{Methods for btmix Objects}
\description{
  Methods for \code{\link{btmix-class}} objects.
}
\usage{
\S4method{worth}{btmix}(object, component = NULL)

% \method{labels}{btmix}(object, ...)
% \method{mscale}{btmix}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{"btmix"}.}
  \item{component}{Indicates which components are returned. Default is
    all components.}
  % \item{...}{Currently not used.}
}
\details{
  \code{worth} returns the worth parameters from the Bradley-Terry model.

  % \code{labels} returns the labels used for the objects which are compared.

  % \code{mscale} returns the measurement scale for the paired comparisons data.
}

\keyword{methods}
\keyword{worth}
