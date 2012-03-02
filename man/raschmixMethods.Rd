\name{raschmix-methods}
\docType{methods}
\alias{show,raschmix-method}
\alias{summary,raschmix-method}
\alias{show,summary.raschmix-method}
\alias{parameters,raschmix-method}
\alias{worth,raschmix-method}
\alias{scoreProbs}

\title{Methods for raschmix Objects}
\description{
  Methods for \code{\link{raschmix-class}} objects.
}
\usage{
\S4method{summary}{raschmix}(object, eps = 1e-4, ...)

\S4method{parameters}{raschmix}(object,
       which = c("model", "item", "score", "concomitant"),
       difficulty = TRUE, component = NULL, simplify = TRUE)

\S4method{worth}{raschmix}(object, difficulty = TRUE, component = NULL)

scoreProbs(object, component = NULL, simplify = TRUE)
}
\arguments{
  \item{object}{An object of class \code{"raschmix"}.}
  \item{eps}{Probabilities below this threshold are treated as zero in the
    summary method.}
  \item{which}{Indicates which type of parameters are used. \code{model}
    refers to both item and score parameters, \code{item} and
    \code{score} to their corresponding parameters separately. The
    parameters of the concomitant model are accessed through
    \code{concomitant}.}
  \item{difficulty}{Indicates whether item difficulty (default) or
    easiness parameters are used.}
  \item{component}{Indicates which components are returned. Default is
    all components.}
  \item{simplify}{Should the result be simplified if possible?}
  \item{...}{Currently not used.}
}
\details{
  \code{worth} transforms the item parameters so that the sum over all
  item parameters (within each component) is zero.

  \code{scoreProbs} does not include any aliased parameters if a certain
  raw score is not present in the data.
}

\keyword{methods}
\keyword{worth}
