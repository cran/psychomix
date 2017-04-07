\name{mptmix-class}
\docType{class}

\alias{mptmix-class}

\title{Class "mptmix"}

\description{A fitted \code{\link{mptmix}} model.}

\section{Slots}{
  \describe{
    \item{\code{model}:}{A \code{FLXMC} object for an MPT mixture model}
    \item{\code{prior}:}{Numeric vector with prior probabilities of classes.}
    \item{\code{posterior}:}{Named list with elements \code{scaled}
      and \code{unscaled}, both matrices with one row per observation
      and one column per class.}
    \item{\code{iter}:}{Number of EM iterations.}
    \item{\code{k}:}{Number of classes after EM.}
    \item{\code{k0}:}{Number of classes at start of EM.}
    \item{\code{cluster}:}{Class assignments of observations.}
    \item{\code{size}:}{Class sizes.}
    \item{\code{logLik}:}{Log-likelihood at EM convergence.}
    \item{\code{df}:}{Total number of parameters of the model.}
    \item{\code{components}:}{List describing
      the fitted components using \code{FLXcomponent} objects.}
    \item{\code{formula}:}{Object of class \code{"formula"}.}
    \item{\code{control}:}{Object of class \code{"FLXcontrol"}.}
    \item{\code{call}:}{The function call used to create the object.}
    \item{\code{group}:}{Object of class \code{"factor"}.}
    \item{\code{converged}:}{Logical, \code{TRUE} if EM algorithm converged.}
    \item{\code{concomitant}:}{Object of class \code{"FLXP"}.}
    \item{\code{weights}:}{Optional weights of the observations.}

    \item{\code{flx.call}:}{Internal call to \code{stepFlexmix}}
    \item{\code{nobs}:}{Number of observations.}
  % \item{\code{labels}:}{Labels of objects compared.}
  }
}
\section{Extends}{
  Class \code{flexmix}, directly.
}

\section{Accessor Functions}{
  The following functions should be used for accessing the corresponding
  slots:
  \describe{
    \item{\code{clusters}:}{Cluster assignments of observations.}
    \item{\code{posterior}:}{A matrix of posterior probabilities for each
      observation.}
  }
}

\keyword{classes}
