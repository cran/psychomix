\name{raschmix-class}
\docType{class}

\alias{raschmix-class}

\title{Class "raschmix"}

\description{A fitted \code{\link{raschmix}} model.}

\section{Slots}{
  \describe{
    \item{\code{model}:}{A \code{FLXMC} object for a Rasch mixture
      model}
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
    \item{\code{df}:}{Total number of
      parameters of the model.}
    \item{\code{components}:}{List describing
      the fitted components using \code{FLXcomponent} objects.}
    \item{\code{formula}:}{Object of class \code{"formula"}.}
    \item{\code{control}:}{Object of class \code{"FLXcontrol"}.}
    \item{\code{call}:}{The function call
      used to create the object.}
    \item{\code{group}:}{Object of class \code{"factor"}.}
    \item{\code{converged}:}{Logical, \code{TRUE} if EM algorithm converged.}
    \item{\code{concomitant}:}{Object of class \code{"FLXP"}..}
    \item{\code{weights}:}{Optional weights of the observations.}

    \item{\code{scores}:}{Type of score model employed.}
    \item{\code{restricted}:}{Logical. Is the score model equal across components?}
    \item{\code{deriv}:}{Type of derivatives used for computing
    gradient and Hessian matrix. Analytical with sum algorithm (\code{"sum"}),
    analytical with difference algorithm (\code{"diff"}, faster but numerically unstable),
    or numerical.}
    \item{\code{extremeScoreProbs}:}{Estimated probability of solving
      either all or no items.}
    \item{\code{rawScoresData}:}{Table of raw scores from the data.}
    \item{\code{flx.call}:}{Internal call to \code{stepFlexmix}}
    \item{\code{nobs}:}{Number of observations without missing values,
    excluding observations with an extreme score.}
    \item{\code{identified.items}:}{Factor indicating which items are identified.}
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
