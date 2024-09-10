\name{mptmix}

\alias{mptmix}
\alias{FLXMCmpt}

\title{Finite Mixtures of Multinomial Processing Tree Models}

\description{
  Fit finite mixtures of multinomial processing tree (MPT) models via maximum
  likelihood with the EM algorithm.
}

\usage{
mptmix(formula, data, k, subset, weights,
          nrep = 3, cluster = NULL, control = NULL,
          verbose = TRUE, drop = TRUE, unique = FALSE, which = NULL,
          spec, treeid = NULL,
          optimargs = list(control = list(reltol =
            .Machine$double.eps^(1/1.2), maxit = 1000)), ...)

FLXMCmpt(formula = . ~ ., spec = NULL, treeid = NULL, optimargs = NULL, ...)
}

\arguments{
  \item{formula}{Symbolic description of the model (of type \code{y ~ 1}
    or \code{y ~ x}).}
  \item{data, subset}{Arguments controlling formula processing.}
  \item{k}{A vector of integers indicating the number of components of
    the finite mixture; passed in turn to the \code{k} argument
    of \code{\link[flexmix]{stepFlexmix}}.}
  \item{weights}{An optional vector of weights to be used in the fitting
    process; passed in turn to the \code{weights} argument of
    \code{\link[flexmix]{flexmix}}.}
  \item{nrep}{Number of runs of the EM algorithm.}
  \item{cluster}{Either a matrix with \code{k} columns of initial
    cluster membership probabilities for each observation; or a factor
    or integer vector with the initial cluster assignments of
    observations at the start of the EM algorithm. Default is random
    assignment into \code{k} clusters.}
  \item{control}{An object of class \code{"FLXcontrol"} or a named list;
    controls the EM algorithm and passed in turn to the \code{control}
    argument of \code{\link[flexmix]{flexmix}}.}
  \item{verbose}{A logical; if \code{TRUE} progress information is shown
    for different starts of the EM algorithm.}
  \item{drop}{A logical; if \code{TRUE} and \code{k} is of length 1,
    then a single \code{mptmix} object is returned instead of a
    \code{stepMPTmix} object.}
  \item{unique}{A logical; if \code{TRUE}, then \code{unique()} is
    called on the result; for details see
    \code{\link[flexmix]{stepFlexmix}}.}
  \item{which}{number of model to get if \code{k} is a vector of
    integers longer than one. If character, interpreted as
    number of components or name of an information criterion.}
  \item{spec, treeid, optimargs}{arguments for the MPT model passed on to
    \code{\link[psychotools]{mptmodel}}.}
  \item{\dots}{Currently not used.}
}

\details{
  Internally \code{\link[flexmix]{stepFlexmix}}
  is called with suitable arguments to fit the finite mixture model with
  the EM algorithm.

  \code{FLXMCmpt} is the \code{\link[flexmix]{flexmix}} driver for
  MPT mixture models.
  
  The interface is designed along the same lines as \code{\link{raschmix}}
  which is introduced in detail in Frick et al. (2012). However, the
  \code{mptmix} function has not yet been fully tested and may change in
  future versions.  

  The latent-class MPT model (Klauer, 2006) is equivalent to an MPT mixture
  model without concomitant variables.

  MPT models are specified using the \code{mptspec} function. See the
  documentation in the \pkg{mpt} package for details.
}

\value{
  Either an object of class \code{"mptmix"} containing the best model
  with respect to the log-likelihood (if \code{k} is a scalar) or the
  one selected according to \code{which} (if specified and \code{k} is a
  vector of integers longer than 1) or an object of class
  \code{"stepMPTmix"} (if \code{which} is not specified and \code{k} is a
  vector of integers longer than 1).
}

\references{
  Frick, H., Strobl, C., Leisch, F., and Zeileis, A. (2012).
  Flexible Rasch Mixture Models with Package psychomix.
  \emph{Journal of Statistical Software}, \bold{48}(7), 1--25.
  \doi{10.18637/jss.v048.i07}

  Klauer, K.C. (2006).
  Hierarchical Multinomial Processing Tree Models: A Latent-Class Approach.
  \emph{Psychometrika}, \bold{71}, 7--31.
  \doi{10.1007/s11336-004-1188-3}
}

\seealso{
  \code{\link[flexmix]{flexmix}}, \code{\link[flexmix]{stepFlexmix}}
}

\examples{

\donttest{
## Data
data("PairClustering", package = "psychotools")
pc <- reshape(PairClustering, timevar = "trial", idvar = "ID",
              direction = "wide")

## Latent-class MPT model (Klauer, 2006)
suppressWarnings(RNGversion("3.5.0"))
set.seed(1)
m <- mptmix(as.matrix(pc[-1]) ~ 1, data = pc, k = 1:3,
            spec = mptspec("SR", .replicates = 2))
m1 <- getModel(m, which = "BIC")

## Inspect results
summary(m1)
parameters(m1)
plot(m1)

library(lattice)
xyplot(m1)
}

}

\keyword{mixture model}
