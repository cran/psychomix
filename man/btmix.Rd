\name{btmix}
\encoding{latin1}

\alias{btmix}
\alias{FLXMCbtreg}

\title{Finite Mixtures of Bradley-Terry Models}

\description{
  Fit finite mixtures of Bradley-Terry models for paired comparisons data via
  maximum likelihood with the EM algorithm.
}

\usage{
btmix(formula, data, k, subset, weights,
         nrep = 3, cluster = NULL, control = NULL,
         verbose = TRUE, drop = TRUE, unique = FALSE, which = NULL,
         type = c("loglin", "logit"), ref = NULL, undecided = NULL,
         position = NULL, ...)

FLXMCbtreg(formula = . ~ ., type = c("loglin", "logit"), ref = NULL,
  undecided = NULL, position = NULL, ...)
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
    then a single \code{raschmix} object is returned instead of a
    \code{stepRaschmix} object.}
  \item{unique}{A logical; if \code{TRUE}, then \code{unique()} is
    called on the result; for details see
    \code{\link[flexmix]{stepFlexmix}}.}
  \item{which}{number of model to get if \code{k} is a vector of
    integers longer than one. If character, interpreted as
    number of components or name of an information criterion.}
  \item{type}{Character. Should an auxiliary log-linear Poisson model or
  logistic binomial be employed for estimation? The latter is only
  available if not undecided effects are estimated.}
  \item{ref}{Character or numeric. Which object parameter should be the
  reference category, i.e., constrained to zero?}
  \item{undecided}{Logical. Should an undecided parameter be estimated?}
  \item{position}{Logical. Should a position effect be estimated?}
  \item{\dots}{Currently not used.}
}

\details{
  Internally \code{\link[flexmix]{stepFlexmix}}
  is called with suitable arguments to fit the finite mixture model with
  the EM algorithm.

  \code{FLXMCbtreg} is the \code{\link[flexmix]{flexmix}}-driver for
  Bradley-Terry mixture models.
  
  The interface is designed along the same lines as \code{\link{raschmix}}
  which is introduced in detail in Frick et al. (2012). However, the
  \code{btmix} function has not yet been fully tested and may change in
  future versions.  
}

\value{
  Either an object of class \code{"btmix"} containing the best model
  with respect to the log-likelihood (if \code{k} is a scalar) or the
  one selected according to \code{which} (if specified and \code{k} is a
  vector of integers longer than 1) or an object of class
  \code{"stepBTmix"} (if \code{which} is not specified and \code{k} is a
  vector of integers longer than 1).
}

\references{
  Bradley, R.A., and Terry, M.E. (1952). Rank Analysis of Incomplete
  Block Designs. I. The Method of Paired Comparisons. \emph{Biometrika},
  \bold{39}(3/4), 324--345.
  
  Dörr, M. (2011). Bradley Terry Mixture Models: Theory, Implementation
  in R and Validation. Diploma Thesis, Ludwig-Maximilians-Universität
  München.
  
  Frick, H., Strobl, C., Leisch, F., and Zeileis, A. (2012).
  Flexible Rasch Mixture Models with Package psychomix.
  \emph{Journal of Statistical Software}, \bold{48}(7), 1--25.
  \url{http://www.jstatsoft.org/v48/i07/}.

  Grün, B., and Leisch, F. (2008). FlexMix Version 2: Finite Mixtures
  with Concomitant Variables and Varying and Constant Parameters.
  \emph{Journal of Statistical Software}, \bold{28}(4), 1--35.
  \url{http://www.jstatsoft.org/v28/i04/}.

  Leisch, F. (2004). FlexMix: A General Framework for Finite Mixture
  Models and Latent Class Regression in R.
  \emph{Journal of Statistical Software}, \bold{11}(8), 1--18.
  \url{http://www.jstatsoft.org/v11/i08/}.
}

\seealso{
  \code{\link[flexmix]{flexmix}}, \code{\link[flexmix]{stepFlexmix}}
}

\examples{

\donttest{
## Data
data("GermanParties2009", package = "psychotools")

## omit single observation with education = 1
gp <- subset(GermanParties2009, education != "1")
gp$education <- factor(gp$education)

## Bradley-Terry mixture models
suppressWarnings(RNGversion("3.5.0"))
set.seed(1)
## fit models for k = 1, ..., 4 with concomitant variables
cm <- btmix(preference ~ gender + education + age + crisis,
  data = gp, k = 1:4, nrep = 3)

## inspect results
plot(cm)

## select model
cm4 <- getModel(cm, which = "4")

## inspect mixture and effects
library("lattice")
xyplot(cm4)
effectsplot(cm4)
effectsplot(cm4, selection = "education")

## vis effects package directly
if(require("effects")) {
eff4 <- allEffects(cm4)
plot(eff4)
}
}

}

\keyword{paired comparisons}
\keyword{Bradley-Terry model}
\keyword{mixture model}
