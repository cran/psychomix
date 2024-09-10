\name{effectsplot}
\alias{effectsplot}
\alias{effectsplot.efflist}
\alias{effectsplot.effpoly}
\alias{allEffects.btmix}
\alias{allEffects.mptmix}
\alias{allEffects.raschmix}
\alias{effect.btmix}
\alias{effect.mptmix}
\alias{effect.raschmix}
\alias{effectsplot,raschmix-method}
\alias{effectsplot,btmix-method}
\alias{effectsplot,mptmix-method}

\title{Effects Displays for Concomitant Variables in Finite Mixture Models}

\description{
  Generic function for visualizing the effects of concomitant variables
  in finite mixture models.
}

\usage{
effectsplot(object, \dots)
}

\arguments{
  \item{object}{Fitted model object.}
  \item{\dots}{Arguments passed to \code{\link[effects]{plot.eff}},
  \code{\link[effects]{plot.efflist}}, or \code{\link[effects]{plot.effpoly}}.}  
}

\details{
  \code{effectsplot} is set up to be both an S3 and S4 generic.
  The idea is that it provides the glue needed to extract the concomitant
  part from a mixture model: First, the concomitant model is
  refitted as a \code{multinom} object or \code{glm} object (in case of
  a mixture with two components). Second, \code{\link[effects]{effect}}
  or \code{\link[effects]{allEffects}} from the \pkg{effects} package
  is called to extract the effects of the concomitants. Third,
  the corresponding \code{plot} methods from the \pkg{effects}
  package create the display.
  
  Currently, this is implemented for \code{\link{raschmix}},
  \code{\link{btmix}}, and \code{\link{mptmix}} objects. The interface is not
  yet fully tested and may change in future versions.
}

\references{
  Fox, J. (2003). Effect Displays in R for Generalised Linear Models.
    \emph{Journal of Statistical Software}, \bold{8}(15), 1--27.
    \doi{10.18637/jss.v008.i15}

  Fox, J., and Hong, J. (2009). Effect Displays in R for Multinomial
    and Proportional-Odds Logit Models: Extensions to the effects Package. 
    \emph{Journal of Statistical Software}, \bold{32}(1), 1--24.
    \doi{10.18637/jss.v032.i01}
}

\seealso{
  \code{\link[effects]{effect}}, \code{\link[effects]{allEffects}},
  \code{\link[nnet]{multinom}}, \code{\link[stats]{glm}}
}

\examples{
\donttest{
## data on party preferences in Germany
## (omit single observation with education = 1)
data("GermanParties2009", package = "psychotools")
gp <- subset(GermanParties2009, education != "1")
gp$education <- factor(gp$education)

## fit Bradley-Terry mixture, see ?btmix for more details
## and a fully-worked example
suppressWarnings(RNGversion("3.5.0"))
set.seed(2)
cm4 <- btmix(preference ~ gender + education + age + crisis, data = gp, k = 4, nrep = 1)

## inspect mixture and effects
library("lattice")
xyplot(cm4)
effectsplot(cm4)

## vis effects package directly
if(require("effects")) {
eff4 <- allEffects(cm4)
plot(eff4)
}
}
}

\keyword{hplot}
\keyword{effects plot}
