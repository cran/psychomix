\name{raschmix}
\encoding{latin1}

\alias{raschmix}
\alias{FLXMCrasch}
\alias{simRaschmix}

\title{Finite Mixtures of Rasch Models}

\description{
  Fit finte mixtures of Rasch models for item response data via
  conditional maximum likelihood with the EM algorithm.
}

\usage{
raschmix(formula, data, k, subset, weights, scores = c("saturated", "meanvar"),
         nrep = 3, cluster = NULL, control = list(minprior = 0),
         verbose = TRUE, drop = TRUE, unique = FALSE, which = NULL,
	 gradtol = 1e-6, deriv = "sum", hessian = FALSE, ...)

FLXMCrasch(formula = . ~ ., scores = c("saturated", "meanvar"),
  nonExtremeProb = 1, gradtol = 1e-6, deriv = "sum", hessian = FALSE, ...)

simRaschmix(nobs = 1800, itemp = NULL, mean = NULL, sd = NULL,
            design = c("rost1", "rost2", "rost3", "cont1",
                       "cont1-2", "cont2", "cont2-2"),
            extremes = FALSE, attributes = TRUE)
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
  \item{scores}{Indicates which model should be fitted
    for the score probabilities: either a \code{saturated} model with a
    separate parameter for each score probability, or, for \code{meanvar},
    a multinomial logit model with a location and a scale parameter.}
  \item{nrep}{Number of runs for the starting values for the EM
    algorithm (if \code{cluster = "mrm"}) or number of runs of the EM
    algorithm itself.}
  \item{cluster}{Either a matrix with \code{k} columns of initial
    cluster membership probabilities for each observation; or a factor
    or integer vector with the initial cluster
    assignments of observations at the start of the EM
    algorithm. If \code{cluster = "mrm"}, the \code{\link[mRm]{mrm}} function is
    used to generate starting values. Default is random assignment into \code{k}
    clusters.}
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
  \item{nonExtremeProb}{A numeric giving the probability of scoring
    either none or all items.}
  \item{gradtol, deriv, hessian}{Control parameters passed to
    \code{\link[psychotools]{RaschModel.fit}} for the M-step.}
  \item{nobs}{Number of observations.}
  \item{itemp}{Optional item parameters for designs \code{"cont1"},
    \code{"cont1-2"}, \code{"cont2"}, and \code{"cont2-2"}. Given as a
    \code{matrix} with each column representing one latent class.}
  \item{mean, sd}{Parameters of the normal distribution used to
    generate person parameters.}
  \item{design}{Type of data generating process.}
  \item{extremes}{Logical. Should observations with none or all items
    solved be included in the data?}
  \item{attributes}{Logical. Should the true group membership as well as
    true item and person parameters be attached to the data as attributes
    \code{"group"}, \code{"item"}, and \code{"person"}?}
  \item{\dots}{Currently not used.}
}

\details{
  Internally \code{\link[flexmix]{stepFlexmix}}
  is called with suitable arguments to fit the finite mixture model with
  the EM algorithm.

  \code{FLXMCrasch} is the \code{\link[flexmix]{flexmix}}-driver for
  the Rasch mixture models with saturated score distribution as
  proposed by Rost (1990), also known as \dQuote{Mixed Rasch Model}.

  For the \code{design} argument of \code{simRaschmix}, \code{"rost1"},
  \code{"rost2"}, and \code{"rost3"} refer to the 3 data generating
  processes (dgps) introduced in Rost (1990). The other arguments refer to
  similar dgps with the same sets of item parameters but the person
  parameters drawn from normal distributions.
  \code{"cont1"} is the counterpart to \code{"rost1"} with the person
  parameters drawn from a standard normal distribution. \code{"cont1-2"}
  is similar, but the person parameters stem from two different normal
  distributions with means 2 and -2.
  \code{"cont2"} and \code{"cont2-2"} are the counterparts to
  \code{"rost2"}. For \code{"cont2"} for all observations with the same
  item parameters, the same normal distribution is used to generate
  person parameters. For \code{"cont2-2"}, in each group regarding the
  item parameters extra heterogenity is introduced via sampling the
  person parameters from two different normal distributions (with means
  2 and -2), similar to dgp \code{"cont1-2"}.
  
}

\value{
  Either an object of class \code{"raschmix"} containing the best model
  with respect to the log-likelihood (if \code{k} is a scalar) or the
  one selected according to \code{which} (if specified and \code{k} is a
  vector of integers longer than 1) or an object of class
  \code{"stepRaschmix"} (if \code{which} is not specified and \code{k} is a
  vector of integers longer than 1).
}

\references{
  Gr�n, B., and Leisch, F. (2008). FlexMix Version 2: Finite Mixtures
  with Concomitant Variables and Varying and Constant Parameters.
  \emph{Journal of Statistical Software}, \bold{28}(4), 1--35.
  \url{http://www.jstatsoft.org/v28/i04/}.

  Leisch, F. (2004). FlexMix: A General Framework for Finite Mixture
  Models and Latent Class Regression in R.
  \emph{Journal of Statistical Software}, \bold{11}(8), 1--18.
  \url{http://www.jstatsoft.org/v11/i08/}.

  Rost, J. (1990). Rasch Models in Latent Classes: An Integration of Two
  Approaches to Item Analysis.
  \emph{Applied Psychological Measurement}, \bold{14}(3), 271--282.

  Rost, J., and von Davier, M. (1995). Mixture Distribution Rasch Models.
  In Fischer, G.H., and Molenaar, I.W. (eds.),
  \emph{Rasch Models: Foundations, Recent Developments, and Applications},
  chapter 14, pp. 257--268. Springer-Verlag, New York.
}

\seealso{
  \code{\link[flexmix]{flexmix}}, \code{\link[flexmix]{stepFlexmix}}
}

\examples{
##########
## Data ##
##########

## simulate response from Rost's scenario 2
set.seed(1)
r2 <- simRaschmix(design = "rost2")

## plus informative and noninformative concomitants
d <- data.frame(
  x1 = rbinom(nrow(r2), prob = c(0.4, 0.6)[attr(r2, "group")], size = 1),
  x2 = rnorm(nrow(r2))
)
d$resp <- r2


####################################################
## Rasch mixture model with saturated score model ##
## (Rost, 1990)                                   ##
####################################################

## fit models for k = 1, 2, 3
m1 <- raschmix(r2, k = 1:3, score = "saturated")
## equivalently: m1 <- raschmix(resp ~ 1, data = d, k = 1:3, score = "saturated")

## inspect results
m1
plot(m1)

## select best BIC model
BIC(m1)
m1b <- getModel(m1, which = "BIC")
summary(m1b)

## compare estimated with true item parameters
parameters(m1b, "item") ##  9 items, item_1 = 0
worth(m1b)              ## 10 items, sum = 0
attr(r2, "item")

## graphical comparison
plot(m1b, pos = "top")
for(i in 1:2) lines(attr(r2, "item")[,i], lty = 2, type = "b")

## extract estimated raw score probabilities
## (approximately equal across components and roughly uniform)
scoreProbs(m1b)

## note: parameters() and worth() take "component" argument
parameters(m1b, "item",  component = 2)
parameters(m1b, "score", component = 1)
worth(m1b, component = 2:1)

## inspect posterior probabilities
histogram(m1b)
head(posterior(m1b)) ## for first observations only

## compare resulting clusters with true groups
table(model = clusters(m1b), true = attr(r2, "group"))

## optionally: leverage mRm package for faster computation of
## starting values
\dontrun{
library("mRm")
## fit 2-component model
m1b_mrm <- raschmix(r2, k = 2, score = "saturated", cluster = "mrm")
## essentially identical to previous solution
table(clusters(m1b), clusters(m1b_mrm))   
worth(m1b) - worth(m1b_mrm)
}

################################################################
##  Rasch mixture model with mean/variance score distribution ##
## (Rost & von Davier, 1995)                                  ##
################################################################

## more parsimonious parametrization,
## fit multinomial logit model for score probabilites

## fit models and select best BIC
m2 <- raschmix(r2, k = 1:3, score = "meanvar")
plot(m2)
m2b <- getModel(m2, which = "BIC")

## compare number of estimated parameters
dim(parameters(m2b)) 
dim(parameters(m1b)) 

## graphical comparison with true parameters
plot(m2b, pos = "top")
for(i in 1:2) lines(attr(r2, "item")[,i], lty = 2, type = "b")

## results from non-parametric and parametric specification
## essentially identical
max(abs(worth(m1b) - worth(m2b, component = 2:1)))


###########################
## Concomitant variables ##
###########################

## employ concomitant variables (x1 = informative, x2 = not)
\dontrun{
## fit model
cm2 <- raschmix(resp ~ x1 + x2, data = d, k = 2:3, score = "meanvar")

## BIC selection
rbind(m2 = BIC(m2), cm2 = c(NA, BIC(cm2)))
cm2b <- getModel(cm2, which = "BIC")

## concomitant coefficients
parameters(cm2b, which = "concomitant")
}


##########
## Misc ##
##########

## note: number of clusters can either be chosen directly
## or directly selected via AIC (or BIC, ICL)
\dontrun{
raschmix(r2, k = 2)
raschmix(r2, k = 1:3, which = "AIC")
}
}

\keyword{item response}
\keyword{Rasch model}
\keyword{mixture model}