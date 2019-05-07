\name{simRaschmix}
\encoding{latin1}

\alias{simRaschmix}

\title{Simulate Data from Rasch Mixture Models}

\description{
  Generate simulated data from mixtures of Rasch models. The latent
  classes of the mixture can differ regarding size as well as item and
  person parameters of the corresponding Rasch model.
}

\usage{
  simRaschmix(design, extremes = FALSE, attributes = TRUE, \dots)
}

\arguments{
  \item{design}{Type of data generating process. Can be provided as a
  character or a named list. See Details.}
  \item{extremes}{Logical. Should observations with none or all items
    solved be included in the data?}
  \item{attributes}{Logical. Should the true group membership as well as
    true item and person parameters be attached to the data as attributes
    \code{"cluster"}, \code{"difficulty"}, and \code{"ability"}?}
  \item{\dots}{Currently not used.}
}

\details{
  The \code{design} of the data generating process (DGP) can be provided
  in essentially three different ways.
  
  If the \code{design} argument is one of \code{"rost1"}, \code{"rost2"}
  or \code{"rost3"}, responses from the three DGPs introduced in Rost
  (1990) will be drawn.

  Alternatively, the \code{design} can be provided as a named list with
  elements \code{nobs}, \code{weights}, \code{ability}, and
  \code{difficulty}.
  The \code{weights} can be provided in three formats: If provided as a
  vector of probabilities (summing to 1), class membership will be drawn
  with these probabilities. If \code{weights} is a vector of integer weights
  (summing to \code{nobs}, or an integer division thereof), Class sizes
  will be either the weights directly or a multiple thereof. As a third
  alternative, the \code{weights} can be provided as a
  function of the number of observations (\code{nobs}).
  The \code{ability} specification can also be provided in three
  formats: If provided as a matrix of dimension 2xk with mean and
  standard deviation for each of the k clusters, the ability parameters
  are drawn from a normal distribution with the corresponding parameters.
  Second, \code{ability} can be an array of dimension (., 2, k) with abilities
  and corresponding weights/probabilities per cluster. Third, it can
  also be provided as a list of k functions which take the number of
  observations as an argument.
  The specification of the item \code{difficulty} can be
  provided either as a matrix with k columns with the item
  difficulties per cluster or as a matrix with \code{nobs} rows with the item
  difficulties per subject.
  
  As a third option, \code{design} may also be a named list containing a
  vector of \code{ability} parameters and a matrix
  \code{difficulty} of dimension (number of observation x number of items).
}

\value{
  A matrix of item responses with dimension (number of observations x
  number of items).
  If \code{attributes = TRUE}, the matrix has attributes \code{cluster},
  \code{ability}, and \code{difficulty}. The class memberships
  \code{cluster} are only returned when not provided implicitly through
  and a vector of abilties and a difficulty matrix with entries for each
  subject.
}

\references{
  Frick, H., Strobl, C., Leisch, F., and Zeileis, A. (2012).
  Flexible Rasch Mixture Models with Package psychomix.
  \emph{Journal of Statistical Software}, \bold{48}(7), 1--25.
  \url{http://www.jstatsoft.org/v48/i07/}.

  Rost, J. (1990). Rasch Models in Latent Classes: An Integration of Two
  Approaches to Item Analysis.
  \emph{Applied Psychological Measurement}, \bold{14}(3), 271--282.
}

\seealso{
  \code{\link{raschmix}}
}

\examples{
#################
## Rost's DGPs ##
#################

suppressWarnings(RNGversion("3.5.0"))
set.seed(1990)

## DGP 1 with just one latent class 
r1 <- simRaschmix(design = "rost1")
## less than 1800 observations because the extreme scorers have been removed
table(attr(r1, "ability"))
table(rowSums(r1))

## DGP 2 with 2 equally large latent classes
r2 <- simRaschmix(design = "rost2", extreme = TRUE)
## exactly 1800 observations including the extreme scorers
table(attr(r2, "ability"))
table(rowSums(r2))

## DGP 3 with 3 latent classes
r3 <- simRaschmix(design = "rost3")
## item parameters in the three latent classes
attr(r3, "difficulty")


####################################
## flexible specification of DGPs ##
####################################

suppressWarnings(RNGversion("3.5.0"))
set.seed(482)

## number of observations
nobs <- 8

## relative weights
weights <- c(1/4, 3/4)
## exact weights: either summing to nobs or an integer division thereof
weights <- c(2, 6)
weights <- c(1, 3)
## weights as function
## here the result is the same as when specifying relative weights
weights <- function(n) sample(size = n, 1:2, prob = c(1/4, 3/4), replace
= TRUE)

## class 1: only ability level 0
## class 2: normally distributed abilities with mean = 2 and sd = 1
ability <- cbind(c(0, 0), c(2, 1))
## class 1: 3 ability levels (-1, 0, 1); class 2: 2 ability levels (-0.5, 0.5)
## with equal probabilities and frequencies, repectively
ability <- array(c(cbind(-1:1, rep(1/3, 3)), cbind(-1:1/2, c(0.5, 0, 0.5))), 
  dim = c(3, 2, 2))
ability <- array(c(cbind(-1:1, rep(1, 3)), cbind(-1:1/2, c(1, 0, 1))),
  dim = c(3, 2, 2))
## ability as function
ability <- list(
  function(n) rnorm(n, mean = 0, sd = 0.5),
  function(n) sample(c(-0.5, 0.5), size = n, replace = TRUE)
)

## difficulty per latent class
difficulty <- cbind(c(-1,1,rep(0,8)), c(rep(0,8),1,-1))

## simulate data
dat <- simRaschmix(design = list(nobs = nobs, weights = weights,
  ability = ability, difficulty = difficulty))

## inspect attributes and raw scores
table(attr(dat, "cluster"))
hist(attr(dat, "ability"))
barplot(table(rowSums(dat)))
attr(dat, "difficulty")


## specification of DGP only via ability and difficulty
## one vector of abilities of all subjects
ability <- c(rnorm(4, mean = 0, sd = 0.5), sample(c(-0.5, 0.5), size = 4, 
  replace = TRUE))
## difficulty per subject
difficulty <- matrix(c(rep(c(-1,1,rep(0,8)), 4), rep(c(rep(0,8),1,-1), 4)),
  nrow = 8, byrow = TRUE)
## simulate data
dat <- simRaschmix(design = list(ability = ability, difficulty = difficulty))

## inspect attributes and raw scores
hist(attr(dat, "ability"))
barplot(table(rowSums(dat)))
attr(dat, "difficulty")
}

\keyword{item response}
\keyword{Rasch model}
\keyword{mixture model}
\keyword{simulated data}
