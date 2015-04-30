\name{raschmix-plot-method}
\docType{methods}
\alias{plot,raschmix-method}
\alias{plot,raschmix,missing-method}
\alias{histogram.raschmix}
\alias{xyplot.raschmix}

\title{Profile Plot of Item Parameters}
\description{
  The \code{plot} method for \code{\link{raschmix-class}} objects gives a
  base plot of the item parameter profiles in each class. A lettice plot
  of the item parameters is returned by \code{xyplot}. A rootogram or
  histogram of the posterior probabilities is plotted via \code{histogram}.
}
\usage{
\S4method{plot}{raschmix,missing}(x, y, component = NULL, difficulty = TRUE,
                   center = TRUE, index = TRUE, names = TRUE,
                   abbreviate = FALSE, ref = TRUE, col =  "black",
                   refcol = "lightgray", linecol = NULL, lty = 2, cex = 1,
                   pch = 19, type = NULL, ylim = NULL, xlab = "Items",
                   ylab = NULL, legend = TRUE, pos = "topright",
                   srt = 45, adj = c(1.1, 1.1), ...)
\method{histogram}{raschmix}(x, data, root = TRUE, ...)
\method{xyplot}{raschmix}(x, data, component = NULL, item = NULL,
  difficulty = TRUE, plot.type = c("multiple", "single"),
  auto.key = NULL, type = "b", lty = NULL, xlab = "Item", ylab = NULL,
  panel = NULL, scales = NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{"raschmix"}.}
  \item{y}{Not used.}
  \item{component}{A vector indicating which components should be
    plotted.}
  \item{difficulty}{Logical. Should item difficulty parameters be used?}
  \item{center}{Logical. Should the item parameters be centered around 0?}
  \item{index}{Logical. Should the index be used for labelling the items?}
  \item{names}{Either logical or an optional vector of names used for
    labeling of the items.}
  \item{abbreviate}{Logical. Should the labels of the items be abbreviated?}
  \item{ref}{Logical. Should a reference line be drawn?}
  \item{col}{Point color. If \code{col} is a vector, it is
    interpreted as the color of the components respectively. Individual
    coloring within components is possible if col is given as a matrix
    with each column representing one component.}
  \item{refcol}{Color of the reference line.}
  \item{linecol}{Line color. Defaults to the point color.}
  \item{lty, cex, pch, type, ylim, xlab, ylab}{Further standard
    graphical parameters.}
  \item{legend}{Logical. Should a legend be included?}
  \item{pos}{Position of the legend.}
  \item{srt, adj}{Passed on to text() if names = TRUE.}
  \item{\dots}{Further graphical parameters.}
  \item{data}{Ignored.}
  \item{root}{Logical. Should a rootogram be drawn?}
  \item{item}{A vector indicating which items should be plotted.}
  \item{plot.type}{Should the item profiles be drawn in multiple panels
    or a single panel?}
  \item{auto.key, panel, scales}{Further graphical parameters for lattice.}
}
\details{
  For a graphical representation of the item parameter in each class use
  \code{plot} (for a base graph) or \code{xyplot} (for a lattice plot).
  
  For a graphical representation of the quality of the mixture use
  \code{histogram}. For details see \code{\link[flexmix]{plot-methods}}.
}
\references{
  Frick, H., Strobl, C., Leisch, F., and Zeileis, A. (2012).
  Flexible Rasch Mixture Models with Package psychomix.
  \emph{Journal of Statistical Software}, \bold{48}(7), 1--25.
  \url{http://www.jstatsoft.org/v48/i07/}.

  Leisch, F. (2004). FlexMix: A General Framework for Finite Mixture
  Models and Latent Class Regression in R.
  \emph{Journal of Statistical Software}, \bold{11}(8), 1--18.
  \url{http://www.jstatsoft.org/v11/i08/}.

  Leisch, F. (2004). Exploring the Structure of Mixture Model
  Components. In Jaromir Antoch, editor, \emph{Compstat 2004 -- Proceedings in
  Computational Statistics}, pages 1405--1412. Physika Verlag, Heidelberg,
  Germany. ISBN 3-7908-1554-3.
}

\keyword{methods}
\keyword{raschmix-plot}
