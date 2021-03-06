## histogram method reusing flexmix's plot() method
histogram.mptmix <- function(x, data, root = TRUE, ...) {
  if(!missing(data)) warning("argument 'data' is ignored")
  ## just calling plot() within an S3 method does not seem to do the right dispatch
  getMethod("plot", c(x = "flexmix", y = "missing"))(as(x, "flexmix"), root = root, ...)
}

## xyplot method
xyplot.mptmix <- function(x, data,
  component = NULL, plot.type = c("multiple", "single"),
  auto.key = NULL, type = "b", lty = NULL, xlab = "MPT model parameter",
  ylab = "Maximum likelihood estimate",
  panel = NULL, scales = NULL, index = NULL, ...)
{
  ## process data
  y <- parameters(x)
  if(!missing(data)) warning("'data' argument is ignored")
  npars <- nrow(y)
  yref <- 0.5  # mean(y[, 1L])
# lab <- labels(x)
  lab <- gsub("coef\\.", "", rownames(y))

  if(is.character(index)) {
    rownames(y) <- lab
    lab <- index     # reorder labels
    y <- y[index, ]
  }

  ## select items/components
  if(is.null(component)) component <- 1:NCOL(y)
  y <- y[, component, drop = FALSE]

  ## set up auxiliary data.frame
  d <- data.frame(
    y = as.vector(y),
    x = rep(1:npars, length(component)),
    z = factor(rep(component, each = npars),
      levels = component, labels = paste("Comp.", component))
  )
  
  ## graphical arguments
  plot.type <- match.arg(plot.type)
  if(plot.type == "single") {
    f <- y ~ x
    groups <- ~ z
    if(is.null(lty)) lty <- trellis.par.get("superpose.line")$lty
  } else {
    f <- y ~ x | z
    groups <- NULL
    if(is.null(lty)) lty <- trellis.par.get("superpose.line")$lty[2]
  }
  if(is.null(auto.key)) auto.key <- plot.type == "single"
  if(is.null(scales)) scales <- list(x = list(at = 1:npars, alternating = 1,
                                     labels = lab))
  if(is.null(panel)) panel <- function(x, y, ...) {
    panel.xyplot(x, y, ...)
    panel.abline(h = yref, reference = TRUE)
  }

  ## call xyplot() formula method
  xyplot(f, groups = groups, data = d,
    type = type, lty = lty, xlab = xlab, ylab = ylab,
    auto.key = auto.key, scales = scales, panel = panel, ...)
}
