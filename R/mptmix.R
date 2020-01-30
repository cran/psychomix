## interface 
mptmix <- function(formula, data, k, subset, weights,
  nrep = 3, cluster = NULL, control = NULL,
  verbose = TRUE, drop = TRUE, unique = FALSE, which = NULL,
  spec, treeid = NULL,
  optimargs = list(control = list(reltol = .Machine$double.eps^(1/1.2),
                                  maxit = 1000)),
  ...)
{
  ## process call
  cl <- match.call()
  has_subset <- !missing(subset)
  has_weights <- !missing(weights)

  ## arrange formula and data arguments correctly
  if(missing(formula)) {
    d <- as.data.frame(matrix(0, nrow = nrow(data), ncol = 0L))
    d$.response <- as.matrix(data)
    conc <- NULL 
    respfreq <- data
  } else if(!inherits(formula, "formula") & missing(data)) {
    d <- as.data.frame(matrix(0, nrow = nrow(formula), ncol = 0L))
    d$.response <- as.matrix(formula)
    conc <- NULL
    respfreq <- formula
  } else {
  
    ## process call
    if(missing(data)) data <- environment(formula)
    aux <- match.call(expand.dots = FALSE)
    aux <- aux[c(1L, match(c("formula", "data", "subset", "weights"), names(aux), 0L))]
    
    ## process formula via Formula
    ff <- Formula(formula)

    ## handle concomitants first
    auxc <- aux
    auxc[[1]] <- as.name("get_all_vars")
    auxc$formula <- formula(ff, lhs = 0, rhs = NULL)
    if(has_subset) names(auxc)[names(auxc) == "subset"] <- ".subset"
    if(has_weights) names(auxc)[names(auxc) == "weights"] <- ".weights"
    d <- eval(auxc, parent.frame())
    if(has_subset) d <- d[d$.subset, , drop = FALSE]
    d$.subset <- NULL
    conc <- ncol(d) > has_weights
    conc <- if(conc) {
      f <- formula(ff, lhs = 0, rhs = NULL, collapse = TRUE)
      FLXPmultinom(f)
    } else {
      NULL
    }
    
    ## handle response
    auxr <- aux
    auxr[[1]] <- as.name("model.frame")
    auxr$formula <- ff
    auxr$lhs <- NULL
    auxr$rhs <- 0
    auxr$na.action <- na.pass
    auxr$weights <- NULL
    mf <- eval(auxr, parent.frame())
    respfreq <- model.part(ff, lhs = NULL, rhs = 0, data = mf)[[1L]]
    d$.response <- as.matrix(respfreq)
  }
  
  ## driver function
  model <- FLXMCmpt(spec = spec, treeid = treeid, optimargs = optimargs)

  ## control parameters
  ctrl <- as(control, "FLXcontrol")  
    
  ## fitting
  z <- if(has_weights) {
    stepFlexmix(.response ~ 1, data = d, k = k, nrep = nrep,
                cluster = cluster, model = model,
                concomitant = conc, control = ctrl, verbose = verbose,
                drop = drop, unique = unique, weights = d$.weights, ...)
  } else {
    stepFlexmix(.response ~ 1, data = d, k = k, nrep = nrep,
                cluster = cluster, model = model,
                concomitant = conc, control = ctrl, verbose = verbose,
                drop = drop, unique = unique, ...)
  }

  ## select model (if desired)
  if (!is.null(which) & inherits(z, "stepFlexmix")){
    z <- getModel(z, which = which)
  }
  
  ## classify
  if (inherits(z, "flexmix")) {
    z <- as(z, "mptmix")
    z@flx.call <- z@call
    z@call <- cl
  # z@labels <- labels(pc)  # possibly parameter names?
  # z@ref <- ref
  # z@type <- match.arg(type, c("loglin", "logit"))
  }
  else {
    z@models <- lapply(z@models, function(model){
      model <- as(model, "mptmix")
      model@flx.call <- model@call
      model@call <- cl
      model@call$k <- model@flx.call$k
    # model@labels <- labels(pc)
    # model@ref <- ref
    # model@type <- match.arg(type, c("loglin", "logit"))
      return(model)
    })
    z <- as(z, "stepMPTmix")
    z@flx.call <- z@call
    z@call <- cl
  # z@labels <- labels(pc)
  }

  return(z)
}


## Flexmix driver for MPT mixture model
FLXMCmpt <- function(formula = . ~ ., spec = NULL, treeid = NULL,
  optimargs = NULL, ...)
{
  retval <- new("FLXMCmpt", weighted = TRUE, formula = formula,
                name = "MPT mixture model")

  retval@defineComponent <- expression({

    logLik <- function(x,y,...) {
      loglikfun_mpt(coef, y, spec, ...)
    }

    new("FLXcomponent", df = df , logLik = logLik,
        parameters = list(coef = coef))
  })

  retval@fit <- function(x,y,w, ...){
    mptmod <- mptmodel(y, weights = w, spec = spec, treeid = treeid,
                       optimargs = optimargs, ...)

    para <- list(coef = coef(mptmod), df = mptmod$npar)

    with(para, eval(retval@defineComponent))
  }
  retval@dist <- "MPT"
  retval
}

## compute individual contributions to log-likelihood function
loglikfun_mpt <- function(par, data, spec, ...)
{
  rval <- drop(data %*% log(spec$par2prob(par)))

  ## return result
  return(rval)
}
