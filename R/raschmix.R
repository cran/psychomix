## interface 
raschmix <- function(formula, data, k, subset, weights,
                     scores = c("saturated", "meanvar"),
                     nrep = 3, cluster = NULL, control = NULL,
                     verbose = TRUE, drop = TRUE, unique = FALSE, which = NULL,
		     gradtol = 1e-6, deriv = "sum", hessian = FALSE,
		     model = NULL, ...){  
  ## process call
  cl <- match.call()
  has_subset <- !missing(subset)
  has_weights <- !missing(weights)

  ## arrange formula and data arguments correctly
  if(missing(formula)) {
    d <- as.data.frame(matrix(0, nrow = nrow(data), ncol = 0L))
    d$.response <- data
    conc <- NULL 
  } else if(!inherits(formula, "formula") & missing(data)) {
    d <- as.data.frame(matrix(0, nrow = nrow(formula), ncol = 0L))
    d$.response <- formula
    conc <- NULL
  } else {
  
    ## process call
    if(missing(data)) data <- environment(formula)
    aux <- match.call(expand.dots = FALSE)
    aux <- aux[c(1L, match(c("formula", "data", "subset", "weights"), names(aux), 0L))]
    
    ## process formula via Formula
    ff <- Formula(formula)

    ## handle conomitants first
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
    d$.response <- as.matrix(model.part(ff, lhs = NULL, rhs = 0, data = mf))
  }
  
  ## data processing: remove observations without any item responses
  ## class(d$.response) <- c("itemresp", "matrix")
  ## missing.obs <- is.na(d$.response)
  missing.obs <- apply(is.na(d$.response), 1, all)
  d <- d[!missing.obs, , drop = FALSE]

  ## data processing: remove observations with any NA in concomitant variables
  #conc.na <- apply(is.na(d[-which(".response" == names(d))]), 2, any)
  #d <- d[-conc.na, , drop = FALSE]
  conc.cc <- complete.cases(d[-which(".response" == names(d))])
  d <- d[conc.cc, , drop = FALSE]

  ## data processing: remove observations with any NA (required for flexmix)
  # (removing only observations with missings in covariates is not enough)
  #cc <- complete.cases(d)
  #d <- d[cc, , drop = FALSE]

  ## data processing: non-identified items
  n.total <- nrow(d$.response)
  cm <- colMeans(d$.response, na.rm = TRUE)
  status <- as.character(cut(cm, c(-Inf, 1/(2 * n.total), 1 - 1/(2 * n.total),
                                   Inf), labels = c("0", "0/1", "1")))
  status[is.na(status)] <- "NA"
  status <- factor(status, levels = c("0/1", "0", "1", "NA"))
  if(is.null(colnames(d$.response)))
    colnames(d$.response) <- paste("Item", gsub(" ", "0",
                                  format(1:ncol(d$.response))), sep = "")
  names(status) <- colnames(d$.response)
  ident <- status == "0/1"
  d$.response <- d$.response[,ident, drop = FALSE]
 
  
  ## data processing: extreme scorers
  ## extreme scorer = subjects who score all non-missing items with either 0 or 1
  score.0 <- rowSums(d$.response, na.rm = TRUE) == 0
  n.0 <- sum(score.0)
  d <- d[!score.0, , drop = FALSE]
  # <FIXME> for observations with NA in item responses: extreme scorer definition
  # of scoring all (non-missing) items needs to be implemented </FIXME>
  ## score.m <- rowSums(d$.response, na.rm = TRUE) == ncol(d$.response)
  score.m <- (rowSums(d$.response, na.rm = TRUE) + rowSums(is.na(d$.response))) == ncol(d$.response)
  n.m <- sum(score.m)
  d <- d[!score.m, , drop = FALSE]

  pi.0 <- n.0 / n.total
  pi.m <- n.m / n.total
  pi.nonex <- nrow(d$.response) / n.total

  ## reference category for saturated score model
  rs <- factor(rowSums(d$.response, na.rm = TRUE), levels = 1:(ncol(d$.response)-1L))
  rs <- table(rs)
  ref <- min(which(rs > 0))

  ## Rasch driver including score model
  scores <- match.arg(head(tolower(scores), 1L), c("saturated", "meanvar", "constant"))
  if(is.null(model)) model <- FLXMCrasch(scores = scores, nonExtremeProb = pi.nonex, ref = ref,
    gradtol = gradtol, deriv = deriv, hessian = hessian)

  ## control parameters
  ctrl <- as(control, "FLXcontrol")  
    
  ## starting values via mrm() from "mRm"
  if(identical(cluster, "mrm")) {
    if(!require("mRm")) stop("starting values from mrm() require package 'mRm'")
    if(length(k) > 1L) warning("starting values from mrm() can only be used with a single 'k', first used")
    k <- k[1L]
    cluster <- if(nrep > 1L) {
      fits <- lapply(1:nrep, function(i) mrm(d$.response, k))
      ix <- which.max(sapply(fits, "[[", "logLik"))
      mrm_clusters(d$.response, k, fits[[ix]])
    } else {
      mrm_clusters(d$.response, k)
    }
    nrep <- 1L
  }

  ## fitting
  z <- if(has_weights) {
    stepFlexmix(.response ~ 1, data = d, k = k, nrep = nrep, cluster = cluster, model = model,
                concomitant = conc, control = ctrl, verbose = verbose,
                drop = drop, unique = unique, weights = d$.weights, ...)
  } else {
    stepFlexmix(.response ~ 1, data = d, k = k, nrep = nrep, cluster = cluster, model = model,
                concomitant = conc, control = ctrl, verbose = verbose,
                drop = drop, unique = unique, ...)
  }

  ## select model (if desired)
  if (!is.null(which) & class(z) == "stepFlexmix"){
    z <- getModel(z, which = which)
  }
  
  ## classify
  ## for Rost model: adjust df, logLik
  if (class(z) == "flexmix") { ## is(z, "flexmix")?
    z@df <- z@df + (pi.0 != 0) + (pi.m != 0)
    if (pi.0 != 0) z@logLik <- z@logLik + n.0*log(pi.0)
    if (pi.m != 0) z@logLik <- z@logLik + n.m*log(pi.m)
    z <- as(z, "raschmix")
    z@scores <- scores
    z@deriv <- deriv
    z@extremeScoreProbs <- c(pi.0, pi.m)
    z@rawScoresData <- rs
    z@flx.call <- z@call
    z@call <- cl
    z@nobs <-  n.total - n.0 - n.m # includes only those observations passed on to stepFlexmix
    z@identified.items <- status
  }
  else {
    if (pi.0 != 0) z@logLiks <- z@logLiks + n.0*log(pi.0)
    if (pi.m != 0) z@logLiks <- z@logLiks + n.m*log(pi.m)
    z@models <- lapply(z@models, function(model){
      model@df <- model@df + (pi.0 != 0) + (pi.m != 0)
      if (pi.0 != 0) model@logLik <- model@logLik + n.0*log(pi.0)
      if (pi.m != 0) model@logLik <- model@logLik + n.m*log(pi.m)
      model <- as(model, "raschmix")
      model@scores <- scores
      model@deriv <- deriv
      model@extremeScoreProbs <- c(pi.0, pi.m)
      model@rawScoresData <- rs
      model@flx.call <- model@call
      model@call <- cl
      model@call$k <- model@flx.call$k
      model@nobs <- n.total - n.0 - n.m # includes only those observations passed on to stepFlexmix
      model@identified.items <- status
      return(model)
    })
    z <- as(z, "stepRaschmix")
    z@flx.call <- z@call
    z@call <- cl
  }

  return(z)
}



## Flexmix driver for Rasch mixture model
FLXMCrasch <- function(formula = . ~ ., scores = "saturated",
  nonExtremeProb = 1, ref = 1, gradtol = 1e-6, deriv = "sum", hessian = FALSE, ...)
{
  scores <- match.arg(scores, c("saturated", "meanvar", "constant"))

  retval <- new("FLXMCrasch", weighted = TRUE, formula = formula, 
                name = sprintf("Rasch mixture model (%s scores)", scores))

  retval@defineComponent <- expression({

    ## include score probabilities pi_rk
    logLik <- function(x,y,...) {
      loglikfun_rasch(item, score, y, deriv = deriv, scores = scores, nonExtremeProb = nonExtremeProb, ...)
    }

    new("FLXcomponent", df = df , logLik = logLik,
        parameters = list(item = item, score = score))
  })

  retval@fit <- function(x,y,w, ...){
    ## set weights to zero if computationally zero
    ## (FIXME: should this be checked by flexmix?)
    w <- ifelse(w < .Machine$double.eps^(1/1.5), 0, w)
    
    ## estimate parameters
    rasch.model <- RaschModel.fit(y, weights = w, gradtol = gradtol,
      deriv = deriv, hessian = hessian, ...)

    # number of items and observations
    m <- ncol(y)
    #nk <- sum(w) ## obsolete

    ## raw scores and associated number of observations
    ## (some of which may be zero)
    rs <- as.integer(rowSums(y))
    nrk <- tapply(w, factor(rs, levels = 1:(m-1L)), sum)
    nrk[is.na(nrk)] <- 0

    ## raw score probabilities via logit model
    switch(scores,    
      "saturated" = {
        #fix <- min(which(nrk > 0))
        xaux <- diag(m - 1L)
        gamma <- log(nrk) - log(nrk[ref])
      },
      "meanvar" = {
        rr <- as.numeric(names(nrk))
        xaux <- cbind(rr / m, 4 * rr * (m - rr) / m^2)
        start <- numeric(2)
        clm <- function(par) {
          eta <- drop(xaux %*% par)
          psi <- exp(eta) / sum(exp(eta))
          -sum(nrk * log(psi))
        }
        gamma <- optim(start, clm)$par
	ref <- NULL
      },      
      "constant" = {
        xaux <- matrix(1, nrow = m - 1L, ncol = 1L)
	gamma <- 0
	ref <- 1
      }
    )

    ## probabilities (colSums instead of %*% to get na.rm = TRUE)
    eta <- colSums(as.vector(gamma) * t(xaux), na.rm = TRUE)
    psi <- exp(eta) / sum(exp(eta))

    ## drop fixed parameters
    if(!is.null(ref)) gamma <- gamma[-ref]
        
    ## any aliased parameters?
    alias <- if(length(gamma) > 0L) !is.finite(gamma) | is.na(gamma) else logical(0L)
    
    ## check for remaining problems
    if(any(is.na(psi))) stop("some parameters in score model not identified")

    # conditional MLE -> r = 0 and r = "number of items" are excluded 
    psi <- psi * nonExtremeProb

    ## df and score parameters
    score <- gamma[!alias]
    df <- rasch.model$df + length(score)
    
    ## collect elements necessary for loglikfun    
    para <- list(item = rasch.model$coef, score = score, df = df)

    with(para, eval(retval@defineComponent))
  }
  retval@dist <- "Rasch"
  retval
}


## compute individual contributions to log-likelihood function
loglikfun_rasch <- function(item, score, y, deriv, scores, nonExtremeProb, ...) {

  ## in case of non-identified parameters return log(0)  
  ##if(any(object$items != "0/1")) return(rep(-Inf, NROW(y)))
  n <- nrow(y)
  #cm <- colMeans(y, na.rm = TRUE)
  #status <- as.character(cut(cm, c(-Inf, 1/(2 * n), 1 - 1/(2 * n), Inf), labels = c("0", "0/1", "1")))
  #status[is.na(status)] <- "NA"
  #status <- factor(status, levels = c("0/1", "0", "1", "NA"))
  #ident <- status == "0/1"
  #names(status) <- colnames(y)
  #if(any(status != "0/1")) return(rep(-Inf, NROW(y)))

  ## missing values?
  any_y_na <- any(is.na(y))
 
  ## in case of NAs
  if(any_y_na) {
    ## compute all NA patterns
    na_patterns <- factor(apply(is.na(y), 1,
                                function(z) paste(which(z), collapse = "\r")))
      
    ## replace NAs
    y[is.na(y)] <- 0
  }
  ## associated elementary symmetric functions (of order 0)
  esf <- if(any_y_na) {
    lgamma <- rep(0, NROW(y))
    for(i in levels(na_patterns)) {
      wi1 <- which(na_patterns == i)
      wi_i <- as.integer(strsplit(i, "\r")[[1]])
      par_i <- if(length(wi_i) < 1) c(0, par)
               else c(0, par)[-wi_i]
      gamma_i <- elementary_symmetric_functions(par_i,
                                                order = 1 - (deriv == "numeric"), 
                                                diff = deriv == "diff")[[1]][-1]
      lgamma[wi1] <- log(gamma_i)[rowSums(y[wi1, , drop = FALSE])]
    }
    lgamma
  } else {
    gamma <- elementary_symmetric_functions(c(0, item),
                                            order = 1 - (deriv == "numeric"), 
                                            diff = deriv == "diff")[[1]][-1]
    log(gamma)[rowSums(y)]
  }

  
  ## estimated parameters (including contrast)
  par <- c(0, item)
   
  ## log-likelihood contribution
  ## = parameter sum on log scale - log(sum of all permutations)
  rval <- -drop(y %*% par) - esf

  ## for the original model by Rost (1990) 
  ## log(probability of row sum r) added to log-likelihood
  ## number of items and observations
  m <- sum(colSums(y) > 0 & colSums(y) < nrow(y))
  rs <- which(tabulate(rowSums(y), ncol(y)-1L) > 0) # FIXME: use ncol(y) - 1L or m?
  nscores <- length(rs)
  #if (nscores != (length(score)+1)) browser()
  switch(scores,
    "saturated" = {
      #xaux <- diag(nscores)
      xaux <- diag(length(score) + 1)
      extra <- 0
    },
    "meanvar" = {
      rr <- 1:(m-1L)
      xaux <- cbind(rr / m, 4 * rr * (m - rr) / m^2)
      extra <- NULL
    },      
    "constant" = {
      xaux <- matrix(1, nrow = m - 1L, ncol = 1L)
      extra <- 0
    }
    )

  ## calculate score probabilities
  eta <- colSums(c(extra, score) * t(xaux), na.rm = TRUE)
  psi <- exp(eta) / sum(exp(eta))
    
  ## conditional MLE -> r = 0 and r = "number of items" are excluded 
  psi <- psi * nonExtremeProb

  ## check for remaining problems
  if(any(is.na(psi))) stop("some parameters in score model not identified")

  ## for saturated model:
  if (scores == "saturated"){
    ## return probability 0 for scores not present in the data
    psi.full <- numeric(m - 1L)
    psi.full[c(rs[1], as.numeric(names(score)))] <- psi
    psi <- psi.full
  }

  ## FIXME: @Achim: stimmt das? Ist es richtig, fuer's constant gar nichts eingeht?
  if (scores == "constant") psi <- NULL
  
  if (!is.null(psi)) rval <- log(psi)[rowSums(y)] + rval
  ## Alternative: employ empirical distribution of the raw scores
  #if (scores == "constant") psi <- tabulate(rowSums(y), 1:(m-1))
  #rval <- log(psi)[rowSums(y)] - drop(y %*% par) - esf
  ## dann auch einfach rval in einem Schritt ausrechnen.

  ## replace -Inf with smallest (machine) value possible
  ## to keep flexmix from stopping
  rval[!is.finite(rval)] <- .Machine$double.xmin
  
  ## return result
  return(rval)
}

## methods for S3 class "itemresp"
is.na.itemresp <- function(x) apply(NextMethod("is.na", x), 1, all)
"[.itemresp" <- function(x, ...) {
  rval <- NextMethod()
  class(rval) <- c("itemresp", class(rval))
  rval
} 

mrm_clusters <- function(y, k, fit = NULL)
{
  ## need mRm package
  stopifnot(require("mRm"))

  ## data
  y <- y[, colSums(y) > 0 & colSums(y) < nrow(y), drop = FALSE]
  m <- ncol(y)
  stopifnot(m > 0)
  y <- na.omit(y)
  y <- y[rowSums(y) > 0 & rowSums(y) < m, ]

  ## raw score frequencies
  rs <- rowSums(y)
  
  ## fit mrm mixture model
  if(is.null(fit)) fit <- mrm(y, k)

  ## extract list of coefficients in RaschModel.fit parametrization
  beta <- lapply(1:k, function(i) {
    b <- fit$beta[,i]
    -b[-1L] + b[1L]
  })
    
  ## extract raw score probabilities
  psi <- lapply(1:k, function(i) fit$pi.r.c[,i])

  ## elementary symmetric functions
  lgamma <- lapply(beta, function(x) log(psychotools:::elementary_symmetric_functions(c(0, x), order = 0)[[1L]][-1L]))

  ## prior probabilties
  prior <- drop(fit$class.size)
  
  ll <- sapply(1:k, function(i) -drop(y %*% c(0, beta[[i]])) - lgamma[[i]][rs] + log(psi[[i]][rs]))
  dens <- t(prior * t(exp(ll)))
  dens / rowSums(dens)
}
