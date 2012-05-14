## raschmix
setMethod("FLXgetParameters", signature(object="raschmix"),
function(object, model) {
  if (missing(model)) model <- seq_along(object@model)
  parms <- callNextMethod(object, model)
  esp <- object@extremeScoreProbs
  c(parms, extremeScoreProbs = esp[esp != 0])
})

setMethod("FLXreplaceParameters", signature(object="raschmix"),
function(object, parms) {
  n_extremeScores <- object@nobs/(1 - sum(object@extremeScoreProbs)) * object@extremeScoreProbs
  comp_names <- names(object@components)
  components <- list()
  for (m in seq_along(object@model)) {
    indices <- grep(paste("^model.", m, sep = ""), names(parms))
    components[[m]] <- FLXreplaceParameters(object@model[[m]], lapply(object@components, "[[", m), parms[indices],
      scores = object@scores, deriv = object@deriv, extremeScoreProbs = object@extremeScoreProbs)
  }
  object@components <- lapply(seq_along(object@components), function(k) lapply(components, "[[", k))
  names(object@components) <- comp_names
  if (object@k > 1) {
    indices <- grep("^concomitant_", names(parms))
    object@concomitant <- FLXreplaceParameters(object@concomitant, parms[indices])
  }
  indices <- grep("^extremeScoreProbs", names(parms))
  object@extremeScoreProbs <- rep(0, 2)
  if (length(indices) > 0) object@extremeScoreProbs[n_extremeScores > 0] <- parms[indices]
  object
})

setMethod("FLXreplaceParameters", signature(object="FLXMCrasch"),
function(object, components, parms, scores, deriv, extremeScoreProbs) {
  Design <- FLXgetDesign(object, components)
  lapply(seq_along(components), function(k) {
    Parameters <- list()
    parms_k <- parms[as.logical(Design[k,])]
    for (i in seq_along(components[[k]]@parameters)) {
      Parameters[[i]] <- parms_k[seq_along(components[[k]]@parameters[[i]])]
      attributes(Parameters[[i]]) <- attributes(components[[k]]@parameters[[i]])
      parms_k <- parms_k[-seq_along(components[[k]]@parameters[[i]])]
    }
    names(Parameters) <- names(components[[k]]@parameters)
    Parameters$df <- components[[k]]@df
    Parameters$scores <- scores
    Parameters$deriv <- deriv
    Parameters$nonExtremeProb <- 1 - extremeScoreProbs
    variables <- c("x", "y")
    for (var in variables) 
      assign(var, slot(object, var))
    with(Parameters, eval(object@defineComponent))
  })
})

setMethod("FLXlogLikfun", signature(object="raschmix"),
function(object, ...) function(parms) {
  n_extremeScores <- object@nobs/(1 - sum(object@extremeScoreProbs)) * object@extremeScoreProbs
  object <- FLXreplaceParameters(object, parms)
  groupfirst <- if (length(object@group) > 1) groupFirst(object@group) else rep(TRUE, FLXgetObs(object@model[[1]]))
  logpostunscaled <- flexmix:::logLikfun_comp(object) +
    log(flexmix:::getPriors(object@concomitant, object@group, groupfirst))
  loglik <- if (is.null(object@weights)) sum(flexmix:::log_row_sums(logpostunscaled[groupfirst,,drop=FALSE]))
            else sum(flexmix:::log_row_sums(logpostunscaled[groupfirst,,drop=FALSE])*object@weights[groupfirst])
  extremeScoreProbs <- parms[grep("^extremeScoreProbs", names(parms))]
  loglik + sum(ifelse(n_extremeScores > 0, n_extremeScores * log(extremeScoreProbs), 0))
})

## btmix
setMethod("FLXreplaceParameters", signature(object="btmix"),
function(object, parms) {
  comp_names <- names(object@components)
  components <- list()
  for (m in seq_along(object@model)) {
    indices <- grep(paste("^model.", m, sep = ""), names(parms))
    components[[m]] <- FLXreplaceParameters(object@model[[m]], lapply(object@components, "[[", m), parms[indices],
                                            labels = object@labels, undecided = object@undecided, ref = object@ref,
                                            type = object@type)
  }
  object@components <- lapply(seq_along(object@components), function(k) lapply(components, "[[", k))
  names(object@components) <- comp_names
  if (object@k > 1) {
    indices <- grep("^concomitant_", names(parms))
    object@concomitant <- FLXreplaceParameters(object@concomitant, parms[indices])
  }
  object
})

setMethod("FLXreplaceParameters", signature(object="FLXMCbt"),
function(object, components, parms, labels, undecided, ref, type) {
  Design <- FLXgetDesign(object, components)
  lapply(seq_along(components), function(k) {
    Parameters <- list()
    parms_k <- parms[as.logical(Design[k,])]
    for (i in seq_along(components[[k]]@parameters)) {
      Parameters[[i]] <- parms_k[seq_along(components[[k]]@parameters[[i]])]
      attributes(Parameters[[i]]) <- attributes(components[[k]]@parameters[[i]])
      parms_k <- parms_k[-seq_along(components[[k]]@parameters[[i]])]
    }
    names(Parameters) <- names(components[[k]]@parameters)
    Parameters$df <- components[[k]]@df
    Parameters$labels <- labels
    Parameters$undecided <- undecided
    Parameters$ref <- ref
    Parameters$type <- type
    variables <- c("x", "y")
    for (var in variables) 
      assign(var, slot(object, var))
    with(Parameters, eval(object@defineComponent))
  })
})
