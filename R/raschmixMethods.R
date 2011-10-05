setGeneric("worth")
setMethod("worth", "raschmix", function(object, difficulty = TRUE,
                                        component = NULL){
  ## get parameters
  coef <- parameters(object, which = "item", difficulty = difficulty,
                     component = component)

  ## apply transformation
  worth <- apply(coef, 2, function(x){c(0,x) - mean(c(0,x))})

  ## include non-identified items (if any)
  rval <- matrix(NA, ncol = ncol(worth), nrow = length(object@identified.items))
  rval[object@identified.items == "0/1", ] <- worth
  rval[object@identified.items == "0", ] <- if(!difficulty) -Inf else Inf
  rval[object@identified.items == "1", ] <- if(!difficulty) Inf else -Inf
  rownames(rval) <- names(object@identified.items)
  colnames(rval) <- colnames(worth)

  return(rval)
})

## turn this into a generic method?
scoreProbs <- function(object, component = NULL, simplify = TRUE,
                        drop = TRUE){

  ## call flexmix method
  para <- parameters(as(object, "flexmix"), component = component, model = NULL,
                         which = "model", simplify = FALSE, drop = drop)

  scores <- lapply(para, function(comp){
    ret <- c(object@extremeScoreProbs[1], comp$scoreProbs,
          object@extremeScoreProbs[2])
    names(ret) <- (1:length(ret))-1
    return(ret)
  })

  if (simplify) {
    scores <- do.call("cbind", scores)
  }

  return(scores)
}


setMethod("parameters", "raschmix", function(object,
       which = c("model", "item", "score", "concomitant"),
       difficulty = TRUE, component = NULL, simplify = TRUE){
  
  which <- match.arg(which)
  which.flx <- if (which %in% c("item", "score")) "model" else which

  if (is.null(component)) component <- 1:object@k
  drop <- length(component) > 1
  
  ## call flexmix method
  para <- callNextMethod(object, component = component, model = NULL,
                         which = which.flx, simplify = FALSE, drop = drop)

  if (which != "concomitant"){

    if (length(component) == 1) para <- para[[1]]

    para <- lapply(para, function(comp){
      comp$scoreProbs <- NULL
      if (!difficulty) comp$item <- -comp$item
      if (object@scores == "meanvar"){
        names(comp$score) <- c("location", "dispersion")
      }
      if (which == "score") comp$item <- NULL
      if (which == "item") comp$score <- NULL
      return(comp)
    })

    if (simplify){
      para <- lapply(para, unlist)
      para <- do.call("cbind", para)
    }
    ## <FIXME> for object@scores == "constant" & which == "score":
    ## return sth "empty" rather than NULL? </FIXME>

  }
  return(para)
})



## <INFO>
## summary-method for raschmix objects:
## it shows:
## call (to raschmix), information about the mixture, item parameters, logLik, AIC, and BIC
## it doesn't show (yet):
## standard errors for item parameters (needs refit function), coefficients for concomitant model
## other notes:
## flexmix separates display of mixture from display of components (summary vs parameters methods) --> stick with that?
## an option to select a certain component is not implemented as it's a summary of the mixture
## there is no summary method for stepFlexmix objects, just a show method -- this gets updated to use
## the raschmix call instead of the flexmix call but is otherwise left untouched
## </INFO>

## modified code from flexmix
setClass("summary.raschmix",
         representation(itemParaTab = "ANY"),
         contains = "summary.flexmix")

setMethod("show", "summary.raschmix",
function(object){
    cat("\nCall:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\n")
    print(object@comptab, digits=3)
    cat("\nItem Parameters:\n")
    print(object@itemParaTab)
    cat("\n")
    print(object@logLik)
    cat("AIC:", object@AIC, "  BIC:", object@BIC, "\n")
    cat("\n")    
})

setMethod("summary", "raschmix",
function(object, eps=1e-4, ...){    
    z <- new("summary.raschmix",
             call = object@call,
             AIC = AIC(object),
             BIC = BIC(object),
             logLik = logLik(object))

    TAB <- data.frame(prior=object@prior,
                      size=object@size)
    rownames(TAB) <- paste("Comp.", seq_len(nrow(TAB)), sep="")
    TAB[["post>0"]] <- colSums(object@posterior$scaled > eps)
    TAB[["ratio"]] <- TAB[["size"]]/TAB[["post>0"]]
    
    z@comptab = TAB
    z@itemParaTab <- worth(object)
    z
})


## <INFO>
## logLik-method for "stepRaschmix" doesn't display df correctly,
## but logLik in not indented for this use anyway
##
## refit
## certainly for raschmix, also for "stepRaschmix"?
##
## ## weights
## ## for raschmix
## setGeneric("weights")
## setMethod("weights", "raschmix", function(object){
##   if(!is.null(object@weights)) object@weights else rep(1, object@nobs)
## })
##
## ## EIC for "stepRaschmix"?
## </INFO>
