## define own class for raschmix and stepRaschmix objects
setClass("raschmix",
         representation(scores = "character",
                        extremeScoreProbs = "numeric",
                        rawScoresData = "table",
                        flx.call = "call",
                        nobs = "numeric",
                        identified.items = "factor"),
         contains = "flexmix")

setClass("stepRaschmix",
         representation(flx.call = "call"),
         contains = "stepFlexmix")


## define own class for btmix and stepBTmix objects
setClass("btmix",
         representation(flx.call = "call",
                        nobs = "numeric",
			labels = "character",
			mscale = "numeric"),
         contains = "flexmix")

setClass("stepBTmix",
         representation(flx.call = "call",
	                labels = "character",
			mscale = "numeric"),
         contains = "stepFlexmix")


