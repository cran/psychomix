### R code from vignette source 'scores.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width = 70, prompt = "R> ", continue = "+  ")
library("psychomix")
set.seed(1090)
cache <- FALSE


###################################################
### code chunk number 2: simPrep
###################################################
## function to generate a design-list for simRaschmix()
generateDesign <- function(nobs = 500, m = 20, weights = NULL,
                           ab = 0, ab.dist = c("fix", "normal"),
                           dif = 2, beta = 1.9, index = 5, coincide = TRUE){
  ## weights
  if (is.null(weights)) weights <- rep(0.25, 4)
  
  ## coincide
  ## can only be set to FALSE if there are differences in both abilities and items
  ## = if either is 0, it has to be TRUE
  if (any(isTRUE(all.equal(c(ab, dif), 0)))) coincide <- TRUE
  
  ## ability
  ab.dist <- match.arg(ab.dist, c("fix", "normal"))
  ab <- c(-ab, ab)
  ab <- if (coincide) rep(ab, each = 2) else rep(ab, times = 2)
  
  ability <- if(ab.dist == "fix"){
    array(c(rbind(ab, 1)), dim = c(1,2,4)) # 1 = rel.frequency in sample()
  } else {
    rbind(ab, 0.5) ## 0.5 = sd for rnorm()
  }

  ## difficulty
  beta <- beta2 <- seq(from = -beta, to = beta, length = m)
  for (i in index){
    beta2[i] <- beta2[i] + dif
    beta2[m-i+1] <- beta2[m-i+1] - dif
  }
  difficulty <- cbind(beta, beta, beta2, beta2)
  
  return(list(nobs = nobs, weights = weights, ability = ability,
              difficulty = difficulty))
  ## NOTE: simRaschmix will return a cluster attribute with 4 different values/classes.
}

stacked_bars <- function(rs, cl = NULL, max = NULL, col = NULL, ...)
{
   if(is.null(max)) max <- max(rs)
   rs <- factor(rs, levels = 0:max)
   if(is.null(cl)) {
     tab <- table(rs)
     names(tab)[diff(-1:max %/% 5) < 1] <- ""
     if(is.null(col)) col <- gray.colors(2)[2]
   } else {
     tab <- table(cl, rs)
     colnames(tab)[diff(-1:max %/% 5) < 1] <- ""
     if(is.null(col)) col <- gray.colors(nrow(tab))
   }
   tab <- prop.table(tab)
   bp <- barplot(tab, space = 0, ...)
}

## colors
mygrays <- gray.colors(2)
myhcl <- psychomix:::qualitative_hcl(3)

## load simulated data
load("scoresim.rda")
scoresim$prop23 <-  1 - scoresim$prop1


###################################################
### code chunk number 3: fig-simD-DIFhomo
###################################################
## frame
par(mfrow = c(1,2))
par(mar = c(2, 4, 2, 2) + 0.1) # c(bottom, left, top, right)

## only DIF
des <- generateDesign(ab = 0, dif = 2)
set.seed(1)
dat <- simRaschmix(des)
rs <- rowSums(dat)
cl <- factor(as.numeric(attr(dat, "cluster") > 2) + 1)
ip <- attr(dat, "difficulty")[,2:3]

plot(ip[,1], type = "n", ylab = "Item difficulty", xlab = "")
points(ip[,2], type = "o", pch = 21, col = 1, bg = mygrays[2], lty = 2)
points(ip[,1], pch = 20, col = mygrays[1])

stacked_bars(rs, cl, max = 20, ylab = "Score frequencies")


###################################################
### code chunk number 4: fig-simD-noDIFhet
###################################################
## frame
par(mfrow = c(1, 2))
par(mar = c(2, 4, 2, 2) + 0.1) 

## only heterogeneous abilities
des <- generateDesign(ab = 1, dif = 0) ## ab = 1 --> impact = 2
set.seed(1)
dat <- simRaschmix(des)
rs <- rowSums(dat)
cl <- factor(as.numeric(attr(dat, "cluster") > 2) + 1)
ip <- attr(dat, "difficulty")[,2:3]

plot(ip[,1], type = "b", pch = 21, bg = mygrays[2], lty = 2,  ylab = "Item difficulty", xlab = "")

stacked_bars(rs, cl = NULL, max = 20, ylab = "Score frequencies", xlab = "")

par(mfrow = c(1, 1))


###################################################
### code chunk number 5: fig-simD-DIFhet
###################################################
## frame
par(mfrow = c(1,2))
par(mar = c(2, 4, 2, 2) + 0.1)

## designs:
## heterogeneous abilities within DIF groups
des <- generateDesign(ab = 1.0, dif = 2, coincide = FALSE) ## ab = 1 --> impact = 2
set.seed(1)
dat <- simRaschmix(des)
rs <- rowSums(dat)
cl <- factor(as.numeric(attr(dat, "cluster") > 2) + 1)
ip <- attr(dat, "difficulty")[,2:3]

stacked_bars(rs, cl, max = 20, ylab = "Score frequencies", xlab = "")

## heterogeneous abilities between DIF groups
des <- generateDesign(ab = 1.0, dif = 2, coincide = TRUE) ## ab = 1 --> impact = 2
set.seed(1)
dat <- simRaschmix(des)
rs <- rowSums(dat)
cl <- factor(as.numeric(attr(dat, "cluster") > 2) + 1)
ip <- attr(dat, "difficulty")[,2:3]

stacked_bars(rs, cl, max = 20, ylab = "Score frequencies", xlab = "")

par(mfrow = c(1,1))


###################################################
### code chunk number 6: fig-simR-DIFhomo
###################################################
par(mar = c(4, 4, 2, 2) + 0.1)

plot(prop23 ~ delta, data = scoresim, subset = theta == 0 & scores == "saturated", 
     ylim = c(0,1), type = "b", 
     xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
     ylab = "Choose more than 1 latent class", col = myhcl[3], lty = 3, pch = 3)
lines(prop23 ~ delta, data = scoresim, subset = theta == 0 & scores == "meanvar", 
     type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(prop23 ~ delta, data = scoresim, subset = theta == 0 & scores == "restricted", 
     type = "b", col = myhcl[1], lty = 2, pch = 6)
legend("topleft", legend = c("saturated", "mean-variance", "restricted"), 
       col = myhcl[3:1], lty = c(3,1,2), pch = c(3,1,6), bty = "n")

par(mar = c(5, 4, 4, 2) + 0.1)


###################################################
### code chunk number 7: fig-simR-noDIFhet
###################################################
par(mar = c(4, 4, 2, 2) + 0.1)

plot(prop23 ~ theta, data = scoresim, subset = delta == 0 & scores == "saturated", 
     ylim = c(0,1), type = "b", 
     xlab = expression(paste("Impact (", Theta, ")", sep = "")),
     ylab = "Choose more than 1 latent class", col = myhcl[3], lty = 3, pch = 3)
lines(prop23 ~ theta, data = scoresim, subset = delta == 0 & scores == "meanvar", 
     type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(prop23 ~ theta, data = scoresim, subset = delta == 0 & scores == "restricted", 
     type = "b", col = myhcl[1], lty = 2, pch = 6)
legend("topleft", legend = c("saturated", "mean-variance", "restricted"), 
       col = myhcl[3:1], lty = c(3,1,2), pch = c(3,1,6), bty = "n")

par(mar = c(5, 4, 4, 2) + 0.1)


###################################################
### code chunk number 8: fig-simR-DIFhetWithin
###################################################
par(mfrow = c(1, 3))
layout(matrix(c(rep(1,5), rep(2,4), rep(3,5)), nrow = 1, byrow = TRUE))

## impact = 0.4
par(mar = c(5, 4, 4, 0.5) + 0.1)
plot(prop23 ~ delta, data = scoresim, 
  subset = theta == 0.4 & scores == "saturated" & (scenario == 4 | delta == 0), 
  main = expression(paste("Impact ", Theta, " = 0.4", sep = "")),
  ylim = c(0,1), type = "b", 
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "Choose more than 1 latent class", col = myhcl[3], lty = 3, pch = 3)
lines(prop23 ~ delta, 
  data = scoresim, subset = theta == 0.4 & scores == "meanvar" & (scenario == 4 | delta == 0), 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(prop23 ~ delta, data = scoresim, 
  subset = theta == 0.4 & scores == "restricted" & (scenario == 4 | delta == 0),
  type = "b", col = myhcl[1], lty = 2, pch = 6)
legend("topleft", legend = c("saturated", "mean-variance", "restricted"), 
  col = myhcl[3:1], lty = c(3,1,2), pch = c(3,1,6), bty = "n")
                
## impact = 2.0
par(mar = c(5, 0.5, 4, 0.5) + 0.1) # c(bottom, left, top, right)
plot(prop23 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "saturated" & (scenario == 4 | delta == 0), 
  main = expression(paste("Impact ", Theta, " = 2.0", sep = "")),
  ylim = c(0,1), type = "b", axes = FALSE,
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "", col = myhcl[3], lty = 3, pch = 3)
box(); axis(1)
lines(prop23 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "meanvar" & (scenario == 4 | delta == 0), 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(prop23 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "restricted" & (scenario == 4 | delta == 0),
  type = "b", col = myhcl[1], lty = 2, pch = 6)

## impact = 3.6
par(mar = c(5, 0.5, 4, 4) + 0.1) # c(bottom, left, top, right)
plot(prop23 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "saturated" & (scenario == 4 | delta == 0), 
  main = expression(paste("Impact ", Theta, " = 3.6", sep = "")),
  ylim = c(0,1), type = "b", axes = FALSE,
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "", col = myhcl[3], lty = 3, pch = 3)
box(); axis(1)
lines(prop23 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "meanvar" & (scenario == 4 | delta == 0), 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(prop23 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "restricted" & (scenario == 4 | delta == 0),
  type = "b", col = myhcl[1], lty = 2, pch = 6)

par(mar = c(5, 4, 4, 2) + 0.1, mfrow = c(1,1)) ## restore default


###################################################
### code chunk number 9: fig-simR-prop23
###################################################
par(mar = c(4, 4, 4, 2) + 0.1)

plot(prop3 ~ delta, data = scoresim,
  subset = theta == 3.6 & scores == "meanvar" & (scenario == 4 | delta == 0),
  type = "b", ylim = c(0, 1), col = myhcl[2], pch = 1,
  main = expression(paste("Impact ", Theta, " = 3.6", sep = "")),
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "Selection proportions")
lines(prop3 ~ delta, data = scoresim,
  subset = theta == 3.6 & scores == "meanvar" & (scenario == 5 | delta == 0),
  type = "b", col = myhcl[2], pch = 16)

lines(prop2 ~ delta, data = scoresim,
  subset = theta == 3.6 & scores == "meanvar" & (scenario == 4 | delta == 0),
  type = "b", col = myhcl[2], pch = 2)
lines(prop2 ~ delta, data = scoresim,
  subset = theta == 3.6 & scores == "meanvar" & (scenario == 5 | delta == 0),
  type = "b", col = myhcl[2], pch = 17)

legend("topleft", 
       legend = c(expression(paste(hat(K), " = 2 - Sc 4", sep = "")),
           expression(paste(hat(K), " = 3 - Sc 4", sep = "")),
           expression(paste(hat(K), " = 2 - Sc 5", sep = "")),
           expression(paste(hat(K), " = 3 - Sc 5", sep = ""))),
           #"2 - Sc 4", "3 - Sc 4", "2 - Sc 5", "3 - Sc 5"), 
  col = myhcl[2], lty = 1, pch = c(2, 1, 17, 16), bty = "n")

par(mar = c(5, 4, 4, 2) + 0.1)


###################################################
### code chunk number 10: fig-simR-DIFhetBetween
###################################################
par(mfrow = c(1, 3))
layout(matrix(c(rep(1,5), rep(2,4), rep(3,5)), nrow = 1, byrow = TRUE))

## impact = 0.4
par(mar = c(5, 4, 4, 0.5) + 0.1)
plot(prop23 ~ delta, data = scoresim, 
  subset = theta == 0.4 & scores == "saturated" & (scenario == 5 | delta == 0), 
  main = expression(paste("Impact ", Theta, " = 0.4", sep = "")),
  ylim = c(0,1), type = "b", 
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "Choose more than 1 latent class", col = myhcl[3], lty = 3, pch = 3)
lines(prop23 ~ delta, 
  data = scoresim, subset = theta == 0.4 & scores == "meanvar" & (scenario == 5 | delta == 0), 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(prop23 ~ delta, data = scoresim, 
  subset = theta == 0.4 & scores == "restricted" & (scenario == 5 | delta == 0),
  type = "b", col = myhcl[1], lty = 2, pch = 6)
legend("topleft", legend = c("saturated", "mean-variance", "restricted"), 
  col = myhcl[3:1], lty = c(3,1,2), pch = c(3,1,6), bty = "n")

## impact = 2.0
par(mar = c(5, 0.5, 4, 0.5) + 0.1) # c(bottom, left, top, right)
plot(prop23 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "saturated" & (scenario == 5 | delta == 0), 
  main = expression(paste("Impact ", Theta, " = 2.0", sep = "")),
  ylim = c(0,1), type = "b", axes = FALSE,
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "", col = myhcl[3], lty = 3, pch = 3)
box(); axis(1)
lines(prop23 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "meanvar" & (scenario == 5 | delta == 0), 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(prop23 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "restricted" & (scenario == 5 | delta == 0),
  type = "b", col = myhcl[1], lty = 2, pch = 6)

## impact = 3.6
par(mar = c(5, 0.5, 4, 4) + 0.1) # c(bottom, left, top, right)
plot(prop23 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "saturated" & (scenario == 5 | delta == 0), 
  main = expression(paste("Impact ", Theta, " = 3.6", sep = "")),
  ylim = c(0,1), type = "b", axes = FALSE,
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "", col = myhcl[3], lty = 3, pch = 3)
box(); axis(1)
lines(prop23 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "meanvar" & (scenario == 5 | delta == 0), 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(prop23 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "restricted" & (scenario == 5 | delta == 0),
  type = "b", col = myhcl[1], lty = 2, pch = 6)

par(mar = c(5, 4, 4, 2) + 0.1, mfrow = c(1,1)) ## restore default


###################################################
### code chunk number 11: fig-simR-DIFhetWithinRand2
###################################################
par(mfrow = c(1, 3))
layout(matrix(c(rep(1,5), rep(2,4), rep(3,5)), nrow = 1, byrow = TRUE))

#par(mar = c(5, 4, 4, 6) + 0.1) # c(bottom, left, top, right)
par(mar = c(5, 4, 4, 0.5) + 0.1)
plot(rand2 ~ delta, data = scoresim, 
  subset = theta == 0.4 & scores == "saturated" & scenario == 4 & delta > 0, 
  main = expression(paste("Impact ", Theta, " = 0.4", sep = "")),
  ylim = c(0.5,1), type = "b", 
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "Rand index", col = myhcl[3], lty = 3, pch = 3)
lines(rand2 ~ delta, 
  data = scoresim, subset = theta == 0.4 & scores == "meanvar" & scenario == 4 & delta > 0, 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(rand2 ~ delta, data = scoresim, 
  subset = theta == 0.4 & scores == "restricted" & scenario == 4 & delta > 0,
  type = "b", col = myhcl[1], lty = 2, pch = 6)
legend("topleft", legend = c("saturated", "mean-variance", "restricted"), 
  col = myhcl[3:1], lty = c(3,1,2), pch = c(3,1,6), bty = "n")
                                        
par(mar = c(5, 0.5, 4, 0.5) + 0.1) # c(bottom, left, top, right)
plot(rand2 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "saturated" & scenario == 4 & delta > 0, 
  main = expression(paste("Impact ", Theta, " = 2.0", sep = "")),
  ylim = c(0.5,1), type = "b", axes = FALSE,
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "", col = myhcl[3], lty = 3, pch = 3)
box(); axis(1)
lines(rand2 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "meanvar" & scenario == 4 & delta > 0, 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(rand2 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "restricted" & scenario == 4 & delta > 0,
  type = "b", col = myhcl[1], lty = 2, pch = 6)

par(mar = c(5, 0.5, 4, 4) + 0.1) # c(bottom, left, top, right)
plot(rand2 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "saturated" & scenario == 4 & delta > 0, 
  main = expression(paste("Impact ", Theta, " = 3.6", sep = "")),
  ylim = c(0.5,1), type = "b", axes = FALSE,
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "", col = myhcl[3], lty = 3, pch = 3)
box(); axis(1)
lines(rand2 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "meanvar" & scenario == 4 & delta > 0, 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(rand2 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "restricted" & scenario == 4 & delta > 0,
  type = "b", col = myhcl[1], lty = 2, pch = 6)

par(mar = c(5, 4, 4, 2) + 0.1, mfrow = c(1,1)) ## restore default


###################################################
### code chunk number 12: fig-simR-DIFhetBetweenRand2
###################################################
par(mfrow = c(1, 3))
layout(matrix(c(rep(1,5), rep(2,4), rep(3,5)), nrow = 1, byrow = TRUE))

#par(mar = c(5, 4, 4, 6) + 0.1) # c(bottom, left, top, right)
par(mar = c(5, 4, 4, 0.5) + 0.1)
plot(rand2 ~ delta, data = scoresim, 
  subset = theta == 0.4 & scores == "saturated" & scenario == 5 & delta > 0, 
  main = expression(paste("Impact ", Theta, " = 0.4", sep = "")),
  ylim = c(0.5,1), type = "b", 
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "Rand index", col = myhcl[3], lty = 3, pch = 3)
lines(rand2 ~ delta, 
  data = scoresim, subset = theta == 0.4 & scores == "meanvar" & scenario == 5 & delta > 0, 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(rand2 ~ delta, data = scoresim, 
  subset = theta == 0.4 & scores == "restricted" & scenario == 5 & delta > 0,
  type = "b", col = myhcl[1], lty = 2, pch = 6)
legend("topleft", legend = c("saturated", "mean-variance", "restricted"), 
  col = myhcl[3:1], lty = c(3,1,2), pch = c(3,1,6), bty = "n")
                                        
par(mar = c(5, 0.5, 4, 0.5) + 0.1) # c(bottom, left, top, right)
plot(rand2 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "saturated" & scenario == 5 & delta > 0, 
  main = expression(paste("Impact ", Theta, " = 2.0", sep = "")),
  ylim = c(0.5,1), type = "b", axes = FALSE,
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "", col = myhcl[3], lty = 3, pch = 3)
box(); axis(1)
lines(rand2 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "meanvar" & scenario == 5 & delta > 0, 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(rand2 ~ delta, data = scoresim, 
  subset = theta == 2.0 & scores == "restricted" & scenario == 5 & delta > 0,
  type = "b", col = myhcl[1], lty = 2, pch = 6)

par(mar = c(5, 0.5, 4, 4) + 0.1) # c(bottom, left, top, right)
plot(rand2 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "saturated" & scenario == 5 & delta > 0, 
  main = expression(paste("Impact ", Theta, " = 3.6", sep = "")),
  ylim = c(0.5,1), type = "b", axes = FALSE,
  xlab = expression(paste("DIF effect size (", Delta, ")", sep = "")),  
  ylab = "", col = myhcl[3], lty = 3, pch = 3)
box(); axis(1)
lines(rand2 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "meanvar" & scenario == 5 & delta > 0, 
  type = "b", col = myhcl[2], lty = 1, pch = 1)
lines(rand2 ~ delta, data = scoresim, 
  subset = theta == 3.6 & scores == "restricted" & scenario == 5 & delta > 0,
  type = "b", col = myhcl[1], lty = 2, pch = 6)

par(mar = c(5, 4, 4, 2) + 0.1, mfrow = c(1,1)) ## restore default


###################################################
### code chunk number 13: VerbalData
###################################################
data("VerbalAggression", package = "psychotools")
VerbalAggression$resp2 <- VerbalAggression$resp2[, 1:12]
va12 <- subset(VerbalAggression, rowSums(resp2) > 0 & rowSums(resp2) < 12)
items <- colnames(va12$resp2)


###################################################
### code chunk number 14: VerbalResFit0 (eval = FALSE)
###################################################
## set.seed(403)
## mvR <- raschmix(resp2 ~ 1, data = va12, k = 1:4, scores = "meanvar",
##   restricted = TRUE)


###################################################
### code chunk number 15: VerbalResFit
###################################################
if(cache & file.exists("va12_mvR.rda")){
    load("va12_mvR.rda")
} else {
set.seed(403)
mvR <- raschmix(resp2 ~ 1, data = va12, k = 1:4, scores = "meanvar",
  restricted = TRUE)
    if(cache){
        save(mvR, file = "va12_mvR.rda")
    } else {
        if(file.exists("va12_mvR.rda")) file.remove("va12_mvR.rda")
    }
}


###################################################
### code chunk number 16: VerbalLC
###################################################
mvRbic <- round(BIC(mvR), digits = 1)
mvR3 <- getModel(mvR, which = "BIC")


###################################################
### code chunk number 17: VerbalScoresFit0 (eval = FALSE)
###################################################
## ## fit all possible score models
## sat3 <- raschmix(resp2 ~ 1, data = va12, k = 3, scores = "saturated")
## satR3 <- raschmix(resp2 ~ 1, data = va12, k = 3, scores = "saturated", 
##                   restricted = TRUE)
## mv3 <- raschmix(resp2 ~ 1, data = va12, k = 3, scores = "meanvar")


###################################################
### code chunk number 18: VerbalScoresFit
###################################################
if(cache & file.exists("va12_m3.rda")){
    load("va12_m3.rda")
} else {
## fit all possible score models
sat3 <- raschmix(resp2 ~ 1, data = va12, k = 3, scores = "saturated")
satR3 <- raschmix(resp2 ~ 1, data = va12, k = 3, scores = "saturated", 
                  restricted = TRUE)
mv3 <- raschmix(resp2 ~ 1, data = va12, k = 3, scores = "meanvar")
    if(cache){
        save(sat3, satR3, mv3, file = "va12_m3.rda")
    } else {
        if(file.exists("va12_m3.rda")) file.remove("va12_m3.rda")
    }
}


###################################################
### code chunk number 19: VerbalLRtest
###################################################
## check which score model fits best
library("lmtest")
lrtest(sat3, mv3, mvR3)
lrtest(sat3, satR3, mvR3)
clsizes <- table(clusters(mvR3))


###################################################
### code chunk number 20: fig-verbalTest
###################################################
par(mar = rep(0, 4))

node <- function(x, y, lab, eps = 0.11, col = mygrays[2]) {
  polygon(x + c(-1, -1, 1, 1, -1) * eps, y + c(-1, 1, 1, -1, -1) * eps, col = col)
  text(x, y, lab)
}
edge <- function(x1, y1, x2, y2, lab, eps = 0.13, length = 0.1, ...) {
  if(x1 == x2) {
    arrows(x1, y1 + eps, x2, y2 - eps, length = length, ...)
    text(mean(c(x1, x2)), mean(c(y1, y2)), lab, pos = 2)
  } else {
    arrows(x1 + eps, y1, x2 - eps, y2, length = length, ...)  
    text(mean(c(x1, x2)), mean(c(y1, y2)), lab, pos = 1)
  }
}

plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")

node(0.2, 0.2, sprintf("saturated\n(BIC: %s)", format(round(BIC(sat3), digits = 1), nsmall = 1)))
node(0.2, 0.8, sprintf("saturated\nrestricted\n(BIC: %s)", format(round(BIC(satR3), digits = 1), nsmall = 1)))
node(0.8, 0.2, sprintf("mean-variance\n(BIC: %s)", format(round(BIC(mv3), digits = 1), nsmall = 1)))
node(0.8, 0.8, sprintf("mean-variance\nrestricted\n(BIC: %s)", format(round(BIC(mvR3), digits = 1), nsmall = 1)))

edge(0.2, 0.2, 0.2, 0.8, format(round(lrtest(sat3, satR3)[2,5], digits = 3), nsmall = 3))
edge(0.2, 0.8, 0.8, 0.8, format(round(lrtest(satR3, mvR3)[2,5], digits = 3), nsmall = 3))
edge(0.2, 0.2, 0.8, 0.2, format(round(lrtest(sat3, mv3  )[2,5], digits = 3), nsmall = 3))
edge(0.8, 0.2, 0.8, 0.8, format(round(lrtest(mv3 , mvR3 )[2,5], digits = 3), nsmall = 3))

par(mar = c(5, 4, 4, 2) + 0.1)


###################################################
### code chunk number 21: VerbalItemPlot
###################################################
trellis.par.set(theme = standard.theme(color = FALSE))
xyplot(mvR3)


