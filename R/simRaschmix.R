## dgp for the 3 designs by Rost (1990)
## and 4 designs with person parameters drawn from normal distributions
simRaschmix <- function(nobs = 1800, itemp = NULL, mean = NULL, sd = NULL,
                        design = c("rost1", "rost2", "rost3", "cont1",
                          "cont1-2", "cont2", "cont2-2"),
                        extremes = FALSE, attributes = TRUE)
{
  design <- match.arg(design)
  if (!is.null(itemp) & is.vector(itemp)) itemp <- as.matrix(itemp, ncol = 1)
  nitems <- 10

  ## warnings / data preparation
  if (design %in% c("rost1", "rost2", "rost3")){
    if(nobs != 1800) {
      warning("For a dataset from Rost (1990) the number of observations is 1800. Other values for nobs will be ignored.")
      nobs <- 1800
    }
    if (!is.null(itemp)){
      warning("For a dataset from Rost(1990) the item parameters cannot be given. Any input will be ignored.")
      itemp <- NULL
    }
    if (!is.null(mean)){
      warning("Mean will be ignored.")
      mean <- NULL
    }
    if (!is.null(sd)){
      warning("SD will be ignored.")
      sd <- NULL
    }
  } else {
    if (design %in% c("cont1", "cont1-2")) {
      if (!is.null(itemp)){
        if (ncol(itemp) > 1)
          warning("Only the first column of itemp will be used.")
        itemp <- itemp[ ,1, drop = FALSE]
      }
      if (!is.null(mean)){
        if (length(mean) > 1) warning("Only the first element of mean will be used.")
        mean <- mean[1]
      }
      if (!is.null(sd)){
        if (length(sd) > 1) warning("Only the first element of sd will be used.")
        sd <- sd[1]
      }
    } else {
      if (!is.null(itemp)){
        if (ncol(itemp) > 2){
          warning("Only the first two columns of itemp will be used.")
          itemp <- itemp[ ,1:2]
        } else if (ncol(itemp) < 2) stop("Not enough item parameters provided for this design.")
      }
      if (!is.null(mean)){
        if (length(mean) > 2){
          warning("Only the first two element of mean will be used.")
          mean <- mean[1:2]
        } else if (length(mean) < 2) stop ("For this design 2 values for mean have to be provided.")
      }
      if (!is.null(sd)){
        if (length(sd) > 2){
          warning("Only the first element of sd will be used.")
          sd <- sd[1:2]
        } else if (length(sd) < 2) stop ("For this design 2 values for sd have to be provided.")
      }
    }
    if (!is.null(itemp)) nitems <- nrow(itemp)
    if (is.null(mean)){
      if (design == "cont1") mean <- 0 else mean <- c(-2, 2)
    }
    if (is.null(sd)) sd <- c(1,1)
  }  
  
  
  ## item (difficulty) parameter
  ## Rost gives item easiness parameter but since the estimation will be done
  ## based on item difficulty parameters are converted to difficulty as well
  beta1 <- -c(-2.7,-2.1,-1.5,-0.9,-0.3,0.3,0.9,1.5,2.1,2.7)
  beta2 <- -c(2.7,2.1,1.5,0.9,0.3,-0.3,-0.9,-1.5,-2.1,-2.7)
  beta3 <- -c(0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5)

  ## response function
  resp <- function(ability, difficulty){
    fsmat <- outer(ability, difficulty, FUN = "-")
    psolve <- exp(fsmat) / (1 + exp(fsmat))
    ret <- matrix(runif(length(ability)*length(difficulty)),
                  nrow = length(ability), ncol = length(difficulty))
    ret <- (ret <= psolve)*1
    return(ret)
  }

  if (design == "rost1"){
    ## 1 group
    ## person parameters:
    ## 4 values, 1/4 of the sample for each value -> 450 for each value

    group <- rep(1, length.out = nobs)
    theta <- sample(rep(c(-2.7, -0.9, 0.9, 2.7), length.out = nobs))
    ret <- resp(theta, beta1)
  }
    
  if (design == "rost2"){
    ## 2 groups, equal size (900 per group)
    ## person parameters:
    ## 4 values, 1/4 of each group = 225 people for each value

    group <- sample(rep(1:2, length.out = nobs))
    theta1 <- sample(rep(c(-2.7, -0.9, 0.9, 2.7), length.out = (nobs/2)))
    theta2 <- sample(rep(c(-2.7, -0.9, 0.9, 2.7), length.out = (nobs/2)))
    theta <- rep(NA, length.out = nobs)
    theta[group == 1] <- theta1
    theta[group == 2] <- theta2    
    ret <- matrix(NA, nrow = nobs, ncol = nitems)
    ret[group == 1, ] <- resp(theta1, beta1)
    ret[group == 2, ] <- resp(theta2, beta2)
  }

  if (design == "rost3"){
    ## 3 groups with sizes: n1 = 800, n2 = 400, n3 = 800
    ## person parameters:
    ## groups 1+2: 4 values, 1/4 of each group for each value
    ## this means: group 1: 200 for each value, group 2: 100
    ## group 3: constant ability with value -0.9 chosen
    
    group <- sample(rep(1:3, times = nobs*c(4/9, 2/9, 3/9)))
    theta1 <- sample(rep(c(-2.7, -0.9, 0.9, 2.7), length.out = (nobs*4/9)))
    theta2 <- sample(rep(c(-2.7, -0.9, 0.9, 2.7), length.out = (nobs*2/9)))
    theta3 <- rep.int(0, (nobs * 3/9))
    theta <- rep(NA, length.out = nobs)
    theta[group == 1] <- theta1
    theta[group == 2] <- theta2
    theta[group == 3] <- theta3
    ret <- matrix(NA, nrow = nobs, ncol = nitems)
    ret[group == 1, ] <- resp(theta1, beta1)
    ret[group == 2, ] <- resp(theta2, beta2)
    ret[group == 3, ] <- resp(theta3, beta3)
  }

  if (design == "cont1"){
    ## 1 group wrt the IP
    ## IP: beta1
    ## PP: from 1 normal distributions with mean 0
    
    group <- rep(1, length.out = nobs)
    theta <- rnorm(nobs, mean = mean[1], sd = sd[1])
    if (!is.null(itemp)) beta1 <- itemp[ ,1]
    ret <- resp(theta, beta1)
  }

  if (design == "cont1-2"){
    ## 1 group wrt the IP, 2 groups wrt the PP
    ## IP: beta1
    ## PP: from 2 normal distributions with mean -2 and 2 (unless given otherwise)

    group <- rep(1, length.out = nobs)
    g <- sample.int(nobs, size = floor(nobs/2))
    theta1 <- rnorm(floor(nobs/2), mean = mean[1], sd = sd[1])
    theta2 <- rnorm(ceiling(nobs/2), mean = mean[2], sd = sd[2])
    theta <- rep(NA, length.out = nobs)
    theta[g] <- theta1
    theta[-g] <- theta2
    if (!is.null(itemp)) beta1 <- itemp[ ,1]
    ret <- matrix(NA, nrow = nobs, ncol = nitems)
    ret[g, ] <- resp(theta1, beta1)
    ret[-g, ] <- resp(theta2, beta1)
  }
  
  if (design == "cont2"){
    ## 2 groups, equal size
    ## IP: beta1 and beta2
    ## PP: from two normal distributions with mean -2 and 2 (unless given otherwise)
    
    group <- sample(rep(1:2, length.out = nobs))
    tab <- table(factor(group, levels = 1:2))
    theta1 <- rnorm(tab[1], mean = mean[1], sd = sd[1])
    theta2 <- rnorm(tab[2], mean = mean[2], sd = sd[2])
    theta <- rep(NA, length.out = nobs)
    theta[group == 1] <- theta1
    theta[group == 2] <- theta2
    if (!is.null(itemp)){ beta1 <- itemp[ ,1]; beta2 <- itemp[ ,2] }
    ret <- matrix(NA, nrow = nobs, ncol = nitems)
    ret[group == 1, ] <- resp(theta1, beta1)
    ret[group == 2, ] <- resp(theta2, beta2)
  }

  if (design == "cont2-2"){
    ## 2 groups wrt the IP, 2 groups wrt the PP, but 4 groups wrt IP*PP
    ## IP: beta1 and beta2
    ## PP: from two normal distributions with mean -2 and 2 (unless given otherwise)
    ## 4 groups: IP1-PP1, IP1-PP2, IP2-PP1, IP2-PP2
    
    g <- sample(rep(1:4, length.out = nobs))
    group <- rep(1, length.out = nobs)
    group[g == 3 | g == 4] <- 2
    tab <- table(factor(g, levels = 1:4))
    theta1 <- rnorm(tab[1], mean = mean[1], sd = sd[1])
    theta2 <- rnorm(tab[2], mean = mean[2], sd = sd[2])
    theta3 <- rnorm(tab[3], mean = mean[1], sd = sd[1])
    theta4 <- rnorm(tab[4], mean = mean[2], sd = sd[2])
    theta <- rep(NA, length.out = nobs)
    theta[g == 1] <- theta1
    theta[g == 2] <- theta2
    theta[g == 3] <- theta3
    theta[g == 4] <- theta4
    if (!is.null(itemp)){ beta1 <- itemp[ ,1]; beta2 <- itemp[ ,2] }
    ret <- matrix(NA, nrow = nobs, ncol = nitems)
    ret[g == 1, ] <- resp(theta1, beta1)
    ret[g == 2, ] <- resp(theta2, beta1)
    ret[g == 3, ] <- resp(theta3, beta2)
    ret[g == 4, ] <- resp(theta4, beta2)    
  }
  
  ## removal of observations with extreme raw scores
  if (!extremes){
    w <- which(rowSums(ret) == 0 | rowSums(ret) == nitems)
    if (length(w) > 0){
      ret <- ret[-w,]
      group <- group[-w]
      theta <- theta[-w]
    }
  }

  ## attach attributes
  if(attributes) {
    attr(ret, "group") <- group
    attr(ret, "person") <- theta
    if (design %in% c("rost1", "cont1", "cont1-2")){
        attr(ret, "item") <- data.frame(beta1 = beta1)
      }
    else {
      if (design %in% c("rost2", "cont2", "cont2-2")){
        attr(ret, "item") <- data.frame(beta1 = beta1, beta2 = beta2)
      }
      else attr(ret, "item") <- data.frame(beta1 = beta1, beta2 = beta2,
                                        beta3 = beta3)
    }
  }

  return(ret)
}
