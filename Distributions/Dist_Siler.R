# Truncated Siler distribution:

## zL = lower truncation point
## zU = upper truncation point

## load libraries
library(lamW)

## R function to run optimiser
optimiseR <- function(interval, a1, a2, b1, b2, c1, zL, zU, p) {
  silerRootFn <- function(x, a1, a2, b1, b2, c1, zL, zU, p) {
    abs(pSiler(x, a1, a2, b1, b2, c1, zL, zU) - p)
  }
  optimise(silerRootFn, interval, a1, a2, b1, b2, c1, zL, zU, p)$minimum
}

## probability density function
dSiler <- function(x, a1, a2, b1, b2, c1, zL = NA, zU = NA,
                   log = FALSE) {
  
  ## all parameters length = 1 or length(x)
  ntot <- length(x)
  if ((length(a1) != 1 & length(a1) != ntot) | 
      (length(a2) != 1 & length(a2) != ntot) | 
      (length(b1) != 1 & length(b1) != ntot) | 
      (length(b2) != 1 & length(b2) != ntot) | 
      (length(c1) != 1 & length(c1) != ntot)) {
    stop("Length of parameters must be = 1 or length(x)")
  }
  if ((length(zL) != 1 & length(zL) != ntot) | 
      (length(zU) != 1 & length(zU) != ntot)) {
    stop("Length of zL/zU must be = 1 or length(x)")
  }
  
  ## expand entries
  x <- cbind(x, zL, zU, a1, a2, b1, b2, c1)
  colnames(x) <- NULL
  
  ## run checks (first column of checkPass corresponds to invalid
  ## inputs [so should return NA]; second column is x outside range
  ## [so should return 0 p.d.f. outside the range])
  checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
  checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 5)
  checkPass <- cbind(checkPass, (x[, 1] >= 0) & (is.na(x[, 2]) | x[, 1] >= x[, 2]) & (is.na(x[, 3]) | x[, 1] <= x[, 3]))
  
  ## extract vectors as needed
  zL <- x[, 2]
  zU <- x[, 3]
  a1 <- x[, 4]
  a2 <- x[, 5]
  b1 <- x[, 6]
  b2 <- x[, 7]
  c1 <- x[, 8]
  x <- x[, 1]
  
  ## No truncation
  logS <- rep(NA, length(x))
  logH <- rep(NA, length(x))
  zN_ind <- which(checkPass[, 1] & checkPass[, 2])
  logS[zN_ind] <- (a1[zN_ind] / b1[zN_ind]) * (exp(-b1[zN_ind] * x[zN_ind]) - 1) - 
      c1[zN_ind] * x[zN_ind] + (a2[zN_ind] / b2[zN_ind]) * (1 - exp(b2[zN_ind] * x[zN_ind]))
  logH[zN_ind] <- log(a1[zN_ind] * exp(-b1[zN_ind] * x[zN_ind]) + c1[zN_ind] + a2[zN_ind] * exp(b2[zN_ind] * x[zN_ind]))
  logProb <- logH + logS
  
  ## Left truncation
  zL_ind <- which(!is.na(zL) & is.na(zU) & checkPass[, 1] & checkPass[, 2])
  if(length(zL_ind) > 0) {
      logS_zL <- (a1[zL_ind] / b1[zL_ind]) * (exp(-b1[zL_ind] * zL[zL_ind]) - 1) - c1[zL_ind] * zL[zL_ind] + (a2[zL_ind] / b2[zL_ind]) * (1 - exp(b2[zL_ind] * zL[zL_ind]))
      logProb[zL_ind] <- logProb[zL_ind] - logS_zL
  }
  
  ## Right truncation
  zU_ind <- which(is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2])
  if(length(zU_ind) > 0) {
      logS_zU <- (a1[zU_ind] / b1[zU_ind]) * (exp(-b1[zU_ind] * zU[zU_ind]) - 1) - c1[zU_ind] * zU[zU_ind] + (a2[zU_ind] / b2[zU_ind]) * (1 - exp(b2[zU_ind] * zU[zU_ind]))
      logProb[zU_ind] <- logProb[zU_ind] - log(1 - exp(logS_zU))
  }
  
  ## Interval truncation
  zI_ind <- which(!is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2])
  if(length(zI_ind) > 0) {
      logS_zL <- (a1[zI_ind] / b1[zI_ind]) * (exp(-b1[zI_ind] * zL[zI_ind]) - 1) - c1[zI_ind] * zL[zI_ind] + (a2[zI_ind] / b2[zI_ind]) * (1 - exp(b2[zI_ind] * zL[zI_ind]))
      logS_zU <- (a1[zI_ind] / b1[zI_ind]) * (exp(-b1[zI_ind] * zU[zI_ind]) - 1) - c1[zI_ind] * zU[zI_ind] + (a2[zI_ind] / b2[zI_ind]) * (1 - exp(b2[zI_ind] * zU[zI_ind]))
      S_zL <- exp(logS_zL)
      S_zU <- exp(logS_zU)
      logProb[zI_ind] <- logProb[zI_ind] - log(S_zL - S_zU)
  }
  
  ## return correctly
  logProb[!checkPass[, 2]] <- -Inf
  logProb[!checkPass[, 1]] <- NA
  
  if(log) {
      return(logProb)
  } else {
      return(exp(logProb))
  }
}

## cumulative distribution function (and survivor function)
pSiler <- function(q, a1, a2, b1, b2, c1, zL = NA, zU = NA,
                   lower.tail = TRUE, log.p = FALSE) {
  
  ## all parameters length = 1 or length(x)
  ntot <- length(q)
  if ((length(a1) != 1 & length(a1) != ntot) | 
      (length(a2) != 1 & length(a2) != ntot) | 
      (length(b1) != 1 & length(b1) != ntot) | 
      (length(b2) != 1 & length(b2) != ntot) | 
      (length(c1) != 1 & length(c1) != ntot)) {
    stop("Length of parameters must be = 1 or length(q)")
  }
  if ((length(zL) != 1 & length(zL) != ntot) | 
      (length(zU) != 1 & length(zU) != ntot)) {
    stop("Length of zL/zU must be = 1 or length(q)")
  }
  
  ## expand entries
  x <- cbind(q, zL, zU, a1, a2, b1, b2, c1)
  colnames(x) <- NULL
  
  ## run checks (first column of checkPass corresponds to invalid
  ## inputs [so should return NA]; second column is x below lower
  ## bound and third column is x above upper bound
  ## [so should return 0 c.d.f. before the lower bound and 1 above 
  ## the upper bound])
  checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
  checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 5)
  checkPass <- cbind(checkPass, (x[, 1] >= 0) & (is.na(x[, 2]) | x[, 1] >= x[, 2]))
  checkPass <- cbind(checkPass, is.na(x[, 3]) | x[, 1] <= x[, 3])
  
  ## extract vectors as needed
  zL <- x[, 2]
  zU <- x[, 3]
  a1 <- x[, 4]
  a2 <- x[, 5]
  b1 <- x[, 6]
  b2 <- x[, 7]
  c1 <- x[, 8]
  q <- x[, 1]
  
  ## No truncation
  logS <- (a1 / b1) * (exp(-b1 * q) - 1) - c1 * q + (a2 / b2) * (1 - exp(b2 * q))
  S <- exp(logS)
  
  ## Left truncation
  zL_ind <- which(!is.na(zL) & is.na(zU) & checkPass[, 1] & checkPass[, 2] & checkPass[, 3])
  if(length(zL_ind) > 0) {
      logS_zL <- (a1[zL_ind] / b1[zL_ind]) * (exp(-b1[zL_ind] * zL[zL_ind]) - 1) - c1[zL_ind] * zL[zL_ind] + (a2[zL_ind] / b2[zL_ind]) * (1 - exp(b2[zL_ind] * zL[zL_ind]))
      logS[zL_ind] <- logS[zL_ind] - logS_zL
  }
  
  ## Right truncation
  zU_ind <- which(is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2] & checkPass[, 3])
  if(length(zU_ind) > 0) {
      logS_zU <- (a1[zU_ind] / b1[zU_ind]) * (exp(-b1[zU_ind] * zU[zU_ind]) - 1) - c1[zU_ind] * zU[zU_ind] + (a2[zU_ind] / b2[zU_ind]) * (1 - exp(b2[zU_ind] * zU[zU_ind]))
      S_zU <- exp(logS_zU)
      logS[zU_ind] <- log(S[zU_ind]-S_zU) - log(1 - S_zU)
  }
  
  ## Interval truncation
  zI_ind <- which(!is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2] & checkPass[, 3])
  if(length(zI_ind) > 0) {
      logS_zLi <- (a1[zI_ind] / b1[zI_ind]) * (exp(-b1[zI_ind] * zL[zI_ind]) - 1) - c1[zI_ind] * zL[zI_ind] + (a2[zI_ind] / b2[zI_ind]) * (1 - exp(b2[zI_ind] * zL[zI_ind]))
      logS_zUi <- (a1[zI_ind] / b1[zI_ind]) * (exp(-b1[zI_ind] * zU[zI_ind]) - 1) - c1[zI_ind] * zU[zI_ind] + (a2[zI_ind] / b2[zI_ind]) * (1 - exp(b2[zI_ind] * zU[zI_ind]))
      S_zLi <- exp(logS_zLi)
      S_zUi <- exp(logS_zUi)
      logS[zI_ind] <- log(S[zI_ind] - S_zUi) - log(S_zLi - S_zUi)
  }

  ## return correctly
  logS[!checkPass[, 2]] <- 0
  logS[!checkPass[, 3]] <- -Inf
  logS[!checkPass[, 1]] <- NA
  
  if(!lower.tail) { 
    if(log.p) return(logS)
    else return(exp(logS))
  } else {
    p <- 1 - exp(logS)
    if(!log.p) return(p)
    else return(log(p))
  }
}

## quantile function
qSiler <- function(p, a1, a2, b1, b2, c1, zL = NA, zU = NA,
                   lower.tail = TRUE, log.p = FALSE) {
  ## all parameters length = 1 or length(x)
  ntot <- length(p)
  if ((length(a1) != 1 & length(a1) != ntot) | 
      (length(a2) != 1 & length(a2) != ntot) | 
      (length(b1) != 1 & length(b1) != ntot) | 
      (length(b2) != 1 & length(b2) != ntot) | 
      (length(c1) != 1 & length(c1) != ntot)) {
    stop("Length of parameters must be = 1 or length(p)")
  }
  if ((length(zL) != 1 & length(zL) != ntot) | 
      (length(zU) != 1 & length(zU) != ntot)) {
    stop("Length of zL/zU must be = 1 or length(p)")
  }
  
  ## check log and lower tail arguments and 
  ## adjust p accordingly
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  
  ## expand entries
  x <- cbind(p, zL, zU, a1, a2, b1, b2, c1)
  colnames(x) <- NULL
  
  ## run checks (first column of checkPass corresponds to invalid
  ## inputs [so should return NA]; second column is p outside range
  ## [so should return NA if p not in (0, 1)])
  checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
  checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 5)
  checkPass <- cbind(checkPass, x[, 1] >= 0 & x[, 1] <= 1)
  
  whichPass <- which(checkPass[, 1] & checkPass[, 2])
  if(length(whichPass) > 0) {
      ## extract valid entries
      x <- x[whichPass, , drop = FALSE]
  
      ## extract vectors as needed
      zL <- x[, 2]
      zU <- x[, 3]
      a1 <- x[, 4]
      a2 <- x[, 5]
      b1 <- x[, 6]
      b2 <- x[, 7]
      c1 <- x[, 8]
      p <- x[, 1]
      
      ## Left truncation
      zL_ind <- which(!is.na(zL) & is.na(zU))
      
      ## Right truncation
      zU_ind <- which(is.na(zL) & !is.na(zU))
      
      ## Interval truncation
      zI_ind <- which(!is.na(zL) & !is.na(zU))
    
      ## No truncation
      zN_ind <- which(is.na(zL) & is.na(zU))
      
      ## set lower bound
      qL <- rep(0, times = length(p))
      qL[zL_ind] <- zL[zL_ind]
      qL[zI_ind] <- zL[zI_ind]
      
      ## set/find upper bound
      qU <- qL + 10
      qU[zU_ind] <- zU[zU_ind]
      qU[zI_ind] <- zU[zI_ind]
      
      ## amend upper bound if necessary
      zNo_ind <- c(zN_ind, zL_ind)
      while(length(zNo_ind) > 0) {
        pTarget <- p[zNo_ind]
        pU <- pSiler(qU[zNo_ind], a1[zNo_ind], a2[zNo_ind], b1[zNo_ind], b2[zNo_ind], c1[zNo_ind], zL[zNo_ind], zU[zNo_ind])
        zNo_ind <- zNo_ind[pU < pTarget]
        pU <- pU[pU < pTarget]
        if(length(zNo_ind) > 0) {
          qU[zNo_ind] <- qU[zNo_ind] + 10
        }
      }
      
      data <- cbind(qL, qU, a1, a2, b1, b2, c1, zL, zU, p)
      out <- rep(NA, ntot)
      out[whichPass] <- apply(data, 1, function(x) {
          optimiseR(x[1:2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10])
      })
  } else {
      out <- rep(NA, ntot)
  }
  return(out)
}

## combined function
rSiler <- function(n, a1, a2, b1, b2, c1, zL = NA, zU = NA) {
  
  ## all parameters length = 1 or length(x)
  if(length(n) > 1) {
    n <- length(n)
    message("Since 'n' is a vector, 'length(n)' is taken to be the number of samples required.")
  }
  if ((length(a1) != 1 & length(a1) != n) | 
      (length(a2) != 1 & length(a2) != n) | 
      (length(b1) != 1 & length(b1) != n) | 
      (length(b2) != 1 & length(b2) != n) | 
      (length(c1) != 1 & length(c1) != n)) {
    stop("Length of parameters must be = 1 or n")
  }
  if ((length(zL) != 1 & length(zL) != n) | 
      (length(zU) != 1 & length(zU) != n)) {
    stop("Length of zL/zU must be = 1 or n")
  }
  
  ## expand entries
  x <- cbind(zL, zU, a1, a2, b1, b2, c1)
  if(n == nrow(x)) {
      x <- cbind(1, x)
  } else {
      if(nrow(x) != 1) {
          stop("Error we're not catching")
      } else {
          x <- matrix(rep(c(1, x[1, ]), n), nrow = n, byrow = TRUE)
      }
  }
  colnames(x) <- NULL
  
  ## run checks (corresponds to invalid
  ## inputs [so should return NA])
  checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
  checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 5)
  
  ## extract vectors as needed
  zL <- x[, 2]
  zU <- x[, 3]
  a1 <- x[, 4]
  a2 <- x[, 5]
  b1 <- x[, 6]
  b2 <- x[, 7]
  c1 <- x[, 8]
  
  ## set up output vector
  rs <- rep(NA, n)
  
  ## no truncation and/or left-truncation
  zN_ind <- which(is.na(zL) & is.na(zU) & checkPass)
  zL_ind <- which(!is.na(zL) & is.na(zU) & checkPass)
  zNL_ind <- c(zN_ind, zL_ind)
  nNL <- length(zNL_ind)
  if(nNL > 0) {
    ## sample from distribution 1
    u1 <- runif(nNL, 0, 1)
    u2 <- runif(nNL, 0, 1)
    
    ## distribution 1
    logS_zNL <- (a1[zNL_ind] / b1[zNL_ind]) * (exp(-b1[zNL_ind] * zL[zNL_ind]) - 1) - c1[zNL_ind] * zL[zNL_ind]
    S_zNL <- exp(logS_zNL)
    
    ## distribution 2
    logS_zNL2 <-  (a2[zNL_ind] / b2[zNL_ind]) * (1 - exp(b2[zNL_ind] * zL[zNL_ind]))
    S_zNL2 <- exp(logS_zNL2)
    
    ## left truncation
    zL_ind1 <- match(zL_ind, zNL_ind)
    if(length(zL_ind1) > 0) {
      u1[zL_ind1] <- u1[zL_ind1] * S_zNL[zL_ind1]
      u2[zL_ind1] <- u2[zL_ind1] * S_zNL2[zL_ind1]
    }
    
    ## sample from distribution 1
    w0 <- (a1[zNL_ind] / c1[zNL_ind]) * exp((log(u1) + a1[zNL_ind] / b1[zNL_ind]) * (b1[zNL_ind] / c1[zNL_ind]))
    w0 <- lambertW0(w0)
    x1 <- (-1 / c1[zNL_ind]) * (log(u1) + a1[zNL_ind] / b1[zNL_ind]) + w0 / b1[zNL_ind]
    ## sample from distribution 2
    x2 <- log(1 - log(u2) * (b2[zNL_ind] / a2[zNL_ind])) / b2[zNL_ind]
    
    rs[zNL_ind] <- ifelse(x1 < x2, x1, x2)
  }
  
  ## Right truncation
  zU_ind <- which(is.na(zL) & !is.na(zU) & checkPass)
  if(length(zU_ind) > 0) {
    rs[zU_ind] <- qSiler(runif(length(zU_ind), 0, 1), a1[zU_ind], a2[zU_ind], b1[zU_ind], b2[zU_ind], c1[zU_ind], zL[zU_ind], zU[zU_ind])
  }
  
  ## Interval truncation
  zI_ind <- which(!is.na(zL) & !is.na(zU) & checkPass)
  if(length(zI_ind) > 0) {
    rs[zI_ind] <- qSiler(runif(length(zI_ind), 0, 1), a1[zI_ind], a2[zI_ind], b1[zI_ind], b2[zI_ind], c1[zI_ind], zL[zI_ind], zU[zI_ind])
  }
  
  ## return correctly
  return(rs)
}
