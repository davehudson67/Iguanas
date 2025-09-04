# Truncated Gompertz Makeham distribution:

## zL = lower truncation point
## zU = upper truncation point

#### Probability density function
dGompzMake <- function(x, a, b, c1, zL = NA, zU = NA, log = FALSE) {
  
    ## all parameters length = 1 or length(x)
    ntot <- length(x)
    if ((length(a) != 1 & length(a) != ntot) | 
        (length(b) != 1 & length(b) != ntot) | 
        (length(c1) != 1 & length(c1) != ntot)) {
      stop("Length of parameters must be = 1 or length(x)")
    }
    if ((length(zL) != 1 & length(zL) != ntot) | 
        (length(zU) != 1 & length(zU) != ntot)) {
      stop("Length of zL/zU must be = 1 or length(x)")
    }
    
    ## expand entries
    x <- cbind(x, zL, zU, a, b, c1)
    colnames(x) <- NULL
    
    ## run checks (first column of checkPass corresponds to invalid
    ## inputs [so should return NA]; second column is x outside range
    ## [so should return 0 p.d.f. outside the range])
    checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
    checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 3)
    checkPass <- cbind(checkPass, (x[, 1] >= 0) & (is.na(x[, 2]) | x[, 1] >= x[, 2]) & (is.na(x[, 3]) | x[, 1] <= x[, 3]))
    
    ## extract vectors as needed
    zL <- x[, 2]
    zU <- x[, 3]
    a <- x[, 4]
    b <- x[, 5]
    c1 <- x[, 6]
    x <- x[, 1]
    
    ## log PDF
    logS <- rep(NA, length(x))
    logH <- rep(NA, length(x))
    zN_ind <- which(checkPass[, 1] & checkPass[, 2])
    logS[zN_ind] <- -c1[zN_ind] * x[zN_ind] - (a[zN_ind] / b[zN_ind]) * (exp(b[zN_ind] * x[zN_ind]) - 1)
    logH[zN_ind] <- log(a[zN_ind] * exp(b[zN_ind] * x[zN_ind]) + c1[zN_ind])
    logProb <- logH + logS
  
    #interval truncation
    zI_ind <- which(!is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2])
    if(length(zI_ind) > 0) {
        logS_zL <- -c1[zI_ind] * zL[zI_ind] - (a[zI_ind] / b[zI_ind]) * (exp(b[zI_ind] * zL[zI_ind]) - 1)
        logS_zU <- -c1[zI_ind] * zU[zI_ind] - (a[zI_ind] / b[zI_ind]) * (exp(b[zI_ind] * zU[zI_ind]) - 1)
        S_zL <-  exp(logS_zL)
        S_zU <- exp(logS_zU)
        logProb[zI_ind] <- logProb[zI_ind] - log(S_zL - S_zU)
    }
    
    #left truncation
    zL_ind <- which(!is.na(zL) & is.na(zU) & checkPass[, 1] & checkPass[, 2])
    if(length(zL_ind) > 0) {
    logS_zL <- -c1[zL_ind] * zL[zL_ind] - (a[zL_ind] / b[zL_ind]) * (exp(b[zL_ind] * zL[zL_ind]) - 1)
    logProb[zL_ind] <- logProb[zL_ind] - logS_zL
    }

    #right truncation
    zU_ind <- which(is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2])
    if(length(zU_ind) > 0) {
        logS_zU <- -c1[zU_ind] * zU[zU_ind] - (a[zU_ind] / b[zU_ind]) * (exp(b[zU_ind] * zU[zU_ind]) - 1)
        logProb[zU_ind] <- logProb[zU_ind] - log(1 - exp(logS_zU))
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

## Cumulative distribution function (and survivor function)
pGompzMake <- function(q, a, b, c1, zL = NA, zU = NA,
                      lower.tail = TRUE, log.p = FALSE) {
  
    ## all parameters length = 1 or length(x)
    ntot <- length(q)
    if ((length(a) != 1 & length(a) != ntot) | 
        (length(b) != 1 & length(b) != ntot) | 
        (length(c1) != 1 & length(c1) != ntot)) {
      stop("Length of parameters must be = 1 or length(q)")
    }
    if ((length(zL) != 1 & length(zL) != ntot) | 
        (length(zU) != 1 & length(zU) != ntot)) {
      stop("Length of zL/zU must be = 1 or length(q)")
    }
    
    ## expand entries
    x <- cbind(q, zL, zU, a, b, c1)
    colnames(x) <- NULL
    
    ## run checks (first column of checkPass corresponds to invalid
    ## inputs [so should return NA]; second column is x below lower
    ## bound and third column is x above upper bound
    ## [so should return 0 c.d.f. before the lower bound and 1 above 
    ## the upper bound])
    checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
    checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 3)
    checkPass <- cbind(checkPass, (x[, 1] >= 0) & (is.na(x[, 2]) | x[, 1] >= x[, 2]))
    checkPass <- cbind(checkPass, is.na(x[, 3]) | x[, 1] <= x[, 3])
    
    ## extract vectors as needed
    zL <- x[, 2]
    zU <- x[, 3]
    a <- x[, 4]
    b <- x[, 5]
    c1 <- x[, 6]
    q <- x[, 1]
    
    ## log survivor function
    logS <- -c1 * q - (a / b) * (exp(b * q) - 1)
    S <- exp(logS)
    
    #interval truncation (same as right truncation for survivor function)
    zI_ind <- which(!is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2] & checkPass[, 3])
    if(length(zI_ind) > 0) {
        logS_zU <- -c1[zI_ind] * zU[zI_ind] - (a[zI_ind] / b[zI_ind]) * (exp(b[zI_ind] * zU[zI_ind]) - 1)
        logS_zL <- -c1[zI_ind] * zL[zI_ind] - (a[zI_ind] / b[zI_ind]) * (exp(b[zI_ind] * zL[zI_ind]) - 1)
        S_zU <- exp(logS_zU)
        S_zL <- exp(logS_zL)
        logS[zI_ind] <- log(S[zI_ind] - S_zU) - log(S_zL - S_zU)
    }
    
    #left truncation (same as no truncation for survivor function)
    zL_ind <- which(!is.na(zL) & is.na(zU) & checkPass[, 1] & checkPass[, 2] & checkPass[, 3])
    if(length(zL_ind) > 0) {
        logS_zL <- -c1[zL_ind] * zL[zL_ind] - (a[zL_ind] / b[zL_ind]) * (exp(b[zL_ind] * zL[zL_ind]) - 1)
        logS[zL_ind] <- logS[zL_ind] - logS_zL
    }

    #right truncation
    zU_ind <- which(is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2] & checkPass[, 3])
    if(length(zU_ind) > 0) {
        logS_zU <- -c1[zU_ind] * zU[zU_ind] - (a[zU_ind] / b[zU_ind]) * (exp(b[zU_ind] * zU[zU_ind]) - 1)
        S_zU <- exp(logS_zU)
        logS[zU_ind] <- log(S[zU_ind] - S_zU) - log(1 - S_zU)
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

## Quantile function
qGompzMake<- function(p, a, b, c1, zL = NA, zU = NA,
                     lower.tail = TRUE, log.p = FALSE) {
    ## all parameters length = 1 or length(x)
    ntot <- length(p)
    if ((length(a) != 1 & length(a) != ntot) | 
        (length(b) != 1 & length(b) != ntot) | 
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
    x <- cbind(p, zL, zU, a, b, c1)
    colnames(x) <- NULL
    
    ## run checks (first column of checkPass corresponds to invalid
    ## inputs [so should return NA]; second column is p outside range
    ## [so should return NA if p not in (0, 1)])
    checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
    checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 3)
    checkPass <- cbind(checkPass, x[, 1] >= 0 & x[, 1] <= 1)
    
    whichPass <- which(checkPass[, 1] & checkPass[, 2])
    if(length(whichPass) > 0) {
        ## extract valid entries
        x <- x[whichPass, , drop = FALSE]
        
        ## extract vectors as needed
        zL <- x[, 2]
        zU <- x[, 3]
        a <- x[, 4]
        b <- x[, 5]
        c1 <- x[, 6]
        p <- x[, 1]
        
        w0 <- (a / c1) * exp(a / c1) * ((1 - p)^(-b / c1))
        w0 <- lambertW0(w0)
        out1 <- (a / (b * c1)) - (1 / c1) * log(1 - p) - (1 / b) *  w0
        
        #interval truncation
        zI_ind <- which(!is.na(zL) & !is.na(zU))
        if(length(zI_ind) > 0) {
            logS_zU <- -c1[zI_ind] * zU[zI_ind] - (a[zI_ind] / b[zI_ind]) * (exp(b[zI_ind] * zU[zI_ind]) - 1)
            S_zU <- exp(logS_zU)
            CDF_zU <- 1 - S_zU
                
            logS_zL <- -c1[zI_ind] * zL[zI_ind] - (a[zI_ind] / b[zI_ind]) * (exp(b[zI_ind] * zL[zI_ind]) - 1)
            S_zL <- exp(logS_zL)
            CDF_zL <- 1 - S_zL
                
            p[zI_ind] <- (p[zI_ind] * (CDF_zU - CDF_zL)) + CDF_zL
            w0[zI_ind] <- (a[zI_ind] / c1[zI_ind]) * exp(a[zI_ind] / c1[zI_ind]) * ((1 - p[zI_ind])^(-b[zI_ind] / c1[zI_ind]))
            w0[zI_ind] <- lambertW0(w0[zI_ind])
            out1[zI_ind] <- (a[zI_ind] / (b[zI_ind] * c1[zI_ind])) - (1 / c1[zI_ind]) * log(1 - p[zI_ind]) - (1 / b[zI_ind]) *  w0[zI_ind]
        }
        
        #left truncation
        zL_ind <- which(!is.na(zL) & is.na(zU))
        if(length(zL_ind) > 0) {
            logS_zL <-  -c1[zL_ind] * zL[zL_ind] - (a[zL_ind] / b[zL_ind]) * (exp(b[zL_ind] * zL[zL_ind]) - 1)
            S_zL <- exp(logS_zL)
            CDF_zL <- 1 - S_zL
                
            p[zL_ind] <- p[zL_ind] * S_zL + CDF_zL
            w0[zL_ind] <- (a[zL_ind] / c1[zL_ind]) * exp(a[zL_ind] / c1[zL_ind]) * ((1 - p[zL_ind])^(-b[zL_ind] / c1[zL_ind]))
            w0[zL_ind] <- lambertW0(w0[zL_ind])
            out1[zL_ind] <- (a[zL_ind] / (b[zL_ind] * c1[zL_ind])) - (1 / c1[zL_ind]) * log(1 - p[zL_ind]) - (1 / b[zL_ind]) *  w0[zL_ind]
        }
        
        #right truncation
        zU_ind <- which(is.na(zL) & !is.na(zU))
        if(length(zU_ind) > 0) {
            logS_zU <-  -c1[zU_ind] * zU[zU_ind] - (a[zU_ind] / b[zU_ind]) * (exp(b[zU_ind] * zU[zU_ind]) - 1)
            S_zU <- exp(logS_zU)
            CDF_zU <- 1 - S_zU
                
            p[zU_ind] <- p[zU_ind] * CDF_zU
            w0[zU_ind] <- (a[zU_ind] / c1[zU_ind]) * exp(a[zU_ind] / c1[zU_ind]) * ((1 - p[zU_ind])^(-b[zU_ind] / c1[zU_ind]))
            w0[zU_ind] <- lambertW0(w0[zU_ind])
            out1[zU_ind] <- (a[zU_ind] / (b[zU_ind] * c1[zU_ind])) - (1 / c1[zU_ind]) * log(1 - p[zU_ind]) - (1 / b[zU_ind]) *  w0[zU_ind]
        }
        out <- rep(NA, ntot)
        out[whichPass] <- out1
    } else {
        out <- rep(NA, ntot)
    }
   return(out)
}

## Random sampler
rGompzMake <- function(n, a, b, c1, zL = NA, zU = NA) {
  
    ## all parameters length = 1 or length(x)
    if(length(n) > 1) {
      n <- length(n)
      message("Since 'n' is a vector, 'length(n)' is taken to be the number of samples required.")
    }
    if ((length(a) != 1 & length(a) != n) | 
        (length(b) != 1 & length(b) != n) | 
        (length(c1) != 1 & length(c1) != n)) {
      stop("Length of parameters must be = 1 or n")
    }
    if ((length(zL) != 1 & length(zL) != n) | 
        (length(zU) != 1 & length(zU) != n)) {
      stop("Length of zL/zU must be = 1 or n")
    }
    
    ## expand entries
    x <- cbind(zL, zU, a, b, c1)
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
    checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 3)
    
    ## extract vectors as needed
    zL <- x[, 2]
    zU <- x[, 3]
    a <- x[, 4]
    b <- x[, 5]
    c1 <- x[, 6]
    
    u <- runif(n, 0, 1)
    w0 <- rep(NA, n)
    rs <- rep(NA, n)
    
    zN_ind <- which(checkPass)
    if(length(zN_ind) > 0) {
        w0[zN_ind] <- (a[zN_ind] / c1[zN_ind]) * exp(a[zN_ind] / c1[zN_ind]) * ((1 - u[zN_ind])^(-b[zN_ind] / c1[zN_ind]))
        w0[zN_ind] <- lambertW0(w0[zN_ind])
        rs[zN_ind] <- (a[zN_ind] / (b[zN_ind] * c1[zN_ind])) - (1 / c1[zN_ind]) * log(1 - u[zN_ind]) - (1 / b[zN_ind]) *  w0[zN_ind]
    } else {
        return(rs)
    }
  
    #interval truncation
    zI_ind <- which(!is.na(zL) & !is.na(zU) & checkPass)
    if(length(zI_ind) > 0) {
        logS_zU <- -c1[zI_ind] * zU[zI_ind] - (a[zI_ind] / b[zI_ind]) * (exp(b[zI_ind] * zU[zI_ind]) - 1)
        S_zU <- exp(logS_zU)
        CDF_zU <- 1 - S_zU
            
        logS_zL <- -c1[zI_ind] * zL[zI_ind] - (a[zI_ind] / b[zI_ind]) * (exp(b[zI_ind] * zL[zI_ind]) - 1)
        S_zL <- exp(logS_zL)
        CDF_zL <- 1 - S_zL
                
        u[zI_ind] <- (u[zI_ind] * (CDF_zU - CDF_zL)) + CDF_zL
        w0[zI_ind] <- (a[zI_ind] / c1[zI_ind]) * exp(a[zI_ind] / c1[zI_ind]) * ((1 - u[zI_ind])^(-b[zI_ind] / c1[zI_ind]))
        w0[zI_ind] <- lambertW0(w0[zI_ind])
        rs[zI_ind] <- (a[zI_ind] / (b[zI_ind] * c1[zI_ind])) - (1 / c1[zI_ind]) * log(1 - u[zI_ind]) - (1 / b[zI_ind]) *  w0[zI_ind]
    }
  
    #left truncation
    zL_ind <- which(!is.na(zL) & is.na(zU) & checkPass)
    if(length(zL_ind) > 0) {
        logS_zL <-  -c1[zL_ind] * zL[zL_ind] - (a[zL_ind] / b[zL_ind]) * (exp(b[zL_ind] * zL[zL_ind]) - 1)
        S_zL <- exp(logS_zL)
        CDF_zL <- 1 - S_zL
              
        u[zL_ind] <- u[zL_ind] * S_zL + CDF_zL
        w0[zL_ind] <- (a[zL_ind] / c1[zL_ind]) * exp(a[zL_ind] / c1[zL_ind]) * ((1 - u[zL_ind])^(-b[zL_ind] / c1[zL_ind]))
        w0[zL_ind] <- lambertW0(w0[zL_ind])
        rs[zL_ind] <- (a[zL_ind] / (b[zL_ind] * c1[zL_ind])) - (1 / c1[zL_ind]) * log(1 - u[zL_ind]) - (1 / b[zL_ind]) *  w0[zL_ind]
    }
            
    #right truncation
    zU_ind <- which(is.na(zL) & !is.na(zU) & checkPass)
    if(length(zU_ind) > 0) {
        logS_zU <-  -c1[zU_ind] * zU[zU_ind] - (a[zU_ind] / b[zU_ind]) * (exp(b[zU_ind] * zU[zU_ind]) - 1)
        S_zU <- exp(logS_zU)
        CDF_zU <- 1 - S_zU
                
        u[zU_ind] <- u[zU_ind] * CDF_zU
        w0[zU_ind] <- (a[zU_ind] / c1[zU_ind]) * exp(a[zU_ind] / c1[zU_ind]) * ((1 - u[zU_ind])^(-b[zU_ind] / c1[zU_ind]))
        w0[zU_ind] <- lambertW0(w0[zU_ind])
        rs[zU_ind] <- (a[zU_ind] / (b[zU_ind] * c1[zU_ind])) - (1 / c1[zU_ind]) * log(1 - u[zU_ind]) - (1 / b[zU_ind]) *  w0[zU_ind]
    }
    
    ## return correctly
    return(rs)
}
