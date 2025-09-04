## Create custom Siler distribution
library(lamW)

RLW0 <- nimbleRcall(function(x = double(0)){}, Rfun = 'lambertW0',
                    returnType = double(0))

## probability density function
dSilerLB <- nimbleFunction(
  run = function(x = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0), lowerBound = double(0, default = 0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    if(a1 < 0 | a2 < 0 | b1 < 0 | b2 < 0 | c < 0) {
      return(NaN)
    }
    x <- x - lowerBound
    logS <- (a1 / b1) * (exp(-b1 * x) - 1) - c * x + (a2 / b2) * (1 - exp(b2 * x))
    logH <- log(a1 * exp(-b1 * x) + c + a2 * exp(b2 * x))
    logProb <- logH + logS
    if(log) return(logProb)
    else return(exp(logProb))
  })

## function to produce random samples
rSilerLB <- nimbleFunction(
  run = function(n = integer(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0), lowerBound = double(0, default = 0)) {
    returnType(double(0))
    if(a1 < 0 | a2 < 0 | b1 < 0 | b2 < 0 | c < 0) {
      return(NaN)
    }
    if(n != 1) print("rsiler only allows n = 1; using n = 1.")
    ## sample from two independent distributions and take minimum
    u <- runif(1, 0, 1)
    w0 <- (a1 / c) * exp((log(u) + a1 / b1) * (b1 / c))
    w0 <- RLW0(w0)
    x1 <- (-1 / c) * (log(u) + a1 / b1) + w0 / b1
    x2 <- log(1 - log(runif(1, 0, 1)) * (b2 / a2)) / b2
    xm <- min(x1, x2)
    xm <- lowerBound + xm
    return(xm)
  })

## cumulative distribution function (and survivor function)
pSilerLB <- nimbleFunction(
  run = function(q = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0), lowerBound = double(0, default = 0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(a1 < 0 | a2 < 0 | b1 < 0 | b2 < 0 | c < 0) {
      return(NaN)
    }
    q <- q - lowerBound
    logS <- (a1 / b1) * (exp(-b1 * q) - 1) - c * q + (a2 / b2) * (1 - exp(b2 * q))
    if(!lower.tail) { 
      if(log.p) return(logS)
      else return(exp(logS))
    } else {
      p <- 1 - exp(logS)
      if(!log.p) return(p)
      else return(log(p))
    }
  })

## quantile function (not yet working)
qSilerLB <- nimbleFunction(
  run = function(p = double(0), a1 = double(0),
                 a2 = double(0), b1 = double(0), b2 = double(0), 
                 c = double(0), lowerBound = double(0, default = 0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(a1 < 0 | a2 < 0 | b1 < 0 | b2 < 0 | c < 0) {
      return(NaN)
    }
    if(log.p) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    print("qsiler() not specified")
    return(NaN)
  })

## register distributions with NIMBLE
registerDistributions(list(
  dSilerLB = list(
    BUGSdist = "dSilerLB(a1, a2, b1, b2, c, lowerBound)",
    Rdist = "dSilerLB(a1, a2, b1, b2, c, lowerBound)",
    pqAvail = TRUE, 
    range = c(0, Inf)
  )
))
