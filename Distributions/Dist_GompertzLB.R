## probability density function
dGompertzLB <- nimbleFunction(
  run = function(x = double(0), a = double(0),
                 b = double(0), lowerBound = double(0, default = 0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    if(a < 0 | b < 0) {
      return(NaN)
    }
    x <- x - lowerBound
    logS <- (a / b) * (1 - exp(b * x))
    logH <- log(a) + (b * x)
    logProb <- logH + logS
    if(log) return(logProb)
    else return(exp(logProb))
  })

## function to produce random samples
rGompertzLB <- nimbleFunction(
  run = function(n = integer(0), a = double(0),
                 b = double(0), lowerBound = double(0, default = 0)) {
    returnType(double(0))  # Return type is now a scalar double
    if (a < 0 | b < 0) {
      return(NaN)
    }
    if (n != 1) nimPrint("rGompertzLB only allows n = 1; using n = 1.")
    
    ## generate a single uniform random number
    u <- runif(1, 0, 1)
    
    ## compute the Gompertz sample
    rs <- (1 / b) * log(1 - (b / a) * log(1 - u))
    rs <- lowerBound + rs
    
    return(rs)  # Return a scalar double
  }
)

## cumulative distribution function (and survivor function)
pGompertzLB <- nimbleFunction(
  run = function(q = double(0), a = double(0),
                 b = double(0), lowerBound = double(0, default = 0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(a < 0 | b < 0) {
      return(NaN)
    }
    q <- q - lowerBound
    logS <- (a / b) * (1 - exp(b * q))
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
qGompertzLB <- nimbleFunction(
  run = function(p = double(0), a = double(0),
                 b = double(0), lowerBound = double(0, default = 0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(a < 0 | b < 0) {
      return(NaN)
    }
    if(log.p) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    print("qGompertz() not specified")
    return(NaN)
  })

## register distributions with NIMBLE
registerDistributions(list(
  dGompertzLB = list(
    BUGSdist = "dGompertzLB(a, b, lowerBound)",
    Rdist = "dGompertzLB(a, b, lowerBound)",
    pqAvail = TRUE, 
    range = c(0, Inf)
  )
))

