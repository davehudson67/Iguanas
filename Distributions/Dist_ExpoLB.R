library(nimble)

## probability density function
dExponentialLB <- nimbleFunction(
  run = function(x = double(0), r = double(0, default = 1),
                 lowerBound = double(0, default = 0),
                 log = integer(0, default = 0)) {
    returnType(double(0))

    x <- x - lowerBound
    logS <- -r * x
    logH <- log(r)
    logProb <- logH + logS
    if(log) return(logProb)
    else return(exp(logProb))
  })


## cumulative distribution function (and survivor function)
pExponentialLB <- nimbleFunction(
  run = function(q = double(0), r = double(0, default = 1),
                 lowerBound = double(0, default = 0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    
    q <- q - lowerBound
    logS <- -r * q
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
qExponentialLB <- nimbleFunction(
  run = function(p = double(0), r = double(0, default = 1), 
                 lowerBound = double(0, default = 0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(log.p) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    print("qExponential() not specified")
    return(NaN)
  })

## function to produce random samples
rExponentialLB <- nimbleFunction(
  run = function(n = integer(0), r = double(0, default = 1),
                 lowerBound = double(0, default = 0)) {
    returnType(double(0))  # Return type is now a scalar double
    if (n != 1) nimPrint("rExponentialLB only allows n = 1; using n = 1.")
    
    ## generate a single uniform random number
    u <- runif(1, 0, 1)
    
    ## compute the Gompertz sample
    rs <- -log(1 - u) / r
    rs <- lowerBound + rs
    
    return(rs)  # Return a scalar double
  }
)

## register distributions with NIMBLE
registerDistributions(list(
  dExponentialLB = list(
    BUGSdist = "dExponentialLB(r, lowerBound)",
    Rdist = "dExponentialLB(r, lowerBound)",
    pqAvail = TRUE, 
    range = c(0, Inf)
  )
))
