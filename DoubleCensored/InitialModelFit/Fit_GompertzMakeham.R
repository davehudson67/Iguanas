library(haven)
library(survival)
library(tidyverse)
library(lubridate)
library(nimble)
library(purrr)
library(mclust)
library(GGally)
library(boot)
library(parallel)
#library(mcmcplots)
rm(list=ls())

source("ModelComparison_FUNCTIONS.R")

#igs <- read_csv("Data/LeafCayThru2019_CleanOct2020.csv")
#load("DoubleCensored/InitialModelFit/LeafClayDataImage_only.RData")
load("DoubleCensored/igs_AllIslands_CleanDH_040925_obsStart.RData")

## load custom distributions
source("Distributions/Dist_GompertzMakehamLB.R")
source("Distributions/Dist_GompzMake.R")

## code for NIMBLE model with censoring
code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood 
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    
    ## Left truncation: must survive to reach time 0 if born before 0.
    L[i] <- max(0, -tB[i])
    
    ## Time-to-death since birth: Exponential, conditional on surviving to L[i].
    tstar[i] ~ dGompertzMakehamLB(a, b, c1, lowerBound = L[i]) 
    
    tD[i] <- tB[i] + tstar[i]
    
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    ## sampling component
    nm[i] <- max(ceiling(tB[i]), 0)
    nM[i] <- min(floor(tD[i]) - nm[i], tMax - nm[i]) + 1
    pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a ~ dexp(1)
  b ~ dexp(1)
  c1 ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

## set up other components of model
consts <- list(nind = nind, tMax = tMax)
data <- list(y = y, cintB = cintB, cintD = cintD,
             censoredD = censoredD, 
             tD = tD, tB = tB, tstar = tstar, dind = dind)

## set initial values
initFn <- function(model, cintB, cintD, censoredD) {
  # browser()
  repeat {
    a <- rexp(1)
    b <- rexp(1)
    c1 <- rexp(1)
    mean.p <- runif(1)
    
    n <- nrow(cintB)
    tBinit <- runif(n, cintB[,1], cintB[,2])
    Linit  <- pmax(0, -tBinit)
    
    tDinit <- numeric(nrow(cintD))
    id_cens <- which(censoredD == 2)
    id_int  <- which(censoredD != 2)
    
    if (length(id_cens)) {
      tDinit[id_cens] <- cintD[id_cens, 2] + rexp(length(id_cens), rate = 1)  # small bump
    }
    if (length(id_int)) {
      tDinit[id_int]  <- runif(length(id_int), cintD[id_int, 1], cintD[id_int, 2])
    }
    
    tiny <- 1e-6
    tDinit <- pmax(tDinit, tBinit + Linit + tiny)
    tstarinit <- tDinit - tBinit
    
    tstarinit <- pmax(tDinit - tBinit, tiny)
    
    inits <- list(
      tB      = tBinit,
      tD      = tDinit,
      tstar   = tstarinit,
      a = a,
      b = b,
      c1 = c1,
      mean.p  = mean.p
    )
    
    model$setInits(inits)
    val <- model$calculate()
    if (is.finite(val)) return(inits)
  }
}



## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data)

## compile the model
cModel <- compileNimble(model)

## find list of valid initial values using compiled model
inits <- list()
for(k in 1:2) {
  inits[[k]] <- initFn(cModel, cintB, cintD, censoredD)
}

## try with default sampler
config <- configureMCMC(cModel, monitors = c("a", "b", "c1", "mean.p"), thin = 1)
#config$removeSamplers(c("a1", "a2", "b1", "b2", "c"))
#config$addSampler(target = c("a1", "a2", "b1", "b2", "c"), type = 'AF_slice')
#config$addSampler(target = c("a2", "c", "b2"), type = 'AF_slice')

#Check monitors and samplers
config$printMonitors()
#config$printSamplers(c("a1", "a2", "b1", "b2", "c"))

#Build the model
built <- buildMCMC(config)
cBuilt <- compileNimble(built)

#Run the model
system.time(run <- runMCMC(cBuilt, 
                           niter = 50000, 
                           nburnin = 9000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run, "DoubleCensored/InitialModelFit/Samples/GM_AllIslands.rds")
run <- readRDS("Samples/GM_dc_LeafClay.rds")
#run$summary

#mcmcplot(run$samples, parms = c("a", "b", "c1", "mean.p"))

## extract birth and death times
samples <- run$samples

#tBpost <- as.matrix(samples) %>%
#  as_tibble() %>%
#  select(starts_with("tB"))
#tDpost <- as.matrix(samples) %>%
#  as_tibble() %>%
#  select(starts_with("tD"))

## examine first few just to check
#tBpost[, 1:20] %>%
#  gather(ind, value) %>%
#  ggplot(aes(x = value)) +
#  geom_density() +
#  facet_wrap(~ ind)

#lifespan_post <- tDpost[,1:20] - tBpost[,1:20]
#mcmcplot(lifespan_post)

################################################################################
#                                                                              #
#        Importance Sampling                                                   #
#                                                                              #
################################################################################
samples <- as.matrix(samples) %>%
  as_tibble() %>%
  select(c("a", "b", "c1", "mean.p"))

## pairs plot
samples <- as.matrix(samples)
samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), ]
samples %>%
  as.data.frame() %>%
  ggpairs()

## fit range of finite mixture models
mod <- densityMclust(samples)

## summary of finite mixture models
summary(mod)
plot(mod, what = "BIC")

## take random samples from mixture
nimp <- 20000
mixind <- rbinom(nimp, size = 1, prob = 0.95)
props <- matrix(NA, nimp, 4)
props[mixind == 1, ] <- sim(mod$modelName, mod$parameters, sum(mixind))[, -1]
colnames(props) <- c("a", "b", "c1", "mean.p")

## take random samples from prior (to create defense mixture)
defense <- runif(sum(mixind == 0), 0, 1)
defense <- cbind(rexp(sum(mixind == 0), 1), rexp(sum(mixind == 0), 1), rexp(sum(mixind == 0), 1), defense)
colnames(defense) <- c("a", "b", "c1", "mean.p")
props[mixind == 0, ] <- defense

## check IS distribution against posterior samples
as.data.frame(props[mixind == 1, ]) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:4)

##########################################
## calculate log likelihood

# Define the log-likelihood function
loglike <- function(cintB, cintD, y, tMax, a, b, c1, p) {
  browser()
  # Sample birth times
  tB <- runif(nrow(cintB), cintB[, 1], cintB[, 2])
  L <- pmax(0, -tB)
  
  # Sample death times from Gompertz-Makeham truncated at upper bound
  tstar <- purrr::map2_dbl(tB, cintD[, 2], function(tB_i, tM_i, a, b, c1) {
    rGompertzMakehamLB(1, a, b, c1, tM_i - tB_i)
  }, a = a, b = b, c1 = c1)
  
  tD <- tB + tstar
  
  # Log PDF for Gompertz-Makeham
  log_pdf <- log(c1 + a * exp(b * tstar)) - c1 * tstar - (a / b) * (exp(b * tstar) - 1)
  
  # Truncation adjustment: log(S(L))
  log_trunc <- - c1 * L - (a / b) * (exp(b * L) - 1)
  
  log_surv <- log_pdf - log_trunc
  
  # Survival function
  S <- function(t) exp(- c1 * t - (a / b) * (exp(b * t) - 1))
  
  # Censoring component
  log_censor <- ifelse(
    is.finite(cintD[, 2]),
    log(pmax(S(cintD[, 1] - tB) - S(cintD[, 2] - tB), .Machine$double.eps)),
    log(pmax(S(cintD[, 1] - tB), .Machine$double.eps))
  )
  
  # Sampling component
  nm <- pmax(ceiling(tB), 0)
  nM <- pmin(floor(tD) - nm, tMax - nm) + 1
  stopifnot(all(nM >= y))
  
  lpd <- y * log(p) + (nM - y) * log(1 - p)
  
  sum(log_surv + log_censor + lpd)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, cintB, cintD, y, tMax) {
                           loglike(cintB, cintD, y, tMax, pars[1], pars[2], pars[3], pars[4])
                         }, cintB = cintB, cintD = cintD,
                         y = y, tMax = tMax)
logimpweight <- reduce(logimpweight, base::c)

## add prior densities
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) +
  dexp(props[, 2], 1, log = TRUE) +
  dexp(props[, 3], 1, log = TRUE) +
  dunif(props[, 4], 0, 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(props, mod$modelName, mod$parameters, FALSE) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) +
                                                                              dexp(props[, 2], 1, log = TRUE) +
                                                                              dexp(props[, 3], 1, log = TRUE) +
                                                                              dunif(props[, 4], 0, 1, log = TRUE)))

saveRDS(logimpweight, "DoubleCensored/LImpWeights/logimpweight_gm.rds")
#logimpweight <- readRDS("DoubleCensored/LImpWeights/logimpweight_GM_KA_ALL.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## bootstrap the importance weights to create 95% intervals
imp.boot <- BootsPlot(logimpweight, 5000, TRUE)



