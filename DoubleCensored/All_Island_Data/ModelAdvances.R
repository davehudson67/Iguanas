code <- nimbleCode({
  for (i in 1:nind) {
    ## Birth time can be before study start (time 0).
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    
    ## Left truncation: must survive to reach time 0 if born before 0.
    L[i] <- max(0, -tB[i])
    
    ## Time-to-death since birth: Gompertz, conditional on surviving to L[i].
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    
    ## Absolute death time
    tD[i] <- tB[i] + tstar[i]
    
    ## Interval/right-censoring for death
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    ## Linear predictors
    log(amult[i]) <- log(a) +
      betaSEX[1]  * sex[i]     * zSEX[1] +
      betaFEED[1] * feeding[i] * zFEED[1] +
      u_a_center[year_index[i]] +
      betaSPECIESa[1, species[i]]  # <- ADDED species random effect for 'a'
    
    log(bmult[i]) <- log(b) +
      betaSEX[2]  * sex[i]     * zSEX[2] +
      betaFEED[2] * feeding[i] * zFEED[2] +
      betaSPECIES[2, species[i]]  # <- ADDED species random effect for 'b'
    
    ## Observation model (binomial kernel via ones trick)
    nm[i] <- max(ceiling(tB[i]), 0)
    nM[i] <- min(floor(tD[i]) - nm[i], tMax - nm[i]) + 1
    pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## Year effect: RW1 on 'a' (cohort effect given current indexing)
  u_a[1] ~ dnorm(0, sd = 10)
  for (t in 2:n_year) {
    u_a[t] ~ dnorm(u_a[t-1], sd = sigma_year_a)
  }
  u_a_bar <- mean(u_a[1:n_year])
  for (t in 1:n_year) {
    u_a_center[t] <- u_a[t] - u_a_bar
  }
  
  ## Priors
  for (k in 1:2) {
    betaSEX[k]  ~ dnorm(0, sd = 1)
    betaFEED[k] ~ dnorm(0, sd = 1)
    zSEX[k]     ~ dbern(0.5)
    zFEED[k]    ~ dbern(0.5)
  }
  
  # ## NEW: Priors for species random effects
  for (j in 1:n_species) {
    betaSPECIES[1, j] ~ dnorm(0, sd = sigma_species_a)
    betaSPECIES[2, j] ~ dnorm(0, sd = sigma_species_b)
  }
  
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
  sigma_year_a ~ dunif(0, 5)
  
  # ## NEW: Priors for species standard deviations
  sigma_species_a ~ dexp(1)
  sigma_species_b ~ dexp(1)
})

censoredD <- rep(2, nind)
censoredD[cintD[,1] > 0] <- 1
n_year <- max(year_index)

## set up other components of model
consts <- list(nind = nind, tMax = ncol(CH), sex = sex, feeding = feeding, species = species, n_species = 2, n_year = n_year, year_index = year_index)
data <- list(y = y, cintB = cintB, cintD = cintD,
             censoredD = censoredD, 
             tD = tD, tB = tB, tstar = tstar, dind = dind)

initFn <- function(model, cintB, cintD, censoredD, n_year) {
  #browser()
  repeat {
    a <- rexp(1)
    b <- rexp(1)
    mean.p <- runif(1)
    tBinit <- runif(nrow(cintB), cintB[,1], cintB[,2])
    Linit  <- pmax(0, -tBinit)
    
    tDinit <- numeric(nrow(cintD))
    id_cens <- which(censoredD == 2)
    id_int  <- which(censoredD != 2)
    
    if (length(id_cens)) {
      tDinit[id_cens] <- cintD[id_cens, 2] + rexp(length(id_cens), rate = 5)  # small bump
    }
    if (length(id_int)) {
      tDinit[id_int]  <- runif(length(id_int), cintD[id_int, 1], cintD[id_int, 2])
    }
    
    tiny <- 1e-6
    tDinit <- pmax(tDinit, tBinit + Linit + tiny)
    tstarinit <- tDinit - tBinit
    
    inits <- list(
      tB = tBinit, tD = tDinit, tstar = tstarinit,
      a = a, b = b, mean.p = mean.p,
      betaSEX = rnorm(2), betaFEED = rnorm(2), betaSPECIES = matrix(rnorm(2 * 2), nrow = 2),
      u_a = rep(0, n_year), sigma_year_a = 0.5, sigma_species_a = rexp(1),
      sigma_species_b = rexp(1), zSEX = c(0,0), zFEED = c(0,0)
    )
    
    model$setInits(inits)
    if (is.finite(model$calculate())) return(inits)
  }
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data)
#model$initializeInfo()

## compile the model
cModel <- compileNimble(model)

## find list of valid initial values using compiled model
inits <- list()
for(k in 1:2) {
  inits[[k]] <- initFn(cModel, cintB, cintD, censoredD, n_year)
}

## configure MCMC
config <- configureMCMC(model, monitors = c("a", "b", "mean.p", "betaSEX", "betaFEED", "betaSPECIES", 
                                            "zSEX", "zFEED", "u_a", "sigma_year_a", "sigma_species_a", "sigma_species_b"))
config$removeSamplers(c("a", "b"))
config$addSampler(target = c("a", "b"), type = 'AF_slice')
#config$addSampler(target = c("betaFEED"), type = 'AF_slice')
#config$addSampler(target = c("b"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
#config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))

## load in custom RJ-MCMC samplers
#source("MCMC_RJ_multi.R")

# Add reversible jump
configureRJ(conf = config,   ## model configuration
            targetNodes = c("betaSEX", "betaFEED"),
            indicatorNodes = c("zSEX", "zFEED"),
            control = list(mean = 0, scale = 0.5))

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

run <- runMCMC(cIndicatorMCMC,
               niter = 250000,
               nburnin = 22000,
               thin = 1,
               nchains = 2, 
               progressBar = TRUE,
               samplesAsCodaMCMC = TRUE)

samples <- as.matrix(run, chains = TRUE)
saveRDS(samples, "DoubleCensored/All_Island_Data/Samples/G_AllIslands_SexFeeding2Species_RETime_a.rds")
out <- readRDS("DoubleCensored/All_Island_Data/Samples/G_AllIslands_SexFeeding2_RETime_a.rds")
plot(as.mcmc(run))
plot(run)
