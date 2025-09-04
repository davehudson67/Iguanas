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
library(mcmcplots)
rm(list=ls())

source("ModelComparison_FUNCTIONS.R")
load("igs_AllIslands_CleanDH_190625.RData")
#-----------------------------------------------------------------------------
# Run in parallel

## set up data and constants
consts <- list(nind = nind, tMax = ncol(CH))

data <- list(y = y, cintB = cintB, cintD = cintD[, 2],
             censoredD = rep(1, nind), 
             tD = tD, tB = tB, tstar = tstar, dind = dind)

siler_path <- normalizePath("Distributions/Dist_SilerLB.R")
silerNim_path <- normalizePath("../NIMBLE_Distributions/Dist_SilerNim.R")

monitors = c("a1", "a2", "b1", "b2", "c1", "mean.p", "tB", "tD")

## --- WORKER FUNCTION --------------------------------------------------------
workflow <- function(seed, dat, consts, monitors, siler_path, silerNim_path, cintB, cintD) {
  library(nimble)
  source(siler_path)
  source(silerNim_path)
  
  ## --- NIMBLE Model Code ---
  code <- nimbleCode({
    for (i in 1:nind) {
      tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
      censoredD[i] ~ dinterval(tD[i], cintD[i])
      tD[i] <- tB[i] + tstar[i]
      tstar[i] ~ dsilerNim(a1, a2, b1, b2, c1)
      
      nm[i] <- max(ceiling(tB[i]), 0)
      nM[i] <- min(floor(tD[i]) - nm[i], tMax - nm[i]) + 1
      pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
      dind[i] ~ dbern(pd[i])
    }
    
    a1 ~ dexp(1)
    a2 ~ dexp(1)
    b1 ~ dexp(1)
    b2 ~ dexp(1)
    c1 ~ dexp(1)
    mean.p ~ dunif(0, 1)
  })
  
  ## Build model
  model <- nimbleModel(code, constants = consts, data = dat)
  cModel <- compileNimble(model)
  
  ## Initial values generator
  initFn <- function(model, cintB, cintD) {
    repeat {
      inits <- list(
        a1 = rexp(1), a2 = rexp(1), b1 = rexp(1),
        b2 = rexp(1), c1 = rexp(1), mean.p = runif(1),
        tB = runif(nrow(cintB), cintB[, 1], cintB[, 2]),
        tstar = cintD[, 2] + rexp(nrow(cintD), 1)
      )
      inits$tD <- inits$tB + inits$tstar
      model$setInits(inits)
      if (is.finite(model$calculate())) return(inits)
    }
  }
  
  inits <- initFn(cModel, cintB, cintD)
  
  ## Configure and build MCMC
  config <- configureMCMC(cModel, monitors = monitors, thin = 1)
  built <- buildMCMC(config)
  cBuilt <- compileNimble(built)
  
  ## Run MCMC
  set.seed(seed)
  run <- runMCMC(cBuilt,
                 niter = 50000,
                 nburnin = 9000,
                 nchains = 1,
                 progressBar = FALSE,
                 summary = FALSE,
                 samplesAsCodaMCMC = TRUE,
                 thin = 1)
  return(run)
}

## Set up cluster
nbcores <- detectCores() - 1
cl <- makeCluster(nbcores)

## Export all needed variables and paths to workers
clusterExport(cl, c("workflow", "data", "consts", "monitors",
                    "siler_path", "silerNim_path", "cintB", "cintD"))

clusterEvalQ(cl, library(nimble))  # load nimble on each core

## Define seed vector for reproducibility
seeds <- sample(1:1e6, nbcores)

## Run in parallel
out <- parLapply(cl, seeds, workflow,
                 dat = data,
                 consts = consts,
                 monitors = monitors,
                 siler_path = siler_path,
                 silerNim_path = silerNim_path,
                 cintB = cintB,
                 cintD = cintD)

stopCluster(cl)

## Combine MCMC output and save
samples <- mcmc.list(out, as.mcmc.list)
saveRDS(samples, "Samples/S_AllIslands.rds")


#run <- readRDS("Samples/S_dc_AllIslandsKA.rds")
run$summary

mcmcplot(run$samples, parms = c("a1", "a2", "b1", "b2", "c1", "mean.p"))

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
  select(c("a1", "a2", "b1", "b2", "c1", "mean.p"))

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
props <- matrix(NA, nimp, 6)
props[mixind == 1, ] <- sim(mod$modelName, mod$parameters, sum(mixind))[, -1]
colnames(props) <- c("a1", "a2", "b1", "b2", "c1", "mean.p")

## take random samples from prior (to create defense mixture)
defense <- runif(sum(mixind == 0), 0, 1)
defense <- cbind(rexp(sum(mixind == 0), 1), rexp(sum(mixind == 0), 1), rexp(sum(mixind == 0), 1), rexp(sum(mixind == 0), 1), rexp(sum(mixind == 0), 1), defense)
colnames(defense) <- c("a1", "a2", "b1", "b2", "c1", "mean.p")
props[mixind == 0, ] <- defense

## check IS distribution against posterior samples
as.data.frame(props[mixind == 1, ]) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:6)

##########################################
## calculate log likelihood

# Define the log-likelihood function
loglike <- function(cintB, cintD, y, tMax, a1, a2, b1, b2, c1, p) {
  
  ## sample birth times across all individuals
  tB <- runif(nrow(cintB), cintB[, 1], cintB[, 2])
  
  ## sample death times across all individuals
  tstar <- map2_dbl(tB, cintD[, 2], function(tB, tM, a1, a2, b1, b2, c1) {
    rSilerLB(1, props[1], props[2], props[3], props[4], props[5], tM - tB)
  }, a1 = a1, a2 = a2, b1 = b1, b2 = b2,  c1 = c1)
  tD <- tstar + tB
  
  ## sampling component
  nm <- pmax(ceiling(tB), 0)
  nM <- pmin(floor(tD) - nm, tMax - nm) + 1
  stopifnot(all(nM >= y))
  lpd <- y * log(p) + (nM - y) * log(1 - p)
  sum(lpd)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, cintB, cintD, y, tMax) {
                           loglike(cintB, cintD, y, tMax, pars[1], pars[2], pars[3], pars[4], pars[5], pars[6])
                         }, cintB = cintB, cintD = cintD,
                         y = y, tMax = tMax, mc.cores = 8)
logimpweight <- reduce(logimpweight, base::c)

## add prior densities
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) +
  dexp(props[, 2], 1, log = TRUE) +
  dexp(props[, 3], 1, log = TRUE) +
  dexp(props[, 4], 1, log = TRUE) +
  dexp(props[, 5], 1, log = TRUE) +
  dunif(props[, 6], 0, 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(props, mod$modelName, mod$parameters, FALSE) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) +
                                                                              dexp(props[, 2], 1, log = TRUE) +
                                                                              dexp(props[, 3], 1, log = TRUE) +
                                                                              dexp(props[, 4], 1, log = TRUE) +
                                                                              dexp(props[, 5], 1, log = TRUE) +
                                                                              dunif(props[, 6], 0, 1, log = TRUE)))

saveRDS(logimpweight, "DoubleCensored/LImpWeights/logimpweight_s.rds")
#logimpweight <- readRDS("DoubleCensored/LImpWeights/logimpweight_S_KA_ALL.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## bootstrap the importance weights to create 95% intervals
imp.boot <- BootsPlot(logimpweight, 5000, TRUE)



