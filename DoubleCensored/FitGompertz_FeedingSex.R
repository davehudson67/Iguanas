##########################################################
##                                                      ###
##       Feeding and Sex effects - Gompertz RJMCMC      ###
##                                                      ###
###########################################################

## load libraries
library(nimble)
library(tidyverse)
library(mvtnorm)
library(boot)
library(lamW)
library(GGally)
library(coda)
library(mclust)
library(parallel)
library(survminer)
library(survival)
library(coda)
library(mcmcplots)
#library(MCMCvis)
library(scales)
library(data.table)

rm(list=ls())

load("igs_ready_AllIslands_KA.RData")

## load custom distributions
source("Distributions/Dist_GompertzLB.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

## set seed
set.seed(11)

sex <- igsAll %>%
  mutate(sex = as.numeric(as.factor(sex))) %>%
  distinct(animal_id, .keep_all = TRUE) %>%
  pull(sex)

feeding <- igsAll %>%
  mutate(feeding = as.factor(igsAll$feeding)) %>%
  distinct(animal_id, .keep_all = TRUE) %>%
  pull(feeding)

summary(feeding)
feeding <- model.matrix(~ feeding, feeding)[, -1]
  
code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood 
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    censoredD[i] ~ dinterval(tD[i], cintD[i])
    tD[i] <- tB[i] + tstar[i]
    tstar[i] ~ dgompzNim(amult[i], bmult[i])
    
    log(amult[i]) <- log(a) + betaSEX[1] * sex[i] * zSEX[1] +
      #betaFEEDmid[1] * feeding[i, 1] * zFEED[1] + 
      #betaFEEDno[1] * feeding[i, 2] * zFEED[1] +
      gammaISLAND[island[i]]
    
    log(bmult[i]) <- log(b) + betaSEX[2] * sex[i] * zSEX[2] +
      #betaFEEDmid[2] * feeding[i, 1]  * zFEED[2] + 
      #betaFEEDno[2] * feeding[i, 2] * zFEED[2] +
      gammaISLAND[island[i]]
    
    ## sampling component
    nm[i] <- max(ceiling(tB[i]), 0)
    nM[i] <- min(floor(tD[i]) - nm[i], tMax - nm[i]) + 1
    pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
    
  }
  
  ## priors
  for (k in 1:2) {
    betaSEX[k] ~ dnorm(0, sd = 1)
    betaFEEDmid[k] ~ dnorm(0, sd = 1)
    betaFEEDno[k] ~ dnorm(0, sd = 1)
    zSEX[k] ~ dbern(0.5)
    zFEED[k] ~ dbern(0.5)
  }  
  
  for (r in 1:nISLANDS){
    gammaISLAND[r] ~ dnorm(0, sd = sigma)
  }
  
  sigma ~ dunif(0, 10)
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

## set up data
consts <- list(nind = nind, tMax = ncol(CH), nISLANDS = 3, sex = (sex - 1), feeding = feeding, island = as.numeric(as.factor(igsAll$island)))

data <- list(y = y, cintB = cintB, cintD = cintD[, 2],
             censoredD = rep(1, nind), 
             tD = tD, tB = tB, tstar = tstar, dind = dind)

## set initial values
initFn <- function(model, cintB, cintD) {
  #browser()
  valid <- 0
  while(valid == 0) {
    a <- rexp(1)
    b <- rexp(1)
    mean.p <- runif(1, 0, 1)
    tBinit <- runif(nrow(cintB), cintB[, 1], cintB[, 2])
    tstarinit <- cintD[, 2] + rexp(nrow(cintD), 1)
    tDinit <- tB + tstarinit
    betaSEX <- rnorm(2, 0, 1)
    betaFEEDmid <- rnorm(2, 0, 1)
    betaFEEDno <- rnorm(2, 0, 1)
    gammaISLAND <- rep(0, 3)
    sigma <- runif(1, 0, 5)
    zSEX <- rep(0, 2)
    zFEED <- rep(0, 2)
    inits <- list(
      tD = tDinit,
      tstar = tstarinit,
      tB = tBinit,
      a = a,
      b = b,
      mean.p = mean.p,
      betaSEX = betaSEX,
      betaFEEDmid = betaFEEDmid,
      betaFEEDno = betaFEEDno,
      gammaISLAND = gammaISLAND,
      sigma = sigma,
      zSEX = zSEX,
      zFEED = zFEED
    )
    model$setInits(inits)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
  }
  return(inits)
}

## build the model without initial values
## (will throw an initialisation warning)
model <- nimbleModel(code, constants = consts, data = data)

## compile the model
cIndicatorModel <- compileNimble(model)

## find list of valid initial values (needs compiled model
## for some reason)
inits <- list()
for(k in 1:2) {
  inits[[k]] <- initFn(cIndicatorModel, cintB, cintD)
}

## configure MCMC
config <- configureMCMC(model, monitors = c("a", "b", "sigma", "mean.p", "betaSEX", "betaFEEDmid", "betaFEEDno", "zSEX", "zFEED"))
#config$removeSamplers(c("a", "b", "c1"))
#config$addSampler(target = c("a"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
#config$addSampler(target = c("b"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
#config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))

## load in custom RJ-MCMC samplers
source("MCMC_RJ_multi.R")

# Add reversible jump
configureRJ_multi(conf = config,   ## model configuration
                  targetNodes = c("betaSEX", "betaFEEDmid", "betaFEEDno"),
                  indicatorNodes = c("zSEX", "zFEED", "zFEED"),
                  control = list(mean = 0, scale = 1))

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 300000, 
                           nburnin = 15000, 
                           nchains = 2,
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

run$summary

## save mcmc ouput
saveRDS(run, "DoubleCensored/G_FeedingSexRun.rds")
run <- readRDS("DoubleCensored/G_FeedingSexRun.rds")

samples <- as.matrix(run$samples)
#saveRDS(samples, "outputs/Sex3Infection_AllParameters_samples.rds")

mcmcplot(run$samples) #, parms = c("betaSEX", "betaFEEDmid", "betaFEEDno", "zSEX", "zFEED"))
plot(run$samples)
MCMCsummary(run$samples)
coda::traceplot(run$samples)
MCMCtrace(run$samples)

samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), ]
samples %>%
  as.data.frame() %>%
  ggpairs()

#dev.off()
#MCMCplot(run$samples)
#MCMCtrace(run$samples, pdf = F)

## save MCMC
#samples <- readRDS("outputs/Sex3Infection_AllParameters_samples.rds")

## Marginal probabilities of inclusion for each variable
zNames <- model$expandNodeNames(c('z', 'zINBR'))
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]
res
res
saveRDS(res, "outputs/FullModel_NoInteractions_z_all_PosteriorModelProbs1.rds")
res <- readRDS("outputs/FullModel_NoInteractions_z_all_PosteriorModelProbs1.rds")

samples <- as.data.frame(samples)

z_indicators <- samples %>%
  select(c(27:36)) %>%
  colSums()

z_indicators <- data.frame(z_indicators/sum(res$N))
#z_indicators$parameter <- c("a1_Inf", "a2_Inf", "b1_Inf", "b2_Inf", "c_Inf",
#                            "a1_Sex", "a2_Sex", "b1_Sex", "b2_Sex", "c_Sex",
#                            "a1_SexInf", "a2SexInf", "b1SexInf", "b2SexInf", "cSexInf")

z_indicators$z <- rownames(z_indicators)
colnames(z_indicators) <- c("Inclusion_Prob", "z")
z_indicators$variable <- rep(c("Inbreeding", "SG"), each = 5)

incl <- ggplot(z_indicators, aes(x = z, y = Inclusion_Prob)) +
  geom_point(aes(colour = variable)) + 
  geom_hline(yintercept = 0.5, colour = "red") +
  scale_y_continuous(limits = c(0,1)) + 
  labs(title = "Main effects only model | Continuous inbreeding", subtitle = "[1] = a1, [2] = a2, [3] = b1, [4] = b2, [5] = c")
incl
