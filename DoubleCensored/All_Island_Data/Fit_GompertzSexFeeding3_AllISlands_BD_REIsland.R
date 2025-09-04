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
#source("DoubleCensored/Prep_All_Data.R")
load("igs_AllIslands_CleanDH_190625.RData")

## load custom distributions
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

#AllData <- AllData[1:200,]
#CH <- CH[1:200,]
#cintB <- cintB[1:200,]
#cintD <- cintD[1:200,]
#censoredB <- censoredB[1:200]
#censoredD <- censoredD[1:200]
#dind <- dind[1:200]
#tB <- tB[1:200]
#tD <- tD[1:200]
#tF <- tF[1:200]
#tKB <- tKB[1:200]
#tKD <- tKD[1:200]
#tL <- tL[1:200]
#tstar <- tstar[1:200]
#y <- y[1:200]
#nind <- 200

sex <- AllData %>%
  mutate(sex = as.numeric(as.factor(sex))) %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(sex) - 1

# sex 0 = female
# sex 1 = male

island <- AllData %>%
  mutate(island = as.numeric(as.factor(Island))) %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(island)

#feeding_counts <- AllData %>%
#  group_by(ID) %>%
#  summarise(n_feeding_levels = n_distinct(feeding)) %>%
#  arrange(desc(n_feeding_levels))

#individuals_multiple_feeding <- feeding_counts %>%
#  filter(n_feeding_levels > 1, .preserve = TRUE)

# Function to calculate the mode (most frequent value)
#get_mode <- function(x) {
#  tab <- table(x)  # Count occurrences
#  mode_value <- names(tab[tab == max(tab)])  # Get the most frequent value(s)
#  
#  if (length(mode_value) > 1) {
#    return(mode_value[1])  # If there's a tie, return the first mode
#  } else {
#    return(mode_value)
#  }
#}

# Add a new column with the mode of feeding per individual
#feedingM <- AllData %>%
#  group_by(ID) %>%
#  summarise(feedingM = get_mode(feeding)) %>%
#  ungroup()

# Merge back into AllData
#AllData <- left_join(AllData, feedingM, by = "ID")

# Extract for the model
#feeding <- AllData %>%
#  #mutate(feeding2 = as.factor(feeding)) %>%
#  mutate(feeding = ifelse(feedingM == "none", 0, 1)) %>%
#  distinct(ID, .keep_all = TRUE) %>%
#  pull(feeding)

# feeding 0 = none
# feeding 1 = high or moderate

summary(as.factor(AllData$feeding))
feeding <- as.factor(AllData$feeding)
feeding <- model.matrix(~ feeding, feeding)[, -1]
colSums((feeding))

## code for NIMBLE model with censoring
code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood 
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    tD[i] <- tB[i] + tstar[i]
    tstar[i] ~ dgompzNim(amult[i], bmult[i])
    
    log(amult[i]) <- log(a) + betaSEX[1] * sex[i] * zSEX[1] +
      betaFEEDmid[1] * feeding[i, 1] * zFEED[1] + 
      betaFEEDno[1] * feeding[i, 2] * zFEED[1] +
      gamma_island_a[island[i]]
    
    log(bmult[i]) <- log(b) + betaSEX[2] * sex[i] * zSEX[2] +
      betaFEEDmid[2] * feeding[i, 1]  * zFEED[2] + 
      betaFEEDno[2] * feeding[i, 2] * zFEED[2] +
      gamma_island_b[island[i]]
    
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
  
  for (j in 1:n_island) {
    gamma_island_a[j] ~ dnorm(0, sd = sigma_island_a)
    gamma_island_b[j] ~ dnorm(0, sd = sigma_island_b)
  }
  
  sigma_island_a ~ dunif(0, 5)
  sigma_island_b ~ dunif(0, 5)
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

censoredD <- rep(2, nind)
censoredD[cintD[,1] > 0] <- 1
n_island <- max(island)


## set up other components of model
consts <- list(nind = nind, tMax = ncol(CH), sex = sex, feeding = feeding, island = island, n_island = n_island)
data <- list(y = y, cintB = cintB, cintD = cintD,
             censoredD = censoredD, 
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
    tDinit <- ifelse(censoredD == 2, 
                     tBinit + (cintD[,2] + rexp(nrow(cintD), 1)), 
                     # For interval-censored individuals: sample tD uniformly within the censoring interval
                     runif(nrow(cintD), cintD[,1], cintD[,2]))
    tstarinit <- tDinit - tBinit
    betaSEX <- rnorm(2, 0, 1)
    betaFEEDmid <- rnorm(2, 0, 1)
    betaFEEDno <- rnorm(2, 0, 1)
    sigmaA <- runif(1, 0, 5)
    sigmaB <- runif(1, 0, 5)
    gamma_island_a <- rnorm(8, 0, 1)
    gamma_island_b <- rnorm(8, 0, 1)
    sigma_island_a <- runif(1, 0, 5)
    sigma_island_b <- runif(1, 0, 5)
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
      gamma_island_a = gamma_island_a,
      gamma_island_b = gamma_island_b,
      sigma_island_a = sigma_island_a,
      sigma_island_b = sigma_island_b,
      zSEX = zSEX,
      zFEED = zFEED
    )
    model$setInits(inits)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
  }
  return(inits)
}


## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data)
#model$initializeInfo()

## compile the model
cModel <- compileNimble(model)

## find list of valid initial values using compiled model
inits <- list()
for(k in 1:2) {
  inits[[k]] <- initFn(cModel, cintB, cintD)
}

## configure MCMC
config <- configureMCMC(model, monitors = c("a", "b", "mean.p", "betaSEX", "betaFEEDmid", "betaFEEDno", "zSEX", "zFEED", "gamma_island_a", "gamma_island_b", "sigma_island_a", "sigma_island_b"))
config$removeSamplers(c("a", "b"))
config$addSampler(target = c("a", "b"), type = 'AF_slice')
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
                           niter = 150000, 
                           nburnin = 19000, 
                           nchains = 2,
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

plot(run$samples)
mcmcplot(run$samples)
run$summary
saveRDS(run, "DoubleCensored/All_Island_Data/Samples/G_AllIslands_SexFeeding3_REIsland.rds")

