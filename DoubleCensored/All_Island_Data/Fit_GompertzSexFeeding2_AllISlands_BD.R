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
source("DoubleCensored/Prep_All_Data.R")

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

#island <- AllData %>%
#  mutate(island = as.numeric(as.factor(Island))) %>%
#  distinct(ID, .keep_all = TRUE) %>%
#  pull(island)

feeding_counts <- AllData %>%
  group_by(ID) %>%
  summarise(n_feeding_levels = n_distinct(feeding)) %>%
  arrange(desc(n_feeding_levels))

individuals_multiple_feeding <- feeding_counts %>%
  filter(n_feeding_levels > 1, .preserve = TRUE)

# Function to calculate the mode (most frequent value)
get_mode <- function(x) {
  tab <- table(x)  # Count occurrences
  mode_value <- names(tab[tab == max(tab)])  # Get the most frequent value(s)
  
  if (length(mode_value) > 1) {
    return(mode_value[1])  # If there's a tie, return the first mode
  } else {
    return(mode_value)
  }
}

# Add a new column with the mode of feeding per individual
feedingM <- AllData %>%
  group_by(ID) %>%
  summarise(feedingM = get_mode(feeding)) %>%
  ungroup()

# Merge back into AllData
AllData <- left_join(AllData, feedingM, by = "ID")

# Extract for the model
feeding <- AllData %>%
  #mutate(feeding2 = as.factor(feeding)) %>%
  mutate(feeding = ifelse(feedingM == "none", 1, 0)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(feeding)

# feeding 0 = none
# feeding 1 = high or moderate

summary(as.factor(AllData$feeding))

#summary(as.factor(AllData$feeding))
#feeding <- model.matrix(~ feeding, feeding)[, -1]

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
      betaFEED[1] * feeding[i] * zFEED[1] 
    log(bmult[i]) <- log(b) + betaSEX[2] * sex[i] * zSEX[2] +
      betaFEED[2] * feeding[i] * zFEED[2] 
 
       ## sampling component
    nm[i] <- max(ceiling(tB[i]), 0)
    nM[i] <- min(floor(tD[i]) - nm[i], tMax - nm[i]) + 1
    pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
    
  }
  
  ## priors
  for (k in 1:2) {
    betaSEX[k] ~ dnorm(0, sd = 1)
    betaFEED[k] ~ dnorm(0, sd = 1)
    zSEX[k] ~ dbern(0.5)
    zFEED[k] ~ dbern(0.5)
  }
  
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

censoredD <- rep(2, nind)
censoredD[cintD[,1] > 0] <- 1

## set up other components of model
consts <- list(nind = nind, tMax = ncol(CH), sex = sex, feeding = feeding)
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
    betaFEED <- rnorm(2, 0, 1)
    sigmaA <- runif(1, 0, 5)
    sigmaB <- runif(1, 0, 5)
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
      betaFEED = betaFEED,
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
config <- configureMCMC(model, monitors = c("a", "b", "mean.p", "betaSEX", "betaFEED", "zSEX", "zFEED"))
#config$removeSamplers(c("a", "b", "c1"))
#config$addSampler(target = c("a"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
#config$addSampler(target = c("b"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
#config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))

## load in custom RJ-MCMC samplers
source("MCMC_RJ_multi.R")

# Add reversible jump
configureRJ_multi(conf = config,   ## model configuration
                  targetNodes = c("betaSEX", "betaFEED"),
                  indicatorNodes = c("zSEX", "zFEED"),
                  control = list(mean = 0, scale = 1))

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

system.time(run <- runMCMC(cIndicatorMCMC, 
                           niter = 100000, 
                           nburnin = 11000, 
                           nchains = 2,
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))


saveRDS(run, "DoubleCensored/All_Island_Data/Samples/G_AllIslands_SexFeeding2.rds")

plot(run$samples)


run <- readRDS("Samples/G_dc_AllIslandsKA.rds")
run$summary

#mcmcplot(run$samples, parms = c("a", "b", "mean.p"))

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
  select(c("a", "b", "mean.p"))

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
props <- matrix(NA, nimp, 3)
props[mixind == 1, ] <- sim(mod$modelName, mod$parameters, sum(mixind))[, -1]
colnames(props) <- c("a", "b", "mean.p")

## take random samples from prior (to create defense mixture)
defense <- runif(sum(mixind == 0), 0, 1)
defense <- cbind(rexp(sum(mixind == 0), 1), rexp(sum(mixind == 0), 1), defense)
colnames(defense) <- c("a", "b", "mean.p")
props[mixind == 0, ] <- defense

## check IS distribution against posterior samples
as.data.frame(props[mixind == 1, ]) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:3)

##########################################
## calculate log likelihood

# Define the log-likelihood function
loglike <- function(cintB, cintD, y, tMax, a, b, p) {
  
  ## sample birth times across all individuals
  tB <- runif(nrow(cintB), cintB[, 1], cintB[, 2])
  
  ## sample death times across all individuals
  tstar <- map2_dbl(tB, cintD[, 2], function(tB, tM, a, b) {
    rGompertzLB(1, a, b, tM - tB)
  }, a = a, b = b)
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
                           loglike(cintB, cintD, y, tMax, pars[1], pars[2], pars[3])
                         }, cintB = cintB, cintD = cintD,
                         y = y, tMax = tMax, mc.cores = 8)
logimpweight <- reduce(logimpweight, base::c)

## add prior densities
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) +
  dexp(props[, 2], 1, log = TRUE) +
  dunif(props[, 3], 0, 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(props, mod$modelName, mod$parameters, FALSE) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) +
                                                                              dexp(props[, 2], 1, log = TRUE) +
                                                                              dunif(props[, 3], 0, 1, log = TRUE)))

saveRDS(logimpweight, "DoubleCensored/LImpWeights/logimpweight_g.rds")
#logimpweight <- readRDS("DoubleCensored/LImpWeights/logimpweight_G_KA_ALL.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## bootstrap the importance weights to create 95% intervals
imp.boot <- BootsPlot(logimpweight, 5000, TRUE)



