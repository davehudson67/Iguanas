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
library(MCMCvis)

rm(list=ls())

source("ModelComparison_FUNCTIONS.R")
#source("DoubleCensored/Prep_All_Data.R")
load("igs_AllIslands_CleanDH_190625.RData")

## load custom distributions
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

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
feeding <- AllData %>%
  mutate(feeding = as.factor(feeding)) %>%
  mutate(feeding = ifelse(feeding == "none", 0, 1)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(feeding)

# feeding 0 = none
# feeding 1 = high or moderate

#summary(as.factor(AllData$feeding))
#feeding <- as.factor(AllData$feeding)
#feeding <- model.matrix(~ feeding, feeding)[, -1]
#colSums((feeding))

#------------------------------------------------------------------------------

n_year <- tMax
t0 <- 0L

# Deterministic year index for each ID from entry-interval midpoint
tB_mid <- rowMeans(cintB)                          # midpoint of [cintB[i,1], cintB[i,2]]
year_index <- floor(tB_mid - t0) + 1L              # map to 1..n_year bins
year_index <- pmin(pmax(year_index, 1L), n_year)   # clamp

#==============================================================================#
# Random temporal effect

## code for NIMBLE model with censoring
code <- nimbleCode({
  for (i in 1:nind) {
    ## Likelihood: continuous-time event model with Gompertz survival
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    tD[i] <- tB[i] + tstar[i]
    tstar[i] ~ dgompzNim(amult[i], bmult[i])
    
    ## Linear predictors
    log(amult[i]) <- log(a) +
      betaSEX[1] * sex[i] * zSEX[1] +
      betaFEED[1] * feeding[i] * zFEED[1] +
      u_a_center[year_index[i]]
    
    log(bmult[i]) <- log(b) +
      betaSEX[2] * sex[i] * zSEX[2] +
      betaFEED[2] * feeding[i] * zFEED[2]
    
    ## Observation model (binomial kernel via ones trick)
    nm[i] <- max(ceiling(tB[i]), 0)
    nM[i] <- min(floor(tD[i]) - nm[i], tMax - nm[i]) + 1
    pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## Year effect: RW1 on calendar year
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
    betaSEX[k] ~ dnorm(0, sd = 1)
    betaFEED[k] ~ dnorm(0, sd = 1)
    zSEX[k] ~ dbern(0.5)
    zFEED[k] ~ dbern(0.5)
  }
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
  sigma_year_a ~ dunif(0, 5)
  
})

censoredD <- rep(2, nind)
censoredD[cintD[,1] > 0] <- 1
n_year <- max(year_index)

## set up other components of model
consts <- list(nind = nind, tMax = ncol(CH), sex = sex, feeding = feeding, n_year = n_year, year_index = year_index)
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
      u_a = rep(0, n_year), 
      sigma_year_a = 0.5,
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
config <- configureMCMC(model, monitors = c("a", "b", "mean.p", "betaSEX", "betaFEED", 
                                            "zSEX", "zFEED", "u_a", "sigma_year_a"))
config$removeSamplers(c("a", "b"))
config$addSampler(target = c("a", "b"), type = 'AF_slice')
#config$addSampler(target = c("b"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
#config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))

## load in custom RJ-MCMC samplers
#source("MCMC_RJ_multi.R")

# Add reversible jump
configureRJ(conf = config,   ## model configuration
                  targetNodes = c("betaSEX", "betaFEED"),
                  indicatorNodes = c("zSEX", "zFEED"),
                  control = list(mean = 0, scale = 1))

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

#-----------------------------------------------------------------------------
# Run in parallel...

## --- PREP ON MASTER ---------------------------------------------------------
# Build constants and data
consts <- list(
  nind = nind,
  tMax = ncol(CH),
  sex = sex,
  feeding = feeding,
  n_year = n_year, 
  year_index = year_index
)

censoredD <- rep(2, consts$nind)
censoredD[cintD[,1] > 0] <- 1

dat <- list(
  y = y, cintB = cintB, cintD = cintD,
  censoredD = censoredD,
  tD = tD, tB = tB, tstar = tstar, dind = dind
)

monitors <- c("a","b","mean.p","betaSEX","betaFEED",
              "zSEX","zFEED","u_a","sigma_year_a")

rj_path <- normalizePath("MCMC_RJ_multi.R")
gompz_path1 <- normalizePath("Distributions/Dist_GompertzLB.R")
gompz_path2 <- normalizePath("Distributions/Dist_Gompertz.R")
gompz_pathN <- normalizePath("../NIMBLE_Distributions/Dist_GompertzNim.R")


## --- WORKER FUNCTION --------------------------------------------------------
workflow <- function(seed, dat, consts, monitors, rj_path, gompz_path1, gompz_path2, gompz_pathN) {
  library(nimble)
  source(gompz_path1)
  source(gompz_path2)
  source(gompz_pathN)
  source(rj_path)
  
  
  code <- nimbleCode({
    for (i in 1:nind) {
      ## Likelihood: continuous-time event model with Gompertz survival
      tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
      censoredD[i] ~ dinterval(tD[i], cintD[i, ])
      tD[i] <- tB[i] + tstar[i]
      tstar[i] ~ dgompzNim(amult[i], bmult[i])
      
      ## Linear predictors
      log(amult[i]) <- log(a) +
        betaSEX[1] * sex[i] * zSEX[1] +
        betaFEED[1] * feeding[i] * zFEED[1] +
        u_a_center[year_index[i]]
      
      log(bmult[i]) <- log(b) +
        betaSEX[2] * sex[i] * zSEX[2] +
        betaFEED[2] * feeding[i] * zFEED[2]
      
      ## Observation model (binomial kernel via ones trick)
      nm[i] <- max(ceiling(tB[i]), 0)
      nM[i] <- min(floor(tD[i]) - nm[i], tMax - nm[i]) + 1
      pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
      dind[i] ~ dbern(pd[i])
    }
    
    ## Year effect: RW1 on calendar year
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
      betaSEX[k] ~ dnorm(0, sd = 1)
      betaFEED[k] ~ dnorm(0, sd = 1)
      zSEX[k] ~ dbern(0.5)
      zFEED[k] ~ dbern(0.5)
    }
    a ~ dexp(1)
    b ~ dexp(1)
    mean.p ~ dunif(0, 1)
    sigma_year_a ~ dunif(0, 5)
    
  })
  
  # Build/compile model with the PASSED constants & data
  model <- nimbleModel(code, constants = consts, data = dat, check = TRUE)
  cModel <- compileNimble(model)
  
  # Inits generator using ONLY dat/consts in scope (no globals)
  initFn <- function(model, dat, consts) {
    repeat {
      inits <- list(
        a = rexp(1), b = rexp(1), mean.p = runif(1),
        tB = runif(consts$nind, dat$cintB[,1], dat$cintB[,2]),
        tD = ifelse(dat$censoredD == 2,
                    NA_real_,  # interval-censored: set below
                    runif(consts$nind, dat$cintD[,1], dat$cintD[,2]))
      )
      # For censoredD==2, ensure tD > tB by some margin
      id_cens <- which(dat$censoredD == 2)
      if (length(id_cens)) inits$tD[id_cens] <- inits$tB[id_cens] + dat$cintD[id_cens,2] + rexp(length(id_cens), 1)
      inits$tstar <- inits$tD - inits$tB
      
      inits$betaSEX <- rnorm(2)
      inits$betaFEED <- rnorm(2)
      inits$zSEX <- c(1,1)
      inits$zFEED <- c(1,1)
      inits$u_a <- rep(0, consts$n_year) 
      inits$sigma_year_a <- 0.5
      
      model$setInits(inits)
      if (is.finite(model$calculate())) return(inits)
    }
  }
  
  inits <- initFn(cModel, dat, consts)
  
  conf <- configureMCMC(model, monitors = monitors)
  conf$removeSamplers(c("a","b", "u_a"))
  conf$addSampler(target = c("a","b"), type = "AF_slice")
  conf$addSampler(target = c("u_a"), type = "AF_slice")
  
  # RJ setup
  configureRJ(
    conf = conf,
    targetNodes   = c("betaSEX", "betaFEED"),
    indicatorNodes= c("zSEX",   "zFEED"),
    control = list(mean = 0, scale = 1)
  )
  
  mcmc  <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  
  set.seed(seed)
  run <- runMCMC(cmcmc,
                 niter = 150000, nburnin = 19000, thin = 1,
                 nchains = 1, progressBar = FALSE,
                 samplesAsCodaMCMC = TRUE)
  run
}

## --- PARALLEL LAUNCH --------------------------------------------------------
showConnections(all = TRUE)
closeAllConnections()

nbcores <- max(1, detectCores() - 1)
cl <- makeCluster(nbcores)

seeds <- c(2022, 653)
out <- parLapply(cl, seeds, workflow,
                 dat = dat,
                 consts = consts,
                 monitors = monitors,
                 rj_path = rj_path,
                 gompz_path1 = gompz_path1,
                 gompz_path2 = gompz_path2,
                 gompz_pathN = gompz_pathN)

stopCluster(cl)

# Combine outputs
samples <- do.call(mcmc.list, lapply(out, `[[`, "samples"))
samples <- mcmc.list(lapply(out, as.mcmc))

MCMCsummary(samples)
plot(samples)
mcmcplot(samples)

saveRDS(samples, "DoubleCensored/All_Island_Data/Samples/G_AllIslands_SexFeeding2_RETime_a.rds")
out <- readRDS("DoubleCensored/All_Island_Data/Samples/G_AllIslands_SexFeeding2_RETime_a.rds")

