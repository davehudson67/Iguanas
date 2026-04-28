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

## load custom distributions
source("Distributions/Dist_ExpoLB.R")
source("Distributions/Dist_Expo.R")

## Load data
igs <- read_csv("Data/LeafCayThru2019_CleanOct2020.csv")

## Set date
igs$date <- mdy(igs$date)
igs$Year <- year(igs$date)
igs$age[igs$age == "2.7?"] <- 2.7

nind <- length(levels(as.factor(igs$animal_id)))
igs$tag <- as.numeric(as.factor(igs$animal_id))
igs <- arrange(igs, tag)
first_year <- 1950
last_year <- 2019
tMax <- last_year - first_year

## sort data
igsAll <- igs %>%
  #filter(!is.na(age), .preserve = TRUE) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(Year = Year - first_year) %>%
  group_by(tag) %>%
  mutate(ageR = round(age)) %>%
  mutate(Birth = min(Year) - min(ageR)) %>%
  mutate(LastSeen = max(Year)) %>%
  distinct(tag, Year, .keep_all = TRUE) %>%
  #mutate(captures = n()) %>%
  #distinct(animal_id, .keep_all = TRUE) %>%
  #mutate(death = NA) %>%
  ungroup() %>%
  #mutate(Max_captures = max(Year)) %>%
  mutate(n = 1)

CH <- matrix(NA, nind, last_year - first_year)

## Fill CH matrix using a loop
for (i in 1:nrow(igsAll)) {
  tag_index <- igsAll$tag[i]
  year_index <- igsAll$Year[i]
  
  # Fill the CH matrix
  CH[tag_index, year_index] <- igsAll$n[i]
}

## Add births
filtered <- filter(igsAll, !is.na(igsAll$age), .preserve = TRUE)
for(i in 1:nrow(filtered)){
  tag <- igsAll$tag[i]
  birth <- igsAll$Birth[i]
  CH[tag, birth] <- 1
}
rm(filtered)

tKD <- rep(NA, nrow(CH))
tKB <- igsAll %>%
  distinct(animal_id, .keep_all = TRUE) %>%
  select(Birth)
tKB <- tKB$Birth

CH[is.na(CH)] <- 0

## extract last alive time
tL <- apply(CH, 1, function(x) max(which(x == 1)))
names(tL) <- NULL

## extract first alive time
tF <- apply(CH, 1, function(x) min(which(x == 1)))
names(tF) <- NULL

## define censoring times for birth time
cintB <- cbind(tF - 1, tF)
cintB[is.na(tKB), 1] <- tF[is.na(tKB)] - 25
cintB[cintB < 0] <- 1
#cintB[is.na(tKB), 2] <- 0
colnames(cintB) <- NULL
censoredB <- rep(1, nrow(cintB))
tB <- rep(NA, length(tKB))

## define censoring matrices for death time
cintD <- cbind(tL, tKD)
cintD[is.na(tKD), 2] <- cintD[is.na(tKD), 1]
cintD[is.na(tKD), 1] <- 0
colnames(cintD) <- NULL
censoredD <- ifelse(!is.na(tKD), 1, 2)
tD <- rep(NA, length(tKD))
tstar <- rep(NA, length(tKD))

## some checks
stopifnot(all(tL >= tF))
stopifnot(all(tKD[!is.na(tKD)] > tL[!is.na(tKD)]))

## set up dummy variables
dind <- rep(1, length(tKD))

## extract number of captures
y <- apply(CH, 1, sum)
names(y) <- NULL

## set up nind
nind <- length(y)

## code for NIMBLE model with censoring
code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood 
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    censoredD[i] ~ dinterval(tD[i], cintD[i])
    tD[i] <- tB[i] + tstar[i]
    tstar[i] ~ dexp(r)
    
    ## sampling component
    nm[i] <- max(ceiling(tB[i]), 0)
    nM[i] <- min(floor(tD[i]) - nm[i], tMax - nm[i]) + 1
    pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  r ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

## set up other components of model
consts <- list(nind = nind, tMax = ncol(CH))
data <- list(y = y, cintB = cintB, cintD = cintD[, 2],
             censoredD = rep(1, nind), 
             tD = tD, tB = tB, tstar = tstar, dind = dind)

## set initial values
initFn <- function(model, cintB, cintD) {
  valid <- 0
  while(valid == 0) {
    r <- rexp(1)
    mean.p <- runif(1, 0, 1)
    tBinit <- runif(nrow(cintB), cintB[, 1], cintB[, 2])
    tstarinit <- cintD[, 2] + rexp(nrow(cintD), r)
    tDinit <- tB + tstarinit
    inits <- list(
      tD = tDinit,
      tstar = tstarinit,
      tB = tBinit,
      r = r,
      mean.p = mean.p
    )
    model$setInits(inits)
    valid <- ifelse(!is.finite(model$calculate()), 0, 1)
  }
  return(inits)
}


## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data)

## compile the model
cModel <- compileNimble(model)

## find list of valid initial values using compiled model
inits <- list()
for(k in 1:2) {
  inits[[k]] <- initFn(cModel, cintB, cintD)
}

## try with default sampler
config <- configureMCMC(cModel, monitors = c("r", "mean.p", "tB", "tD"), thin = 1)
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

saveRDS(run, "Samples/E_dc_LeafCay.rds")
run <- readRDS("Samples/E_dc_LeafCay.rds")
run$summary

mcmcplot(run$samples, parms = c("r", "mean.p"))

## extract birth and death times
samples <- run$samples

tBpost <- as.matrix(samples) %>%
  as_tibble() %>%
  select(starts_with("tB"))
tDpost <- as.matrix(samples) %>%
  as_tibble() %>%
  select(starts_with("tD"))

## examine first few just to check
tBpost[, 1:20] %>%
  gather(ind, value) %>%
  ggplot(aes(x = value)) +
  geom_density() +
  facet_wrap(~ ind)

lifespan_post <- tDpost[,1:20] - tBpost[,1:20]
mcmcplot(lifespan_post)

################################################################################
#                                                                              #
#        Importance Sampling                                                   #
#                                                                              #
################################################################################
samples <- as.matrix(samples) %>%
  as_tibble() %>%
  select(c("r", "mean.p"))

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
nimp <- 10000
mixind <- rbinom(nimp, size = 1, prob = 0.95)
props <- matrix(NA, nimp, 2)
props[mixind == 1, ] <- sim(mod$modelName, mod$parameters, sum(mixind))[, -1]
colnames(props) <- c("r", "mean.p")

## take random samples from prior (to create defense mixture)
defense <- runif(sum(mixind == 0), 0, 1)
defense <- cbind(rexp(sum(mixind == 0), 1), defense)
colnames(defense) <- c("r", "mean.p")
props[mixind == 0, ] <- defense

## check IS distribution against posterior samples
as.data.frame(props[mixind == 1, ]) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:2)

##########################################
## calculate log likelihood

# Define the log-likelihood function
loglike <- function(cintB, cintD, y, tMax, r, p) {
  
  ## sample birth times across all individuals
  tB <- runif(nrow(cintB), cintB[, 1], cintB[, 2])
  
  ## sample death times across all individuals
  tstar <- map2_dbl(tB, cintD[, 2], function(tB, tM, r) {
    rExponentialLB(1, r, tM - tB)
  }, r = r)
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
                           loglike(cintB, cintD, y, tMax, pars[1], pars[2])
                         }, cintB = cintB, cintD = cintD,
                         y = y, tMax = tMax, mc.cores = 8)
logimpweight <- reduce(logimpweight, base::c)

## add prior densities
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) +
                               dunif(props[, 2], 0, 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(props, mod$modelName, mod$parameters, FALSE) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) +
                                                                              dunif(props[, 2], 0, 1, log = TRUE)))

saveRDS(logimpweight, "LImpWeights/logimpweight_e1.rds")
logimpweight <- readRDS("LImpWeights/logimpweight_e1.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## bootstrap the importance weights to create 95% intervals
imp.boot <- BootsPlot(logimpweight, 5000, TRUE)

logmarg_e <- log_sum_exp_marg(logimpweight)




