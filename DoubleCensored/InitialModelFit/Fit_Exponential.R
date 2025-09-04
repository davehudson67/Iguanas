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

## load custom distributions
source("Distributions/Dist_ExpoLB.R")
source("Distributions/Dist_Expo.R")

## Load data
#igs <- read_csv("Data/LeafCayThru2019_CleanOct2020.csv")
load("DoubleCensored/igs_AllIslands_CleanDH_040925_obsStart.RData")

## Set date
#igs$date <- mdy(igs$date)
#igs$Year <- year(igs$date)
#igs$age[igs$age == "2.7?"] <- 2.7 # Adjust presumed typo

#nind <- length(levels(as.factor(igs$animal_id)))
#igs$tag <- as.numeric(as.factor(igs$animal_id))
#igs <- arrange(igs, tag)

#min(igs$Year) # First record is 1980, iguana lifespan approx 30yrs so presume earliest iguana in the dataset has earliest birth of 1950
#max(igs$Year)
#first_year <- 1950
#last_year <- 2019 # Final year of recording on db
#tMax <- last_year - first_year # Potential yearly captures

## sort data
#igsAll <- igs %>%
#  #filter(!is.na(age), .preserve = TRUE) %>%
#  mutate(age = as.numeric(age)) %>%
#  mutate(Year = Year - first_year) %>%
#  group_by(tag) %>%
#  mutate(ageR = round(age)) %>%
#  mutate(Birth = min(Year) - min(ageR)) %>%
#  mutate(LastSeen = max(Year)) %>%
#  distinct(tag, Year, .keep_all = TRUE) %>%
#  #mutate(captures = n()) %>%
#  #distinct(animal_id, .keep_all = TRUE) %>%
#  #mutate(death = NA) %>%
#  ungroup() %>%
#  #mutate(Max_captures = max(Year)) %>%
#  mutate(n = 1)

#CH <- matrix(NA, nind, last_year - first_year)

## Fill CH matrix using a loop
#for (i in 1:nrow(igsAll)) {
#  tag_index <- igsAll$tag[i]
#  year_index <- igsAll$Year[i]
#  
  # Fill the CH matrix
#  CH[tag_index, year_index] <- igsAll$n[i]
#}

## Add births
#filtered <- filter(igsAll, !is.na(igsAll$age), .preserve = TRUE)
#for(i in 1:nrow(filtered)){
#  tag <- igsAll$tag[i]
#  birth <- igsAll$Birth[i]
#  CH[tag, birth] <- 1
#}
#rm(filtered)

#tKD <- rep(NA, nrow(CH))
#tKB <- igsAll %>%
#  distinct(animal_id, .keep_all = TRUE) %>%
#  select(Birth)
#tKB <- tKB$Birth

#CH[is.na(CH)] <- 0

## extract last alive time
#tL <- apply(CH, 1, function(x) max(which(x == 1)))
#names(tL) <- NULL

## extract first alive time
#tF <- apply(CH, 1, function(x) min(which(x == 1)))
#names(tF) <- NULL

## define censoring times for birth time
#cintB <- cbind(tF - 1, tF)
#cintB[is.na(tKB), 1] <- tF[is.na(tKB)] - 25
#cintB[cintB < 0] <- 1
##cintB[is.na(tKB), 2] <- 0
#colnames(cintB) <- NULL
#censoredB <- rep(1, nrow(cintB))
#tB <- rep(NA, length(tKB))

## define censoring matrices for death time
#cintD <- cbind(tL, tKD)
#cintD[is.na(tKD), 2] <- cintD[is.na(tKD), 1]
#cintD[is.na(tKD), 1] <- 0
#colnames(cintD) <- NULL
#censoredD <- ifelse(!is.na(tKD), 1, 2)
#tD <- rep(NA, length(tKD))
#tstar <- rep(NA, length(tKD))

## some checks
stopifnot(all(tL >= tF))
stopifnot(all(tKD[!is.na(tKD)] > tL[!is.na(tKD)]))

## set up dummy variables
#dind <- rep(1, length(tKD))

## extract number of captures
#y <- apply(CH, 1, sum)
#names(y) <- NULL

## set up nind
#nind <- length(y)

#save.image("LeafClayDataImage_only.RData")

## code for NIMBLE model with censoring
code <- nimbleCode({
  
  ## survival components for dead badgers
  for (i in 1:nind) {
    
    ## likelihood 
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])

    ## Left truncation: must survive to reach time 0 if born before 0.
    L[i] <- max(0, -tB[i])

    ## Time-to-death since birth: Exponential, conditional on surviving to L[i].
    tstar[i] ~ dExponentialLB(r, lowerBound = L[i]) 

    tD[i] <- tB[i] + tstar[i]
    
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
   
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
consts <- list(nind = nind, tMax = tMax)
data <- list(y = y, cintB = cintB, cintD = cintD,
             censoredD = censoredD, 
             tD = tD, tB = tB, tstar = tstar, dind = dind)

## set initial values
initFn <- function(model, cintB, cintD, censoredD) {
  browser()
  repeat {
    r      <- rexp(1)
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
      r       = r,
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
for (k in 1:2) {
  inits[[k]] <- initFn(cModel, cintB, cintD, censoredD)
}

## try with default sampler
config <- configureMCMC(cModel, monitors = c("r", "mean.p"), thin = 1)
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

saveRDS(run, "DoubleCensored/InitialModelFit/Samples/E_AllIslands.rds")
run <- readRDS("Samples/E_dc_LeafCay.rds")
run$summary

plot(run$samples, parms = c("r", "mean.p"))

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
  #browser()
  # Sample birth times
  tB <- runif(nrow(cintB), cintB[, 1], cintB[, 2])
  L <- pmax(0, -tB)
  
  # Sample death times
  tstar <- purrr::map2_dbl(tB, cintD[, 2], function(tB_i, tM_i, r) {
    rExponentialLB(1, r, tM_i - tB_i)
  }, r = r)
  
  tD <- tB + tstar
  
  # Survival component with left truncation
  log_surv <- log(r) - r * tstar - log(exp(-r * L))
  
  # Censoring component: interval or right-censored
  log_censor <- ifelse(
    is.finite(cintD[, 2]),
    log(pmax(pexp(cintD[, 2] - tB, rate = r) - pexp(cintD[, 1] - tB, rate = r), .Machine$double.eps)),
    log(pmax(1 - pexp(cintD[, 1] - tB, rate = r), .Machine$double.eps))
  )
  
  # Sampling component
  nm <- pmax(ceiling(tB), 0)
  nM <- pmin(floor(tD) - nm, tMax - nm) + 1
  stopifnot(all(nM >= y))
  
  lpd <- y * log(p) + (nM - y) * log(1 - p)
  
  # Total log-likelihood
  sum(log_surv + log_censor + lpd)
}


## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, cintB, cintD, y, tMax) {
                           loglike(cintB, cintD, y, tMax, pars[1], pars[2])
                         }, cintB = cintB, cintD = cintD,
                         y = y, tMax = tMax)
logimpweight <- reduce(logimpweight, base::c)

## add prior densities
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) +
                               dunif(props[, 2], 0, 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(props, mod$modelName, mod$parameters, FALSE) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) +
                                                                              dunif(props[, 2], 0, 1, log = TRUE)))

saveRDS(logimpweight, "DoubleCensored/InitialModelFit/LImpWeights/logimpweight_e.rds")
logimpweight <- readRDS("LImpWeights/logimpweight_e1.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## bootstrap the importance weights to create 95% intervals
imp.boot <- BootsPlot(logimpweight, 5000, TRUE)

logmarg_e <- log_sum_exp_marg(logimpweight)




