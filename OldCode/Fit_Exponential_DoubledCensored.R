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

set.seed(11)

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
    
    ## likelihood for interval-truncated Siler
    censoredB[i] ~ dinterval(tB[i], cintB[i, ])
    tB[i] ~ dflat()
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    tD[i] ~ dExponentialLB(r, tB[i])
    
    ## sampling component
    nm[i] <- max(ceiling(tB[i]), 0)
    nM[i] <- min(floor(tD[i]) - nm[i], tMax - nm[i])
    pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  r ~ dexp(1)
  mean.p ~ dunif(0, 1)
  
})

## set up other components of model
consts <- list(nind = nind, tMax = ncol(CH))
data <- list(y = y, cintB = cintB, cintD = cintD,
             censoredB = censoredB, censoredD = censoredD, 
             tD = tD, tB = tB, dind = dind)

## set initial values
tBinit <- apply(cbind(cintB, censoredB), 1, function(x) {
  if(x[3] == 0) {
    y <- x[1] - rexp(1, 0.1)
  } else {
    y <- runif(1, x[1], x[2])
  }
  y
})
tDinit <- apply(cbind(cintD, censoredD), 1, function(x) {
  if(x[3] == 2) {
    y <- x[2] + rexp(1, 0.1)
  } else {
    y <- runif(1, x[1], x[2])
  }
  y
})

initFn <- function(model) {
  valid <- 0
    while(valid == 0) {
  ## output initial values
  inits <- list(
    tD = tDinit,
    tB = tBinit,
    r = rexp(1), 
    mean.p = runif(1, 0, 1))
  
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
  inits[[k]] <- initFn(cModel)
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
nmix <- rbinom(1, size = nimp, prob = 0.95)
props <- sim(mod$modelName, mod$parameters, nmix)
props <- props[, -1, drop = FALSE]
colnames(props) <- c("r", "mean.p")

## take random samples from prior (to create defense mixture)
dmp <- runif((nimp - nmix), 0, 1)
defense <- matrix(rexp(nimp- nmix, 1))
defense <- cbind(defense, dmp)
colnames(defense) <- c("r", "mean.p")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:2)

## combine defense and importance samples
props <- rbind(props, defense)

##########################################
## calculate log likelihood

bL <- cintB[,1]
bU <- cintB[,2]
dU <- cintD[,2]

# Define the log-likelihood function
# Define the log-likelihood function
loglike <- function(censoredB, bL, bU, censoredD, tL, dU, y, tMax, r, p) {
  
  ## loop over individuals
  fsurv1 <- pmap_dbl(list(censoredB, bL, bU, censoredD, tL, dU, y, tMax), function(censoredB, bL, bU, censoredD, tL, dU, y, tMax, r, p) {
    
    ll <- y * log(p)
    
    # Combined Birth, Death, and Sampling components
    if (censoredB == 1) {
      t <- seq(1:(bU - bL))
      birth_ll <- log(1 / length(t))
      birth_ll <- log_sum_exp_marg(birth_ll, mn = FALSE)
    } else {
      birth_ll <- 0
    }
    
    if (censoredD == 1) { #interval censored
      t1 <- tL:1 + dU
      temp <- sapply(bL:bU, function(tB_val) {
        pExpo(t1 - tB_val, r) - pExpo(t1 - 1 - tB_val, r)
      })
      death_ll <- log(temp)
      death_ll <- log_sum_exp_marg(death_ll, mn = FALSE)
      sampling_ll <- sapply(bL:bU, function(tB_val) {
        (t1 - tB_val - y) * log(1 - p)
      })
      ll <- ll + birth_ll + death_ll + log_sum_exp_marg(sampling_ll, mn = FALSE)
    } else { #right censored
      if (tL < tMax) {
        t1 <- (tL + 1):tMax
        temp <- sapply(bL:bU, function(tB_val) {
          pExpo(t1 - tB_val, r) - pExpo(t1 - 1 - tB_val, r)
        })
        temp1 <- log(temp)
        temp_tail <- sapply(bL:bU, function(tB_val) {
          pExpo(tMax - tB_val, r, lower.tail = FALSE)
        })
        temp_tail <- log(temp_tail)
        death_ll <- c(temp1, temp_tail)
        death_ll <- log_sum_exp_marg(death_ll, mn = FALSE)
        sampling_ll <- sapply(bL:bU, function(tB_val) {
          (t1 - tB_val - y) * log(1 - p)
        })
        ll <- ll + birth_ll + death_ll + log_sum_exp_marg(sampling_ll, mn = FALSE)
      } else { #tL == tMax
        temp <- sapply(bL:bU, function(tB_val) {
          pExpo(tMax - tB_val, r, lower.tail = FALSE)
        })
        temp <- log(temp)
        death_ll <- temp
        death_ll <- log_sum_exp_marg(death_ll, mn = FALSE)
        sampling_ll <- sapply(bL:bU, function(tB_val) {
          (tL - tB_val - y) * log(1 - p)
        })
        ll <- ll + birth_ll + death_ll + log_sum_exp_marg(sampling_ll, mn = FALSE)
      }
    }
    ll
  }, r = r, p = p)
  
  sum(fsurv1)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, censoredB, bL, bU, censoredD, tL, dU, y, tMax) {
                           loglike(censoredB, bL, bU, censoredD, tL, dU, y, tMax, pars[1], pars[2])
                         }, censoredB = censoredB, bL = bL, bU = bU, censoredD = censoredD, tL = tL, dU = dU,
                         y = y, tMax = tMax, mc.cores = 24)
logimpweight <- reduce(logimpweight, base::c)

## add prior densities
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) +
                               dunif(props[, 2], 0, 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(props, mod$modelName, mod$parameters, FALSE) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) +
                                                                              dunif(props[, 2], 0, 1, log = TRUE)))

logimpweight1 <- logimpweight

saveRDS(logimpweight, "LImpWeights/logimpweight_e.rds")
logimpweight <- readRDS("LImpWeights/logimpweight_e.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## bootstrap the importance weights to create 95% intervals
imp.boot <- BootsPlot(logimpweight, 5000, TRUE)

logmarg_e <- log_sum_exp_marg(logimpweight)




