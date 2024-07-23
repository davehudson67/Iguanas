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

## load custom distributions
source("Dist_Gompertz.R")

## Load data
igs <- read_csv("Data/LeafCayThru2019_CleanOct2020.csv")

## Set date
igs$date <- mdy(igs$date)
igs$Year <- year(igs$date)
igs$age[igs$age == "2.7?"] <- 2.7

nind <- length(levels(as.factor(igs$animal_id)))
igs$tag <- as.numeric(as.factor(igs$animal_id))
igs <- arrange(igs, tag)
first_year <- 1970
last_year <- 2019

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
cintB[is.na(tKB), 1] <- tF[is.na(tKB)] - 20
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
    tD[i] ~ dGompertzLB(a, b, tB[i])
    
    ## sampling component
    nm[i] <- max(ceiling(tB[i]), 0)
    nM[i] <- min(floor(tD[i]) - nm[i], tMax - nm[i])
    pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a ~ dexp(1)
  b ~ dexp(1)
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
inits <- list(
  tD = tDinit,
  tB = tBinit,
  a = 0.1, 
  b = 0.1,
  mean.p = runif(1, 0, 1)
)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cModel <- compileNimble(model)

## try with adaptive slice sampler
config <- configureMCMC(cModel, monitors = c("a", "b", "mean.p", "tB", "tD"), thin = 1)
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
system.time(runAF <- runMCMC(cBuilt, 
                             niter = 50000, 
                             nburnin = 9000, 
                             nchains = 2, 
                             progressBar = TRUE, 
                             summary = TRUE, 
                             samplesAsCodaMCMC = TRUE, 
                             thin = 1))

runAF$summary
mcmcplot(runAF$samples, parms = c("a", "b", "mean.p"))

## extract birth and death times
samples <- runAF$samples
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


lifespan_post <- tDpost - tBpost
mcmcplot(lifespan_post)













#######
tBpost_df <- as.data.frame(tBpost[, 1:20])
colnames(tBpost_df) <- gsub("\\[|\\]", "", colnames(tBpost_df))
tBpost_long <- tBpost_df %>%
  gather(key = "ind", value = "value", 1:20)

# Add true birth times to the long format data frame for plotting
true_births <- rep(births_Kn[1:20], each = 320000/20)

# Merge true values with the posterior data
tBpost_long <- tBpost_long %>%
  mutate(Birth = true_births)

# Plot density of posterior samples with vertical lines for true values
ggplot(tBpost_long, aes(x = value)) +
  geom_density() +
  facet_wrap(~ ind, scales = "free") +
  geom_vline(aes(xintercept = Birth), color = "red") +
  labs(x = "Posterior Samples", y = "Density", title = "Posterior Distributions vs True Birth Times") +
  theme_minimal()

births_Kn[13]
cintB[15,]

