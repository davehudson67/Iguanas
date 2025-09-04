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

## Load data
igs <- read_csv("Data/LeafCayThru2019_CleanOct2020.csv")

## Set date
igs$date <- mdy(igs$date)
igs$Year <- year(igs$date)

igs$age[igs$age == "2.7?"] <- 2.7

## sort data
igsKA <- igs %>%
  filter(!is.na(age), .preserve = TRUE) %>%
  mutate(age = as.numeric(age)) %>%
  group_by(animal_id) %>%
  mutate(ageR = round(age)) %>%
  mutate(Birth = min(Year) - min(ageR)) %>%
  mutate(LastSeen = max(Year)) %>%
  distinct(animal_id, Year, .keep_all = TRUE) %>%
  mutate(captures = n()) %>%
  distinct(animal_id, .keep_all = TRUE) %>%
  mutate(death = NA) %>%
  ungroup() %>%
  mutate(Max_captures = max(Year))

## read in data
tKD <- igsKA$death
tB <- igsKA$Birth

## extract max possible capture time
tM <- igsKA$Max_captures

## extract last alive time
tL <- igsKA$LastSeen

## some checks
stopifnot(all(tL >= tB))
stopifnot(all(tKD[!is.na(tKD)] >= tL[!is.na(tKD)]))

## normalise to survival times
## (necessary at the moment due to censoring
## constraints)
#tKD <- tKD - tB
tL <- tL - tB
tM <- tM - tB
  
## define censoring matrices
cint <- cbind(tL, tKD)
cint[is.na(tKD), 2] <- cint[is.na(tKD), 1]
cint[is.na(tKD), 1] <- 0
colnames(cint) <- NULL
censored <- ifelse(!is.na(tKD), 1, 2)
igsKA$censored <- censored
tD <- rep(NA, length(tKD))
dind <- rep(1, length(tKD))

## extract number of captures
y <- igsKA$captures
names(y) <- NULL

## some checks
stopifnot(all(tM >= y))

## set up nind
nind <- length(y)

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated exponential
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dexp(r)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  r ~ dexp(1)
  mean.p ~ dunif(0,1)
})

## set up other components of model
consts <- list(nind = nind, tM = tM)
data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind)

## find overdispersed initial values
tinitFn <- function(cint, censored) {
  apply(cbind(cint, censored), 1, function(x) {
    if(x[3] == 2) {
      y <- x[2] + rexp(1, 1)
    } else {
      y <- runif(1, x[1], x[2])
    }
    y
  })
}
initFn <- function(cint, censored) {
  ## get ML estimates as initial values
  optFn <- function(pars, t) {
    if(pars < 0) {
      return(NA)
    }
    sum(dExpo(t, r = pars, log = TRUE))
  }
  ## sample missing values
  tD <- tinitFn(cint, censored)
  pars <- optimise(optFn, c(0, 100), t = tD, maximum = TRUE)
  pars <- pars$maximum
  list(
    tD = tD,
    r = pars,
    mean.p = runif(1, 0, 1)
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored))

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = c("r", "mean.p", "tD"), thin = 1)

## check monitors and samplers
config$printMonitors()
config$printSamplers("r")

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run_e <- runMCMC(cbuilt, 
                             niter = 500000, 
                             nburnin = 19000, 
                             nchains = 2, 
                             progressBar = TRUE, 
                             summary = TRUE, 
                             samplesAsCodaMCMC = TRUE, 
                             thin = 1))
run_e$summary

mcmcplot(run_e$samples, parms = "tD")

## plot mcmcm
plot(run_e$samples)
samples <- run_e$samples




## load data
load("Data/badgerSexInb_AdultInfvCubInfvUninf.RData")


CH$death[is.na(CH$death)] <- CH$last_seen[is.na(CH$death)] + rnorm(1, 2, 1)
CH$dur <- CH$death - CH$birth
CH$censored[CH$censored == 2] <- 0

badgers <- CH %>%
  filter(censored > 0)


badgers <- CH

sex <- badgers$sex
#inbreeding <- ifelse(badgers$hom > median(badgers$hom), 1, 0)
infection <- rep(0, times = nrow(badgers))
infection[badgers$infected_as_cub == 1] <- 1
infection[badgers$infected_lifetime > badgers$infected_as_cub] <- 2


inbr_level0 <- quantile(badgers$hom, probs = 0.333)
inbr_level1 <- quantile(badgers$hom, probs = 0.666)
inbr <- rep(0, times = nrow(badgers))

inbr[badgers$hom >= inbr_level0] <- 1
inbr[badgers$hom >= inbr_level1] <- 2


badgersS <- Surv(badgers$dur, badgers$censored)


fit <- survfit(badgersS ~ inbr + sex)

ggsurvplot(fit, data = badgersS, risk.table = TRUE,
           conf.int = FALSE, pval = TRUE, title = "Kaplan-Meier plot of survival probability",
           legend = "right")
