# Load necessary library
library(nimble)
library(tidyverse)
library(lubridate)
source("Rotger_et_al_PaperCode/functions.R")

# Load data
igs <- read_csv("Data/LeafCayThru2019_CleanOct2020.csv")
igs$age <- as.numeric(igs$age)

# Randomly select 50 individuals
set.seed(42)  # for reproducibility
#sampled_ids <- sample(unique(igs$animal_id), 5)
#igs <- igs %>%
#  filter(animal_id %in% sampled_ids)

# Select individuals with svl
igsKA <- igs %>%
  filter(!is.na(svl)) %>%
  mutate(animal_id = as.numeric(as.factor(animal_id))) %>%
  mutate(date = mdy(date)) %>%
  mutate(year = year(date)) %>%
  group_by(animal_id) %>%
  distinct(year, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(Year = year - 1979) %>%
  group_by(animal_id) %>%
  mutate(last = max(Year)) %>%
  arrange(animal_id) %>%
  ungroup()

nind <- length(unique(igsKA$animal_id))
nT_ints <- 40
deltaT <- rep(1, nT_ints - 1)

# Data matrices
svl_matrix <- matrix(NA, nrow = nind, ncol = nT_ints)
y <- matrix(0, nrow = nind, ncol = nT_ints)
z <- matrix(NA, nrow = nind, ncol = nT_ints)
#age_matrix <- matrix(NA, nrow = N, ncol = nT_ints)

# Fill the matrix with SVL measurements
for (i in 1:nrow(igsKA)) {
  id <- igsKA$animal_id[i]
  year <- igsKA$Year[i]
  svl <- igsKA$svl[i]
  svl_matrix[id, year] <- svl
  y[id, year] <- 1
}

# last capture occasion
last <- igsKA %>%
  distinct(animal_id, .keep_all = TRUE) %>%
  pull(last)

# Determine the first capture times
f <- apply(!is.na(svl_matrix), 1, which.max)

# Add z known data
for (i in 1:nind){
  z[i, ((f[i]):last[i])] <- 1
  z[i, f[i]] <- NA
}

code <- nimbleCode({    
  #####################
  ####Growth model####
  ##################### 
  for(i in 1:nind){
    #SVL at hatching size
    Hs[i] <- mu_Hs + eps_Hs[i]
    eps_Hs[i] ~ dnorm(0, sd = sd_Hs)
    #maximum size
    y2[i] <- mu_y2 + maleEff_y2 * sex[i] + eps_y2[i]
    eps_y2[i] ~ dnorm(0, sd = sd_y2)
    #characteristic grwth rate
    log(k[i]) <- mu_k + maleEff_k * sex[i] + eps_k[i]
    eps_k[i] ~ dnorm(0, sd = sd_k)
    
    # first capture
    Lr[i, f[i]] ~ dnorm(Lr_hat[i, f[i]], sd = sd_Lr)         # SVL - observed
    Lr_hat[i, f[i]] ~ dunif(4, 50)   # expected size of animal at unknown age at first capture
    
    #sex[i] ~ dbern(psi)
  }
  
  #Priors Hsize
  mu_Hs ~ dunif(5, 20)
  sd_Hs ~ dunif(0, 10) 
  
  #Priors for y2
  mu_y2 ~ dunif(20, 70)
  maleEff_y2 ~ dnorm(0, 1)
  sd_y2 ~ dunif(0, 5)
  
  #Priors for K
  mu_k ~ dunif(-10, 10)
  maleEff_k ~ dnorm(0, 1)
  sd_k ~ dunif(0, 2)
  
  #psi ~ dunif(0, 1) #p of being male
  sd_Lr ~ dunif(0, 5)
  
  for(i in 1:nind){
    for(t in (f[i]+1):Ti){
      #Schnute growth equation
      Lr[i, t] ~ dnorm(Lr_hat[i, t], sd = sd_Lr)   
      Lr_hat[i, t] <- Lr_hat[i, t - 1] * exp(-k[i] * deltaT[t - 1]) + (y2[i] - Hs[i] * exp(-k[i] * (55 - 0))) * (1 - exp(-k[i] * deltaT[t - 1])) / (1 - exp(-k[i] * (55 - 0)))
    }
    Linf[i] <- (y2[i] - Hs[i] * exp(-k[i] * (55 - 0))) / 1 - exp(-k[i] * (55 - 0))
    t0[i] <- 1 / k[i] * log((y2[i] - Hs[i]) / (y2[i] - Hs[i] * exp(-k[i] * (55 - 0))))
    
    # standardize SVL
    for(t in (f[i]):Ti){
      SVL_st[i, t] <- (Lr_hat[i, t] - mean.svl) / sd.svl
    }
  }
  
  #####################
  ####CJS model####
  ##################### 
  
  #------------------------------------------------------------------#
  #                            Survival                     #
  #------------------------------------------------------------------#
  
  #priors
  mu_phi <- log(mean.phi / (1-mean.phi))
  mean.phi ~ dunif(0, 1)
  
  betaS.sex ~ dnorm(0, sd = 1)
  betaS.svl ~ dnorm(0, sd = 1)
  
  #Survival
  for (i in 1:nind){
    for (t in f[i]:(Ti - 1)){
      logit(phi[i, t]) <- mu_phi +  betaS.sex * sex[i] + betaS.svl * SVL_st[i, t] + eps_phi[t] 
    }
  }
  
  #------------------------------------------------------------------#
  #                            Recapture                             #
  #------------------------------------------------------------------#
  
  #priors
  mu_p <- log(mean.p / (1 - mean.p))
  mean.p ~ dunif(0, 1)
  
  betaP.sex ~ dnorm(0, sd = 1)
  betaP.svl ~ dnorm(0, sd = 1)
  
  #Detectability
  for (i in 1:nind){
    for (t in f[i]:(Ti - 1)){
      logit(p[i, t]) <-  mu_p + betaP.sex * sex[i] + betaP.svl * SVL_st[i, t] + eps_p[t]
    }
  }
  
  for (t in 1:(Ti - 1)){
    eps_phi[t] ~ dnorm(0, sd = sigma_phi)
    eps_p[t] ~ dnorm(0, sd = sigma_p)
  }
  
  sigma_phi ~ dunif(0, 5)
  sigma_p ~ dunif(0, 5)
  
  #likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i, f[i]] <- 1
    for (t in (f[i] + 1):Ti){
      # State process
      z[i, t] ~ dbern(mu1[i, t])
      mu1[i, t] <- phi[i, t - 1] * z[i, t - 1]
      # Observation process
      y[i, t] ~ dbern(mu2[i, t])
      mu2[i, t] <- p[i, t - 1] * z[i, t]
    } 
  }
})

# remove f = Ti individuals
rem <- which(f == nT_ints)

sex <- as.numeric(as.factor(igsKA$sex)) - 1

# Create lists for NIMBLE
data <- list(y = y[-rem,], 
             z = z[-rem, ],
             sex = sex[-rem],
             #meshpoints = meshpoints,
             Lr = svl_matrix[-rem,])

#data <- list(y = y, 
#             z = z,
#             sex = sex,
#             #meshpoints = meshpoints,
#             Lr = svl_matrix)

mean_svl <- mean(igs$svl, na.rm = TRUE)
sd_svl <- sd(igs$svl, na.rm = TRUE)

consts <- list(f = f[-rem], 
               nind = nind - length(rem),
               mean.svl = mean_svl,
               sd.svl = sd_svl,
               deltaT = deltaT,
               Ti = nT_ints)

#consts <- list(f = f, 
#               nind = nind,
#               mean.svl = mean_svl,
#               sd.svl = sd_svl,
#               deltaT = deltaT,
#               Ti = nT_ints)

#Initial values
zinit <- cjs.init.z(data$y, consts$f)

inits <- list(mu_Hs = 6, sd_Hs = 3,
              mu_y2 = runif(1, 20, 60),
              maleEff_y2 = rnorm(1, 10), sd_y2 = runif(1, 1, 2),
              mu_k = runif(1, -7, -6.5), maleEff_k = rnorm(1, 10), sd_k = runif(1, 0.5, 0.6),
              sd_Lr = runif(1, 0.1, 0.5),
              z = zinit,
              betaS.sex = runif(1, -5, 5), betaS.svl = runif(1, -5, 5),
              betaP.sex = runif(1, -5, 5), betaP.svl = runif(1, -5, 5),
              sigma_phi = runif(1, 0.25, 0.30), sigma_p = runif(1, 0.4, 0.5))

#Parameters
params<- c('mu_Hs','sd_Hs','mu_y2','sd_y2','maleEff_y2',
           'mu_k','maleEff_k','sd_k','sd_Lr',
           'mean.p','mean.phi',
           'mu_phi','betaS.sex','betaS.svl',
           'mu_p','betaP.sex','betaP.svl',
           'sigma_phi','sigma_p')

params_growth<- c('Hs', 'Lr_hat', 'y2', 'k','sex') #Individual Growth curve parameters

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data)#, inits = inits)
model$initializeInfo()

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = params, monitors2 = params_growth, thin = 10, thin2 = 10)

## check monitors and samplers
config$printMonitors()

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run <- runMCMC(cbuilt, 
                             niter = 50000, 
                             nburnin = 19000, 
                             nchains = 2, 
                             progressBar = TRUE, 
                             summary = TRUE, 
                             samplesAsCodaMCMC = TRUE, 
                             thin = 1))

saveRDS(run, "ig_modelrun_out.rds")
out <- readRDS("ig_modelrun_out.rds")

library(MCMCvis)
MCMCsummary(out$samples)
run$summary

plot(out$samples)
run_e$samples




