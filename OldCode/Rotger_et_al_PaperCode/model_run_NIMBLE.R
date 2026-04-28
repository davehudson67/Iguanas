library(nimble)
library(coda)
library(tidyverse)

#load data
load("data.Rdata")
source("functions.R")

code <- nimbleCode({    
  #####################
  ####Growth model####
  ##################### 
  for(i in 1:nind){
    #SVL at hatching size
    Hs[i] <- mu_Hs + eps_Hs[i]
    eps_Hs[i] ~ dnorm(0, tau_Hs)
    #maximum size
    y2[i] <- mu_y2 + maleEff_y2 * sex[i] + eps_y2[i]
    eps_y2[i] ~ dnorm(0, tau_y2)
    #characteristic grwth rate
    log(k[i]) <- mu_k + maleEff_k * sex[i] + eps_k[i]
    eps_k[i] ~ dnorm(0, tau_k)
    
    # first capture
    Lr[i,f[i]] ~ dnorm(Lr_hat[i,f[i]], tau_Lr)         # SVL - observed
    Lr_hat[i,f[i]] ~ dunif(30,80)   # expected size of animal at unknown age at first capture
    
    sex[i] ~ dbern(psi)
  }
  
  #Priors Hsize
  mu_Hs ~ dunif(25, 35)
  tau_Hs <- 1 / sd_Hs^2
  sd_Hs ~ dunif(2, 3) 
  
  #Priors for y2
  mu_y2 ~ dunif(65, 70)
  maleEff_y2 ~ dunif(0, 15)
  tau_y2 <- 1 / sd_y2^2
  sd_y2 ~ dunif(0, 5)
  
  #Priors for K
  mu_k ~ dunif(-8,-5)
  maleEff_k ~ dunif(0,1.5)
  tau_k <- 1 / sd_k^2
  sd_k ~ dunif(0,2)
  
  psi ~ dunif(0, 1)
  tau_Lr <- 1 / sd_Lr^2
  sd_Lr ~ dunif(0, 5)
  
  for(i in 1:nind){
    for(t in (f[i]+1):Ti){
      #Schnute growth equation
      Lr[i, t] ~ dnorm(Lr_hat[i, t], tau_Lr)   
      Lr_hat[i, t] <- Lr_hat[i, t - 1] * exp(-k[i] * deltaT[t - 1]) + (y2[i] - Hs[i] * exp(-k[i] * (5110 - 0))) * (1 - exp(-k[i] * deltaT[t - 1])) / (1 - exp(-k[i] * (5110 - 0)))
    }
    Linf[i] <- (y2[i] - Hs[i] * exp(-k[i] * (5110 - 0))) / 1 - exp(-k[i] * (5110 - 0))
    t0[i] <- 1 / k[i] * log((y2[i] - Hs[i]) / (y2[i] - Hs[i] * exp(-k[i] * (5110-0))))
    
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
  
  betaS.sex~dunif(-5, 5)
  betaS.svl~dunif(-5, 5)
  
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
  
  betaP.sex~dunif(-5, 5)
  betaP.svl~dunif(-5, 5)
  
  #Detectability
  for (i in 1:nind){
    for (t in f[i]:(Ti - 1)){
      logit(p[i, t]) <-  mu_p + betaP.sex * sex[i] + betaP.svl * SVL_st[i, t] + eps_p[t]
    }
  }
  
  for (t in 1:(Ti - 1)){
    eps_phi[t] ~ dnorm(0, tau_phi)
    eps_p[t] ~ dnorm(0, tau_p)
  }
  
  tau_phi <- 1 / sigma_phi^2
  sigma_phi ~ dunif(0, 5)
  
  tau_p <- 1 / sigma_p^2
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

# Unlist JAGS data
list2env(jags.data, envir = .GlobalEnv)
Ti = 10

# remove f = Ti individuals
rem <- which(f == Ti)

# Create lists for NIMBLE
data <- list(y = y[-rem,], 
             z = z[-rem, ],
             sex = sex[-rem],
             #meshpoints = meshpoints,
             Lr = Lr[-rem,])

consts <- list(f = f[-rem], 
               nind = nind - length(rem),
               mean.svl = mean.svl,
               sd.svl = sd.svl,
               deltaT = deltaT,
               Ti = 10
               )

#Initial values
zinit <- cjs.init.z(data$y, consts$f)

inits <- list(mu_Hs=29.3,sd_Hs=2.57,
              mu_y2=runif(1,65,68),
              maleEff_y2=runif(1,7,8),sd_y2=runif(1,1,2),
              mu_k=runif(1,-7,-6.5),maleEff_k=runif(1,0.45,0.50),sd_k=runif(1,0.5,0.6),
              sd_Lr=runif(1,0.1,0.5),psi=runif(1,0.45,0.55),
              z=zinit,
              betaS.sex=runif(1,-5,5),betaS.svl=runif(1,-5,5),
              betaP.sex=runif(1,-5,5),betaP.svl=runif(1,-5,5),
              sigma_phi=runif(1,0.25,0.30),sigma_p=runif(1,0.4,0.5))

#Parameters
params<- c('mu_Hs','sd_Hs','mu_y2','sd_y2','maleEff_y2',
           'mu_k','maleEff_k','sd_k','sd_Lr','psi',
           'mean.p','mean.phi',
           'mu_phi','betaS.sex','betaS.svl',
           'mu_p','betaP.sex','betaP.svl',
           'sigma_phi','sigma_p')

params_growth<- c('Hs', 'Lr_hat', 'y2', 'k','sex') #Individual Growth curve parameters

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = inits)

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = params, monitors2 = params_growth, thin = 1, thin2 = 10)

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


plot(run$samples)

library(MCMCvis)
MCMCsummary(run$samples)



#####################################################################################################################################
####################################################################################################################################
set.seed(123)


n.adapt <- 20000
n.burnin <- 50000
n.iter <- 150000
n.thin <- 10  
n.chains <- 5

start.time = Sys.time()
out <- jags(data = jags.data,
            inits=inits,
            params,
            "Best_model_code.txt",
            n.chains=n.chains,
            n.adapt=n.adapt,
            n.thin=n.thin,
            n.iter=n.iter,
            n.burnin=n.burnin)
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'),2)
cat(paste(paste('Posterior computed in ', elapsed.time, sep=' '),'minutes\n', sep=' '))


library(MCMCvis)
MCMCsummary(out)
MCMCplot(out, params=c("betaS.sex","betaS.svl","betaP.sex","betaP.svl"))
MCMCtrace(out)

save(out, file = 'Plilfordi_MCMCsamples.RData')

############################################################################################################################################