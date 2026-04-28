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

setwd("~/Iguanas")
source("ModelComparison_FUNCTIONS.R")
#source("DoubleCensored/Prep_All_Data.R")
load("igs_AllIslands_CleanDH_061025_obsStart.RData")
rm(inorn)
rm(sep_ids)

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

species <- AllData %>%
  mutate(species = as.numeric(as.factor(Species))) %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(species) - 1

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

summary(as.factor(AllData$feeding))
#feeding <- as.factor(AllData$feeding)
#feeding <- model.matrix(~ feeding, feeding)[, -1]
#colSums((feeding))
n_year <- ncol(CH)
t0 <- 0L

# Deterministic year index for each ID from entry-interval midpoint
tB_mid <- rowMeans(cintB)                          # midpoint of cint interval
year_index <- floor(tB_mid - t0) + 1L              # map to 1..n_year bins
year_index <- pmin(pmax(year_index, 1L), n_year)   # clamp
b_year <- AllData %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(Birth)

#-------------------------------------------------------------------------------

code <- nimbleCode({
  for (i in 1:nind) {
    # Latent Birth Year
    # tB is the continuous time of birth within that year.
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    b_year[i] <- floor(tB[i])
    
    # Left Truncation
    # The animal must survive from birth until surveys started on its island.
    # IslandStart is now a continuous time point on the global scale (e.g., 24.0).
    L[i] <- max(0, IslandStart[i] - tB[i])
    
    # Survival Process
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    # Linear Predictors (Cohort Model)
    # Use the integer b_year[i] as the cohort index.
    log(amult[i]) <- log(a) +
      betaSEX[1] * sex[i] * zSEX[1] +
      betaFEED[1] * feeding[i] * zFEED[1] +
      betaSPECIES[1] * species[i] * zSPECIES[1]
      #betaDENSITY[1] * density[b_year[i], island_index[i]] * zDENSITY[1] +
      #u_a_center[b_year[i]]
    
    log(bmult[i]) <- log(b) +
      betaSEX[2] * sex[i] * zSEX[2] +
      betaFEED[2] * feeding[i] * zFEED[2] +
      betaSPECIES[2] * species[i] * zSPECIES[2]
    #  betaDENSITY[2] * density[b_year[i], island_index[i]] * zDENSITY[2]
    
    # Observation Model
    # The first possible occasion to SEE the animal is the later of:
    #   a) The year it was born (as an integer index).
    #   b) The year surveys started on its island (as an integer index).
    
    # The integer year index of birth is simply b_year[i].
    # We need to provide the integer start year for the island.
    nm[i] <- max(b_year[i], IslandStart[i])
    
    # The number of observation opportunities
    # tMax is now the total number of years (e.g., 65).
    nM[i] <- min(floor(tD[i]), tMax) - nm[i] + 1
    pd[i] <- exp(y[i] * log(mean.p) + (nM[i] - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  # --- Priors and Random Walk
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
    betaSEX[k]  ~ dnorm(0, sd = 1)
    betaFEED[k] ~ dnorm(0, sd = 1)
    betaSPECIES[k] ~ dnorm(0, sd = 1)
    zSEX[k] ~ dbern(0.5)
    zFEED[k] ~ dbern(0.5)
    zSPECIES[k] ~ dbern(0.5)
  }
  
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
  sigma_year_a ~ dunif(0, 5)

  })

censoredD <- rep(2, nind)
censoredD[cintD[,1] > 0] <- 1
n_year <- max(year_index)
IslandStart <- AllData %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(IslandStart)

consts <- list(
  nind = nind,
  tMax = ncol(CH),
  n_year = n_year,
  sex = sex,
  feeding = feeding,
  species = species,
  island_index = island,
  IslandStart = IslandStart
)

data <- list(
  y = rowSums(CH),
  cintB = cintB,
  cintD = cintD,
  censoredD = censoredD,
  b_year = b_year,
  density = scaled_density_matrix
)

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data)
#model$initializeInfo()

## compile the model
cModel <- compileNimble(model)

## find list of valid initial values using compiled model
inits <- list()
for(k in 1:2) {
  inits[[k]] <- initFn(cModel, cintB, cintD, censoredD, n_year)
}

## configure MCMC
config <- configureMCMC(model, monitors = c("a", "b", "mean.p", "betaSEX", "betaFEED", "betaSPECIES",
                                            "zSEX", "zFEED", "zSPECIES", "sigma_year_a"))
config$removeSamplers(c("a", "b"))
config$addSampler(target = c("a", "b"), type = 'AF_slice')
#config$addSampler(target = c("betaFEED"), type = 'AF_slice')
#config$addSampler(target = c("b"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))
#config$addSampler(target = c("c1"), type = 'slice', control = list(sliceWidth = 1.5, adaptInterval = 50))

## load in custom RJ-MCMC samplers
#source("MCMC_RJ_multi.R")

# Add reversible jump
configureRJ(conf = config,   ## model configuration
            targetNodes = c("betaSEX", "betaFEED", "betaSPECIES"),
            indicatorNodes = c("zSEX", "zFEED", "zSPECIES"),
            control = list(mean = 0, scale = 0.5))

rIndicatorMCMC <- buildMCMC(config)
cIndicatorMCMC <- compileNimble(rIndicatorMCMC, project = model)

run <- runMCMC(cIndicatorMCMC,
               niter = 100000,
               nburnin = 22000,
               thin = 1,
               nchains = 2, 
               progressBar = TRUE,
               samplesAsCodaMCMC = TRUE)

samples <- as.matrix(run, chains = TRUE)
saveRDS(samples, "DoubleCensored/All_Island_Data/Samples/G_AllIslands_SexFeeding2_RETime_a.rds")
out <- readRDS("DoubleCensored/All_Island_Data/Samples/G_AllIslands_SexFeeding2Species_RETime_a.rds")
plot(as.mcmc(out))
out <- as.mcmc(out)
mcmcplot(run, parms = c("a", "b","betaSEX", "betaFEED", "betaSPECIES"))
MCMCsummary(run)
