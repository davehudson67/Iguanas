## ============================================================
## COMPLETE PIPELINE: Decadal Dual-Effect Gompertz Model
## ============================================================
library(nimble)
library(MCMCvis)
library(tidyverse)
library(data.table)
library(lubridate)
library(coda)

# --- 1. DATA LOADING & CLEANING ---
# Load base data (assuming AllData, CH, cintB, cintD, etc. are in this RData)
load("igs_AllIslands_CleanDH_280426_obsStart.RData")
hab_area <- fread("Data/Iguana_islands_habitat.csv") %>% select(islands, Useable)

# Custom Distributions (Ensure these files are in your working directory)
source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

# --- 2. TIMELINE & INDIVIDUAL INDEXING ---
YearZ <- 1950
tMax  <- max(AllData$Year)
nind  <- length(unique(AllData$ID))

# One-row-per-individual dataframe (ordered by ID to match CH)
id_df <- AllData %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  arrange(ID) 

ordered_islands <- levels(as.factor(id_df$Island))
n_islands <- length(ordered_islands)
island_idx <- as.numeric(factor(id_df$Island, levels = ordered_islands))

# --- 3. DECADAL INDEXING ---
# Create a lookup: Year 1-10 -> Decade 1, etc.
decade_lookup <- floor((1:tMax - 1) / 10) + 1
n_decades <- max(decade_lookup)

# --- 4. ISLAND-LEVEL COVARIATES (AREA) ---
island_area_ref <- hab_area %>%
  rename(Island = islands) %>%
  filter(Island %in% ordered_islands) %>%
  mutate(idx = as.numeric(factor(Island, levels = ordered_islands))) %>%
  arrange(idx)

# Standardize log-area for model stability
area_raw <- island_area_ref$Useable[island_idx]
area_z   <- as.numeric(scale(log(area_raw))) 

# --- 5. BIRTH/DEATH BOUNDS & EFFORT ---
# Birth bounds (cintB): Max age 30, birth >= 1950 (Year 1)
cintB <- matrix(NA, nrow = nind, ncol = 2)
known_idx <- !is.na(tKB)
cintB[known_idx, 1] <- pmax(1, tKB[known_idx] - 1)
cintB[known_idx, 2] <- pmax(2, tKB[known_idx])

unknown_idx <- is.na(tKB)
cintB[unknown_idx, 1] <- pmax(1, tF[unknown_idx] - 30) 
cintB[unknown_idx, 2] <- tF[unknown_idx]

# Survey Effort Matrix (1..tMax)
effort_matrix <- matrix(0, nrow = n_islands, ncol = tMax)
actual_surveys <- AllData %>% 
  distinct(Island, Year) %>% 
  mutate(idx = as.numeric(factor(Island, levels = ordered_islands)))

for(r in 1:nrow(actual_surveys)) {
  effort_matrix[actual_surveys$idx[r], actual_surveys$Year[r]] <- 1
}
cum_effort <- t(apply(effort_matrix, 1, cumsum))
cum_effort <- cbind(rep(0, n_islands), cum_effort) # Column 1 is Year 0

# Island Start Year Vector
island_starts <- AllData %>% group_by(Island) %>% summarise(Start = min(IslandStart), .groups = "drop")
island_start_lookup <- island_starts %>% 
  mutate(idx = as.numeric(factor(Island, levels = ordered_islands))) %>%
  arrange(idx) %>%
  pull(Start)

IslandStart_ind <- island_start_lookup[island_idx]
is_start_vec <- pmax(1, pmin(tMax, floor(IslandStart_ind) + 1))

# Covariates
sex <- as.numeric(as.factor(id_df$sex)) - 1
feeding <- ifelse(as.character(id_df$feeding) == "none", 0, 1)

# Last seen (for inits)
last_occasion <- apply(CH, 1, function(r) {
  w <- which(r > 0); if (length(w) == 0) 0L else max(w)
})

## ------------------------------------------------------------
## 6) NIMBLE MODEL: DECADAL DUAL-EFFECT
## ------------------------------------------------------------
code_Decadal <- nimbleCode({
  for (i in 1:nind) {
    ## Latent Times
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    cohort_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    L[i] <- max(0, IslandStart_ind[i] - tB[i])
    
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    ## log(a): Baseline Hazard (Initial Mortality)
    log(amult[i]) <- log_a + 
      betaSEX[1]  * sex[i] + 
      betaFEED[1] * feeding[i] +
      betaAREA * area_z[i]
    
    ## log(b): Senescence Rate (Rate of Aging)
    log(bmult[i]) <- log_b + 
      betaSEX[2]  * sex[i] + 
      betaFEED[2] * feeding[i]
    
    ## Observation Model
    nm_start[i] <- max(cohort_idx[i], is_start_vec[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    
    nMpos_raw[i] <- max(0, cum_effort[ island[i], end_year[i] + 1 ] - 
                          cum_effort[ island[i], nm_start[i] ])
    nMpos_safe[i] <- max(y_val[i], nMpos_raw[i])
    y[i] ~ dbin(mean.p, nMpos_safe[i])
  }
  
  ## Decadal Random Walk (u_decade[1] is anchored to 0)
  #u_decade[1] <- 0
  #for (d in 2:n_decades) { 
  #  u_decade[d] ~ dnorm(u_decade[d-1], sd = sigma_decade) 
  #}
  #sigma_decade ~ dexp(1)
  
  ## Priors
  for(k in 1:2) {
    betaSEX[k]  ~ dnorm(0, sd = 1.5)
    betaFEED[k] ~ dnorm(0, sd = 1.5)
  }
  log_a ~ dnorm(-8, sd = 5) 
  log_b ~ dnorm(-2, sd = 5) 
  mean.p ~ dunif(0, 1)
  betaAREA ~ dnorm(0, sd = 1.5)
})

## ------------------------------------------------------------
## 7) EXECUTION
## ------------------------------------------------------------
consts <- list(
  nind = nind, 
  tMax = tMax, 
  n_decades = n_decades,
  sex = sex, 
  feeding = feeding, 
  island = island_idx,
  y_val = as.numeric(y), 
  is_start_vec = as.numeric(is_start_vec), 
  IslandStart_ind = IslandStart_ind, 
  area_z = area_z
  # decade_lookup REMOVED from here
)

data <- list(
  y = as.numeric(y), 
  cintB = cintB_ultra, 
  cintD = cintD, 
  censoredD = censoredD,
  cum_effort = cum_effort,
  # decade_lookup ADDED here to allow dynamic indexing
  decade_lookup = decade_lookup 
)

init_Decadal <- function() {
  tBinit <- runif(nind, cintB[,1], cintB[,2])
  tDinit <- pmax(tBinit + 0.1, last_occasion + 0.1)
  list(tB = tBinit, tstar = tDinit - tBinit, 
       u_decade = c(NA, rnorm(n_decades-1, 0, 0.1)),
       sigma_decade = 0.1, log_a = -8, log_b = -2, mean.p = 0.2, 
       betaSEX = c(0,0), betaFEED = c(0,0), betaAREA = 0)
}

# Build and Compile Model
model <- nimbleModel(code_Decadal, constants = consts, data = data, inits = init_Decadal())
cModel <- compileNimble(model)

# Configure MCMC
conf <- configureMCMC(model, monitors = c("log_a", "log_b", "mean.p", "betaSEX", 
                                          "betaFEED", "betaAREA"))

# --- SAMPLER OPTIMIZATION ---
# 1. Block individuals
for (i in 1:nind) {
  conf$removeSamplers(paste0('tB[', i, ']'))
  conf$removeSamplers(paste0('tstar[', i, ']'))
  conf$addSampler(target = c(paste0('tB[', i, ']'), paste0('tstar[', i, ']')), type = "AF_slice")
}

# 2. Block Intercept with the Decadal Random Walk
#conf$removeSamplers(c("log_a", "u_decade"))
#conf$addSampler(target = c("log_a", "u_decade[2:n_decades]"), type = "AF_slice")

# 3. Use AF_slice for global parameters
conf$removeSamplers(c("log_a", "betaSEX[1]", "betaFEED[1]", "betaAREA[1]"))
conf$addSampler(target = c("log_a", "betaSEX[1]", "betaFEED[1]", "betaAREA"), type = "AF_slice")
conf$removeSamplers(c("log_b", "betaSEX[2]", "betaFEED[2]", "betaAREA[2]"))
conf$addSampler(target = c("log_b", "betaSEX[2]", "betaFEED[2]"), type = "AF_slice")

# Compile and Run
cMCMC <- compileNimble(buildMCMC(conf), project = model)
results <- runMCMC(cMCMC, niter = 60000, nburnin = 15000, nchains = 2, summary = TRUE)
results$summary
saveRDS(results, "Model_A4_Samples.rds")

# --- 8) OUTPUT ---
MCMCtrace(results$samples)
MCMCsummary(results$samples)

AllData <- group_by(AllData, ID)
