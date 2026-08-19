## ============================================================
## FINAL INTEGRATED PIPELINE: Dual-Effect Gompertz Model
## ============================================================
library(nimble)
library(MCMCvis)
library(tidyverse)
library(data.table)

# Load data
load("igs_AllIslands_CleanDH_280426_obsStart.RData")
hab_area <- fread("Data/Iguana_islands_habitat.csv") %>% select(islands, Useable)

source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

# --- 1. TIMELINE & INDIVIDUAL INDEXING ---
YearZ <- 1950
tMax  <- max(AllData$Year)
nind  <- length(unique(AllData$ID))

id_df <- AllData %>% distinct(ID, .keep_all = TRUE) %>% arrange(ID) 

# Explicitly define island levels to ensure consistency across all matrices
ordered_islands <- levels(as.factor(id_df$Island))
n_islands <- length(ordered_islands)
island_idx <- as.numeric(factor(id_df$Island, levels = ordered_islands))

# --- 2. ISLAND-LEVEL COVARIATE (AREA) ---
island_area_ref <- hab_area %>%
  rename(Island = islands) %>%
  filter(Island %in% ordered_islands) %>%
  mutate(island_idx = as.numeric(factor(Island, levels = ordered_islands))) %>%
  arrange(island_idx) # Ensure order matches 1..n_islands

area_raw <- island_area_ref$Useable[island_idx]
area_z   <- as.numeric(scale(log(area_raw))) 

# --- 3. BIRTH BOUNDS & EFFORT ---
cintB <- matrix(NA, nrow = nind, ncol = 2)
known_idx <- !is.na(tKB)
cintB[known_idx, 1] <- pmax(1, tKB[known_idx] - 1)
cintB[known_idx, 2] <- pmax(2, tKB[known_idx])

unknown_idx <- is.na(tKB)
cintB[unknown_idx, 1] <- pmax(1, tF[unknown_idx] - 30) 
cintB[unknown_idx, 2] <- tF[unknown_idx]

# --- FIX: Survey Effort Matrix Indexing ---
effort_matrix <- matrix(0, nrow = n_islands, ncol = tMax)
actual_surveys <- AllData %>% 
  distinct(Island, Year) %>% 
  # Use the SAME factor levels as island_idx
  mutate(idx = as.numeric(factor(Island, levels = ordered_islands))) %>%
  filter(!is.na(idx))

for(r in 1:nrow(actual_surveys)) {
  effort_matrix[actual_surveys$idx[r], actual_surveys$Year[r]] <- 1
}
cum_effort <- t(apply(effort_matrix, 1, cumsum))
cum_effort <- cbind(rep(0, n_islands), cum_effort)

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
## 4) NIMBLE MODEL
## ------------------------------------------------------------
code_DualEffect <- nimbleCode({
  for (i in 1:nind) {
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    cohort_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    L[i] <- max(0, IslandStart_ind[i] - tB[i])
    
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    log(amult[i]) <- log_a + 
      betaSEX[1]  * sex[i] + 
      betaFEED[1] * feeding[i] +
      betaAREA[1] * area_z[i] +
      u_global[ cohort_idx[i] ]
    
    log(bmult[i]) <- log_b + 
      betaSEX[2]  * sex[i] + 
      betaFEED[2] * feeding[i] +
      betaAREA[2] * area_z[i]
    
    nm_start[i] <- max(cohort_idx[i], is_start_vec[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    
    # Accessing cum_effort[island, column]
    # Column 1 is Year 0, Column tMax+1 is Year tMax
    nMpos_raw[i] <- max(0, cum_effort[ island[i], end_year[i] + 1 ] - 
                          cum_effort[ island[i], nm_start[i] ])
    nMpos_safe[i] <- max(y_val[i], nMpos_raw[i])
    y[i] ~ dbin(mean.p, nMpos_safe[i])
  }
  
  u_global[1] <- 0
  for (t in 2:tMax) { 
    u_global[t] ~ dnorm(u_global[t-1], sd = sigma_drift) 
  }
  sigma_drift ~ dexp(1)
  
  for(k in 1:2) {
    betaSEX[k]  ~ dnorm(0, sd = 1.5)
    betaFEED[k] ~ dnorm(0, sd = 1.5)
    betaAREA[k] ~ dnorm(0, sd = 1.5)
  }
  
  log_a ~ dnorm(-8, sd = 5) 
  log_b ~ dnorm(-2, sd = 5) 
  mean.p ~ dunif(0, 1)
})

## ------------------------------------------------------------
## 5) EXECUTION
## ------------------------------------------------------------
consts <- list(nind = nind, tMax = tMax, sex = sex, feeding = feeding, island = island_idx,
               y_val = as.numeric(y), is_start_vec = as.numeric(is_start_vec), 
               IslandStart_ind = IslandStart_ind, area_z = area_z)

data <- list(y = as.numeric(y), cintB = cintB, cintD = cintD, censoredD = censoredD,
             cum_effort = cum_effort)

init_Final <- function() {
  tBinit <- runif(nind, cintB[,1], cintB[,2])
  tDinit <- pmax(tBinit + 0.1, last_occasion + 0.1)
  list(tB = tBinit, tstar = tDinit - tBinit, u_global = c(NA, rnorm(tMax-1, 0, 0.1)),
       sigma_drift = 0.1, log_a = -8, log_b = -2, mean.p = 0.2, 
       betaSEX = c(0,0), betaFEED = c(0,0), betaAREA = c(0,0))
}

model <- nimbleModel(code_DualEffect, constants = consts, data = data, inits = init_Final())
cModel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c("log_a", "log_b", "mean.p", "betaSEX", 
                                          "betaFEED", "betaAREA", "sigma_drift", "u_global"))

# Sampler Optimization
for (i in 1:nind) {
  conf$removeSamplers(paste0('tB[', i, ']'))
  conf$removeSamplers(paste0('tstar[', i, ']'))
  conf$addSampler(target = c(paste0('tB[', i, ']'), paste0('tstar[', i, ']')), type = "RW_block")
}

#conf$removeSamplers(c("log_a", "u_global"))
#conf$addSampler(target = c("log_a", "u_global[2:tMax]"), type = "RW_block")

conf$removeSamplers(c("log_a", "log_b"))
conf$addSampler(target = c("log_a", "log_b"), type = "AF_slice")

cMCMC <- compileNimble(buildMCMC(conf), project = model)
results <- runMCMC(cMCMC, niter = 200000, nburnin = 75000, nchains = 2, summary = TRUE)