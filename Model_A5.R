## ============================================================
## UNIFIED PIPELINE: SVL-Tightened Decadal Gompertz Model
## ============================================================
library(nimble)
library(MCMCvis)
library(tidyverse)
library(data.table)
library(quantreg)

# --- 1. DATA LOADING & CLEANING ---
# Load base data (assuming AllData, CH, cintB, cintD, etc. are in this RData)
load("igs_AllIslands_CleanDH_280426_obsStart.RData")
hab_area <- fread("Data/Iguana_islands_habitat.csv") %>% select(islands, Useable)

# Custom Distributions (Ensure these files are in your working directory)
source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

# --- 1. DATA PREPARATION & CLEANING ---
# Ensure SVL and Age are numeric
AllData <- AllData %>% 
  mutate(svl = as.numeric(as.character(svl)), 
         age = as.numeric(as.character(age)))

# Create individual-level metadata (one row per ID, ordered by ID)
id_df <- AllData %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  arrange(ID) %>%
  mutate(sex_num = as.numeric(as.factor(sex)) - 1,
         feed_num = ifelse(as.character(feeding) == "none", 0, 1))

nind <- nrow(id_df)
ordered_islands <- levels(as.factor(id_df$Island))
n_islands <- length(ordered_islands)
island_idx <- as.numeric(factor(id_df$Island, levels = ordered_islands))

# --- 2. GROWTH CALIBRATION (SVL -> Age) ---
calibration_data <- AllData %>%
  filter(!is.na(age), !is.na(svl), age > 0, !is.na(sex)) %>%
  mutate(log_age = log(age), sex = as.factor(sex))

# 95% Biological Age Envelope
fit_min_log_age <- rq(log_age ~ svl + I(svl^2) + sex, tau = 0.025, data = calibration_data)
fit_max_log_age <- rq(log_age ~ svl + I(svl^2) + sex, tau = 0.975, data = calibration_data)

# --- 3. GROWTH INFO & BIRTH WINDOWS (The "Squeeze") ---

# Safe summary: handles individuals with NO SVL data without crashing
growth_info <- AllData %>%
  group_by(ID) %>%
  summarise(
    # Returns the first non-NA SVL or NA if none exist
    first_svl = svl[!is.na(svl)][1],
    first_svl_year = Year[!is.na(svl)][1],
    # Total growth during study
    total_growth = if(any(!is.na(svl))) max(svl, na.rm=TRUE) - min(svl, na.rm=TRUE) else 0,
    tF = min(Year, na.rm = TRUE), # First capture year
    .groups = "drop"
  ) %>% 
  arrange(ID)

# Predict birth windows using the Global Quantile Regression from previous step
growth_info <- growth_info %>%
  mutate(sex = as.factor(id_df$sex), svl = first_svl) %>%
  mutate(
    # Predict with fallback for individuals with NO SVL data
    p_min = if_else(!is.na(svl), exp(predict(fit_min_log_age, newdata = .)), 0),
    p_max = if_else(!is.na(svl), exp(predict(fit_max_log_age, newdata = .)), 35)
  ) %>%
  mutate(
    # Logic Gate A: Juvenile Cap (Small SVL < 25cm cannot be old)
    p_max = if_else(!is.na(svl) & svl < 25, pmin(p_max, 4), p_max),
    # Logic Gate B: Growth Squeeze (If they grew > 5cm, they were young at start)
    p_max = if_else(total_growth > 5, pmin(p_max, 15), p_max)
  ) %>%
  mutate(pred_min_age = pmax(0, p_min), 
         pred_max_age = pmax(pred_min_age + 1, p_max))

# --- 4. TIMELINE ALIGNMENT (The "Origin" Fix) ---

# Calculate birth bounds relative to your current Year scale (since 1950)
cintB_raw <- matrix(NA, nrow = nind, ncol = 2)
for(i in 1:nind) {
  if(!is.na(tKB[i])) {
    cintB_raw[i, ] <- tKB[i] + c(-1, 0) # Known birth
  } else {
    # Use SVL year if available, otherwise first capture year
    ref_yr <- if(!is.na(growth_info$first_svl_year[i])) growth_info$first_svl_year[i] else growth_info$tF[i]
    cintB_raw[i, 1] <- floor(ref_yr - growth_info$pred_max_age[i])
    cintB_raw[i, 2] <- floor(ref_yr - growth_info$pred_min_age[i])
  }
}

# Define Year 1 as the earliest possible birth year in the dataset
Origin <- min(cintB_raw[, 1]) 

# Shift everything so the model starts at Year 1
cintB_adj <- cintB_raw - Origin + 1
cintD_adj <- cintD - Origin + 1
last_occ_adj <- apply(CH, 1, function(x) max(which(x == 1))) - Origin + 1
IslandStart_adj <- id_df$IslandStart - Origin + 1
is_start_vec_adj <- pmax(1, floor(IslandStart_adj) + 1)
tMax_adj <- max(AllData$Year) - Origin + 1

# Decadal lookup (10-year chunks)
decade_lookup <- floor((1:tMax_adj - 1) / 10) + 1
n_decades <- max(decade_lookup)

# Final logical check for cintB (Birth cannot be after first capture)
for(i in 1:nind) {
  tF_adj <- growth_info$tF[i] - Origin + 1
  cintB_adj[i, 2] <- min(cintB_adj[i, 2], tF_adj)
  if(cintB_adj[i, 1] >= cintB_adj[i, 2]) cintB_adj[i, 1] <- cintB_adj[i, 2] - 1
  cintB_adj[i, 1] <- pmax(1, cintB_adj[i, 1])
}

# --- 5. COVARIATES & EFFORT ---

# Standardize Log-Area
area_raw <- hab_area$Useable[match(ordered_islands, hab_area$islands)][island_idx]
area_z <- as.numeric(scale(log(area_raw)))

# Build Effort Matrix on the adjusted timeline
effort_matrix <- matrix(0, nrow = n_islands, ncol = tMax_adj)
actual_surveys <- AllData %>% 
  distinct(Island, Year) %>% 
  mutate(idx = as.numeric(factor(Island, levels = ordered_islands)),
         yr_adj = Year - Origin + 1)
for(r in 1:nrow(actual_surveys)) { effort_matrix[actual_surveys$idx[r], actual_surveys$yr_adj[r]] <- 1 }
cum_effort_adj <- t(apply(effort_matrix, 1, cumsum))
cum_effort_adj <- cbind(rep(0, n_islands), cum_effort_adj)

# --- 6. NIMBLE MODEL CODE ---

code_Final <- nimbleCode({
  for (i in 1:nind) {
    ## Latent Times
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    cohort_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    L[i] <- max(0, IslandStart_ind[i] - tB[i])
    
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    ## log(a): Baseline Hazard
    log(amult[i]) <- log_a + betaSEX[1]*sex[i] + betaFEED[1]*feeding[i] + 
      betaAREA[1]*area_z[i] + u_decade[decade_lookup[cohort_idx[i]]]
    
    ## log(b): Senescence Rate
    log(bmult[i]) <- log_b + betaSEX[2]*sex[i] + betaFEED[2]*feeding[i] + 
      betaAREA[2]*area_z[i]
    
    ## Observation Model
    nm_start[i] <- max(cohort_idx[i], is_start_vec[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    nMpos_raw[i] <- max(0, cum_effort[island[i], end_year[i] + 1] - cum_effort[island[i], nm_start[i]])
    nMpos_safe[i] <- max(y_val[i], nMpos_raw[i])
    y[i] ~ dbin(mean.p, nMpos_safe[i])
  }
  
  u_decade[1] <- 0 
  for (d in 2:n_decades) { u_decade[d] ~ dnorm(u_decade[d-1], sd = sigma_decade) }
  sigma_decade ~ dexp(1)
  
  for(k in 1:2) {
    betaSEX[k] ~ dnorm(0, sd = 1.5); betaFEED[k] ~ dnorm(0, sd = 1.5); betaAREA[k] ~ dnorm(0, sd = 1.5)
  }
  log_a ~ dnorm(-8, sd = 4); log_b ~ dnorm(-2, sd = 2); mean.p ~ dunif(0, 1)
})

# --- 7. MCMC EXECUTION ---

init_Robust <- function() {
  tBinit <- runif(nind, cintB_adj[, 1], cintB_adj[, 2])
  tDinit <- pmax(tBinit + 0.1, last_occ_adj + 0.1)
  list(tB = tBinit, tstar = tDinit - tBinit, u_decade = c(NA, rnorm(n_decades-1, 0, 0.1)),
       sigma_decade = 0.1, log_a = -8, log_b = -2, mean.p = 0.2, 
       betaSEX = c(0,0), betaFEED = c(0,0), betaAREA = c(0,0))
}

consts <- list(nind = nind, tMax = tMax_adj, n_decades = n_decades, sex = id_df$sex_num, 
               feeding = id_df$feed_num, island = island_idx, area_z = area_z,
               y_val = as.numeric(y), is_start_vec = as.numeric(is_start_vec_adj), 
               IslandStart_ind = IslandStart_adj)

data <- list(y = as.numeric(y), cintB = cintB_adj, cintD = cintD_adj, 
             censoredD = censoredD, cum_effort = cum_effort_adj, decade_lookup = decade_lookup)

model <- nimbleModel(code_Final, constants = consts, data = data, inits = init_Robust())
cModel <- compileNimble(model)
conf <- configureMCMC(model, monitors = c("log_a", "log_b", "mean.p", "betaSEX", "betaFEED", "betaAREA", "sigma_decade"))

# Optimized Samplers
for (i in 1:nind) {
  conf$removeSamplers(paste0('tB[', i, ']'))
  conf$removeSamplers(paste0('tstar[', i, ']'))
  conf$addSampler(target = paste0('tB[', i, ']'), type = "slice")
  conf$addSampler(target = paste0('tstar[', i, ']'), type = "slice")
}

#conf$removeSamplers(c("log_a", "u_decade"))
#conf$addSampler(target = c("log_a", "u_decade[2:n_decades]"), type = "AF_slice")
#conf$removeSamplers(c("log_b", "betaSEX", "betaFEED", "betaAREA"))
#conf$addSampler(target = c("log_b", "betaSEX", "betaFEED", "betaAREA"), type = "AF_slice")

# 3. Use AF_slice for global parameters
conf$removeSamplers(c("log_a", "betaSEX[1]", "betaFEED[1]", "betaAREA[1]"))
conf$addSampler(target = c("log_a", "betaSEX[1]", "betaFEED[1]", "betaAREA"), type = "AF_slice")
conf$removeSamplers(c("log_b", "betaSEX[2]", "betaFEED[2]", "betaAREA[2]"))
conf$addSampler(target = c("log_b", "betaSEX[2]", "betaFEED[2]", "betaAREA[2]"), type = "AF_slice")

conf$getUnsampledNodes()

cMCMC <- compileNimble(buildMCMC(conf), project = model)
results <- runMCMC(cMCMC, niter = 100000, nburnin = 40000, nchains = 2, summary = TRUE)

# --- 8. OUTPUT ---
MCMCsummary(results$samples, params = c("log_a", "log_b", "betaSEX", "betaFEED", "betaAREA"))
MCMCtrace(results$samples, params = c("log_a", "log_b", "betaSEX", "betaFEED", "betaAREA"))
