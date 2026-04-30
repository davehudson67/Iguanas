## ============================================================
## FINAL INTEGRATED PIPELINE: Data Refinement + Isolated Feeding Model
## ============================================================
library(nimble)
library(MCMCvis)
library(tidyverse)
library(data.table)

# Load the base data
load("igs_AllIslands_CleanDH_280426_obsStart.RData")
hab_area <- fread("Data/Iguana_islands_habitat.csv") %>% select(islands, Useable)

## ---- Custom functions/distributions ----
source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

# Ensure YearZ is set and tMax is calculated correctly
YearZ <- 1950
tMax  <- max(AllData$Year)
nind  <- length(unique(AllData$ID))
ordered_islands <- levels(as.factor(AllData$Island))
n_islands <- length(ordered_islands)

## ------------------------------------------------------------
## 1) REFINED DENSITY MATRIX (Neutral Filling)
## ------------------------------------------------------------
# 1. Calculate MNA per island per year
mna_dens <- AllData %>%
  filter(Status == 1) %>%
  group_by(ID, Island) %>%
  summarise(First = min(Year), Last = max(Year), .groups = "drop") %>%
  mutate(Year = map2(First, Last, seq)) %>%
  unnest(Year) %>%
  group_by(Island, Year) %>%
  summarise(MNA = n(), .groups = "drop") %>%
  left_join(hab_area %>% rename(Island = islands), by = "Island") %>%
  mutate(Density = MNA / Useable)

# 2. Map to matrix and fill buffer years (1950 -> IslandStart)
dens_matrix <- matrix(NA, nrow = n_islands, ncol = tMax)
island_starts <- AllData %>% group_by(Island) %>% summarise(Start = min(IslandStart))

for(j in 1:n_islands) {
  island_nm <- ordered_islands[j]
  obs <- mna_dens %>% filter(Island == island_nm)
  
  # Fill observed years
  dens_matrix[j, obs$Year] <- obs$Density
  
  # Fill buffer years with the mean observed density (neutral prior)
  mean_d <- mean(obs$Density, na.rm = TRUE)
  start_yr <- island_starts$Start[island_starts$Island == island_nm]
  dens_matrix[j, 1:(start_yr-1)] <- mean_d
}

# 3. Standardize for model stability
dens_matrix_scaled <- (dens_matrix - mean(dens_matrix, na.rm=T)) / sd(dens_matrix, na.rm=T)
dens_matrix_scaled[is.na(dens_matrix_scaled)] <- 0 

## ------------------------------------------------------------
## 2) REFINED BIRTH BOUNDS (cintB) & EFFORT
## ------------------------------------------------------------
cintB <- matrix(NA, nrow = nind, ncol = 2)

# Known Ages (using tKB from your RData)
known_idx <- !is.na(tKB)
cintB[known_idx, 1] <- pmax(1, tKB[known_idx] - 1)
cintB[known_idx, 2] <- pmax(2, tKB[known_idx])

# Unknown Ages (Adults): Assume max age 30, birth >= 1950 (Year 1)
unknown_idx <- is.na(tKB)
cintB[unknown_idx, 1] <- pmax(1, tF[unknown_idx] - 30)
cintB[unknown_idx, 2] <- tF[unknown_idx]

# Effort Matrix (1..tMax)
effort_matrix <- matrix(0, nrow = n_islands, ncol = tMax)
actual_surveys <- AllData %>% distinct(Island, Year) %>% 
  mutate(idx = as.numeric(as.factor(Island)))

for(r in 1:nrow(actual_surveys)) {
  effort_matrix[actual_surveys$idx[r], actual_surveys$Year[r]] <- 1
}
cum_effort <- t(apply(effort_matrix, 1, cumsum))
cum_effort <- cbind(rep(0, n_islands), cum_effort)

# Create a data frame with one row per individual in the same order as CH
id_df <- AllData %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  arrange(ID) # CH is ordered by ID

# Define the unique island names and their numeric index
ordered_islands <- levels(as.factor(id_df$Island))
n_islands <- length(ordered_islands)

# Create island_idx (a vector of length nind, values 1 to 8)
island_idx <- as.numeric(factor(id_df$Island, levels = ordered_islands))

# Create the start year lookup (ordered by our island index)
island_start_lookup <- island_starts %>%
  arrange(match(Island, ordered_islands)) %>%
  pull(Start)

# Map the start year to every individual (length nind)
IslandStart_ind <- island_start_lookup[island_idx]

# Calculate the vector for NIMBLE (length nind)
is_start_vec <- pmax(1, pmin(tMax, floor(IslandStart_ind) + 1))

sex <- as.numeric(as.factor(id_df$sex)) - 1
feeding <- ifelse(as.character(id_df$feeding) == "none", 0, 1)

# Identify the last occasion each individual was seen (used for inits)
last_occasion <- apply(CH, 1, function(r) {
  w <- which(r > 0)
  if (length(w) == 0) 0L else max(w)
})

# --- SANITY CHECKS ---
message("nind: ", nind)
message("Length of island_idx: ", length(island_idx))
message("Length of IslandStart_ind: ", length(IslandStart_ind))
message("Length of is_start_vec: ", length(is_start_vec))

# Check if any are NA
if(any(is.na(island_idx))) stop("NA found in island_idx!")
if(any(is.na(is_start_vec))) stop("NA found in is_start_vec!")

## ------------------------------------------------------------
## 3) NIMBLE MODEL: ISOLATED FEEDING
## ------------------------------------------------------------
code_Final <- nimbleCode({
  for (i in 1:nind) {
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    cohort_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    L[i] <- max(0, IslandStart_ind[i] - tB[i])
    
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    log(amult[i]) <- log_a + betaSEX * sex[i] + betaFEED * feeding[i] +
      betaDENS * dens_matrix[island[i], cohort_idx[i]] +
      u_global[ cohort_idx[i] ]
    log(bmult[i]) <- log_b
    
    nm_start[i] <- max(cohort_idx[i], is_start_vec[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    nMpos_raw[i] <- max(0, cum_effort[ island[i], end_year[i] + 1 ] - 
                          cum_effort[ island[i], nm_start[i] ])
    nMpos_safe[i] <- max(y_val[i], nMpos_raw[i])
    y[i] ~ dbin(mean.p, nMpos_safe[i])
  }
  
  u_global[1] <- 0
  for (t in 2:tMax) { u_global[t] ~ dnorm(u_global[t-1], sd = sigma_drift) }
  sigma_drift ~ dexp(1)
  
  betaSEX ~ dnorm(0, sd = 1.5); betaFEED ~ dnorm(0, sd = 1.5); betaDENS ~ dnorm(0, sd = 1.5)
  log_a ~ dnorm(-5, sd = 5); log_b ~ dnorm(-2, sd = 5); mean.p ~ dunif(0, 1)
})

## ------------------------------------------------------------
## 4) EXECUTION
## ------------------------------------------------------------
consts <- list(nind = nind, tMax = tMax, sex = sex, feeding = feeding, island = island_idx,
               y_val = as.numeric(y), is_start_vec = as.numeric(is_start_vec), 
               IslandStart_ind = IslandStart_ind)

data <- list(y = as.numeric(y), cintB = cintB, cintD = cintD, censoredD = censoredD,
             dens_matrix = dens_matrix_scaled, cum_effort = cum_effort)

init_Final <- function() {
  tBinit <- runif(nind, cintB[,1], cintB[,2])
  tDinit <- pmax(tBinit + 0.1, last_occasion + 0.1)
  list(tB = tBinit, tstar = tDinit - tBinit, u_global = c(NA, rnorm(tMax-1, 0, 0.1)),
       sigma_drift = 0.1, log_a = -7, log_b = -2, mean.p = 0.2, 
       betaSEX = 0, betaFEED = 0, betaDENS = 0)
}

model <- nimbleModel(code_Final, constants = consts, data = data, inits = init_Final())
cModel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c("log_a", "log_b", "mean.p", "betaSEX", "betaFEED", "betaDENS", "sigma_drift", "u_global"))

# Sampler Optimization
for (i in 1:nind) {
  conf$removeSamplers(paste0('tB[', i, ']'))
  conf$removeSamplers(paste0('tstar[', i, ']'))
  conf$addSampler(target = c(paste0('tB[', i, ']'), paste0('tstar[', i, ']')), type = "RW_block")
}
#conf$removeSamplers("u_global")
#conf$addSampler(target = "u_global[2:tMax]", type = "RW_block")
conf$removeSamplers(c("log_a", "log_b", "betaSEX", "betaFEED", "betaDENS"))
conf$addSampler(target = c("log_a", "log_b", "betaSEX", "betaFEED", "betaDENS"), type = "AF_slice")

# Compile and Run
cMCMC <- compileNimble(buildMCMC(conf), project = model)
results <- runMCMC(cMCMC, niter = 150000, nburnin = 50000, nchains = 2, summary = TRUE)

MCMCsummary(results$samples, params = c("log_a", "log_b", "betaSEX", "betaFEED", "betaDENS", "sigma_drift"))