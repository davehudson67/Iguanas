## ============================================================
## IMPROVED MODEL A: Gompertz LB + Cohort RW1 + Blocked Sampling
## ============================================================
library(data.table)
library(tidyverse)
library(nimble)
library(mcmcplots)
library(MCMCvis)

setwd("~/Iguanas")

## ---- Load data ----
load("igs_AllIslands_CleanDH_280426_obsStart.RData")
hab_area <- fread("Data/Iguana_islands_habitat.csv") %>% select(islands, Useable)

## ---- Custom functions/distributions ----
source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

## ============================================================
## 1) DATA CORRECTION: FFRC 'figginsi' -> 'inornata'
## ============================================================
problem_ids <- AllData %>%
  filter(Island == "FFRC", Species == "figginsi") %>%
  pull(ID)

if (length(problem_ids) > 0) {
  message("Correcting Species for FFRC IDs: ", paste(problem_ids, collapse = ", "))
  AllData <- AllData %>%
    mutate(Species = ifelse(ID %in% problem_ids & Island == "FFRC", "inornata", Species))
} else {
  message("No FFRC 'figginsi' records found to correct.")
}

## ============================================================
## 2) ONE-ROW-PER-ID covariates
## ============================================================
id_df <- AllData %>% distinct(ID, .keep_all = TRUE)
nind  <- nrow(id_df)
tMax  <- ncol(CH)

stopifnot(nrow(CH) == nind, nrow(cintB) == nind, nrow(cintD) == nind)

ordered_islands <- levels(as.factor(id_df$Island))
n_islands <- length(ordered_islands)

sex <- as.numeric(as.factor(id_df$sex)) - 1
feeding <- ifelse(as.character(id_df$feeding) == "none", 0, 1)
island_idx <- as.numeric(factor(id_df$Island, levels = ordered_islands))

island_start_by_island <- id_df %>%
  distinct(Island, .keep_all = TRUE) %>%
  mutate(Island = as.character(Island)) %>%
  arrange(factor(Island, levels = ordered_islands)) %>%
  pull(IslandStart)

stopifnot(length(island_start_by_island) == n_islands)
IslandStart_ind <- island_start_by_island[island_idx]

censoredD <- rep(2, nind)
censoredD[cintD[,1] > 0] <- 1

y <- rowSums(CH)
dind <- rep(1, nind)

last_occasion <- apply(CH, 1, function(r) {
  w <- which(r > 0)
  if (length(w) == 0) 0L else max(w)
})

## ============================================================
## 3) TIME-VARYING DENSITY MATRIX
## ============================================================
area_vec <- hab_area %>%
  mutate(islands = as.character(islands)) %>%
  right_join(tibble(islands = ordered_islands), by = "islands") %>%
  pull(Useable)

if (any(is.na(area_vec))) {
  stop("Missing Useable area for islands: ", paste(ordered_islands[is.na(area_vec)], collapse = ", "))
}
if (any(area_vec <= 0, na.rm = TRUE)) stop("Non-positive island area(s) in hab_area.")

id_df <- id_df %>% mutate(row_idx = 1:nind)

first_last <- as.data.frame(CH) %>%
  mutate(row_idx = 1:n()) %>%
  pivot_longer(-row_idx, names_to = "Year_Idx", values_to = "Detected") %>%
  mutate(Year_Idx = as.integer(str_remove(Year_Idx, "V"))) %>%
  left_join(id_df %>% select(row_idx, Island), by = "row_idx") %>%
  group_by(row_idx, Island) %>%
  summarise(
    First = suppressWarnings(min(Year_Idx[Detected == 1], na.rm = TRUE)),
    Last  = suppressWarnings(max(Year_Idx[Detected == 1], na.rm = TRUE)),
    anyDet = any(Detected == 1),
    .groups = "drop"
  ) %>%
  filter(anyDet)

mna_df <- first_last %>%
  mutate(Year_Idx = map2(First, Last, seq)) %>%
  unnest(Year_Idx) %>%
  group_by(Island, Year_Idx) %>%
  summarise(MNA = n(), .groups = "drop") %>%
  mutate(island_idx = as.numeric(factor(Island, levels = ordered_islands))) %>%
  left_join(tibble(island_idx = 1:n_islands, Useable = area_vec), by = "island_idx") %>%
  mutate(Density = MNA / Useable) %>%
  filter(Year_Idx >= 1, Year_Idx <= tMax)

dens_matrix <- matrix(0, nrow = n_islands, ncol = tMax)

mna_scaled <- mna_df %>%
  group_by(island_idx) %>%
  mutate(dens_z = as.numeric(scale(log(Density)))) %>%
  ungroup()

mna_scaled$dens_z[!is.finite(mna_scaled$dens_z)] <- 0

for (r in 1:nrow(mna_scaled)) {
  dens_matrix[mna_scaled$island_idx[r], mna_scaled$Year_Idx[r]] <- mna_scaled$dens_z[r]
}

## ============================================================
## 4) SURVEY EFFORT MATRIX
## ============================================================
## sanity: if AllData$Year is not 1..tMax, you must remap it
# range(AllData$Year, na.rm=TRUE)

actual_surveys <- AllData %>%
  filter(Status == 1) %>%
  distinct(Island, Year) %>%
  mutate(
    island_idx = as.numeric(factor(Island, levels = ordered_islands)),
    year_idx   = as.integer(Year)
  ) %>%
  filter(year_idx >= 1, year_idx <= tMax)

effort_matrix <- matrix(0, nrow = n_islands, ncol = tMax)
effort_matrix[as.matrix(actual_surveys[, c("island_idx", "year_idx")])] <- 1

cum_effort <- t(apply(effort_matrix, 1, cumsum))
cum_effort <- cbind(rep(0, nrow(cum_effort)), cum_effort)
cum_effort <- matrix(as.numeric(cum_effort), nrow = n_islands, ncol = tMax + 1)

## ============================================================
## 5) COHORT INDEX RANGE
## ============================================================
n_cohort <- tMax


## ============================================================
## FINAL INTEGRATED MODEL A: Gompertz LB + Cohort RW1
## ============================================================
library(nimble)
library(MCMCvis)
library(tidyverse)
library(data.table)

# --- 1. PRE-CALCULATIONS (R Environment) ---

# Ensure IslandStart is a vector of length nind
is_start_vec <- pmax(1, pmin(tMax, floor(IslandStart_ind) + 1))

# Identify the last occasion each individual was seen (used for inits)
last_occasion <- apply(CH, 1, function(r) {
  w <- which(r > 0)
  if (length(w) == 0) 0L else max(w)
})

# --- 2. NIMBLE MODEL CODE ---

code_A <- nimbleCode({
  
  for (i in 1:nind) {
    
    ## Latent birth/entry time (occasion units)
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    
    ## Discrete cohort index (1..tMax)
    cohort_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    
    ## Left truncation Age (Age at which surveys started on that island)
    L[i] <- max(0, IslandStart_ind[i] - tB[i])
    
    ## Latent lifespan (Age at death)
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    
    ## Death time
    tD[i] <- tB[i] + tstar[i]
    
    ## Interval / right censoring constraint
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    ## Density-at-entry lookup (Dynamic indexing)
    dens_entry[i] <- dens_matrix[island[i], cohort_idx[i]]
    
    ## Gompertz parameters (Anchored RW1: u_cohort[1]=0)
    log(amult[i]) <- log(a) +
      betaSEX[1]  * sex[i]        * zSEX[1] +
      betaFEED[1] * feeding[i]    * zFEED[1] +
      betaDENS[1] * dens_entry[i] * zDENS[1] +
      u_cohort[ cohort_idx[i] ]
    
    log(bmult[i]) <- log(b) +
      betaSEX[2]  * sex[i]        * zSEX[2] +
      betaFEED[2] * feeding[i]    * zFEED[2] +
      betaDENS[2] * dens_entry[i] * zDENS[2]
    
    ## Observation model indices
    # When surveys could first see this specific animal
    nm_start[i] <- max(cohort_idx[i], is_start_vec[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    
    # Map to cum_effort indices (which are 1..tMax+1)
    effort_idx_end[i]   <- max(1, min(tMax + 1, end_year[i] + 1))
    effort_idx_start[i] <- max(1, min(tMax + 1, nm_start[i]))
    
    # Calculate number of opportunities
    nMpos_raw[i] <- max(0, cum_effort[ island[i], effort_idx_end[i] ] - 
                          cum_effort[ island[i], effort_idx_start[i] ])
    
    # Safety check (Ensures trials >= successes for dbin)
    # y_val is a constant copy of y used to break the graph cycle
    nMpos_safe[i] <- max(y_val[i], nMpos_raw[i])
    
    ## Likelihood
    y[i] ~ dbin(mean.p, nMpos_safe[i])
  }
  
  ## --- Cohort RW1: Anchored ---
  u_cohort[1] <- 0
  for (t in 2:tMax) {
    u_cohort[t] ~ dnorm(u_cohort[t-1], sd = sigma_year_a)
  }
  sigma_year_a ~ T(dnorm(0, sd = 0.5), 0, )
  
  ## --- Priors ---
  for (k in 1:2) {
    betaSEX[k]  ~ dnorm(0, sd = 1.5)
    betaFEED[k] ~ dnorm(0, sd = 1.5)
    betaDENS[k] ~ dnorm(0, sd = 1.5)
    
    zSEX[k]  ~ dbern(0.5)
    zFEED[k] ~ dbern(0.5)
    zDENS[k] ~ dbern(0.5)
  }
  
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
})

# --- 3. DATA AND CONSTANTS ---

consts_A <- list(
  nind = nind,
  tMax = tMax,
  sex = sex,
  feeding = feeding,
  island = island_idx,
  y_val = as.numeric(y), # Constant for cycle-breaking
  is_start_vec = as.numeric(is_start_vec),
  IslandStart_ind = IslandStart_ind
)

data_A <- list(
  y = as.numeric(y),
  cintB = as.matrix(cintB),
  cintD = as.matrix(cintD),
  censoredD = as.numeric(censoredD),
  dens_matrix = dens_matrix,
  cum_effort = cum_effort
)

# --- 4. INITIALIZATION FUNCTION ---

init_A <- function() {
  tBinit <- runif(nind, cintB[,1], cintB[,2])
  tDinit <- numeric(nind)
  for(i in 1:nind) {
    # Ensure death is after birth and after the animal was last seen alive
    low_limit <- max(tBinit[i], last_occasion[i]) + 0.1
    if (censoredD[i] == 1) { 
      tDinit[i] <- runif(1, max(low_limit, cintD[i,1]), cintD[i,2])
    } else { 
      tDinit[i] <- low_limit + rexp(1, 0.5)
    }
  }
  list(
    tB = tBinit,
    tstar = tDinit - tBinit,
    u_cohort = c(NA, rnorm(tMax-1, 0, 0.1)),
    sigma_year_a = 0.1,
    a = 0.5,
    b = 0.05,
    mean.p = 0.2,
    betaSEX  = rnorm(2, 0, 0.1),
    betaFEED = rnorm(2, 0, 0.1),
    betaDENS = rnorm(2, 0, 0.1),
    zSEX  = c(0,0),
    zFEED = c(0,0),
    zDENS = c(0,0)
  )
}

# --- 5. MCMC CONFIGURATION ---

# Build model
model_A <- nimbleModel(code_A, constants = consts_A, data = data_A, inits = init_A())
cModel_A <- compileNimble(model_A)

conf_A <- configureMCMC(model_A,
                        monitors = c("a","b","mean.p",
                                     "betaSEX","betaFEED","betaDENS",
                                     "zSEX","zFEED","zDENS",
                                     "sigma_year_a","u_cohort"),
                        enableWAIC = TRUE)

## --- CRITICAL: BLOCK SAMPLING ---
# This is the most important part for mixing birth/death.
cat("Blocking tB and tstar nodes...\n")
for (i in 1:nind) {
  conf_A$removeSamplers(paste0('tB[', i, ']'))
  conf_A$removeSamplers(paste0('tstar[', i, ']'))
  conf_A$addSampler(target = c(paste0('tB[', i, ']'), paste0('tstar[', i, ']')),
                    type = "RW_block")
}

# Add Slice sampling for Gompertz parameters
conf_A$removeSamplers(c("a","b"))
conf_A$addSampler(target = c("a","b"), type = "AF_slice")

# Reversible Jump for variable selection
configureRJ(conf = conf_A,
            targetNodes    = c("betaSEX","betaFEED","betaDENS"),
            indicatorNodes = c("zSEX","zFEED","zDENS"),
            control = list(mean = 0, scale = 0.5))

mcmc_A  <- buildMCMC(conf_A)
cMCMC_A <- compileNimble(mcmc_A, project = model_A)

# --- 6. RUN MCMC ---

run_A <- runMCMC(cMCMC_A,
                 niter = 200000, nburnin = 50000, nchains = 2,
                 inits = list(init_A(), init_A()),
                 setSeed = TRUE,
                 samplesAsCodaMCMC = TRUE,
                 summary = TRUE,
                 WAIC = TRUE,
                 progressBar = TRUE)



# --- 7. RESULTS ---

MCMCsummary(run_A$samples, params = c("a","b","mean.p","sigma_year_a","zSEX","zFEED","zDENS"))