## ============================================================
## IGUANA SURVIVAL MODEL (CLEAN + CORRECTED)
## ============================================================

library(data.table)
library(tidyverse)
library(nimble)

setwd("~/Iguanas")

## ---- Load data ----
load("igs_AllIslands_CleanDH_160226_obsStart.RData")
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
species <- ifelse(as.character(id_df$Species) == "inornata", 1, 0)  
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

last_occasion <- apply(CH, 1, function(r) { w <- which(r > 0); if (length(w) == 0) 0L else max(w) })

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

# Modern tidyverse approach replacing deprecated do()
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
## 6) NIMBLE MODEL
## ============================================================
code <- nimbleCode({
  
  for (i in 1:nind) {
    
    ## Latent birth/entry time in OCCASION units (continuous)
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    
    ## Cohort index for RW1: floor(tB)+1 clamped to 1..n_cohort
    cohort_idx[i] <- max(1, min(n_cohort, floor(tB[i]) + 1))
    
    ## Island-specific start (occasion units)
    L[i] <- max(0, IslandStart_ind[i] - tB[i])
    
    ## Left-truncated Gompertz survival time (age since birth)
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    
    ## Death time (occasion units)
    tD[i] <- tB[i] + tstar[i]
    
    ## Interval / right censoring constraint
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    ## Density-at-entry lookup: use the cohort entry year (floor(tB)+1)
    dens_entry[i] <- dens_matrix[ island[i], cohort_idx[i] ]
    
    ## Linear predictors: log(a_i) and log(b_i)
    ## REMOVED: u_island_a and betaSPEC_FEED to prevent confounding
    log(amult[i]) <- log(a) +
      betaSEX[1]       * sex[i]            * zSEX[1] +
      betaFEED[1]      * feeding[i]        * zFEED[1] +
      betaSPEC[1]      * species[i]        * zSPEC[1] +
      betaDENS[1]      * dens_entry[i]     * zDENS[1] +
      u_cohort_center[ cohort_idx[i] ]
    
    log(bmult[i]) <- log(b) +
      betaSEX[2]       * sex[i]            * zSEX[2] +
      betaFEED[2]      * feeding[i]        * zFEED[2] +
      betaSPEC[2]      * species[i]        * zSPEC[2] +
      betaDENS[2]      * dens_entry[i]     * zDENS[2]
    
    ## Observation model: count of detections y[i] out of actual survey opportunities
    tB_year_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    island_start_year_idx[i] <- max(1, min(tMax, floor(IslandStart_ind[i]) + 1))
    
    nm_start[i] <- max(tB_year_idx[i], island_start_year_idx[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    
    nMpos[i] <- cum_effort[ island[i], end_year[i] + 1 ] - cum_effort[ island[i], nm_start[i] ]
    
    ## FIX 1: Ensure survey opportunities are at least the number of actual detections
    nMpos_safe[i] <- max(y[i], nMpos[i])
    
    ## ones-trick likelihood kernel
    pd_raw[i] <- exp(y[i] * log(mean.p + 1e-10) + (nMpos_safe[i] - y[i]) * log(1 - mean.p + 1e-10))
    
    ## FIX 2: Cap probability at 0.999999 to prevent dbern() NaNs
    pd[i] <- min(0.999999, pd_raw[i])
    dind[i] ~ dbern(pd[i])
  }
  
  ## Cohort RW1 on log(a) random effect (length n_cohort)
  u_cohort[1] ~ dnorm(0, sd = 10)
  for (t in 2:n_cohort) {
    u_cohort[t] ~ dnorm(u_cohort[t-1], sd = sigma_year_a)
  }
  u_bar <- mean(u_cohort[1:n_cohort])
  for (t in 1:n_cohort) {
    u_cohort_center[t] <- u_cohort[t] - u_bar
  }
  sigma_year_a ~ dunif(0, 5)
  
  ## Fixed effects + RJ indicators (REMOVED: betaSPEC_FEED and zSPEC_FEED)
  for (k in 1:2) {
    betaSEX[k] ~ dnorm(0, sd = 1.5)
    betaFEED[k] ~ dnorm(0, sd = 1.5)
    betaSPEC[k] ~ dnorm(0, sd = 1.5)
    betaDENS[k] ~ dnorm(0, sd = 1.5)
    
    zSEX[k] ~ dbern(0.5)
    zFEED[k] ~ dbern(0.5)
    zSPEC[k] ~ dbern(0.5)
    zDENS[k] ~ dbern(0.5)
  }
  
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
})

## ============================================================
## 7) CONSTANTS + DATA
## ============================================================
consts <- list(
  nind = nind,
  tMax = tMax,
  n_cohort = n_cohort,
  sex = sex,
  feeding = feeding,
  species = species,
  island = island_idx,
  IslandStart_ind = IslandStart_ind
)

data_list <- list(
  y = as.numeric(y),
  cintB = as.matrix(cintB),
  cintD = as.matrix(cintD),
  censoredD = as.numeric(censoredD),
  dind = as.numeric(dind),
  dens_matrix = dens_matrix,
  cum_effort = cum_effort
)

## ============================================================
## 8) FAST INIT FUNCTION 
## ============================================================
initFn_fast <- function(cintB, cintD, censoredD,
                        IslandStart_ind, island_idx,
                        last_occasion,
                        n_cohort) {
  
  n <- nrow(cintB)
  tiny <- 1e-6
  
  tBinit <- runif(n, cintB[,1], cintB[,2])
  
  min_death <- pmax(tBinit, IslandStart_ind) + tiny
  min_death <- pmax(min_death, last_occasion + tiny)
  
  tDinit <- numeric(n)
  id_right <- which(censoredD == 2)
  id_int   <- which(censoredD != 2)
  
  if (length(id_int)) {
    lo <- pmax(cintD[id_int,1], min_death[id_int])
    hi <- cintD[id_int,2]
    bad <- which(lo > hi | is.na(lo) | is.na(hi))
    if (length(bad)) {
      stop("Init impossible for some interval-censored individuals. Check time scales/cintD.")
    }
    tDinit[id_int] <- runif(length(id_int), lo, hi)
  }
  
  if (length(id_right)) {
    tDprop <- cintD[id_right,2] + rexp(length(id_right), rate = 5)
    tDinit[id_right] <- pmax(tDprop, min_death[id_right])
  }
  
  list(
    tB = tBinit,
    tD = tDinit,
    tstar = tDinit - tBinit,
    
    ## Explicitly initialize cohort_idx to prevent NA indexing during model build
    cohort_idx = pmax(1, pmin(n_cohort, floor(tBinit) + 1)), 
    
    a = rexp(1, rate = 10),
    b = rexp(1, rate = 10),
    mean.p = runif(1, 0.2, 0.8),
    
    betaSEX = rnorm(2, 0, 0.1),
    betaFEED = rnorm(2, 0, 0.1),
    betaSPEC = rnorm(2, 0, 0.1),
    betaDENS = rnorm(2, 0, 0.1),
    
    zSEX = c(0,0),
    zFEED = c(0,0),
    zSPEC = c(0,0),
    zDENS = c(0,0),
    
    u_cohort = rnorm(n_cohort, 0, 0.1),
    sigma_year_a = runif(1, 0.1, 1)
  )
}

inits_list <- list(
  initFn_fast(cintB, cintD, censoredD, IslandStart_ind, island_idx, last_occasion, n_cohort),
  initFn_fast(cintB, cintD, censoredD, IslandStart_ind, island_idx, last_occasion, n_cohort)
)

## ============================================================
## 9) BUILD / COMPILE MODEL + MCMC
## ============================================================
## Pass inits directly to prevent NA evaluation warnings during build
model <- nimbleModel(code, constants = consts, data = data_list, inits = inits_list[[1]])
cModel <- compileNimble(model)

conf <- configureMCMC(model,
                      monitors = c("a","b","mean.p",
                                   "betaSEX","betaFEED","betaSPEC","betaDENS",
                                   "zSEX","zFEED","zSPEC","zDENS",
                                   "sigma_year_a","u_cohort_center"),
                      enableWAIC = TRUE)

## Better samplers for positive a/b
conf$removeSamplers(c("a","b"))
conf$addSampler(target = c("a","b"), type = "AF_slice")

## RJ configuration (REMOVED: betaSPEC_FEED and zSPEC_FEED)
configureRJ(conf = conf,
            targetNodes = c("betaSEX","betaFEED","betaSPEC","betaDENS"),
            indicatorNodes = c("zSEX","zFEED","zSPEC","zDENS"),
            control = list(mean = 0, scale = 0.5))

mcmc  <- buildMCMC(conf)
cMCMC <- compileNimble(mcmc, project = model)

## ============================================================
## 10) RUN
## ============================================================
niter_val   <- 250000
nburnin_val <- 75000
nchains_val <- 2

run <- runMCMC(cMCMC,
               niter = niter_val,
               nburnin = nburnin_val,
               nchains = nchains_val,
               inits = inits_list,
               setSeed = TRUE,
               samplesAsCodaMCMC = TRUE,
               summary = TRUE,
               WAIC = TRUE,
               progressBar = TRUE)

saveRDS(run, "iguana_hierarchical_density_noIsland_samples.rds")
run <- readRDS("iguana_hierarchical_density_noIsland_samples.rds")

library(mcmcplots)
library(MCMCvis)
mcmcplots::mcmcplot(run$samples)
run$summary
MCMCtrace(run$samples, params = "betaFEED")
plot(run$samples, pars = "betaFEED")

#------------------------------------------------------------------------------#
## ============================================================
## 6) NIMBLE MODEL (No Island RE, No Cohort RE, No Interaction)
## ============================================================
code <- nimbleCode({
  
  for (i in 1:nind) {
    
    ## Latent birth/entry time in OCCASION units (continuous)
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    
    ## Cohort index: floor(tB)+1 clamped to 1..n_cohort (Used for density lookup)
    cohort_idx[i] <- max(1, min(n_cohort, floor(tB[i]) + 1))
    
    ## Island-specific start (occasion units)
    L[i] <- max(0, IslandStart_ind[i] - tB[i])
    
    ## Left-truncated Gompertz survival time (age since birth)
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    
    ## Death time (occasion units)
    tD[i] <- tB[i] + tstar[i]
    
    ## Interval / right censoring constraint
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    ## Density-at-entry lookup: use the cohort entry year
    dens_entry[i] <- dens_matrix[ island[i], cohort_idx[i] ]
    
    ## Linear predictors: log(a_i) and log(b_i)
    ## REMOVED: u_cohort_center (Cohort RE) and u_island_a (Island RE)
    log(amult[i]) <- log(a) +
      betaSEX[1]       * sex[i]            * zSEX[1] +
      betaFEED[1]      * feeding[i]        * zFEED[1] +
      betaSPEC[1]      * species[i]        * zSPEC[1] +
      betaDENS[1]      * dens_entry[i]     * zDENS[1]
    
    log(bmult[i]) <- log(b) +
      betaSEX[2]       * sex[i]            * zSEX[2] +
      betaFEED[2]      * feeding[i]        * zFEED[2] +
      betaSPEC[2]      * species[i]        * zSPEC[2] +
      betaDENS[2]      * dens_entry[i]     * zDENS[2]
    
    ## Observation model: count of detections y[i] out of actual survey opportunities
    tB_year_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    island_start_year_idx[i] <- max(1, min(tMax, floor(IslandStart_ind[i]) + 1))
    
    nm_start[i] <- max(tB_year_idx[i], island_start_year_idx[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    
    nMpos[i] <- cum_effort[ island[i], end_year[i] + 1 ] - cum_effort[ island[i], nm_start[i] ]
    
    ## FIX 1: Ensure survey opportunities are at least the number of actual detections
    nMpos_safe[i] <- max(y[i], nMpos[i])
    
    ## ones-trick likelihood kernel
    pd_raw[i] <- exp(y[i] * log(mean.p + 1e-10) + (nMpos_safe[i] - y[i]) * log(1 - mean.p + 1e-10))
    
    ## FIX 2: Cap probability at 0.999999 to prevent dbern() NaNs
    pd[i] <- min(0.999999, pd_raw[i])
    dind[i] ~ dbern(pd[i])
  }
  
  ## Fixed effects + RJ indicators
  for (k in 1:2) {
    betaSEX[k] ~ dnorm(0, sd = 1.5)
    betaFEED[k] ~ dnorm(0, sd = 1.5)
    betaSPEC[k] ~ dnorm(0, sd = 1.5)
    betaDENS[k] ~ dnorm(0, sd = 1.5)
    
    zSEX[k] ~ dbern(0.5)
    zFEED[k] ~ dbern(0.5)
    zSPEC[k] ~ dbern(0.5)
    zDENS[k] ~ dbern(0.5)
  }
  
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
})

## ============================================================
## 7) CONSTANTS + DATA
## ============================================================
consts <- list(
  nind = nind,
  tMax = tMax,
  n_cohort = n_cohort,
  sex = sex,
  feeding = feeding,
  species = species,
  island = island_idx,
  IslandStart_ind = IslandStart_ind
)

data_list <- list(
  y = as.numeric(y),
  cintB = as.matrix(cintB),
  cintD = as.matrix(cintD),
  censoredD = as.numeric(censoredD),
  dind = as.numeric(dind),
  dens_matrix = dens_matrix,
  cum_effort = cum_effort
)

## ============================================================
## 8) FAST INIT FUNCTION 
## ============================================================
initFn_fast <- function(cintB, cintD, censoredD,
                        IslandStart_ind, island_idx,
                        last_occasion,
                        n_cohort) {
  
  n <- nrow(cintB)
  tiny <- 1e-6
  
  tBinit <- runif(n, cintB[,1], cintB[,2])
  
  min_death <- pmax(tBinit, IslandStart_ind) + tiny
  min_death <- pmax(min_death, last_occasion + tiny)
  
  tDinit <- numeric(n)
  id_right <- which(censoredD == 2)
  id_int   <- which(censoredD != 2)
  
  if (length(id_int)) {
    lo <- pmax(cintD[id_int,1], min_death[id_int])
    hi <- cintD[id_int,2]
    bad <- which(lo > hi | is.na(lo) | is.na(hi))
    if (length(bad)) {
      stop("Init impossible for some interval-censored individuals. Check time scales/cintD.")
    }
    tDinit[id_int] <- runif(length(id_int), lo, hi)
  }
  
  if (length(id_right)) {
    tDprop <- cintD[id_right,2] + rexp(length(id_right), rate = 5)
    tDinit[id_right] <- pmax(tDprop, min_death[id_right])
  }
  
  list(
    tB = tBinit,
    tD = tDinit,
    tstar = tDinit - tBinit,
    
    ## Explicitly initialize cohort_idx to prevent NA indexing during model build
    cohort_idx = pmax(1, pmin(n_cohort, floor(tBinit) + 1)), 
    
    a = rexp(1, rate = 10),
    b = rexp(1, rate = 10),
    mean.p = runif(1, 0.2, 0.8),
    
    betaSEX = rnorm(2, 0, 0.1),
    betaFEED = rnorm(2, 0, 0.1),
    betaSPEC = rnorm(2, 0, 0.1),
    betaDENS = rnorm(2, 0, 0.1),
    
    zSEX = c(0,0),
    zFEED = c(0,0),
    zSPEC = c(0,0),
    zDENS = c(0,0)
    
    ## REMOVED: u_cohort and sigma_year_a
  )
}

inits_list <- list(
  initFn_fast(cintB, cintD, censoredD, IslandStart_ind, island_idx, last_occasion, n_cohort),
  initFn_fast(cintB, cintD, censoredD, IslandStart_ind, island_idx, last_occasion, n_cohort)
)

## ============================================================
## 9) BUILD / COMPILE MODEL + MCMC
## ============================================================
model <- nimbleModel(code, constants = consts, data = data_list, inits = inits_list[[1]])
cModel <- compileNimble(model)

conf <- configureMCMC(model,
                      ## REMOVED: sigma_year_a and u_cohort_center
                      monitors = c("a","b","mean.p",
                                   "betaSEX","betaFEED","betaSPEC","betaDENS",
                                   "zSEX","zFEED","zSPEC","zDENS"),
                      enableWAIC = TRUE)

conf$removeSamplers(c("a","b"))
conf$addSampler(target = c("a","b"), type = "AF_slice")

configureRJ(conf = conf,
            targetNodes = c("betaSEX","betaFEED","betaSPEC","betaDENS"),
            indicatorNodes = c("zSEX","zFEED","zSPEC","zDENS"),
            control = list(mean = 0, scale = 0.5))

mcmc  <- buildMCMC(conf)
cMCMC <- compileNimble(mcmc, project = model)

## ============================================================
## 10) RUN
## ============================================================
niter_val   <- 250000
nburnin_val <- 75000
nchains_val <- 2

run <- runMCMC(cMCMC,
               niter = niter_val,
               nburnin = nburnin_val,
               nchains = nchains_val,
               inits = inits_list,
               setSeed = TRUE,
               samplesAsCodaMCMC = TRUE,
               summary = TRUE,
               WAIC = TRUE,
               progressBar = TRUE)

## Renamed output file to reflect the dropped cohort effect
saveRDS(run, "iguana_hierarchical_density_noIsland_noCohort_samples.rds")

run$summary

#-------------------------------------------------------------------------------

## ============================================================
## SPECIES-SPECIFIC IGUANA SURVIVAL MODELS
## (No Island RE, No Cohort RE, No Species Covariate)
## ============================================================

library(data.table)
library(tidyverse)
library(nimble)

setwd("~/Iguanas")

## ---- Load data ----
load("igs_AllIslands_CleanDH_160226_obsStart.RData")
hab_area <- fread("Data/Iguana_islands_habitat.csv") %>% select(islands, Useable)

## ---- Custom functions/distributions ----
source("ModelComparison_FUNCTIONS.R")     
source("Distributions/Dist_GompertzLB.R") 
source("Distributions/Dist_Gompertz.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

## ============================================================
## 1) DATA CORRECTION
## ============================================================
problem_ids <- AllData %>%
  filter(Island == "FFRC", Species == "figginsi") %>%
  pull(ID)

if (length(problem_ids) > 0) {
  message("Correcting Species for FFRC IDs: ", paste(problem_ids, collapse = ", "))
  AllData <- AllData %>%
    mutate(Species = ifelse(ID %in% problem_ids & Island == "FFRC", "inornata", Species))
}

## Master ID list to track original row indices for CH, cintB, cintD
id_df_master <- AllData %>% 
  distinct(ID, .keep_all = TRUE) %>%
  mutate(master_row_idx = 1:n())

tMax <- ncol(CH)

## ============================================================
## 2) DEFINE THE NIMBLE MODEL (Species removed)
## ============================================================
code <- nimbleCode({
  for (i in 1:nind) {
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    cohort_idx[i] <- max(1, min(n_cohort, floor(tB[i]) + 1))
    L[i] <- max(0, IslandStart_ind[i] - tB[i])
    
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    dens_entry[i] <- dens_matrix[ island[i], cohort_idx[i] ]
    
    ## Linear predictors (Species removed)
    log(amult[i]) <- log(a) +
      betaSEX[1]  * sex[i]        * zSEX[1] +
      betaFEED[1] * feeding[i]    * zFEED[1] +
      betaDENS[1] * dens_entry[i] * zDENS[1]
    
    log(bmult[i]) <- log(b) +
      betaSEX[2]  * sex[i]        * zSEX[2] +
      betaFEED[2] * feeding[i]    * zFEED[2] +
      betaDENS[2] * dens_entry[i] * zDENS[2]
    
    tB_year_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    island_start_year_idx[i] <- max(1, min(tMax, floor(IslandStart_ind[i]) + 1))
    
    nm_start[i] <- max(tB_year_idx[i], island_start_year_idx[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    
    nMpos[i] <- cum_effort[ island[i], end_year[i] + 1 ] - cum_effort[ island[i], nm_start[i] ]
    nMpos_safe[i] <- max(y[i], nMpos[i])
    
    pd_raw[i] <- exp(y[i] * log(mean.p + 1e-10) + (nMpos_safe[i] - y[i]) * log(1 - mean.p + 1e-10))
    pd[i] <- min(0.999999, pd_raw[i])
    dind[i] ~ dbern(pd[i])
  }
  
  for (k in 1:2) {
    betaSEX[k] ~ dnorm(0, sd = 1.5)
    betaFEED[k] ~ dnorm(0, sd = 1.5)
    betaDENS[k] ~ dnorm(0, sd = 1.5)
    
    zSEX[k] ~ dbern(0.5)
    zFEED[k] ~ dbern(0.5)
    zDENS[k] ~ dbern(0.5)
  }
  
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
})

## ============================================================
## 3) INIT FUNCTION (Species removed)
## ============================================================
initFn_fast <- function(cintB, cintD, censoredD, IslandStart_ind, last_occasion, n_cohort) {
  n <- nrow(cintB)
  tiny <- 1e-6
  tBinit <- runif(n, cintB[,1], cintB[,2])
  
  min_death <- pmax(tBinit, IslandStart_ind) + tiny
  min_death <- pmax(min_death, last_occasion + tiny)
  
  tDinit <- numeric(n)
  id_right <- which(censoredD == 2)
  id_int   <- which(censoredD != 2)
  
  if (length(id_int)) {
    lo <- pmax(cintD[id_int,1], min_death[id_int])
    hi <- cintD[id_int,2]
    tDinit[id_int] <- runif(length(id_int), lo, hi)
  }
  if (length(id_right)) {
    tDprop <- cintD[id_right,2] + rexp(length(id_right), rate = 5)
    tDinit[id_right] <- pmax(tDprop, min_death[id_right])
  }
  
  list(
    tB = tBinit, tD = tDinit, tstar = tDinit - tBinit,
    cohort_idx = pmax(1, pmin(n_cohort, floor(tBinit) + 1)), 
    a = rexp(1, rate = 10), b = rexp(1, rate = 10), mean.p = runif(1, 0.2, 0.8),
    betaSEX = rnorm(2, 0, 0.1), betaFEED = rnorm(2, 0, 0.1), betaDENS = rnorm(2, 0, 0.1),
    zSEX = c(0,0), zFEED = c(0,0), zDENS = c(0,0)
  )
}

## ============================================================
## 4) LOOP OVER SPECIES
## ============================================================
target_species_list <- c("figginsi", "inornata")

for (sp in target_species_list) {
  message("\n==================================================")
  message("PREPPING AND RUNNING MODEL FOR: ", sp)
  message("==================================================")
  
  ## ---- Subset Data ----
  id_df <- id_df_master %>% filter(Species == sp)
  keep_idx <- id_df$master_row_idx
  
  nind <- nrow(id_df)
  CH_sub <- CH[keep_idx, , drop = FALSE]
  cintB_sub <- cintB[keep_idx, , drop = FALSE]
  cintD_sub <- cintD[keep_idx, , drop = FALSE]
  
  ordered_islands <- levels(as.factor(id_df$Island))
  n_islands <- length(ordered_islands)
  
  sex <- as.numeric(as.factor(id_df$sex)) - 1
  feeding <- ifelse(as.character(id_df$feeding) == "none", 0, 1)
  island_idx <- as.numeric(factor(id_df$Island, levels = ordered_islands))  
  
  island_start_by_island <- id_df %>%
    distinct(Island, .keep_all = TRUE) %>%
    arrange(factor(Island, levels = ordered_islands)) %>%
    pull(IslandStart)
  
  IslandStart_ind <- island_start_by_island[island_idx]
  
  censoredD <- rep(2, nind)
  censoredD[cintD_sub[,1] > 0] <- 1
  y <- rowSums(CH_sub)
  dind <- rep(1, nind)
  last_occasion <- apply(CH_sub, 1, function(r) { w <- which(r > 0); if (length(w) == 0) 0L else max(w) })
  
  ## ---- Density Matrix (Subset) ----
  area_vec <- hab_area %>%
    mutate(islands = as.character(islands)) %>%
    right_join(tibble(islands = ordered_islands), by = "islands") %>%
    pull(Useable)
  
  id_df <- id_df %>% mutate(row_idx = 1:nind)
  
  first_last <- as.data.frame(CH_sub) %>%
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
  
  ## ---- Effort Matrix (Subset) ----
  actual_surveys <- AllData %>%
    filter(Status == 1, Species == sp) %>%
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
  
  n_cohort <- tMax
  
  ## ---- Constants & Data ----
  consts <- list(
    nind = nind, tMax = tMax, n_cohort = n_cohort,
    sex = sex, feeding = feeding, island = island_idx, IslandStart_ind = IslandStart_ind
  )
  
  data_list <- list(
    y = as.numeric(y), cintB = as.matrix(cintB_sub), cintD = as.matrix(cintD_sub),
    censoredD = as.numeric(censoredD), dind = as.numeric(dind),
    dens_matrix = dens_matrix, cum_effort = cum_effort
  )
  
  inits_list <- list(
    initFn_fast(cintB_sub, cintD_sub, censoredD, IslandStart_ind, last_occasion, n_cohort),
    initFn_fast(cintB_sub, cintD_sub, censoredD, IslandStart_ind, last_occasion, n_cohort)
  )
  
  ## ---- Build & Run Model ----
  model <- nimbleModel(code, constants = consts, data = data_list, inits = inits_list[[1]])
  cModel <- compileNimble(model)
  
  conf <- configureMCMC(model,
                        monitors = c("a","b","mean.p", "betaSEX","betaFEED","betaDENS", "zSEX","zFEED","zDENS"),
                        enableWAIC = TRUE)
  
  conf$removeSamplers(c("a","b"))
  conf$addSampler(target = c("a","b"), type = "AF_slice")
  
  configureRJ(conf = conf,
              targetNodes = c("betaSEX","betaFEED","betaDENS"),
              indicatorNodes = c("zSEX","zFEED","zDENS"),
              control = list(mean = 0, scale = 0.5))
  
  mcmc  <- buildMCMC(conf)
  cMCMC <- compileNimble(mcmc, project = model)
  
  run <- runMCMC(cMCMC,
                 niter = 250000,
                 nburnin = 75000,
                 nchains = 2,
                 inits = inits_list,
                 setSeed = TRUE,
                 samplesAsCodaMCMC = TRUE,
                 summary = TRUE,
                 WAIC = TRUE,
                 progressBar = TRUE)
  
  ## Save with species name in the file
  save_name <- paste0("iguana_survival_", sp, "_samples.rds")
  saveRDS(run, save_name)
  message("Finished ", sp, "! Saved to ", save_name)
  
  ## Optional: Clear Nimble memory before next loop iteration
  rm(model, cModel, conf, mcmc, cMCMC, run)
  gc()
}


run_f <- readRDS("iguana_survival_figginsi_samples.rds")
run_f$summary
run_i <- readRDS("iguana_survival_inornata_samples.rds")
run_i$summary
