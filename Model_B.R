## ============================================================
## MODEL B: Gompertz LB + Density-at-entry + Cohort RW1 on log(a) AND log(b)
## - species removed
## - NO persistent island RE (to avoid competing with feeding)
## - Feeding interpreted as between-island contrast conditional on:
##   (i) density-at-entry and (ii) smooth cohort time drift (RW1)
## - RJ on: Sex, Feeding, Density
## ============================================================

library(data.table)
library(tidyverse)
library(nimble)

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
## 1b) OPTIONAL: REMOVE the 6 Bitter "none" individuals (if not already removed)
##     - You said you removed them: keep this OFF if already done upstream.
##     - If you *do* use it, you must know which IDs they are.
## ============================================================
# bitter_none_ids <- AllData %>%
#   filter(Island == "Bitter") %>%
#   group_by(ID) %>%
#   summarise(any_none = any(as.character(feeding) == "none"), .groups="drop") %>%
#   filter(any_none) %>%
#   pull(ID)
# message("Dropping Bitter mixed-regime IDs: ", paste(bitter_none_ids, collapse=", "))
# AllData <- AllData %>% filter(!(ID %in% bitter_none_ids))
# CH    <- CH[!(rownames(CH) %in% bitter_none_ids), , drop=FALSE]  # only if CH rownames are IDs
# cintB <- cintB[!(rownames(cintB) %in% bitter_none_ids), , drop=FALSE]
# cintD <- cintD[!(rownames(cintD) %in% bitter_none_ids), , drop=FALSE]

## ============================================================
## 2) ONE-ROW-PER-ID covariates
## ============================================================
id_df <- AllData %>% distinct(ID, .keep_all = TRUE)
nind  <- nrow(id_df)
tMax  <- ncol(CH)

stopifnot(nrow(CH) == nind, nrow(cintB) == nind, nrow(cintD) == nind)

ordered_islands <- levels(as.factor(id_df$Island))
n_islands <- length(ordered_islands)

sex     <- as.numeric(as.factor(id_df$sex)) - 1
feeding <- ifelse(as.character(id_df$feeding) == "none", 0, 1)
island_idx <- as.numeric(factor(id_df$Island, levels = ordered_islands))

## IslandStart is already in OCCASION units, but varies by island
island_start_by_island <- id_df %>%
  distinct(Island, .keep_all = TRUE) %>%
  mutate(Island = as.character(Island)) %>%
  arrange(factor(Island, levels = ordered_islands)) %>%
  pull(IslandStart)

stopifnot(length(island_start_by_island) == n_islands)
IslandStart_ind <- island_start_by_island[island_idx]

## Censoring indicator
censoredD <- rep(2, nind)
censoredD[cintD[,1] > 0] <- 1

## Detections
y <- rowSums(CH)
dind <- rep(1, nind)

last_occasion <- apply(CH, 1, function(r) {
  w <- which(r > 0)
  if (length(w) == 0) 0L else max(w)
})

## ============================================================
## 3) TIME-VARYING DENSITY MATRIX by Island × YearIdx
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
## 4) SURVEY EFFORT MATRIX (Island × YearIdx)
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
## 6) MODEL B: RW1 on log(a) and RW1 on log(b)
## ============================================================
code_B <- nimbleCode({
  
  for (i in 1:nind) {
    
    ## Latent birth/entry time (occasion units)
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    
    ## Entry year index (1..tMax) for density lookup and RW1 indexing
    cohort_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    
    ## Left truncation lower bound in AGE units
    L[i] <- max(0, IslandStart_ind[i] - tB[i])
    
    ## Left-truncated Gompertz age-at-death
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    
    ## Death time
    tD[i] <- tB[i] + tstar[i]
    
    ## Interval / right censoring constraint
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    ## Density-at-entry (island x entry-year)
    dens_entry[i] <- dens_matrix[ island[i], cohort_idx[i] ]
    
    ## log(a_i): baseline + covariates + cohort RW1(a) (centered)
    log(amult[i]) <- log(a) +
      betaSEX[1]  * sex[i]        * zSEX[1] +
      betaFEED[1] * feeding[i]    * zFEED[1] +
      betaDENS[1] * dens_entry[i] * zDENS[1] +
      uA_center[ cohort_idx[i] ]
    
    ## log(b_i): baseline + covariates + cohort RW1(b) (centered)
    log(bmult[i]) <- log(b) +
      betaSEX[2]  * sex[i]        * zSEX[2] +
      betaFEED[2] * feeding[i]    * zFEED[2] +
      betaDENS[2] * dens_entry[i] * zDENS[2] +
      uB_center[ cohort_idx[i] ]
    
    ## Observation model: count of detections out of actual survey opportunities
    tB_year_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    island_start_year_idx[i] <- max(1, min(tMax, floor(IslandStart_ind[i]) + 1))
    
    nm_start[i] <- max(tB_year_idx[i], island_start_year_idx[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    
    ## effort-based opportunities in [nm_start, end_year]
    nMpos[i] <- cum_effort[ island[i], end_year[i] + 1 ] - cum_effort[ island[i], nm_start[i] ]
    
    ## ensure trials >= successes
    nMpos_safe[i] <- max(y[i], nMpos[i])
    
    ## ones-trick kernel
    pd_raw[i] <- exp(y[i] * log(mean.p + 1e-10) +
                       (nMpos_safe[i] - y[i]) * log(1 - mean.p + 1e-10))
    pd[i] <- min(0.999999, pd_raw[i])
    dind[i] ~ dbern(pd[i])
  }
  
  ## -------------------------
  ## RW1 on log(a): uA[1..tMax], centered
  ## -------------------------
  uA[1] ~ dnorm(0, sd = 2)
  for (t in 2:tMax) {
    uA[t] ~ dnorm(uA[t-1], sd = sigmaA)
  }
  uA_bar <- mean(uA[1:tMax])
  for (t in 1:tMax) {
    uA_center[t] <- uA[t] - uA_bar
  }
  sigmaA ~ T(dnorm(0, sd = 0.5), 0, )
  
  ## -------------------------
  ## RW1 on log(b): uB[1..tMax], centered
  ## -------------------------
  uB[1] ~ dnorm(0, sd = 2)
  for (t in 2:tMax) {
    uB[t] ~ dnorm(uB[t-1], sd = sigmaB)
  }
  uB_bar <- mean(uB[1:tMax])
  for (t in 1:tMax) {
    uB_center[t] <- uB[t] - uB_bar
  }
  sigmaB ~ T(dnorm(0, sd = 0.5), 0, )
  
  ## Fixed effects + RJ indicators (two equations: a and b)
  for (k in 1:2) {
    betaSEX[k]  ~ dnorm(0, sd = 1.5)
    betaFEED[k] ~ dnorm(0, sd = 1.5)
    betaDENS[k] ~ dnorm(0, sd = 1.5)
    
    zSEX[k]  ~ dbern(0.5)
    zFEED[k] ~ dbern(0.5)
    zDENS[k] ~ dbern(0.5)
  }
  
  ## Baseline Gompertz parameters + detection
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
})

## ---- constants + data ----
consts_B <- list(
  nind = nind,
  tMax = tMax,
  sex = sex,
  feeding = feeding,
  island = island_idx,
  IslandStart_ind = IslandStart_ind
)

data_B <- list(
  y = as.numeric(y),
  cintB = as.matrix(cintB),
  cintD = as.matrix(cintD),
  censoredD = as.numeric(censoredD),
  dind = rep(1, nind),
  dens_matrix = dens_matrix,
  cum_effort = cum_effort
)

## ---- inits ----
init_B <- function() {
  
  tiny <- 1e-6
  tBinit <- runif(nind, cintB[,1], cintB[,2])
  
  min_death <- pmax(tBinit, IslandStart_ind) + tiny
  min_death <- pmax(min_death, last_occasion + tiny)
  
  tDinit <- numeric(nind)
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
    tB = tBinit,
    tD = tDinit,
    tstar = tDinit - tBinit,
    
    ## helps avoid NA during build + used for RW1 indexing
    cohort_idx = pmax(1, pmin(tMax, floor(tBinit) + 1)),
    
    ## RW1 terms
    uA = rnorm(tMax, 0, 0.05),
    uB = rnorm(tMax, 0, 0.05),
    sigmaA = abs(rnorm(1, 0, 0.2)),
    sigmaB = abs(rnorm(1, 0, 0.2)),
    
    a = rexp(1, 10),
    b = rexp(1, 10),
    mean.p = runif(1, 0.2, 0.8),
    
    betaSEX  = rnorm(2, 0, 0.1),
    betaFEED = rnorm(2, 0, 0.1),
    betaDENS = rnorm(2, 0, 0.1),
    
    zSEX  = c(0,0),
    zFEED = c(0,0),
    zDENS = c(0,0)
  )
}

inits_list_B <- list(init_B(), init_B())

## ---- build + MCMC ----
model_B <- nimbleModel(code_B, constants = consts_B, data = data_B, inits = inits_list_B[[1]])
cModel_B <- compileNimble(model_B)

conf_B <- configureMCMC(model_B,
                        monitors = c("a","b","mean.p",
                                     "betaSEX","betaFEED","betaDENS",
                                     "zSEX","zFEED","zDENS",
                                     "sigmaA","sigmaB",
                                     "uA_center","uB_center"),
                        enableWAIC = TRUE)

## Better samplers for positive a/b
conf_B$removeSamplers(c("a","b"))
conf_B$addSampler(target = c("a","b"), type = "AF_slice")

## RJ configuration
configureRJ(conf = conf_B,
            targetNodes    = c("betaSEX","betaFEED","betaDENS"),
            indicatorNodes = c("zSEX","zFEED","zDENS"),
            control = list(mean = 0, scale = 0.5))

mcmc_B  <- buildMCMC(conf_B)
cMCMC_B <- compileNimble(mcmc_B, project = model_B)

run_B <- runMCMC(cMCMC_B,
                 niter = 250000, nburnin = 75000, nchains = 2,
                 inits = inits_list_B,
                 setSeed = TRUE,
                 samplesAsCodaMCMC = TRUE,
                 summary = TRUE,
                 WAIC = TRUE,
                 progressBar = TRUE)

saveRDS(run_B, "iguana_modelB_RW1cohort_on_a_and_b_density_feeding_noIslandRE.rds")

run_B <- readRDS("iguana_modelB_RW1cohort_on_a_and_b_density_feeding_noIslandRE.rds")
run_B$summary
library(mcmcplots)
mcmcplot(run_B$samples)
