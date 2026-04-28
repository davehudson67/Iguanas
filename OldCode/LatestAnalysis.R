## ============================================================
## COMPLETE IGUANA SURVIVAL MODEL
## - Left-truncated Gompertz (Island-specific start)
## - RJMC Variable Selection: Sex, Feeding, Species, Density
## - Time-Varying Density (Cohort-specific at entry)
## ============================================================

library(data.table)
library(tidyverse)
library(nimble)
library(mcmcplots)

# Set working directory and load data
setwd("~/Iguanas")
load("igs_AllIslands_CleanDH_160226_obsStart.RData")
hab_area <- fread("Data/Iguana_islands_habitat.csv") %>% select(islands, Useable)

## Load custom distributions and functions
source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

## ------------------------------------------------------------
## 1) Build individual-level covariates
## ------------------------------------------------------------
id_df <- AllData %>% distinct(ID, .keep_all = TRUE)

# SEX: 0/1
sex <- as.numeric(as.factor(id_df$sex)) - 1
# FEEDING: 0 = none, 1 = moderate/high
feeding <- ifelse(id_df$feeding == "none", 0, 1)
# SPECIES: 0 = figginsi, 1 = inornata
species <- ifelse(id_df$Species == "inornata", 1, 0)
# ISLAND index: 1..n_island
island_idx <- as.numeric(as.factor(id_df$Island))
# Island-specific start time
IslandStart <- id_df$IslandStart

## ------------------------------------------------------------
## 2) Time-varying density calculation (MNA)
## ------------------------------------------------------------
# Calculate Minimum Number Alive (MNA) per Island per Year
mna_df <- as.data.frame(CH) %>%
  mutate(ID = 1:n()) %>%
  pivot_longer(-ID, names_to = "Year_Idx", values_to = "Detected") %>%
  mutate(Year_Idx = as.numeric(str_remove(Year_Idx, "V"))) %>%
  left_join(id_df %>% select(ID, Island), by = "ID") %>%
  group_by(ID, Island) %>%
  filter(any(Detected == 1)) %>%
  summarise(First = min(Year_Idx[Detected == 1]),
            Last = max(Year_Idx[Detected == 1]), .groups = "drop") %>%
  rowwise() %>%
  do(data.frame(ID = .$ID, Island = .$Island, Year_Idx = .$First:.$Last)) %>%
  group_by(Island, Year_Idx) %>%
  summarise(MNA = n(), .groups = "drop") %>%
  left_join(hab_area, by = c("Island" = "islands")) %>%
  mutate(Density = MNA / Useable)

# Create strictly numeric Density Matrix [Island x Year]
dens_matrix_numeric <- mna_df %>%
  group_by(Island) %>%
  mutate(dens_z = as.numeric(scale(log(Density)))) %>%
  ungroup() %>%
  select(Island, Year_Idx, dens_z) %>%
  pivot_wider(names_from = Year_Idx, values_from = dens_z) %>%
  arrange(match(Island, levels(as.factor(id_df$Island)))) %>%
  select(-Island) %>% # Remove character column
  as.matrix()

dens_matrix_numeric[is.na(dens_matrix_numeric)] <- 0 # Fill gaps with mean


# 1. Get the dimensions right
n_islands <- length(levels(as.factor(id_df$Island)))
tMax <- ncol(CH)

# 2. Create a blank matrix filled with 0 (the mean for Z-scores)
# Rows = Islands, Cols = Years
dens_matrix_fixed <- matrix(0, nrow = n_islands, ncol = tMax)

# 3. Fill it with your calculated MNA data
# Ensure Year_Idx is within 1:tMax
mna_clean <- mna_df %>%
  group_by(Island) %>%
  mutate(dens_z = as.numeric(scale(log(Density)))) %>%
  ungroup() %>%
  filter(Year_Idx >= 1, Year_Idx <= tMax) %>%
  mutate(island_idx = as.numeric(as.factor(Island)))

for(i in 1:nrow(mna_clean)) {
  dens_matrix_fixed[mna_clean$island_idx[i], mna_clean$Year_Idx[i]] <- mna_clean$dens_z[i]
}

# 4. Final check: Should be [n_islands x tMax] with no NAs
stopifnot(nrow(dens_matrix_fixed) == n_islands)
stopifnot(ncol(dens_matrix_fixed) == tMax)
stopifnot(!any(is.na(dens_matrix_fixed)))


## ------------------------------------------------------------
## 3) NIMBLE Model Code
## ------------------------------------------------------------
code <- nimbleCode({
  for (i in 1:nind) {
    # Latent birth time
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    
    # Left truncation
    L[i] <- max(0, IslandStart[i] - tB[i])
    
    # Survival
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    
    # Censoring
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    # --- FIXED INDEXING ---
    # Use trunc() and clamp strictly to matrix dimensions
    # This prevents the "Out of Bounds" warning
    entry_yr[i] <- max(1, min(tMax, trunc(tB[i])))
    idx_dens[i] <- dens_matrix[ island[i], entry_yr[i] ]
    
    # Linear Predictors
    log(amult[i]) <- log(a) +
      betaSEX[1]     * sex[i]     * zSEX[1] +
      betaFEED[1]    * feeding[i] * zFEED[1] +
      betaSPEC[1]    * species[i] * zSPEC[1] +
      betaDENS[1]    * idx_dens[i] * zDENS[1]
    
    log(bmult[i]) <- log(b) +
      betaSEX[2]     * sex[i]     * zSEX[2] +
      betaFEED[2]    * feeding[i] * zFEED[2] +
      betaSPEC[2]    * species[i] * zSPEC[2] +
      betaDENS[2]    * idx_dens[i] * zDENS[2]
    
    # Observation Model
    nm[i] <- max(max(ceiling(tB[i]), IslandStart[i]), 1)
    nM[i] <- min(floor(tD[i]), tMax) - nm[i] + 1
    nMpos[i] <- max(nM[i], 0)
    
    pd[i] <- exp(y[i] * log(mean.p + 1e-10) + (nMpos[i] - y[i]) * log(1 - mean.p + 1e-10))
    dind[i] ~ dbern(pd[i])
  }
  
  # Priors
  for (k in 1:2) {
    betaSEX[k] ~ dnorm(0, sd = 1.5); betaFEED[k] ~ dnorm(0, sd = 1.5)
    betaSPEC[k] ~ dnorm(0, sd = 1.5); betaDENS[k] ~ dnorm(0, sd = 1.5)
    zSEX[k] ~ dbern(0.5); zFEED[k] ~ dbern(0.5)
    zSPEC[k] ~ dbern(0.5); zDENS[k] ~ dbern(0.5)
  }
  a ~ dexp(1); b ~ dexp(1); mean.p ~ dunif(0, 1)
})

## ------------------------------------------------------------
## 4) Setup: Constants, Data, and Inits
## ------------------------------------------------------------
y <- rowSums(CH)
tMax <- ncol(CH)

# Sanitize Data (Must be strictly numeric arrays/vectors)
data_list <- list(
  y           = as.numeric(unname(y)),
  cintB       = as.matrix(unname(cintB)),
  cintD       = as.matrix(unname(cintD)),
  censoredD   = as.numeric(unname(censoredD)),
  dind        = rep(1, nind),
  dens_matrix = dens_matrix_fixed # The pre-allocated numeric matrix
)

# Clean constants
consts <- list(
  nind        = nind,
  tMax        = tMax,
  sex         = as.numeric(unname(sex)),
  feeding     = as.numeric(unname(feeding)),
  species     = as.numeric(unname(species)),
  island      = as.numeric(unname(island_idx)),
  IslandStart = as.numeric(unname(IslandStart))
)

# Robust Initial Values
last_occasion <- apply(CH, 1, function(r) {
  w <- which(r > 0); if (length(w) == 0) 0L else max(w)
})

initFn <- function(cintB, cintD, censoredD, IslandStart, last_occasion) {
  n <- nrow(cintB)
  tBinit <- runif(n, cintB[,1], cintB[,2])
  min_death <- pmax(tBinit + pmax(0, IslandStart - tBinit) + 1e-6, last_occasion + 1e-6)
  tDinit <- numeric(n)
  id_right <- which(censoredD == 2); id_int <- which(censoredD != 2)
  if (length(id_int)) tDinit[id_int] <- runif(length(id_int), pmax(cintD[id_int,1], min_death[id_int]), cintD[id_int,2])
  if (length(id_right)) tDinit[id_right] <- pmax(cintD[id_right,2] + rexp(length(id_right), 5), min_death[id_right])
  
  list(
    tB = tBinit, tD = tDinit, tstar = tDinit - tBinit,
    a = rexp(1, 10), b = rexp(1, 10), mean.p = runif(1, 0.2, 0.8),
    betaSEX = rnorm(2, 0, 0.1), betaFEED = rnorm(2, 0, 0.1),
    betaSPEC = rnorm(2, 0, 0.1), betaDENS = rnorm(2, 0, 0.1),
    zSEX = c(0,0), zFEED = c(0,0), zSPEC = c(0,0), zDENS = c(0,0)
  )
}

inits_list <- list(initFn(cintB, cintD, censoredD, IslandStart, last_occasion),
                   initFn(cintB, cintD, censoredD, IslandStart, last_occasion))

## ------------------------------------------------------------
## 5) Compilation and Execution
## ------------------------------------------------------------
model <- nimbleModel(code, constants = consts, data = data_list, inits = inits_list)
cModel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c("a", "b", "mean.p", "betaSEX", "betaFEED", 
                                          "betaSPEC", "betaDENS", "zSEX", "zFEED", 
                                          "zSPEC", "zDENS"))

# Specialized samplers
conf$removeSamplers(c("a", "b"))
conf$addSampler(target = c("a", "b"), type = "AF_slice")

configureRJ(conf = conf,
            targetNodes = c("betaSEX", "betaFEED", "betaSPEC", "betaDENS"),
            indicatorNodes = c("zSEX", "zFEED", "zSPEC", "zDENS"),
            control = list(mean = 0, scale = 0.5))

mcmc <- buildMCMC(conf)
cMCMC <- compileNimble(mcmc, project = model)

# Run
run <- runMCMC(cMCMC, niter = 150000, nburnin = 20000, nchains = 2,
               samplesAsCodaMCMC = TRUE, summary = TRUE, progressBar = TRUE)

saveRDS(run, "iguana_final_results.rds")


run$summary
