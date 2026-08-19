## ============================================================
## MODEL_PaceOfLife.R
##
## Gompertz pace-of-life model:
##   a: sex + cohort density + temporal cohort RW1 + island RE
##   b: sex + feeding + species + island RE
##
## Feeding is an island/block-level treatment.
## Missing historical cohort densities use island mean density.
## ============================================================

source("DataPrep.R")

source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

library(nimble)
library(MCMCvis)
library(mcmcplots)
library(tidyverse)


## ============================================================
## 1. SPECIES + FEEDING
## ============================================================

## Species coding: 0 = figginsi, 1 = inornata
species_island <- as.integer(ordered_islands %in% c("FFRC", "UCay", "U Cay"))
names(species_island) <- ordered_islands
species <- species_island[island_idx]

## Feeding coding: 0 = none, 1 = any feeding
feeding_by_island <- growth_info %>%
  mutate(feeding = tolower(trimws(feeding))) %>%
  group_by(island) %>%
  summarise(feed = as.integer(any(feeding != "none")), .groups = "drop")

feed_island <- feeding_by_island$feed[match(ordered_islands, feeding_by_island$island)]
feed <- feed_island[island_idx]

island_table <- tibble(
  island_idx = 1:n_islands,
  island = ordered_islands,
  species = ifelse(species_island == 1, "inornata", "figginsi"),
  feeding = ifelse(feed_island == 1, "fed", "unfed")
)

print(island_table)
stopifnot(!anyNA(species), !anyNA(feed))

## ============================================================
## 2. COHORT DENSITY
##
## Density = MNA / usable habitat.
## Observed island-year density where available.
## Missing historical years = mean observed density for island.
## ============================================================

habitat <- read_csv("Data/Iguana_islands_habitat.csv", show_col_types = FALSE) %>%
  transmute(island = islands, usable_m2 = Useable)

## Observed span of each individual
first_last <- growth_info %>%
  transmute(ID, island,
            first_year = tF_cal - Origin + 1,
            last_year  = tL_cal - Origin + 1)

## Minimum number alive (MNA) during observed life-history span
mna_long <- first_last %>%
  rowwise() %>%
  mutate(year_idx = list(seq(first_year, last_year))) %>%
  unnest(year_idx) %>%
  ungroup() %>%
  count(island, year_idx, name = "MNA") %>%
  left_join(habitat, by = "island") %>%
  mutate(density_raw = MNA / (usable_m2 / 10000))       # animals / usable ha

## Only island-years in which surveys occurred
survey_long <- actual_surveys %>%
  transmute(island = Island, year_idx = yr_adj)

density_observed <- survey_long %>%
  left_join(mna_long %>% select(island, year_idx, MNA, density_raw),
            by = c("island", "year_idx")) %>%
  mutate(MNA = replace_na(MNA, 0L), density_raw = replace_na(density_raw, 0))

## Mean observed density for each island
island_mean_density <- density_observed %>%
  group_by(island) %>%
  summarise(mean_density_raw = mean(density_raw, na.rm = TRUE), .groups = "drop")

## Complete island x year series
density_complete <- expand_grid(island = ordered_islands, year_idx = 1:tMax) %>%
  left_join(density_observed %>% select(island, year_idx, density_raw),
            by = c("island", "year_idx")) %>%
  left_join(island_mean_density, by = "island") %>%
  mutate(
    density_source = if_else(is.na(density_raw), "Island mean", "Observed cohort"),
    density_filled_raw = coalesce(density_raw, mean_density_raw)
  )

## Global standardisation after filling
dens_mu <- mean(log1p(density_complete$density_filled_raw), na.rm = TRUE)
dens_sd <- sd(log1p(density_complete$density_filled_raw), na.rm = TRUE)
if(!is.finite(dens_sd) || dens_sd <= 0) stop("Invalid density SD.")

density_complete <- density_complete %>%
  mutate(density_z = (log1p(density_filled_raw) - dens_mu) / dens_sd)

## Island x cohort/year matrix
density_matrix <- matrix(NA_real_, n_islands, tMax,
                         dimnames = list(ordered_islands, 1:tMax))

for(r in 1:nrow(density_complete)) {
  j <- match(density_complete$island[r], ordered_islands)
  density_matrix[j, density_complete$year_idx[r]] <- density_complete$density_z[r]
}

stopifnot(!anyNA(density_matrix))

cat("\nDensity source by island:\n")
print(density_complete %>% count(island, density_source))

cat("\nRaw island mean densities:\n")
print(island_mean_density)


## ============================================================
## 3. MODEL
## ============================================================

n_cohort <- tMax

code_POL <- nimbleCode({
  
  for(i in 1:nind) {
    
    ## Latent birth and death
    tB[i] ~ dunif(cintB[i,1], cintB[i,2])
    L[i] <- max(0, island_start_time[i] - tB[i])
    
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    censoredD[i] ~ dinterval(tD[i], cintD[i,1:2])
    
    ## Latent birth cohort and associated density
    cohort_idx[i] <- max(1, min(n_cohort, floor(tB[i]) + 1))
    dens_cohort[i] <- density_matrix[island_idx[i], cohort_idx[i]]
    
    ## --------------------------------------------------------
    ## Gompertz a: background mortality / pace
    ## --------------------------------------------------------
    log(amult[i]) <- log(a) +
      betaSEX[1] * sex[i] +
      betaDENS_a * dens_cohort[i] +
      u_cohort[cohort_idx[i]] +
      island_a[island_idx[i]]
    
    ## --------------------------------------------------------
    ## Gompertz b: actuarial ageing / senescence
    ## --------------------------------------------------------
    log(bmult[i]) <- log(b) +
      betaSEX[2] * sex[i] +
      betaFEED_b * feed[i] +
      betaSPEC_b * species[i] +
      island_b[island_idx[i]]
    
    ## Observation / detection model
    tB_year_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    nm_start[i] <- max(tB_year_idx[i], island_start_idx[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    
    nMpos[i] <- cum_effort[island_idx[i], end_year[i] + 1] -
      cum_effort[island_idx[i], nm_start[i]]
    
    nMpos_safe[i] <- max(y[i], nMpos[i])
    
    pd_raw[i] <- exp(
      y[i] * log(mean.p + 1e-10) +
        (nMpos_safe[i] - y[i]) * log(1 - mean.p + 1e-10)
    )
    
    pd[i] <- min(0.999999, pd_raw[i])
    dind[i] ~ dbern(pd[i])
  }
  
  ## Individual fixed effects
  for(k in 1:2) betaSEX[k] ~ dnorm(0, sd = 1)
  
  betaDENS_a ~ dnorm(0, sd = 1)
  betaFEED_b ~ dnorm(0, sd = 1)
  betaSPEC_b ~ dnorm(0, sd = 1)
  
  ## Residual island variation: partial pooling
  sigma_island_a ~ dunif(0, 2)
  sigma_island_b ~ dunif(0, 2)
  
  for(r in 1:n_islands) {
    island_a_raw[r] ~ dnorm(0, sd = 1)
    island_b_raw[r] ~ dnorm(0, sd = 1)
    
    island_a[r] <- sigma_island_a * island_a_raw[r]
    island_b[r] <- sigma_island_b * island_b_raw[r]
  }
  
  ## Temporal birth-cohort effect: RW1
  sigma_cohort ~ dunif(0, 1)
  
  u_cohort[1] <- 0
  for(tt in 2:n_cohort) {
    u_cohort[tt] ~ dnorm(u_cohort[tt-1], sd = sigma_cohort)
  }
  
  ## Baselines
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0, 1)
})


## ============================================================
## 4. CONSTANTS + DATA
## ============================================================

consts_POL <- list(
  nind = nind, tMax = tMax, n_cohort = n_cohort, n_islands = n_islands,
  cintB = as.matrix(cintB_adj), cintD = as.matrix(cintD_adj),
  sex = as.numeric(sex), feed = as.numeric(feed), species = as.numeric(species),
  island_idx = as.numeric(island_idx),
  island_start_time = as.numeric(island_start_time),
  island_start_idx = as.numeric(island_start_idx),
  y = as.numeric(y)
)

data_POL <- list(
  censoredD = as.numeric(censoredD_vec),
  dind = rep(1, nind),
  cum_effort = as.matrix(cum_effort_for_model),
  density_matrix = as.matrix(density_matrix)
)


## ============================================================
## 5. INITIAL VALUES
## ============================================================

make_inits_POL <- function(chain_id = 1) {
  
  set.seed(1000 + chain_id)
  
  tB_start <- runif(nind, cintB_adj[,1], cintB_adj[,2])
  tD_start <- numeric(nind)
  
  for(i in 1:nind) {
    if(censoredD_vec[i] == 1) {
      tD_start[i] <- runif(1, cintD_adj[i,1] + 0.001, cintD_adj[i,2] - 0.001)
    } else {
      tD_start[i] <- cintD_adj[i,2] + runif(1, 1, 5)
    }
    
    if(tD_start[i] <= tB_start[i]) tD_start[i] <- tB_start[i] + runif(1,1,5)
  }
  
  u_start <- numeric(n_cohort)
  if(n_cohort > 1) u_start[2:n_cohort] <- cumsum(rnorm(n_cohort - 1, 0, 0.02))
  
  list(
    tB = tB_start,
    tstar = pmax(0.001, tD_start - tB_start),
    
    a = exp(rnorm(1, log(0.002), 0.25)),
    b = exp(rnorm(1, log(0.15), 0.20)),
    mean.p = runif(1, 0.20, 0.35),
    
    betaSEX = rnorm(2, 0, 0.10),
    betaDENS_a = rnorm(1, 0, 0.10),
    betaFEED_b = rnorm(1, 0, 0.10),
    betaSPEC_b = rnorm(1, 0, 0.10),
    
    sigma_island_a = runif(1, 0.05, 0.4),
    sigma_island_b = runif(1, 0.03, 0.3),
    island_a_raw = rnorm(n_islands, 0, 0.2),
    island_b_raw = rnorm(n_islands, 0, 0.2),
    
    sigma_cohort = runif(1, 0.01, 0.12),
    u_cohort = c(NA, u_start[2:n_cohort])
  )
}

inits_POL <- lapply(1:3, make_inits_POL)

## ============================================================
## 6. BUILD + CHECK MODEL
## ============================================================

model_POL <- nimbleModel(
  code = code_POL,
  constants = consts_POL,
  data = data_POL,
  inits = inits_POL[[1]],
  calculate = FALSE
)

lp <- model_POL$calculate()
cat("\nInitial log probability:", lp, "\n")
if(!is.finite(lp)) stop("Initial log probability is not finite.")

## ============================================================
## 7. CONFIGURE MCMC
## ============================================================

conf_POL <- configureMCMC(
  model_POL,
  monitors = c(
    "a","b","mean.p",
    "betaSEX","betaDENS_a","betaFEED_b","betaSPEC_b",
    "sigma_island_a","sigma_island_b","island_a","island_b",
    "sigma_cohort","u_cohort"
  ),
  enableWAIC = TRUE
)

## Baseline Gompertz
conf_POL$removeSamplers(c("a","b"))
conf_POL$addSampler(target = c("a","b"), type = "AF_slice")

## Variance parameters
for(v in c("sigma_island_a","sigma_island_b","sigma_cohort")) {
  conf_POL$removeSamplers(v)
  conf_POL$addSampler(target = v, type = "slice")
}

## Fixed effects
for(v in c("betaSEX[1]","betaSEX[2]","betaDENS_a","betaFEED_b","betaSPEC_b")) {
  conf_POL$removeSamplers(v)
  conf_POL$addSampler(target = v, type = "slice")
}

## Latent birth/lifespan
for(i in 1:nind) {
  conf_POL$removeSamplers(c(paste0("tB[",i,"]"), paste0("tstar[",i,"]")))
  conf_POL$addSampler(
    target = c(paste0("tB[",i,"]"), paste0("tstar[",i,"]")),
    type = "AF_slice"
  )
}

mcmc_POL <- buildMCMC(conf_POL)

## ============================================================
## 8. COMPILE
## ============================================================

cModel_POL <- compileNimble(model_POL, resetFunctions = TRUE)
cMCMC_POL <- compileNimble(mcmc_POL, project = model_POL)


## ============================================================
## 9. TEST RUN
## ============================================================

run_POL_test <- runMCMC(
  cMCMC_POL,
  niter = 100000,
  nburnin = 30000,
  nchains = 3,
  thin = 10,
  inits = inits_POL,
  setSeed = TRUE,
  samplesAsCodaMCMC = TRUE,
  summary = TRUE,
  WAIC = TRUE,
  progressBar = TRUE
)

saveRDS(run_POL_test, "Model_PaceOfLife_TEST_samples.rds")


## ============================================================
## 10. DIAGNOSTICS
## ============================================================

main_pars <- c(
  "a","b","betaSEX",
  "betaDENS_a","betaFEED_b","betaSPEC_b",
  "sigma_island_a","sigma_island_b","sigma_cohort","mean.p"
)

print(MCMCsummary(run_POL_test$samples, params = main_pars))
print(run_POL_test$WAIC)

cat("\nResidual island effects:\n")
print(MCMCsummary(run_POL_test$samples, params = c("island_a","island_b")))

cat("\nTemporal cohort effects:\n")
print(MCMCsummary(run_POL_test$samples, params = "u_cohort"))

mcmcplot(run_POL_test$samples, parms = main_pars)


## ============================================================
## 11. FINAL RUN - ONLY AFTER TEST LOOKS GOOD
## ============================================================

# run_POL <- runMCMC(
#   cMCMC_POL,
#   niter = 500000,
#   nburnin = 100000,
#   nchains = 3,
#   thin = 10,
#   inits = inits_POL,
#   setSeed = TRUE,
#   samplesAsCodaMCMC = TRUE,
#   summary = TRUE,
#   WAIC = TRUE,
#   progressBar = TRUE
# )
#
# saveRDS(run_POL, "Model_PaceOfLife_samples.rds")
# print(MCMCsummary(run_POL$samples))
# print(run_POL$WAIC)