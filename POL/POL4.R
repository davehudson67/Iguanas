## ============================================================
## MODEL_PaceOfLife.R
##
## REPARAMETERISED VERSION
##
## log(a): centred sex + cohort density + centred 5-year temporal effect
##         + centred residual island effect
##
## log(b): centred sex + centred feeding + centred species
##
## Key changes:
##   - model log_a and log_b directly
##   - centre sex, feeding and species
##   - exact zero-centred island/time effects
##   - no island effect on b
##   - sigma_time prior widened to 0-5
##   - jointly sample highly correlated b-side parameters
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
## 1. CENTRED SEX, SPECIES + FEEDING
## ============================================================

## Sex: original coding 0 = female, 1 = male
sex_mean <- mean(sex)
sex_c <- sex - sex_mean

## Species: 0 = figginsi, 1 = inornata
species_island <- as.integer(ordered_islands %in% c("FFRC", "UCay", "U Cay"))
names(species_island) <- ordered_islands

## Feeding: 0 = none, 1 = any feeding
feeding_by_island <- growth_info %>%
  mutate(feeding = tolower(trimws(feeding))) %>%
  group_by(island) %>%
  summarise(feed = as.integer(any(feeding != "none")), .groups = "drop")

feed_island <- feeding_by_island$feed[match(ordered_islands, feeding_by_island$island)]

## Centre these at the ISLAND level because treatment/species vary by island
feed_island_mean <- mean(feed_island)
species_island_mean <- mean(species_island)

feed_c_island <- feed_island - feed_island_mean
species_c_island <- species_island - species_island_mean

## Individual vectors simply inherit their island value
feed_c <- feed_c_island[island_idx]
species_c <- species_c_island[island_idx]

island_table <- tibble(
  island_idx = seq_len(n_islands),
  island = ordered_islands,
  species = ifelse(species_island == 1, "inornata", "figginsi"),
  feeding = ifelse(feed_island == 1, "fed", "unfed"),
  feed_c = feed_c_island,
  species_c = species_c_island
)

cat("\nIsland structure:\n")
print(island_table)

cat("\nCentred covariate means:\n")
cat("sex_c =", mean(sex_c), "\n")
cat("feed_c island mean =", mean(feed_c_island), "\n")
cat("species_c island mean =", mean(species_c_island), "\n")

stopifnot(
  length(sex_c) == nind,
  length(feed_c) == nind,
  length(species_c) == nind,
  !anyNA(sex_c),
  !anyNA(feed_c),
  !anyNA(species_c)
)

## ============================================================
## 2. COHORT DENSITY
## ============================================================

habitat <- read_csv("Data/Iguana_islands_habitat.csv", show_col_types = FALSE) %>%
  transmute(island = islands, usable_m2 = Useable)

first_last <- growth_info %>%
  transmute(
    ID,
    island,
    first_year = tF_cal - Origin + 1,
    last_year = tL_cal - Origin + 1
  )

mna_long <- first_last %>%
  rowwise() %>%
  mutate(year_idx = list(seq(first_year, last_year))) %>%
  unnest(year_idx) %>%
  ungroup() %>%
  count(island, year_idx, name = "MNA") %>%
  left_join(habitat, by = "island") %>%
  mutate(density_raw = MNA / (usable_m2 / 10000))

survey_long <- actual_surveys %>%
  transmute(island = Island, year_idx = yr_adj)

density_observed <- survey_long %>%
  left_join(
    mna_long %>%
      select(island, year_idx, MNA, density_raw),
    by = c("island", "year_idx")
  ) %>%
  mutate(
    MNA = replace_na(MNA, 0L),
    density_raw = replace_na(density_raw, 0)
  )

island_mean_density <- density_observed %>%
  group_by(island) %>%
  summarise(mean_density_raw = mean(density_raw, na.rm = TRUE), .groups = "drop")

cat("\nMean observed density by island:\n")
print(island_mean_density)

## Missing historical cohorts receive their island mean density
density_complete <- expand_grid(island = ordered_islands, year_idx = seq_len(tMax)) %>%
  left_join(
    density_observed %>%
      select(island, year_idx, density_raw),
    by = c("island", "year_idx")
  ) %>%
  left_join(island_mean_density, by = "island") %>%
  mutate(
    density_source = if_else(is.na(density_raw), "Island mean", "Observed cohort"),
    density_filled_raw = coalesce(density_raw, mean_density_raw)
  )

cat("\nDensity source by island:\n")

print(
  density_complete %>%
    count(island, density_source) %>%
    arrange(island, density_source)
)

## Global standardisation after filling
dens_mu <- mean(log1p(density_complete$density_filled_raw), na.rm = TRUE)
dens_sd <- sd(log1p(density_complete$density_filled_raw), na.rm = TRUE)

if(!is.finite(dens_sd) || dens_sd <= 0) stop("Invalid density standard deviation.")

density_complete <- density_complete %>%
  mutate(density_z = (log1p(density_filled_raw) - dens_mu) / dens_sd)

density_matrix <- matrix(
  NA_real_,
  nrow = n_islands,
  ncol = tMax,
  dimnames = list(ordered_islands, seq_len(tMax))
)

for(r in seq_len(nrow(density_complete))) {
  j <- match(density_complete$island[r], ordered_islands)
  tt <- density_complete$year_idx[r]
  density_matrix[j, tt] <- density_complete$density_z[r]
}

stopifnot(
  !anyNA(density_matrix),
  nrow(density_matrix) == n_islands,
  ncol(density_matrix) == tMax
)

cat("\nDensity matrix summary:\n")
print(summary(as.vector(density_matrix)))

## ============================================================
## 3. FIVE-YEAR TEMPORAL COHORT BLOCKS
## ============================================================

block_width <- 5
n_time_blocks <- ceiling(tMax / block_width)

cat("\nTemporal structure:", n_time_blocks, "blocks of", block_width, "years\n")

## ============================================================
## 4. MODEL
## ============================================================

code_POL <- nimbleCode({
  
  for(i in 1:nind) {
    
    ## --------------------------------------------------------
    ## LATENT LIFE HISTORY
    ## --------------------------------------------------------
    
    tB[i] ~ dunif(cintB[i,1], cintB[i,2])
    L[i] <- max(0, island_start_time[i] - tB[i])
    
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    
    censoredD[i] ~ dinterval(tD[i], cintD[i,1:2])
    
    ## --------------------------------------------------------
    ## LATENT BIRTH COHORT
    ## --------------------------------------------------------
    
    cohort_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    
    time_block[i] <- max(
      1,
      min(n_time_blocks, floor((cohort_idx[i] - 1) / block_width) + 1)
    )
    
    dens_cohort[i] <- density_matrix[island_idx[i], cohort_idx[i]]
    
    ## ========================================================
    ## GOMPERTZ a
    ##
    ## log_a = population-average log mortality intercept
    ##
    ## sex + cohort density + temporal block + island deviation
    ## ========================================================
    
    log(amult[i]) <- log_a +
      betaSEX[1] * sex_c[i] +
      betaDENS_a * dens_cohort[i] +
      u_time[time_block[i]] +
      island_a[island_idx[i]]
    
    ## ========================================================
    ## GOMPERTZ b
    ##
    ## log_b = population-average log senescence parameter
    ##
    ## sex + feeding + species
    ## ========================================================
    
    log(bmult[i]) <- log_b +
      betaSEX[2] * sex_c[i] +
      betaFEED_b * feed_c[i] +
      betaSPEC_b * species_c[i]
    
    ## --------------------------------------------------------
    ## OBSERVATION / DETECTION MODEL
    ## --------------------------------------------------------
    
    tB_year_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    nm_start[i] <- max(tB_year_idx[i], island_start_idx[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    
    nMpos[i] <- cum_effort[island_idx[i], end_year[i] + 1] -
      cum_effort[island_idx[i], nm_start[i]]
    
    nMpos_safe[i] <- max(y[i], nMpos[i])
    
    pd_raw[i] <- exp(
      y[i] * log(mean.p + 1.0E-10) +
        (nMpos_safe[i] - y[i]) * log(1 - mean.p + 1.0E-10)
    )
    
    pd[i] <- min(0.999999, pd_raw[i])
    
    dind[i] ~ dbern(pd[i])
  }
  
  ## ==========================================================
  ## FIXED EFFECTS
  ## ==========================================================
  
  for(k in 1:2) betaSEX[k] ~ dnorm(0, sd = 1)
  
  betaDENS_a ~ dnorm(0, sd = 1)
  betaFEED_b ~ dnorm(0, sd = 1)
  betaSPEC_b ~ dnorm(0, sd = 1)
  
  ## ==========================================================
  ## RESIDUAL ISLAND EFFECT ON a
  ##
  ## Raw effects are centred exactly at zero before scaling.
  ## This prevents the island effects shifting the global log_a.
  ## ==========================================================
  
  sigma_island_a ~ dunif(0, 2)
  
  for(r in 1:n_islands) {
    island_a_raw[r] ~ dnorm(0, sd = 1)
  }
  
  mean_island_a_raw <- sum(island_a_raw[1:n_islands]) / n_islands
  
  for(r in 1:n_islands) {
    island_a[r] <- sigma_island_a * (island_a_raw[r] - mean_island_a_raw)
  }
  
  ## ==========================================================
  ## FIVE-YEAR TEMPORAL EFFECT
  ##
  ## Also exactly centred around zero.
  ## sigma_time now allowed up to 5.
  ## ==========================================================
  
  sigma_time ~ dunif(0, 5)
  
  for(tt in 1:n_time_blocks) {
    u_time_raw[tt] ~ dnorm(0, sd = 1)
  }
  
  mean_u_time_raw <- sum(u_time_raw[1:n_time_blocks]) / n_time_blocks
  
  for(tt in 1:n_time_blocks) {
    u_time[tt] <- sigma_time * (u_time_raw[tt] - mean_u_time_raw)
  }
  
  ## ==========================================================
  ## POPULATION-AVERAGE GOMPERTZ PARAMETERS
  ##
  ## Directly estimate on log scale.
  ## ==========================================================
  
  log_a ~ dnorm(log(0.002), sd = 2)
  log_b ~ dnorm(log(0.20), sd = 1)
  
  ## Deterministic natural-scale values for output/interpretation
  a <- exp(log_a)
  b <- exp(log_b)
  
  mean.p ~ dunif(0, 1)
})

## ============================================================
## 5. CONSTANTS + DATA
## ============================================================

consts_POL <- list(
  nind = nind,
  tMax = tMax,
  n_time_blocks = n_time_blocks,
  block_width = block_width,
  n_islands = n_islands,
  cintB = as.matrix(cintB_adj),
  cintD = as.matrix(cintD_adj),
  sex_c = as.numeric(sex_c),
  feed_c = as.numeric(feed_c),
  species_c = as.numeric(species_c),
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
## 6. INITIAL VALUES
## ============================================================

make_inits_POL <- function(chain_id = 1) {
  
  set.seed(2000 + chain_id)
  
  tB_start <- runif(nind, min = cintB_adj[,1], max = cintB_adj[,2])
  tD_start <- numeric(nind)
  
  for(i in seq_len(nind)) {
    
    if(censoredD_vec[i] == 1) {
      tD_start[i] <- runif(
        1,
        min = cintD_adj[i,1] + 0.001,
        max = cintD_adj[i,2] - 0.001
      )
    }
    
    if(censoredD_vec[i] == 2) {
      tD_start[i] <- cintD_adj[i,2] + runif(1, 1, 5)
    }
    
    if(tD_start[i] <= tB_start[i]) {
      tD_start[i] <- tB_start[i] + runif(1, 1, 5)
    }
  }
  
  list(
    tB = tB_start,
    tstar = pmax(0.001, tD_start - tB_start),
    
    log_a = rnorm(1, log(0.002), 0.20),
    log_b = rnorm(1, log(0.20), 0.15),
    
    betaSEX = rnorm(2, 0, 0.10),
    betaDENS_a = rnorm(1, 0, 0.10),
    betaFEED_b = rnorm(1, 0, 0.10),
    betaSPEC_b = rnorm(1, 0, 0.10),
    
    sigma_island_a = runif(1, 0.5, 1.5),
    island_a_raw = rnorm(n_islands, 0, 0.5),
    
    sigma_time = runif(1, 0.8, 2.0),
    u_time_raw = rnorm(n_time_blocks, 0, 0.5),
    
    mean.p = runif(1, 0.20, 0.35)
  )
}

inits_POL <- lapply(seq_len(2), make_inits_POL)

## ============================================================
## 7. BUILD MODEL
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
## 8. CONFIGURE MCMC
## ============================================================

conf_POL <- configureMCMC(
  model_POL,
  monitors = c(
    "a",
    "b",
    "log_a",
    "log_b",
    "mean.p",
    "betaSEX",
    "betaDENS_a",
    "betaFEED_b",
    "betaSPEC_b",
    "sigma_island_a",
    "island_a",
    "sigma_time",
    "u_time"
  ),
  enableWAIC = TRUE
)

## ============================================================
## 9. CUSTOM SAMPLERS
## ============================================================

## ------------------------------------------------------------
## a-side:
## jointly update global log(a), sex effect on a, density effect
## ------------------------------------------------------------

conf_POL$removeSamplers(c("log_a", "betaSEX[1]", "betaDENS_a"))

conf_POL$addSampler(
  target = c("log_a", "betaSEX[1]", "betaDENS_a"),
  type = "AF_slice"
)

## ------------------------------------------------------------
## b-side:
## jointly update global log(b), sex, feeding and species
##
## This is the key change intended to improve the strong
## correlation between b and the island-level predictors.
## ------------------------------------------------------------

conf_POL$removeSamplers(c("log_b", "betaSEX[2]", "betaFEED_b", "betaSPEC_b"))

conf_POL$addSampler(
  target = c("log_b", "betaSEX[2]", "betaFEED_b", "betaSPEC_b"),
  type = "AF_slice"
)

## ------------------------------------------------------------
## Scale parameters
## ------------------------------------------------------------

for(v in c("sigma_island_a", "sigma_time")) {
  
  conf_POL$removeSamplers(v)
  
  conf_POL$addSampler(
    target = v,
    type = "slice"
  )
}

## ------------------------------------------------------------
## Latent birth + lifespan
## ------------------------------------------------------------

for(i in seq_len(nind)) {
  
  conf_POL$removeSamplers(
    c(
      paste0("tB[", i, "]"),
      paste0("tstar[", i, "]")
    )
  )
  
  conf_POL$addSampler(
    target = c(
      paste0("tB[", i, "]"),
      paste0("tstar[", i, "]")
    ),
    type = "AF_slice"
  )
}

## ============================================================
## 10. BUILD + COMPILE
## ============================================================

mcmc_POL <- buildMCMC(conf_POL)

cModel_POL <- compileNimble(
  model_POL,
  resetFunctions = TRUE
)

cMCMC_POL <- compileNimble(
  mcmc_POL,
  project = model_POL
)

## ============================================================
## 11. TEST RUN
## ============================================================

run_POL_test4 <- runMCMC(
  cMCMC_POL,
  niter = 50000,
  nburnin = 12000,
  nchains = 2,
  thin = 10,
  inits = inits_POL,
  setSeed = TRUE,
  samplesAsCodaMCMC = TRUE,
  summary = TRUE,
  WAIC = TRUE,
  progressBar = TRUE
)

## ============================================================
## 12. SAVE
## ============================================================

saveRDS(
  run_POL_test4,
  "Model_PaceOfLife_5yr_centred_TEST_samples.rds"
)

## ============================================================
## 13. MAIN PARAMETER SUMMARY
## ============================================================

main_pars <- c(
  "a",
  "b",
  "log_a",
  "log_b",
  "betaSEX",
  "betaDENS_a",
  "betaFEED_b",
  "betaSPEC_b",
  "sigma_island_a",
  "sigma_time",
  "mean.p"
)

cat("\nMAIN PARAMETER SUMMARY\n")

print(
  MCMCsummary(
    run_POL_test4$samples,
    params = main_pars
  )
)

## ============================================================
## 14. TEMPORAL EFFECTS
## ============================================================

cat("\nTEMPORAL 5-YEAR EFFECTS\n")

print(
  MCMCsummary(
    run_POL_test4$samples,
    params = "u_time"
  )
)

## ============================================================
## 15. RESIDUAL ISLAND EFFECT ON a
## ============================================================

cat("\nRESIDUAL ISLAND EFFECT ON a\n")

print(
  MCMCsummary(
    run_POL_test4$samples,
    params = "island_a"
  )
)

## ============================================================
## 16. WAIC
## ============================================================

cat("\nWAIC\n")
print(run_POL_test4$WAIC)

## ============================================================
## 17. TRACE PLOTS
## ============================================================

mcmcplot(
  run_POL_test4$samples,
  parms = main_pars
)

## ============================================================
## 18. USEFUL EFFECT-SCALE SUMMARIES
##
## Because predictors are on log(a) / log(b):
##
## exp(beta) = multiplicative change in parameter.
##
## For feeding:
## exp(betaFEED_b) = ratio of b in fed vs unfed islands.
##
## For sex:
## exp(betaSEX[2]) = male/female ratio in b.
##
## Density:
## exp(betaDENS_a) = multiplicative change in a per +1 SD density.
## ============================================================

samps <- as.matrix(run_POL_test4$samples)

effect_summary <- function(x) {
  c(
    median = median(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975)
  )
}

cat("\nMultiplicative effect on a per +1 SD density:\n")
print(effect_summary(exp(samps[, "betaDENS_a"])))

cat("\nMale / female ratio in b:\n")
print(effect_summary(exp(samps[, "betaSEX[2]"])))

cat("\nFed / unfed ratio in b:\n")
print(effect_summary(exp(samps[, "betaFEED_b"])))

cat("\nInornata / figginsi ratio in b:\n")
print(effect_summary(exp(samps[, "betaSPEC_b"])))

## ============================================================
## 19. FINAL LONG RUN
##
## Only launch once this model shows satisfactory convergence.
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
# saveRDS(
#   run_POL,
#   "Model_PaceOfLife_5yr_centred_samples.rds"
# )
#
# print(
#   MCMCsummary(
#     run_POL$samples
#   )
# )