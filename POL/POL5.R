## ============================================================
## MODEL_PaceOfLife_FINAL.R
##
## FINAL CANDIDATE MODEL
##
## log(a):
##   sex + cohort density + feeding + 5-year temporal effect
##   + residual island effect
##
## log(b):
##   sex + feeding
##
## Species excluded following negligible species effect.
## Final two temporal blocks pooled.
## All continuous/intercept components reparameterised for mixing.
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
## 1. CENTRED SEX + FEEDING
## ============================================================

## Sex: 0 = female, 1 = male
sex_mean <- mean(sex)
sex_c <- sex - sex_mean

## Feeding: 0 = none, 1 = any feeding
feeding_by_island <- growth_info %>%
  mutate(feeding = tolower(trimws(feeding))) %>%
  group_by(island) %>%
  summarise(feed = as.integer(any(feeding != "none")), .groups = "drop")

feed_island <- feeding_by_island$feed[match(ordered_islands, feeding_by_island$island)]

## Centre feeding at island/block level
feed_island_mean <- mean(feed_island)
feed_c_island <- feed_island - feed_island_mean
feed_c <- feed_c_island[island_idx]

island_table <- tibble(
  island_idx = seq_len(n_islands),
  island = ordered_islands,
  feeding = ifelse(feed_island == 1, "fed", "unfed"),
  feed_c = feed_c_island
)

cat("\nIsland feeding structure:\n")
print(island_table)

cat("\nCentred means:\n")
cat("sex_c =", mean(sex_c), "\n")
cat("feed_c island mean =", mean(feed_c_island), "\n")

stopifnot(
  length(sex_c) == nind,
  length(feed_c) == nind,
  !anyNA(sex_c),
  !anyNA(feed_c)
)

## ============================================================
## 2. COHORT DENSITY
##
## MNA / usable habitat area.
## Missing historical cohorts use island mean observed density.
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

## Standardise globally AFTER filling
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
##
## Pool final two raw blocks due sparse final cohort.
## ============================================================

block_width <- 5
n_time_blocks_raw <- ceiling(tMax / block_width)
n_time_blocks <- n_time_blocks_raw - 1

cat(
  "\nTemporal structure:",
  n_time_blocks_raw,
  "raw blocks ->",
  n_time_blocks,
  "modelled blocks after pooling final two\n"
)

birth_mid <- rowMeans(cintB_adj)
birth_block_raw <- floor((birth_mid - 1) / block_width) + 1
birth_block_raw <- pmax(1, pmin(n_time_blocks_raw, birth_block_raw))
birth_block <- pmin(birth_block_raw, n_time_blocks)

cat("\nApproximate birth cohort counts:\n")
print(table(birth_block))

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
    
    time_block_raw[i] <- max(
      1,
      min(n_time_blocks_raw, floor((cohort_idx[i] - 1) / block_width) + 1)
    )
    
    ## Pool final two temporal blocks
    time_block[i] <- min(time_block_raw[i], n_time_blocks)
    
    ## Density associated with latent birth cohort
    dens_cohort[i] <- density_matrix[island_idx[i], cohort_idx[i]]
    
    ## ========================================================
    ## GOMPERTZ a
    ##
    ## General/background mortality level
    ##
    ## feeding here tests whether supplementation shifts
    ## the overall mortality schedule.
    ## ========================================================
    
    log(amult[i]) <- log_a +
      betaSEX[1] * sex_c[i] +
      betaDENS_a * dens_cohort[i] +
      betaFEED_a * feed_c[i] +
      u_time[time_block[i]] +
      island_a[island_idx[i]]
    
    ## ========================================================
    ## GOMPERTZ b
    ##
    ## Rate of actuarial senescence
    ##
    ## feeding here tests whether supplementation changes
    ## the age-related acceleration in mortality.
    ## ========================================================
    
    log(bmult[i]) <- log_b +
      betaSEX[2] * sex_c[i] +
      betaFEED_b * feed_c[i]
    
    ## --------------------------------------------------------
    ## OBSERVATION / DETECTION COMPONENT
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
  betaFEED_a ~ dnorm(0, sd = 1)
  betaFEED_b ~ dnorm(0, sd = 1)
  
  ## ==========================================================
  ## RESIDUAL ISLAND EFFECT ON a
  ##
  ## Exact zero-centering prevents drift against log_a.
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
  ## Exact zero-centering prevents drift against log_a.
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
  ## ==========================================================
  
  log_a ~ dnorm(log(0.002), sd = 2)
  log_b ~ dnorm(log(0.20), sd = 1)
  
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
  n_time_blocks_raw = n_time_blocks_raw,
  n_time_blocks = n_time_blocks,
  block_width = block_width,
  n_islands = n_islands,
  cintB = as.matrix(cintB_adj),
  cintD = as.matrix(cintD_adj),
  sex_c = as.numeric(sex_c),
  feed_c = as.numeric(feed_c),
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
  
  set.seed(4000 + chain_id)
  
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
    betaFEED_a = rnorm(1, 0, 0.10),
    betaFEED_b = rnorm(1, 0, 0.10),
    
    sigma_island_a = runif(1, 0.5, 1.5),
    island_a_raw = rnorm(n_islands, 0, 0.5),
    
    sigma_time = runif(1, 1.0, 2.5),
    u_time_raw = rnorm(n_time_blocks, 0, 0.5),
    
    mean.p = runif(1, 0.20, 0.35)
  )
}

inits_POL <- lapply(seq_len(3), make_inits_POL)

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
    "betaFEED_a",
    "betaFEED_b",
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

## a-side joint update
conf_POL$removeSamplers(
  c(
    "log_a",
    "betaSEX[1]",
    "betaDENS_a",
    "betaFEED_a"
  )
)

conf_POL$addSampler(
  target = c(
    "log_a",
    "betaSEX[1]",
    "betaDENS_a",
    "betaFEED_a"
  ),
  type = "AF_slice"
)

## b-side joint update
conf_POL$removeSamplers(
  c(
    "log_b",
    "betaSEX[2]",
    "betaFEED_b"
  )
)

conf_POL$addSampler(
  target = c(
    "log_b",
    "betaSEX[2]",
    "betaFEED_b"
  ),
  type = "AF_slice"
)

## Scale parameters
for(v in c("sigma_island_a", "sigma_time")) {
  conf_POL$removeSamplers(v)
  conf_POL$addSampler(target = v, type = "slice")
}

## Latent birth + lifespan
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
## 11. FINAL LONG RUN
## ============================================================

run_POL <- runMCMC(
  cMCMC_POL,
  niter = 500000,
  nburnin = 100000,
  nchains = 3,
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
  run_POL,
  "Model_PaceOfLife_FINAL_samples.rds"
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
  "betaFEED_a",
  "betaFEED_b",
  "sigma_island_a",
  "sigma_time",
  "mean.p"
)

cat("\nMAIN PARAMETER SUMMARY\n")

main_summary <- MCMCsummary(
  run_POL$samples,
  params = main_pars
)

print(main_summary)

## ============================================================
## 14. TEMPORAL EFFECTS
## ============================================================

cat("\nTEMPORAL 5-YEAR EFFECTS\n")

time_summary <- MCMCsummary(
  run_POL$samples,
  params = "u_time"
)

print(time_summary)

## ============================================================
## 15. RESIDUAL ISLAND EFFECT ON a
## ============================================================

cat("\nRESIDUAL ISLAND EFFECT ON a\n")

island_summary <- MCMCsummary(
  run_POL$samples,
  params = "island_a"
)

print(island_summary)

## ============================================================
## 16. WAIC
## ============================================================

cat("\nWAIC\n")
print(run_POL$WAIC)

## ============================================================
## 17. EFFECT-SCALE SUMMARIES
## ============================================================

samps <- as.matrix(run_POL$samples)

effect_summary <- function(x) {
  c(
    median = median(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975),
    Pr_gt_1 = mean(x > 1)
  )
}

cat("\nMale / female ratio in a:\n")
print(effect_summary(exp(samps[, "betaSEX[1]"])))

cat("\nMale / female ratio in b:\n")
print(effect_summary(exp(samps[, "betaSEX[2]"])))

cat("\nEffect on a per +1 SD cohort density:\n")
print(effect_summary(exp(samps[, "betaDENS_a"])))

cat("\nFed / unfed ratio in a:\n")
print(effect_summary(exp(samps[, "betaFEED_a"])))

cat("\nFed / unfed ratio in b:\n")
print(effect_summary(exp(samps[, "betaFEED_b"])))

## Also report posterior probabilities on coefficient scale
cat("\nPosterior P(betaDENS_a < 0):", mean(samps[, "betaDENS_a"] < 0), "\n")
cat("Posterior P(betaFEED_a < 0):", mean(samps[, "betaFEED_a"] < 0), "\n")
cat("Posterior P(betaFEED_b < 0):", mean(samps[, "betaFEED_b"] < 0), "\n")
cat("Posterior P(betaSEX[2] > 0):", mean(samps[, "betaSEX[2]"] > 0), "\n")

## ============================================================
## 18. TRACE PLOTS
## ============================================================

mcmcplot(
  run_POL$samples,
  parms = main_pars
)

## ============================================================
## 19. SAVE SUMMARIES
## ============================================================

saveRDS(
  list(
    main = main_summary,
    temporal = time_summary,
    island = island_summary,
    WAIC = run_POL$WAIC
  ),
  "Model_PaceOfLife_FINAL_summaries.rds"
)