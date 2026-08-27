## ============================================================
## MODEL_PaceOfLife_LOO_Island_LONG.R
##
## LEAVE-ONE-ISLAND-OUT SENSITIVITY ANALYSIS
##
## LONG RUN:
##   500,000 iterations
##   100,000 burn-in
##   thin = 10
##   3 chains
##
## Final model:
##
## log(a):
##   sex + cohort density + feeding + temporal cohort
##   + residual island effect
##
## log(b):
##   sex + feeding
##
## Each model removes one complete island.
##
## Primary question:
##   Is the estimated feeding effect on b robust to removal
##   of individual island treatment blocks?
## ============================================================

source("DataPrep.R")

source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

library(nimble)
library(MCMCvis)
library(tidyverse)
library(coda)

## ============================================================
## 1. SAVE FULL-DATA OBJECTS
## ============================================================

growth_info_FULL <- growth_info
cintB_FULL <- cintB_adj
cintD_FULL <- cintD_adj
censoredD_FULL <- censoredD_vec
sex_FULL <- sex
y_FULL <- y
island_idx_FULL <- island_idx

ordered_islands_FULL <- ordered_islands
n_islands_FULL <- n_islands
nind_FULL <- nind
tMax_FULL <- tMax

cum_effort_FULL <- cum_effort_for_model
island_start_time_FULL <- island_start_time
island_start_idx_FULL <- island_start_idx
actual_surveys_FULL <- actual_surveys

island_ind_FULL <- ordered_islands_FULL[island_idx_FULL]

stopifnot(
  length(island_ind_FULL) == nind_FULL,
  length(sex_FULL) == nind_FULL,
  length(y_FULL) == nind_FULL,
  length(island_start_time_FULL) == nind_FULL,
  length(island_start_idx_FULL) == nind_FULL,
  nrow(cintB_FULL) == nind_FULL,
  nrow(cintD_FULL) == nind_FULL,
  nrow(cum_effort_FULL) == n_islands_FULL,
  ncol(cum_effort_FULL) >= tMax_FULL + 1
)

cat("\nFull dataset islands:\n")
print(table(island_ind_FULL))

## ============================================================
## 2. FULL-DATA FEEDING STRUCTURE
## ============================================================

feeding_by_island_FULL <- growth_info_FULL %>%
  mutate(feeding = tolower(trimws(feeding))) %>%
  group_by(island) %>%
  summarise(feed = as.integer(any(feeding != "none")), .groups = "drop")

feed_island_FULL <- feeding_by_island_FULL$feed[
  match(ordered_islands_FULL, feeding_by_island_FULL$island)
]

stopifnot(!anyNA(feed_island_FULL))

cat("\nFull-data feeding structure:\n")

print(
  tibble(
    island = ordered_islands_FULL,
    feeding = ifelse(feed_island_FULL == 1, "fed", "unfed")
  )
)

## ============================================================
## 3. BUILD FULL-DATA DENSITY MATRIX
##
## Same density construction/scaling as final model.
## LOO models remove the omitted island row only.
## ============================================================

habitat_FULL <- read_csv("Data/Iguana_islands_habitat.csv", show_col_types = FALSE) %>%
  transmute(island = islands, usable_m2 = Useable)

first_last_FULL <- growth_info_FULL %>%
  transmute(ID, island, first_year = tF_cal - Origin + 1, last_year = tL_cal - Origin + 1)

mna_long_FULL <- first_last_FULL %>%
  rowwise() %>%
  mutate(year_idx = list(seq(first_year, last_year))) %>%
  unnest(year_idx) %>%
  ungroup() %>%
  count(island, year_idx, name = "MNA") %>%
  left_join(habitat_FULL, by = "island") %>%
  mutate(density_raw = MNA / (usable_m2 / 10000))

survey_long_FULL <- actual_surveys_FULL %>%
  transmute(island = Island, year_idx = yr_adj)

density_observed_FULL <- survey_long_FULL %>%
  left_join(
    mna_long_FULL %>%
      select(island, year_idx, MNA, density_raw),
    by = c("island", "year_idx")
  ) %>%
  mutate(MNA = replace_na(MNA, 0L), density_raw = replace_na(density_raw, 0))

island_mean_density_FULL <- density_observed_FULL %>%
  group_by(island) %>%
  summarise(mean_density_raw = mean(density_raw, na.rm = TRUE), .groups = "drop")

density_complete_FULL <- expand_grid(island = ordered_islands_FULL, year_idx = seq_len(tMax_FULL)) %>%
  left_join(
    density_observed_FULL %>%
      select(island, year_idx, density_raw),
    by = c("island", "year_idx")
  ) %>%
  left_join(island_mean_density_FULL, by = "island") %>%
  mutate(
    density_source = if_else(is.na(density_raw), "Island mean", "Observed cohort"),
    density_filled_raw = coalesce(density_raw, mean_density_raw)
  )

if(anyNA(density_complete_FULL$density_filled_raw)) {
  stop("Missing density values remain after full-data density filling.")
}

dens_mu_FULL <- mean(log1p(density_complete_FULL$density_filled_raw))
dens_sd_FULL <- sd(log1p(density_complete_FULL$density_filled_raw))

if(!is.finite(dens_sd_FULL) || dens_sd_FULL <= 0) {
  stop("Invalid full-data density standard deviation.")
}

density_complete_FULL <- density_complete_FULL %>%
  mutate(density_z = (log1p(density_filled_raw) - dens_mu_FULL) / dens_sd_FULL)

density_matrix_FULL <- matrix(
  NA_real_,
  nrow = n_islands_FULL,
  ncol = tMax_FULL,
  dimnames = list(ordered_islands_FULL, seq_len(tMax_FULL))
)

for(r in seq_len(nrow(density_complete_FULL))) {
  jj <- match(density_complete_FULL$island[r], ordered_islands_FULL)
  tt <- density_complete_FULL$year_idx[r]
  density_matrix_FULL[jj,tt] <- density_complete_FULL$density_z[r]
}

stopifnot(!anyNA(density_matrix_FULL))

cat("\nFull density matrix summary:\n")
print(summary(as.vector(density_matrix_FULL)))

## ============================================================
## 4. TEMPORAL STRUCTURE
##
## 5-year blocks; final two raw blocks pooled.
## ============================================================

block_width <- 5
n_time_blocks_raw <- ceiling(tMax_FULL / block_width)
n_time_blocks <- n_time_blocks_raw - 1

cat(
  "\nTemporal structure:",
  n_time_blocks_raw,
  "raw blocks ->",
  n_time_blocks,
  "modelled blocks after pooling final two\n"
)

## ============================================================
## 5. MODEL CODE
## ============================================================

code_LOO <- nimbleCode({
  
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
    
    time_block[i] <- min(time_block_raw[i], n_time_blocks)
    
    dens_cohort[i] <- density_matrix[island_idx[i], cohort_idx[i]]
    
    ## ========================================================
    ## GOMPERTZ a
    ## ========================================================
    
    log(amult[i]) <- log_a +
      betaSEX[1] * sex_c[i] +
      betaDENS_a * dens_cohort[i] +
      betaFEED_a * feed_c[i] +
      u_time[time_block[i]] +
      island_a[island_idx[i]]
    
    ## ========================================================
    ## GOMPERTZ b
    ## ========================================================
    
    log(bmult[i]) <- log_b +
      betaSEX[2] * sex_c[i] +
      betaFEED_b * feed_c[i]
    
    ## --------------------------------------------------------
    ## DETECTION COMPONENT
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
  ## TEMPORAL COHORT EFFECT ON a
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
## 6. FUNCTION TO FIT ONE LOO MODEL
## ============================================================

fit_LOO_island <- function(
    omit_island,
    niter = 500000,
    nburnin = 100000,
    thin = 10,
    nchains = 3,
    print_samplers = FALSE
) {
  
  cat("\n")
  cat("============================================================\n")
  cat("LEAVING OUT:", omit_island, "\n")
  cat("============================================================\n")
  
  ## ----------------------------------------------------------
  ## INDIVIDUAL SUBSETTING
  ## ----------------------------------------------------------
  
  keep <- island_ind_FULL != omit_island
  
  growth_info_loo <- growth_info_FULL[keep,]
  cintB_loo <- cintB_FULL[keep,,drop = FALSE]
  cintD_loo <- cintD_FULL[keep,,drop = FALSE]
  censoredD_loo <- censoredD_FULL[keep]
  sex_loo <- sex_FULL[keep]
  y_loo <- y_FULL[keep]
  
  island_start_time_loo <- island_start_time_FULL[keep]
  island_start_idx_loo <- island_start_idx_FULL[keep]
  
  island_ind_loo <- island_ind_FULL[keep]
  
  ordered_islands_loo <- ordered_islands_FULL[
    ordered_islands_FULL != omit_island
  ]
  
  n_islands_loo <- length(ordered_islands_loo)
  nind_loo <- sum(keep)
  
  island_idx_loo <- match(island_ind_loo, ordered_islands_loo)
  
  cat("\nIndividuals retained:", nind_loo, "\n")
  cat("\nRemaining islands:\n")
  print(table(island_ind_loo))
  
  ## ----------------------------------------------------------
  ## HARD DIMENSION / INDEX CHECKS
  ## ----------------------------------------------------------
  
  stopifnot(
    length(keep) == nind_FULL,
    nrow(growth_info_loo) == nind_loo,
    nrow(cintB_loo) == nind_loo,
    nrow(cintD_loo) == nind_loo,
    length(censoredD_loo) == nind_loo,
    length(sex_loo) == nind_loo,
    length(y_loo) == nind_loo,
    length(island_idx_loo) == nind_loo,
    length(island_start_time_loo) == nind_loo,
    length(island_start_idx_loo) == nind_loo,
    !anyNA(island_idx_loo),
    !anyNA(island_start_time_loo),
    !anyNA(island_start_idx_loo),
    all(island_idx_loo >= 1),
    all(island_idx_loo <= n_islands_loo),
    all(island_start_idx_loo >= 1),
    all(island_start_idx_loo <= tMax_FULL)
  )
  
  ## ----------------------------------------------------------
  ## CENTRED SEX
  ## ----------------------------------------------------------
  
  sex_mean_loo <- mean(sex_loo)
  sex_c_loo <- sex_loo - sex_mean_loo
  
  ## ----------------------------------------------------------
  ## FEEDING
  ## ----------------------------------------------------------
  
  feed_island_loo <- feed_island_FULL[
    match(ordered_islands_loo, ordered_islands_FULL)
  ]
  
  if(anyNA(feed_island_loo)) {
    stop("Missing feeding classification when omitting ", omit_island)
  }
  
  feed_island_mean_loo <- mean(feed_island_loo)
  feed_c_island_loo <- feed_island_loo - feed_island_mean_loo
  feed_c_loo <- feed_c_island_loo[island_idx_loo]
  
  cat("\nFeeding structure:\n")
  
  print(
    tibble(
      island = ordered_islands_loo,
      feeding = ifelse(feed_island_loo == 1, "fed", "unfed"),
      feed_c = feed_c_island_loo
    )
  )
  
  stopifnot(
    length(feed_c_loo) == nind_loo,
    !anyNA(feed_c_loo),
    length(unique(feed_island_loo)) == 2
  )
  
  ## ----------------------------------------------------------
  ## ISLAND-LEVEL EFFORT
  ## ----------------------------------------------------------
  
  island_rows <- match(ordered_islands_loo, ordered_islands_FULL)
  
  cum_effort_loo <- cum_effort_FULL[
    island_rows,
    ,
    drop = FALSE
  ]
  
  stopifnot(
    !anyNA(island_rows),
    nrow(cum_effort_loo) == n_islands_loo,
    ncol(cum_effort_loo) >= tMax_FULL + 1,
    !anyNA(cum_effort_loo)
  )
  
  ## ----------------------------------------------------------
  ## DENSITY MATRIX
  ## ----------------------------------------------------------
  
  density_matrix_loo <- density_matrix_FULL[
    island_rows,
    ,
    drop = FALSE
  ]
  
  stopifnot(
    nrow(density_matrix_loo) == n_islands_loo,
    ncol(density_matrix_loo) == tMax_FULL,
    !anyNA(density_matrix_loo)
  )
  
  ## ----------------------------------------------------------
  ## EXTRA INDEX AUDIT BEFORE NIMBLE
  ## ----------------------------------------------------------
  
  cat("\nIndex audit:\n")
  cat("nind =", nind_loo, "\n")
  cat("n_islands =", n_islands_loo, "\n")
  cat("island_idx range =", range(island_idx_loo), "\n")
  cat("island_start_idx range =", range(island_start_idx_loo), "\n")
  cat("cum_effort dimensions =", dim(cum_effort_loo), "\n")
  cat("density_matrix dimensions =", dim(density_matrix_loo), "\n")
  
  ## ----------------------------------------------------------
  ## CONSTANTS + DATA
  ## ----------------------------------------------------------
  
  consts_loo <- list(
    nind = nind_loo,
    tMax = tMax_FULL,
    n_time_blocks_raw = n_time_blocks_raw,
    n_time_blocks = n_time_blocks,
    block_width = block_width,
    n_islands = n_islands_loo,
    cintB = as.matrix(cintB_loo),
    cintD = as.matrix(cintD_loo),
    sex_c = as.numeric(sex_c_loo),
    feed_c = as.numeric(feed_c_loo),
    island_idx = as.numeric(island_idx_loo),
    island_start_time = as.numeric(island_start_time_loo),
    island_start_idx = as.numeric(island_start_idx_loo),
    y = as.numeric(y_loo)
  )
  
  data_loo <- list(
    censoredD = as.numeric(censoredD_loo),
    dind = rep(1, nind_loo),
    cum_effort = as.matrix(cum_effort_loo),
    density_matrix = as.matrix(density_matrix_loo)
  )
  
  ## ----------------------------------------------------------
  ## INITIAL VALUES
  ## ----------------------------------------------------------
  
  make_inits_loo <- function(chain_id = 1) {
    
    set.seed(
      5000 +
        match(omit_island, ordered_islands_FULL) * 100 +
        chain_id
    )
    
    tB_start <- runif(
      nind_loo,
      min = cintB_loo[,1],
      max = cintB_loo[,2]
    )
    
    tD_start <- numeric(nind_loo)
    
    for(i in seq_len(nind_loo)) {
      
      if(censoredD_loo[i] == 1) {
        tD_start[i] <- runif(
          1,
          min = cintD_loo[i,1] + 0.001,
          max = cintD_loo[i,2] - 0.001
        )
      }
      
      if(censoredD_loo[i] == 2) {
        tD_start[i] <- cintD_loo[i,2] + runif(1, 1, 5)
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
      island_a_raw = rnorm(n_islands_loo, 0, 0.5),
      
      sigma_time = runif(1, 1.0, 2.5),
      u_time_raw = rnorm(n_time_blocks, 0, 0.5),
      
      mean.p = runif(1, 0.20, 0.35)
    )
  }
  
  inits_loo <- lapply(seq_len(nchains), make_inits_loo)
  
  ## ----------------------------------------------------------
  ## BUILD MODEL
  ## ----------------------------------------------------------
  
  model_loo <- nimbleModel(
    code = code_LOO,
    constants = consts_loo,
    data = data_loo,
    inits = inits_loo[[1]],
    calculate = FALSE
  )
  
  lp <- model_loo$calculate()
  
  cat("\nInitial log probability:", lp, "\n")
  
  if(!is.finite(lp)) {
    stop("Initial log probability not finite for ", omit_island)
  }
  
  ## ----------------------------------------------------------
  ## MCMC CONFIGURATION
  ## ----------------------------------------------------------
  
  conf_loo <- configureMCMC(
    model_loo,
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
      "sigma_time"
    ),
    enableWAIC = FALSE
  )
  
  ## ----------------------------------------------------------
  ## a-SIDE JOINT SAMPLER
  ## ----------------------------------------------------------
  
  conf_loo$removeSamplers(
    c(
      "log_a",
      "betaSEX[1]",
      "betaDENS_a",
      "betaFEED_a"
    )
  )
  
  conf_loo$addSampler(
    target = c(
      "log_a",
      "betaSEX[1]",
      "betaDENS_a",
      "betaFEED_a"
    ),
    type = "AF_slice"
  )
  
  ## ----------------------------------------------------------
  ## b-SIDE JOINT SAMPLER
  ## ----------------------------------------------------------
  
  conf_loo$removeSamplers(
    c(
      "log_b",
      "betaSEX[2]",
      "betaFEED_b"
    )
  )
  
  conf_loo$addSampler(
    target = c(
      "log_b",
      "betaSEX[2]",
      "betaFEED_b"
    ),
    type = "AF_slice"
  )
  
  ## ----------------------------------------------------------
  ## SCALE PARAMETERS
  ## ----------------------------------------------------------
  
  for(v in c("sigma_island_a", "sigma_time")) {
    conf_loo$removeSamplers(v)
    conf_loo$addSampler(target = v, type = "slice")
  }
  
  ## ----------------------------------------------------------
  ## LATENT BIRTH + LIFESPAN
  ## ----------------------------------------------------------
  
  for(i in seq_len(nind_loo)) {
    
    conf_loo$removeSamplers(
      c(
        paste0("tB[", i, "]"),
        paste0("tstar[", i, "]")
      )
    )
    
    conf_loo$addSampler(
      target = c(
        paste0("tB[", i, "]"),
        paste0("tstar[", i, "]")
      ),
      type = "AF_slice"
    )
  }
  
  ## ----------------------------------------------------------
  ## OPTIONAL FINAL SAMPLER CHECK
  ##
  ## IMPORTANT:
  ## configureMCMC prints the ORIGINAL default RW configuration
  ## before we replace samplers above.
  ##
  ## This check prints the FINAL sampler configuration.
  ##
  ## Set print_samplers = TRUE for the first island only.
  ## ----------------------------------------------------------
  
  if(print_samplers) {
    cat("\nFINAL SAMPLER CONFIGURATION AFTER REPLACEMENTS:\n")
    conf_loo$printSamplers()
  }
  
  ## ----------------------------------------------------------
  ## BUILD + COMPILE
  ## ----------------------------------------------------------
  
  mcmc_loo <- buildMCMC(conf_loo)
  
  cModel_loo <- compileNimble(
    model_loo,
    resetFunctions = TRUE
  )
  
  cMCMC_loo <- compileNimble(
    mcmc_loo,
    project = model_loo
  )
  
  ## ----------------------------------------------------------
  ## RUN
  ## ----------------------------------------------------------
  
  cat(
    "\nRunning LONG LOO model excluding",
    omit_island,
    "with",
    nchains,
    "chains x",
    niter,
    "iterations...\n"
  )
  
  fit <- runMCMC(
    cMCMC_loo,
    niter = niter,
    nburnin = nburnin,
    nchains = nchains,
    thin = thin,
    inits = inits_loo,
    setSeed = TRUE,
    samplesAsCodaMCMC = TRUE,
    summary = TRUE,
    WAIC = FALSE,
    progressBar = TRUE
  )
  
  ## ----------------------------------------------------------
  ## SAVE FULL FIT IMMEDIATELY
  ## ----------------------------------------------------------
  
  filename <- paste0(
    "Model_PaceOfLife_LOO_LONG_",
    omit_island,
    ".rds"
  )
  
  saveRDS(
    fit,
    filename
  )
  
  cat("\nSaved:", filename, "\n")
  
  ## ----------------------------------------------------------
  ## SUMMARISE PRIMARY EFFECTS
  ## ----------------------------------------------------------
  
  s <- as.matrix(fit$samples)
  
  sm <- MCMCsummary(
    fit$samples,
    params = c(
      "betaFEED_a",
      "betaFEED_b",
      "betaDENS_a",
      "betaSEX[2]"
    )
  )
  
  feed_a <- s[, "betaFEED_a"]
  feed_b <- s[, "betaFEED_b"]
  
  result <- tibble(
    omitted_island = omit_island,
    
    omitted_treatment = ifelse(
      feed_island_FULL[
        match(omit_island, ordered_islands_FULL)
      ] == 1,
      "fed",
      "unfed"
    ),
    
    nind = nind_loo,
    
    betaFEED_a_median = median(feed_a),
    betaFEED_a_low = unname(quantile(feed_a, 0.025)),
    betaFEED_a_high = unname(quantile(feed_a, 0.975)),
    Pr_FEED_a_negative = mean(feed_a < 0),
    
    feed_ratio_a_median = median(exp(feed_a)),
    feed_ratio_a_low = unname(quantile(exp(feed_a), 0.025)),
    feed_ratio_a_high = unname(quantile(exp(feed_a), 0.975)),
    
    betaFEED_b_median = median(feed_b),
    betaFEED_b_low = unname(quantile(feed_b, 0.025)),
    betaFEED_b_high = unname(quantile(feed_b, 0.975)),
    Pr_FEED_b_negative = mean(feed_b < 0),
    
    feed_ratio_b_median = median(exp(feed_b)),
    feed_ratio_b_low = unname(quantile(exp(feed_b), 0.025)),
    feed_ratio_b_high = unname(quantile(exp(feed_b), 0.975)),
    
    Rhat_FEED_a = sm["betaFEED_a", "Rhat"],
    Rhat_FEED_b = sm["betaFEED_b", "Rhat"],
    ESS_FEED_a = sm["betaFEED_a", "n.eff"],
    ESS_FEED_b = sm["betaFEED_b", "n.eff"],
    
    Rhat_DENS_a = sm["betaDENS_a", "Rhat"],
    Rhat_SEX_b = sm["betaSEX[2]", "Rhat"]
  )
  
  cat("\nLONG LOO feeding summary excluding", omit_island, ":\n")
  print(result)
  
  ## ----------------------------------------------------------
  ## CLEAN LARGE COMPILED OBJECTS
  ## ----------------------------------------------------------
  
  rm(
    model_loo,
    conf_loo,
    mcmc_loo,
    cModel_loo,
    cMCMC_loo
  )
  
  gc()
  
  return(result)
}

## ============================================================
## 7. OPTIONAL SAMPLER-CHECK RUN
##
## We already know the indexing works from the previous LOO
## analysis, so no 1,000-iteration test is needed.
##
## If you want to inspect the FINAL sampler configuration before
## launching all 8 fits, run this manually:
##
## sampler_check <- fit_LOO_island(
##   omit_island = "Bitter",
##   niter = 1000,
##   nburnin = 200,
##   thin = 1,
##   nchains = 1,
##   print_samplers = TRUE
## )
##
## Do NOT run this section unless you want the sampler printout.
## ============================================================


## ============================================================
## 8. RUN ALL EIGHT LONG LOO MODELS
##
## 500k iterations
## 100k burn-in
## 3 chains
## thin 10
##
## print_samplers = TRUE for Bitter only so we can verify the
## final sampler configuration without printing thousands of
## sampler lines eight times.
## ============================================================

loo_results_LONG <- vector(
  "list",
  length(ordered_islands_FULL)
)

names(loo_results_LONG) <- ordered_islands_FULL

for(j in seq_along(ordered_islands_FULL)) {
  
  omit_island <- ordered_islands_FULL[j]
  
  cat("\n\n")
  cat("############################################################\n")
  cat("LONG LOO MODEL", j, "OF", length(ordered_islands_FULL), "\n")
  cat("OMITTING:", omit_island, "\n")
  cat("############################################################\n")
  
  loo_results_LONG[[j]] <- tryCatch(
    
    fit_LOO_island(
      omit_island = omit_island,
      niter = 500000,
      nburnin = 100000,
      thin = 10,
      nchains = 3,
      print_samplers = j == 1
    ),
    
    error = function(e) {
      
      cat(
        "\nERROR FOR",
        omit_island,
        ":\n",
        conditionMessage(e),
        "\n"
      )
      
      tibble(
        omitted_island = omit_island,
        error = conditionMessage(e)
      )
    }
  )
  
  ## Save progress after every island
  saveRDS(
    loo_results_LONG,
    "Model_PaceOfLife_LOO_LONG_PROGRESS.rds"
  )
}

## ============================================================
## 9. COMBINE LONG LOO RESULTS
## ============================================================

loo_summary_LONG <- bind_rows(
  loo_results_LONG
)

cat("\n")
cat("============================================================\n")
cat("FINAL LONG LOO SUMMARY\n")
cat("============================================================\n")

print(
  loo_summary_LONG,
  n = Inf
)

write_csv(
  loo_summary_LONG,
  "Model_PaceOfLife_LOO_LONG_summary.csv"
)

saveRDS(
  loo_summary_LONG,
  "Model_PaceOfLife_LOO_LONG_summary.rds"
)

## ============================================================
## 10. ADD FULL MODEL FOR COMPARISON
## ============================================================

if(!exists("run_POL") && file.exists("Model_PaceOfLife_FINAL_samples.rds")) {
  run_POL <- readRDS(
    "Model_PaceOfLife_FINAL_samples.rds"
  )
}

if(exists("run_POL")) {
  
  full_samps <- as.matrix(
    run_POL$samples
  )
  
  full_feed_a <- full_samps[, "betaFEED_a"]
  full_feed_b <- full_samps[, "betaFEED_b"]
  
  full_sm <- MCMCsummary(
    run_POL$samples,
    params = c(
      "betaFEED_a",
      "betaFEED_b"
    )
  )
  
  full_result <- tibble(
    omitted_island = "FULL MODEL",
    omitted_treatment = "-",
    nind = nind_FULL,
    
    betaFEED_a_median = median(full_feed_a),
    betaFEED_a_low = unname(quantile(full_feed_a, 0.025)),
    betaFEED_a_high = unname(quantile(full_feed_a, 0.975)),
    Pr_FEED_a_negative = mean(full_feed_a < 0),
    
    feed_ratio_a_median = median(exp(full_feed_a)),
    feed_ratio_a_low = unname(quantile(exp(full_feed_a), 0.025)),
    feed_ratio_a_high = unname(quantile(exp(full_feed_a), 0.975)),
    
    betaFEED_b_median = median(full_feed_b),
    betaFEED_b_low = unname(quantile(full_feed_b, 0.025)),
    betaFEED_b_high = unname(quantile(full_feed_b, 0.975)),
    Pr_FEED_b_negative = mean(full_feed_b < 0),
    
    feed_ratio_b_median = median(exp(full_feed_b)),
    feed_ratio_b_low = unname(quantile(exp(full_feed_b), 0.025)),
    feed_ratio_b_high = unname(quantile(exp(full_feed_b), 0.975)),
    
    Rhat_FEED_a = full_sm["betaFEED_a", "Rhat"],
    Rhat_FEED_b = full_sm["betaFEED_b", "Rhat"],
    ESS_FEED_a = full_sm["betaFEED_a", "n.eff"],
    ESS_FEED_b = full_sm["betaFEED_b", "n.eff"],
    
    Rhat_DENS_a = NA_real_,
    Rhat_SEX_b = NA_real_
  )
  
  loo_with_full_LONG <- bind_rows(
    full_result,
    loo_summary_LONG
  )
  
  cat("\nLONG LOO + FULL MODEL COMPARISON:\n")
  
  print(
    loo_with_full_LONG,
    n = Inf
  )
  
  write_csv(
    loo_with_full_LONG,
    "Model_PaceOfLife_LOO_LONG_with_FULL_summary.csv"
  )
  
  saveRDS(
    loo_with_full_LONG,
    "Model_PaceOfLife_LOO_LONG_with_FULL_summary.rds"
  )
}

## ============================================================
## 11. FEEDING-b ROBUSTNESS TABLE
## ============================================================

cat("\n")
cat("============================================================\n")
cat("FEEDING EFFECT ON b ACROSS LONG LOO MODELS\n")
cat("============================================================\n")

feeding_b_summary_LONG <- loo_summary_LONG %>%
  select(
    omitted_island,
    omitted_treatment,
    nind,
    betaFEED_b_median,
    betaFEED_b_low,
    betaFEED_b_high,
    Pr_FEED_b_negative,
    feed_ratio_b_median,
    feed_ratio_b_low,
    feed_ratio_b_high,
    Rhat_FEED_b,
    ESS_FEED_b
  ) %>%
  arrange(betaFEED_b_median)

print(
  feeding_b_summary_LONG,
  n = Inf
)

## ============================================================
## 12. CONVERGENCE CHECK
## ============================================================

cat("\n")
cat("============================================================\n")
cat("LOO betaFEED_b CONVERGENCE\n")
cat("============================================================\n")

convergence_b_LONG <- loo_summary_LONG %>%
  select(
    omitted_island,
    Rhat_FEED_b,
    ESS_FEED_b
  ) %>%
  arrange(desc(Rhat_FEED_b))

print(
  convergence_b_LONG,
  n = Inf
)

cat(
  "\nMaximum betaFEED_b Rhat:",
  max(loo_summary_LONG$Rhat_FEED_b, na.rm = TRUE),
  "\n"
)

cat(
  "Minimum betaFEED_b ESS:",
  min(loo_summary_LONG$ESS_FEED_b, na.rm = TRUE),
  "\n"
)

## ============================================================
## 13. ROBUSTNESS SUMMARY
## ============================================================

cat(
  "\nNumber of LOO fits with median betaFEED_b < 0:",
  sum(loo_summary_LONG$betaFEED_b_median < 0, na.rm = TRUE),
  "of",
  nrow(loo_summary_LONG),
  "\n"
)

cat(
  "Range of median fed/unfed b ratios:",
  paste(
    round(
      range(
        loo_summary_LONG$feed_ratio_b_median,
        na.rm = TRUE
      ),
      3
    ),
    collapse = " - "
  ),
  "\n"
)

cat(
  "Range of P(betaFEED_b < 0):",
  paste(
    round(
      range(
        loo_summary_LONG$Pr_FEED_b_negative,
        na.rm = TRUE
      ),
      3
    ),
    collapse = " - "
  ),
  "\n"
)

## ============================================================
## 14. SAVE FINAL ROBUSTNESS TABLE
## ============================================================

write_csv(
  feeding_b_summary_LONG,
  "Model_PaceOfLife_LOO_LONG_feeding_b_robustness.csv"
)

saveRDS(
  feeding_b_summary_LONG,
  "Model_PaceOfLife_LOO_LONG_feeding_b_robustness.rds"
)