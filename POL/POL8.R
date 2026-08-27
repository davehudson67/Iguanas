## ============================================================
## MODEL_PaceOfLife_REFAGE10_FEEDREF_100k.R
##
## COMPLETE SETUP + MODEL + RUN
##
## h(t) = a10 * exp[b * (t - 10)]
##
## a10 = mortality hazard at age 10
## b   = actuarial senescence rate
##
## FOCAL BIOLOGICAL EFFECTS ON BOTH:
##   sex
##   cohort density
##   feeding
##
## Feeding reference:
##   fed   =  0
##   unfed = -1
##
## Therefore:
##   betaFEED_a = log(a10_fed) - log(a10_unfed)
##   betaFEED_b = log(b_fed)   - log(b_unfed)
##
## Hierarchical controls on a10:
##   island random effect
##   5-year cohort temporal effect
##
## TEST RUN:
##   niter   = 100000
##   burnin  = 30000
##   chains  = 3
##   thin    = 10
## ============================================================


## ============================================================
## 1. CANONICAL DATA PREP
## ============================================================

source("DataPrep.R")

source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

library(tidyverse)
library(nimble)
library(MCMCvis)
library(coda)


cat("\n")
cat("============================================================\n")
cat("REBUILDING CLEAN BIRTH-TIME SETUP\n")
cat("============================================================\n")


## ============================================================
## 2. SAVE ORIGINAL BIRTH WINDOWS
##
## DataPrep.R currently contains the old blanket
## birth/death-censor adjustment.
##
## We rebuild the intended biological birth windows below and
## overwrite cintB_adj with the corrected version.
## ============================================================

cintB_adj_original <- cintB_adj

stopifnot(
  nrow(growth_info) == nind,
  nrow(cintB_adj_original) == nind,
  nrow(cintD_adj) == nind
)


## ============================================================
## 3. REBUILD CALENDAR-SCALE BIRTH WINDOWS
## ============================================================

cintB_cal_clean <- matrix(
  NA_real_,
  nrow = nind,
  ncol = 2
)

birth_source <- character(nind)
birth_svl_capture_conflict <- logical(nind)


for(i in seq_len(nind)) {
  
  ## ----------------------------------------------------------
  ## KNOWN AGE
  ## ----------------------------------------------------------
  
  if(!is.na(growth_info$known_birth_cal[i])) {
    
    cintB_cal_clean[i, ] <-
      growth_info$known_birth_cal[i] +
      c(-1, 0)
    
    birth_source[i] <- "Known birth"
    
  } else {
    
    ## --------------------------------------------------------
    ## UNKNOWN AGE
    ## --------------------------------------------------------
    
    has_svl <-
      is.finite(growth_info$first_svl[i]) &&
      is.finite(growth_info$first_svl_year[i])
    
    if(has_svl) {
      
      ref_yr <- growth_info$first_svl_year[i]
      
      raw_low <- floor(
        ref_yr -
          growth_info$pred_max_age[i]
      )
      
      raw_high <- floor(
        ref_yr -
          growth_info$pred_min_age[i]
      )
      
      birth_svl_capture_conflict[i] <-
        raw_high >
        growth_info$tF_cal[i]
      
      birth_source[i] <- if(
        birth_svl_capture_conflict[i]
      ) {
        "SVL informed - constrained by first capture"
      } else {
        "SVL informed"
      }
      
    } else {
      
      ## No known age and no SVL:
      ## retain broad 0-40 year age uncertainty.
      
      ref_yr <- growth_info$tF_cal[i]
      
      raw_low <- floor(
        ref_yr -
          growth_info$pred_max_age[i]
      )
      
      raw_high <- floor(
        ref_yr -
          growth_info$pred_min_age[i]
      )
      
      birth_source[i] <-
        "Unknown birth - no SVL"
    }
    
    
    ## --------------------------------------------------------
    ## FIRST CAPTURE IS HARD EVIDENCE OF BEING ALIVE
    ## --------------------------------------------------------
    
    final_high <- min(
      raw_high,
      growth_info$tF_cal[i]
    )
    
    final_low <- raw_low
    
    ## If the SVL-derived interval lies entirely after first
    ## capture, use a one-year interval ending at first capture.
    if(final_low >= final_high) {
      final_low <- final_high - 1
    }
    
    cintB_cal_clean[i, ] <- c(
      final_low,
      final_high
    )
  }
  
  
  if(any(!is.finite(cintB_cal_clean[i, ]))) {
    stop(
      "Invalid rebuilt birth interval for ",
      growth_info$ID_raw[i]
    )
  }
  
  if(
    cintB_cal_clean[i, 1] >=
    cintB_cal_clean[i, 2]
  ) {
    stop(
      "Non-positive rebuilt birth interval for ",
      growth_info$ID_raw[i]
    )
  }
}


## ============================================================
## 4. CHECK TIME ORIGIN
## ============================================================

Origin_rebuilt <- min(
  cintB_cal_clean[, 1],
  na.rm = TRUE
)

cat(
  "\nExisting Origin:",
  Origin,
  "\n"
)

cat(
  "Rebuilt Origin:",
  Origin_rebuilt,
  "\n"
)

if(Origin_rebuilt != Origin) {
  stop(
    "Rebuilt birth windows alter Origin. Inspect before proceeding."
  )
}


## ============================================================
## 5. CONVERT CLEAN BIRTH WINDOWS TO MODEL TIME
## ============================================================

cintB_adj_clean <-
  cintB_cal_clean -
  Origin +
  1


for(i in seq_len(nind)) {
  
  tF_adj_i <-
    growth_info$tF_cal[i] -
    Origin +
    1
  
  cintB_adj_clean[i, 2] <- min(
    cintB_adj_clean[i, 2],
    tF_adj_i
  )
  
  cintB_adj_clean[i, 1] <- max(
    1,
    cintB_adj_clean[i, 1]
  )
  
  if(
    cintB_adj_clean[i, 1] >=
    cintB_adj_clean[i, 2]
  ) {
    
    cintB_adj_clean[i, 1] <- max(
      1,
      cintB_adj_clean[i, 2] - 1
    )
  }
  
  if(
    cintB_adj_clean[i, 1] >=
    cintB_adj_clean[i, 2]
  ) {
    stop(
      "Invalid final birth interval for ",
      growth_info$ID_raw[i]
    )
  }
}


## ============================================================
## 6. AUDIT DIFFERENCE FROM OLD DataPrep.R WINDOWS
## ============================================================

birth_changed <- which(
  rowSums(
    abs(
      cintB_adj_clean -
        cintB_adj_original
    ) > 1e-8
  ) > 0
)


birth_change_audit <- tibble()

if(length(birth_changed) > 0) {
  
  birth_change_audit <- tibble(
    i = birth_changed,
    ID = growth_info$ID_raw[birth_changed],
    island = growth_info$island[birth_changed],
    censor = censoredD_vec[birth_changed],
    
    old_low =
      cintB_adj_original[birth_changed, 1] +
      Origin - 1,
    
    old_high =
      cintB_adj_original[birth_changed, 2] +
      Origin - 1,
    
    clean_low =
      cintB_adj_clean[birth_changed, 1] +
      Origin - 1,
    
    clean_high =
      cintB_adj_clean[birth_changed, 2] +
      Origin - 1
  )
}


cat(
  "\nBirth windows changed after removing old correction:",
  length(birth_changed),
  "\n"
)

if(nrow(birth_change_audit) > 0) {
  
  cat("Changes by censoring category:\n")
  
  birth_change_audit %>%
    count(censor) %>%
    print()
}


## ============================================================
## 7. USE CLEAN BIRTH WINDOWS
## ============================================================

cintB_adj <- cintB_adj_clean


growth_info <- growth_info %>%
  mutate(
    birth_source = birth_source,
    
    birth_svl_capture_conflict =
      birth_svl_capture_conflict,
    
    birth_no_svl =
      is.na(known_birth_cal) &
      !is.finite(first_svl),
    
    birth_weak_or_adjusted =
      birth_no_svl |
      birth_svl_capture_conflict,
    
    birth_lower_cal =
      cintB_adj[, 1] +
      Origin - 1,
    
    birth_upper_cal =
      cintB_adj[, 2] +
      Origin - 1,
    
    birth_window_width =
      cintB_adj[, 2] -
      cintB_adj[, 1]
  )


cat("\nBirth-source summary:\n")

growth_info %>%
  count(birth_source) %>%
  print()


cat(
  "\nWeak/adjusted birth histories:",
  sum(growth_info$birth_weak_or_adjusted),
  "\n"
)


## ============================================================
## 8. BIRTH-WINDOW VALIDITY
## ============================================================

stopifnot(
  all(cintB_adj[, 1] < cintB_adj[, 2]),
  all(cintB_adj[, 1] >= 1),
  
  all(
    growth_info$birth_upper_cal <=
      growth_info$tF_cal
  )
)


## ============================================================
## 9. CHECK TRUE BIRTH/DEATH IMPOSSIBILITIES
##
## We do NOT alter birth windows for right-censored individuals.
##
## For interval-censored deaths, however, an entire birth
## interval occurring after death would be impossible.
## ============================================================

impossible_interval_birth <- which(
  censoredD_vec == 1 &
    cintB_adj[, 1] >=
    cintD_adj[, 2]
)

if(length(impossible_interval_birth) > 0) {
  
  print(
    tibble(
      ID =
        growth_info$ID_raw[
          impossible_interval_birth
        ],
      
      island =
        growth_info$island[
          impossible_interval_birth
        ],
      
      birth_low =
        cintB_adj[
          impossible_interval_birth,
          1
        ],
      
      birth_high =
        cintB_adj[
          impossible_interval_birth,
          2
        ],
      
      death_low =
        cintD_adj[
          impossible_interval_birth,
          1
        ],
      
      death_high =
        cintD_adj[
          impossible_interval_birth,
          2
        ]
    ),
    n = Inf
  )
  
  stop(
    "Impossible interval-censored birth/death history detected."
  )
}


birth_censor_overlap <- which(
  cintB_adj[, 2] >=
    cintD_adj[, 2]
)

cat(
  "\nBirth upper >= death/censor upper:",
  length(birth_censor_overlap),
  "\n"
)

if(length(birth_censor_overlap) > 0) {
  
  tibble(
    censor =
      censoredD_vec[
        birth_censor_overlap
      ]
  ) %>%
    count(censor) %>%
    print()
}


## ============================================================
## 10. REFERENCE AGE
## ============================================================

age_ref <- 10

cat(
  "\nReference age:",
  age_ref,
  "years\n"
)


## ============================================================
## 11. CENTRED SEX
##
## 0 = female
## 1 = male
## ============================================================

stopifnot(
  length(sex) == nind,
  !anyNA(sex),
  all(sex %in% c(0, 1))
)

sex_mean <- mean(sex)

sex_c <-
  sex -
  sex_mean


cat(
  "Sex mean:",
  sex_mean,
  "\n"
)


## ============================================================
## 12. FEEDING BY ISLAND
##
## Fed islands are REFERENCE:
##
## fed   =  0
## unfed = -1
##
## betaFEED therefore represents:
##
##   log(fed parameter) - log(unfed parameter)
## ============================================================

feeding_by_island <- growth_info %>%
  mutate(
    feeding =
      tolower(
        trimws(feeding)
      )
  ) %>%
  group_by(island) %>%
  summarise(
    feed =
      as.integer(
        any(
          feeding != "none",
          na.rm = TRUE
        )
      ),
    .groups = "drop"
  )


feed_island <-
  feeding_by_island$feed[
    match(
      ordered_islands,
      feeding_by_island$island
    )
  ]

stopifnot(
  !anyNA(feed_island),
  all(feed_island %in% c(0, 1))
)


feed_ref_island <- ifelse(
  feed_island == 1,
  0,
  -1
)

feed_ref <-
  feed_ref_island[
    island_idx
  ]


stopifnot(
  length(feed_ref) == nind,
  !anyNA(feed_ref),
  all(feed_ref %in% c(-1, 0))
)


cat("\nFeeding reference coding:\n")

tibble(
  island = ordered_islands,
  
  feeding = ifelse(
    feed_island == 1,
    "fed",
    "unfed"
  ),
  
  feed_ref =
    feed_ref_island
) %>%
  print(n = Inf)


## ============================================================
## 13. COHORT DENSITY
##
## MNA per usable hectare.
##
## Observed survey years:
##   annual MNA-derived density.
##
## Pre-monitoring/missing cohort years:
##   island mean observed density.
##
## Transform:
##   log1p()
##   global z-standardisation.
## ============================================================

habitat <- read_csv(
  "Data/Iguana_islands_habitat.csv",
  show_col_types = FALSE
) %>%
  transmute(
    island = islands,
    usable_m2 = Useable
  )


first_last <- growth_info %>%
  transmute(
    ID,
    island,
    
    first_year =
      tF_cal -
      Origin +
      1,
    
    last_year =
      tL_cal -
      Origin +
      1
  )


mna_long <- first_last %>%
  rowwise() %>%
  mutate(
    year_idx = list(
      seq(
        first_year,
        last_year
      )
    )
  ) %>%
  unnest(year_idx) %>%
  ungroup() %>%
  count(
    island,
    year_idx,
    name = "MNA"
  ) %>%
  left_join(
    habitat,
    by = "island"
  ) %>%
  mutate(
    density_raw =
      MNA /
      (usable_m2 / 10000)
  )


survey_long <- actual_surveys %>%
  transmute(
    island = Island,
    year_idx = yr_adj
  ) %>%
  distinct(
    island,
    year_idx
  )


density_observed <- survey_long %>%
  left_join(
    mna_long %>%
      select(
        island,
        year_idx,
        MNA,
        density_raw
      ),
    by = c(
      "island",
      "year_idx"
    )
  ) %>%
  mutate(
    MNA =
      replace_na(
        MNA,
        0L
      ),
    
    density_raw =
      replace_na(
        density_raw,
        0
      )
  )


island_mean_density <- density_observed %>%
  group_by(island) %>%
  summarise(
    mean_density_raw =
      mean(
        density_raw,
        na.rm = TRUE
      ),
    .groups = "drop"
  )


density_complete <- expand_grid(
  island = ordered_islands,
  year_idx = seq_len(tMax)
) %>%
  left_join(
    density_observed %>%
      select(
        island,
        year_idx,
        density_raw
      ),
    by = c(
      "island",
      "year_idx"
    )
  ) %>%
  left_join(
    island_mean_density,
    by = "island"
  ) %>%
  mutate(
    density_source =
      if_else(
        is.na(density_raw),
        "Island mean",
        "Observed cohort"
      ),
    
    density_filled_raw =
      coalesce(
        density_raw,
        mean_density_raw
      )
  )


if(anyNA(density_complete$density_filled_raw)) {
  stop(
    "Missing density values remain."
  )
}


dens_mu <- mean(
  log1p(
    density_complete$density_filled_raw
  )
)

dens_sd <- sd(
  log1p(
    density_complete$density_filled_raw
  )
)

if(
  !is.finite(dens_sd) ||
  dens_sd <= 0
) {
  stop(
    "Invalid density standard deviation."
  )
}


density_complete <- density_complete %>%
  mutate(
    density_z =
      (
        log1p(
          density_filled_raw
        ) -
          dens_mu
      ) /
      dens_sd
  )


density_matrix <- matrix(
  NA_real_,
  nrow = n_islands,
  ncol = tMax,
  dimnames = list(
    ordered_islands,
    seq_len(tMax)
  )
)


for(r in seq_len(nrow(density_complete))) {
  
  jj <- match(
    density_complete$island[r],
    ordered_islands
  )
  
  tt <-
    density_complete$year_idx[r]
  
  density_matrix[jj, tt] <-
    density_complete$density_z[r]
}


stopifnot(
  !anyNA(density_matrix)
)


cat("\nDensity matrix summary:\n")

print(
  summary(
    as.vector(
      density_matrix
    )
  )
)


## ============================================================
## 14. TEMPORAL COHORT BLOCKS
##
## 5-year independent cohort blocks.
##
## Raw final two sparse blocks are pooled.
## ============================================================

block_width <- 5

n_time_blocks_raw <- ceiling(
  tMax /
    block_width
)

n_time_blocks <-
  n_time_blocks_raw -
  1


cat(
  "\nTemporal structure:",
  n_time_blocks_raw,
  "raw blocks ->",
  n_time_blocks,
  "model blocks\n"
)


## ============================================================
## 15. FINAL SETUP CHECK
## ============================================================

stopifnot(
  nrow(growth_info) == nind,
  
  nrow(cintB_adj) == nind,
  ncol(cintB_adj) == 2,
  
  nrow(cintD_adj) == nind,
  ncol(cintD_adj) == 2,
  
  length(censoredD_vec) == nind,
  length(sex_c) == nind,
  length(feed_ref) == nind,
  
  length(island_idx) == nind,
  length(island_start_time) == nind,
  length(island_start_idx) == nind,
  
  nrow(density_matrix) == n_islands,
  ncol(density_matrix) == tMax,
  
  nrow(cum_effort_for_model) == n_islands,
  ncol(cum_effort_for_model) >= tMax + 1,
  
  all(cintB_adj[, 1] < cintB_adj[, 2]),
  all(cintB_adj[, 1] >= 1),
  
  all(island_idx >= 1),
  all(island_idx <= n_islands)
)


cat("\n")
cat("============================================================\n")
cat("SETUP COMPLETE\n")
cat("============================================================\n")
cat("Individuals:", nind, "\n")
cat("Islands:", n_islands, "\n")
cat("Reference age:", age_ref, "\n")
cat(
  "Weak/adjusted birth histories:",
  sum(growth_info$birth_weak_or_adjusted),
  "\n"
)
cat("============================================================\n")


## ============================================================
## 16. NIMBLE MODEL
## ============================================================

code_REF10_FR <- nimbleCode({
  
  for(i in 1:nind) {
    
    ## --------------------------------------------------------
    ## LATENT BIRTH
    ## --------------------------------------------------------
    
    tB[i] ~ dunif(
      cintB[i, 1],
      cintB[i, 2]
    )
    
    
    ## --------------------------------------------------------
    ## LEFT TRUNCATION
    ##
    ## Age survived before monitoring began.
    ## --------------------------------------------------------
    
    L[i] <- max(
      0,
      island_start_time[i] -
        tB[i]
    )
    
    
    ## --------------------------------------------------------
    ## LATENT BIRTH COHORT
    ## --------------------------------------------------------
    
    cohort_idx[i] <- max(
      1,
      min(
        tMax,
        floor(tB[i]) + 1
      )
    )
    
    
    time_block_raw[i] <- max(
      1,
      min(
        n_time_blocks_raw,
        floor(
          (cohort_idx[i] - 1) /
            block_width
        ) + 1
      )
    )
    
    
    ## Pool final sparse raw block with previous block.
    time_block[i] <- min(
      time_block_raw[i],
      n_time_blocks
    )
    
    
    dens_cohort[i] <-
      density_matrix[
        island_idx[i],
        cohort_idx[i]
      ]
    
    
    ## ========================================================
    ## MORTALITY LEVEL AT AGE 10
    ##
    ## log_a10_fed:
    ## mortality at age 10 for FED reference group,
    ## centred sex/density and zero random effects.
    ##
    ## betaFEED_a:
    ## log(a10_fed) - log(a10_unfed)
    ## ========================================================
    
    log_a10_ind[i] <-
      log_a10_fed +
      betaSEX[1] * sex_c[i] +
      betaDENS_a * dens_cohort[i] +
      betaFEED_a * feed_ref[i] +
      u_time[time_block[i]] +
      island_a[island_idx[i]]
    
    
    a10_ind[i] <- exp(
      log_a10_ind[i]
    )
    
    
    ## ========================================================
    ## ACTUARIAL SENESCENCE
    ##
    ## log_b_fed:
    ## senescence rate for FED reference group,
    ## centred sex/density.
    ##
    ## betaFEED_b:
    ## log(b_fed) - log(b_unfed)
    ## ========================================================
    
    log_b_ind[i] <-
      log_b_fed +
      betaSEX[2] * sex_c[i] +
      betaDENS_b * dens_cohort[i] +
      betaFEED_b * feed_ref[i]
    
    
    bmult[i] <- exp(
      log_b_ind[i]
    )
    
    
    ## ========================================================
    ## CONVERT AGE-10 PARAMETER TO CONVENTIONAL GOMPERTZ a0
    ##
    ## h(t) = a10 * exp[b(t - 10)]
    ##
    ## therefore:
    ##
    ## a0 = a10 * exp(-b*10)
    ## ========================================================
    
    log_amult[i] <-
      log_a10_ind[i] -
      bmult[i] * age_ref
    
    
    amult[i] <- exp(
      log_amult[i]
    )
    
    
    ## --------------------------------------------------------
    ## LATENT LIFESPAN
    ## --------------------------------------------------------
    
    tstar[i] ~ dGompertzLB(
      amult[i],
      bmult[i],
      lowerBound = L[i]
    )
    
    
    tD[i] <-
      tB[i] +
      tstar[i]
    
    
    censoredD[i] ~ dinterval(
      tD[i],
      cintD[i, 1:2]
    )
    
    
    ## ========================================================
    ## DETECTION
    ## ========================================================
    
    tB_year_idx[i] <- max(
      1,
      min(
        tMax,
        floor(tB[i]) + 1
      )
    )
    
    
    nm_start[i] <- max(
      tB_year_idx[i],
      island_start_idx[i]
    )
    
    
    end_year[i] <- max(
      1,
      min(
        tMax,
        floor(tD[i]) + 1
      )
    )
    
    
    nMpos[i] <-
      cum_effort[
        island_idx[i],
        end_year[i] + 1
      ] -
      cum_effort[
        island_idx[i],
        nm_start[i]
      ]
    
    
    nMpos_safe[i] <- max(
      y[i],
      nMpos[i]
    )
    
    
    pd_raw[i] <- exp(
      y[i] *
        log(
          mean.p +
            1.0E-10
        ) +
        (
          nMpos_safe[i] -
            y[i]
        ) *
        log(
          1 -
            mean.p +
            1.0E-10
        )
    )
    
    
    pd[i] <- min(
      0.999999,
      pd_raw[i]
    )
    
    
    dind[i] ~ dbern(
      pd[i]
    )
  }
  
  
  ## ==========================================================
  ## 17. FOCAL BIOLOGICAL EFFECTS
  ## ==========================================================
  
  ## Sex
  ## [1] mortality at age 10
  ## [2] senescence
  
  for(k in 1:2) {
    betaSEX[k] ~ dnorm(
      0,
      sd = 1
    )
  }
  
  
  ## Density
  
  betaDENS_a ~ dnorm(
    0,
    sd = 1
  )
  
  betaDENS_b ~ dnorm(
    0,
    sd = 1
  )
  
  
  ## Feeding
  
  betaFEED_a ~ dnorm(
    0,
    sd = 1
  )
  
  betaFEED_b ~ dnorm(
    0,
    sd = 1
  )
  
  
  ## ==========================================================
  ## 18. ISLAND HETEROGENEITY IN MORTALITY AT AGE 10
  ##
  ## Non-centred + exact zero centring.
  ## ==========================================================
  
  sigma_island_a ~ dunif(
    0,
    2
  )
  
  
  for(r in 1:n_islands) {
    
    island_a_raw[r] ~ dnorm(
      0,
      sd = 1
    )
  }
  
  
  mean_island_a_raw <-
    sum(
      island_a_raw[1:n_islands]
    ) /
    n_islands
  
  
  for(r in 1:n_islands) {
    
    island_a[r] <-
      sigma_island_a *
      (
        island_a_raw[r] -
          mean_island_a_raw
      )
  }
  
  
  ## ==========================================================
  ## 19. TEMPORAL COHORT HETEROGENEITY IN a10
  ## ==========================================================
  
  sigma_time ~ dunif(
    0,
    5
  )
  
  
  for(tt in 1:n_time_blocks) {
    
    u_time_raw[tt] ~ dnorm(
      0,
      sd = 1
    )
  }
  
  
  mean_u_time_raw <-
    sum(
      u_time_raw[1:n_time_blocks]
    ) /
    n_time_blocks
  
  
  for(tt in 1:n_time_blocks) {
    
    u_time[tt] <-
      sigma_time *
      (
        u_time_raw[tt] -
          mean_u_time_raw
      )
  }
  
  
  ## ==========================================================
  ## 20. FED-REFERENCE PARAMETERS
  ## ==========================================================
  
  log_a10_fed ~ dnorm(
    log(0.015),
    sd = 2
  )
  
  
  log_b_fed ~ dnorm(
    log(0.20),
    sd = 1
  )
  
  
  a10_fed <- exp(
    log_a10_fed
  )
  
  
  b_fed <- exp(
    log_b_fed
  )
  
  
  ## ----------------------------------------------------------
  ## DERIVED UNFED PARAMETERS
  ##
  ## Because feed_ref = -1 for unfed:
  ##
  ## log(parameter_unfed) =
  ## log(parameter_fed) - betaFEED
  ## ----------------------------------------------------------
  
  log_a10_unfed <-
    log_a10_fed -
    betaFEED_a
  
  
  log_b_unfed <-
    log_b_fed -
    betaFEED_b
  
  
  a10_unfed <- exp(
    log_a10_unfed
  )
  
  
  b_unfed <- exp(
    log_b_unfed
  )
  
  
  ## ==========================================================
  ## 21. DETECTION
  ## ==========================================================
  
  mean.p ~ dunif(
    0,
    1
  )
})


## ============================================================
## 22. CONSTANTS
## ============================================================

consts_REF10_FR <- list(
  nind = nind,
  tMax = tMax,
  
  n_time_blocks_raw =
    n_time_blocks_raw,
  
  n_time_blocks =
    n_time_blocks,
  
  block_width =
    block_width,
  
  n_islands =
    n_islands,
  
  age_ref =
    age_ref,
  
  cintB =
    as.matrix(
      cintB_adj
    ),
  
  cintD =
    as.matrix(
      cintD_adj
    ),
  
  sex_c =
    as.numeric(
      sex_c
    ),
  
  feed_ref =
    as.numeric(
      feed_ref
    ),
  
  island_idx =
    as.numeric(
      island_idx
    ),
  
  island_start_time =
    as.numeric(
      island_start_time
    ),
  
  island_start_idx =
    as.numeric(
      island_start_idx
    ),
  
  y =
    as.numeric(
      y
    )
)


## ============================================================
## 23. DATA
##
## density_matrix and cum_effort stay DATA because of dynamic
## indexing in NIMBLE.
## ============================================================

data_REF10_FR <- list(
  censoredD =
    as.numeric(
      censoredD_vec
    ),
  
  dind =
    rep(
      1,
      nind
    ),
  
  cum_effort =
    as.matrix(
      cum_effort_for_model
    ),
  
  density_matrix =
    as.matrix(
      density_matrix
    )
)


## ============================================================
## 24. CHAIN-SPECIFIC INITIAL VALUES
## ============================================================

make_inits_REF10_FR <- function(chain_id = 1) {
  
  set.seed(
    9100 +
      chain_id
  )
  
  
  ## ----------------------------------------------------------
  ## BIRTH
  ## ----------------------------------------------------------
  
  tB_start <- runif(
    nind,
    min = cintB_adj[, 1],
    max = cintB_adj[, 2]
  )
  
  
  tD_start <- numeric(
    nind
  )
  
  
  ## ----------------------------------------------------------
  ## DEATH
  ## ----------------------------------------------------------
  
  for(i in seq_len(nind)) {
    
    ## --------------------------------------------------------
    ## INTERVAL-CENSORED DEATH
    ## --------------------------------------------------------
    
    if(censoredD_vec[i] == 1) {
      
      loD <-
        cintD_adj[i, 1] +
        0.001
      
      hiD <-
        cintD_adj[i, 2] -
        0.001
      
      
      if(hiD <= loD) {
        stop(
          "Invalid death interval for ",
          growth_info$ID_raw[i]
        )
      }
      
      
      if(tB_start[i] >= hiD) {
        
        possible_birth_high <- min(
          cintB_adj[i, 2],
          hiD - 0.01
        )
        
        
        if(
          possible_birth_high <=
          cintB_adj[i, 1]
        ) {
          stop(
            "Cannot initialise compatible birth/death for ",
            growth_info$ID_raw[i]
          )
        }
        
        
        tB_start[i] <- runif(
          1,
          min = cintB_adj[i, 1],
          max = possible_birth_high
        )
      }
      
      
      loD_use <- max(
        loD,
        tB_start[i] + 0.01
      )
      
      
      if(loD_use >= hiD) {
        stop(
          "No valid interval death initial value for ",
          growth_info$ID_raw[i]
        )
      }
      
      
      tD_start[i] <- runif(
        1,
        min = loD_use,
        max = hiD
      )
    }
    
    
    ## --------------------------------------------------------
    ## RIGHT-CENSORED DEATH
    ## --------------------------------------------------------
    
    if(censoredD_vec[i] == 2) {
      
      tD_start[i] <- max(
        cintD_adj[i, 2] +
          runif(
            1,
            1,
            5
          ),
        
        tB_start[i] +
          runif(
            1,
            1,
            5
          )
      )
    }
  }
  
  
  list(
    tB =
      tB_start,
    
    tstar =
      pmax(
        0.001,
        tD_start -
          tB_start
      ),
    
    log_a10_fed =
      rnorm(
        1,
        log(0.015),
        0.20
      ),
    
    log_b_fed =
      rnorm(
        1,
        log(0.205),
        0.10
      ),
    
    betaSEX =
      rnorm(
        2,
        0,
        0.10
      ),
    
    betaDENS_a =
      rnorm(
        1,
        0,
        0.10
      ),
    
    betaDENS_b =
      rnorm(
        1,
        -0.15,
        0.05
      ),
    
    betaFEED_a =
      rnorm(
        1,
        0,
        0.10
      ),
    
    betaFEED_b =
      rnorm(
        1,
        0,
        0.15
      ),
    
    sigma_island_a =
      runif(
        1,
        0.5,
        1.5
      ),
    
    island_a_raw =
      rnorm(
        n_islands,
        0,
        0.5
      ),
    
    sigma_time =
      runif(
        1,
        1,
        2.5
      ),
    
    u_time_raw =
      rnorm(
        n_time_blocks,
        0,
        0.5
      ),
    
    mean.p =
      runif(
        1,
        0.20,
        0.35
      )
  )
}


inits_REF10_FR <- lapply(
  1:3,
  make_inits_REF10_FR
)


## ============================================================
## 25. BUILD MODEL
## ============================================================

model_REF10_FR <- nimbleModel(
  code = code_REF10_FR,
  constants = consts_REF10_FR,
  data = data_REF10_FR,
  inits = inits_REF10_FR[[1]],
  calculate = FALSE
)


lp_REF10_FR <- model_REF10_FR$calculate()


cat(
  "\nInitial log probability:",
  lp_REF10_FR,
  "\n"
)


if(!is.finite(lp_REF10_FR)) {
  stop(
    "Initial model log probability is not finite."
  )
}


## ============================================================
## 26. CONFIGURE MCMC
## ============================================================

conf_REF10_FR <- configureMCMC(
  model_REF10_FR,
  
  monitors = c(
    "a10_fed",
    "a10_unfed",
    
    "b_fed",
    "b_unfed",
    
    "log_a10_fed",
    "log_a10_unfed",
    
    "log_b_fed",
    "log_b_unfed",
    
    "betaSEX",
    
    "betaDENS_a",
    "betaDENS_b",
    
    "betaFEED_a",
    "betaFEED_b",
    
    "sigma_island_a",
    "sigma_time",
    
    "island_a",
    "u_time",
    
    "mean.p"
  ),
  
  enableWAIC = FALSE
)


## ============================================================
## 27. JOINT a10-SIDE SAMPLER
## ============================================================

a_pars <- c(
  "log_a10_fed",
  "betaSEX[1]",
  "betaDENS_a",
  "betaFEED_a"
)


conf_REF10_FR$removeSamplers(
  a_pars
)


conf_REF10_FR$addSampler(
  target = a_pars,
  type = "AF_slice"
)


## ============================================================
## 28. JOINT b-SIDE SAMPLER
## ============================================================

b_pars <- c(
  "log_b_fed",
  "betaSEX[2]",
  "betaDENS_b",
  "betaFEED_b"
)


conf_REF10_FR$removeSamplers(
  b_pars
)


conf_REF10_FR$addSampler(
  target = b_pars,
  type = "AF_slice"
)


## ============================================================
## 29. SCALE PARAMETERS
## ============================================================

for(v in c(
  "sigma_island_a",
  "sigma_time"
)) {
  
  conf_REF10_FR$removeSamplers(
    v
  )
  
  conf_REF10_FR$addSampler(
    target = v,
    type = "slice"
  )
}


## ============================================================
## 30. JOINT LATENT BIRTH + LIFESPAN SAMPLERS
## ============================================================

for(i in seq_len(nind)) {
  
  life_nodes <- c(
    paste0(
      "tB[",
      i,
      "]"
    ),
    
    paste0(
      "tstar[",
      i,
      "]"
    )
  )
  
  
  conf_REF10_FR$removeSamplers(
    life_nodes
  )
  
  
  conf_REF10_FR$addSampler(
    target = life_nodes,
    type = "AF_slice"
  )
}


cat("\n")
cat("Custom samplers assigned:\n")
cat(
  "a10 block:",
  paste(a_pars, collapse = ", "),
  "\n"
)

cat(
  "b block:",
  paste(b_pars, collapse = ", "),
  "\n"
)

cat(
  "latent life-history blocks:",
  nind,
  "\n"
)


## ============================================================
## 31. BUILD + COMPILE
## ============================================================

mcmc_REF10_FR <- buildMCMC(
  conf_REF10_FR
)


cModel_REF10_FR <- compileNimble(
  model_REF10_FR,
  resetFunctions = TRUE
)


cMCMC_REF10_FR <- compileNimble(
  mcmc_REF10_FR,
  project = model_REF10_FR
)


## ============================================================
## 32. RUN 100k TEST
## ============================================================

cat("\n")
cat("============================================================\n")
cat("RUNNING REFAGE10 FED-REFERENCE 100k MODEL\n")
cat("============================================================\n")


run_REF10_FR <- runMCMC(
  cMCMC_REF10_FR,
  
  niter = 100000,
  nburnin = 30000,
  nchains = 3,
  thin = 10,
  
  inits = inits_REF10_FR,
  
  setSeed = TRUE,
  
  samplesAsCodaMCMC = TRUE,
  
  summary = TRUE,
  WAIC = FALSE,
  progressBar = TRUE
)


## ============================================================
## 33. SAVE RAW MODEL IMMEDIATELY
## ============================================================

saveRDS(
  run_REF10_FR,
  "Model_PaceOfLife_REFAGE10_FEEDREF_100k.rds"
)


cat(
  "\nSaved: Model_PaceOfLife_REFAGE10_FEEDREF_100k.rds\n"
)


## ============================================================
## 34. PRIMARY SUMMARY
## ============================================================

sm_REF10_FR <- MCMCsummary(
  run_REF10_FR$samples,
  
  params = c(
    "a10_fed",
    "a10_unfed",
    
    "b_fed",
    "b_unfed",
    
    "log_a10_fed",
    "log_a10_unfed",
    
    "log_b_fed",
    "log_b_unfed",
    
    "betaSEX",
    
    "betaDENS_a",
    "betaDENS_b",
    
    "betaFEED_a",
    "betaFEED_b",
    
    "sigma_island_a",
    "sigma_time",
    
    "mean.p"
  )
)


cat("\n")
cat("============================================================\n")
cat("PRIMARY POSTERIOR SUMMARY\n")
cat("============================================================\n")


print(
  sm_REF10_FR
)


## ============================================================
## 35. POSTERIOR MATRIX
## ============================================================

s_REF10_FR <- as.matrix(
  run_REF10_FR$samples
)


## ============================================================
## 36. BIOLOGICAL EFFECT SUMMARY
##
## exp(beta) interpretation:
##
## sex:
##   male / female
##
## density:
##   +1 SD density / reference density
##
## feeding:
##   fed / unfed
## ============================================================

effect_REF10_FR <- tibble(
  effect = c(
    "Sex -> mortality at age 10",
    "Sex -> senescence",
    
    "Density -> mortality at age 10",
    "Density -> senescence",
    
    "Feeding -> mortality at age 10",
    "Feeding -> senescence"
  ),
  
  parameter = c(
    "betaSEX[1]",
    "betaSEX[2]",
    
    "betaDENS_a",
    "betaDENS_b",
    
    "betaFEED_a",
    "betaFEED_b"
  )
) %>%
  rowwise() %>%
  mutate(
    median_beta =
      median(
        s_REF10_FR[, parameter]
      ),
    
    low_beta =
      unname(
        quantile(
          s_REF10_FR[, parameter],
          0.025
        )
      ),
    
    high_beta =
      unname(
        quantile(
          s_REF10_FR[, parameter],
          0.975
        )
      ),
    
    median_ratio =
      median(
        exp(
          s_REF10_FR[, parameter]
        )
      ),
    
    low_ratio =
      unname(
        quantile(
          exp(
            s_REF10_FR[, parameter]
          ),
          0.025
        )
      ),
    
    high_ratio =
      unname(
        quantile(
          exp(
            s_REF10_FR[, parameter]
          ),
          0.975
        )
      ),
    
    Pr_positive =
      mean(
        s_REF10_FR[, parameter] > 0
      ),
    
    Pr_negative =
      mean(
        s_REF10_FR[, parameter] < 0
      )
  ) %>%
  ungroup()


cat("\n")
cat("============================================================\n")
cat("BIOLOGICAL EFFECT SUMMARY\n")
cat("============================================================\n")


print(
  effect_REF10_FR,
  n = Inf
)


## ============================================================
## 37. KEY CORRELATIONS
##
## Main diagnostic:
##
## Has reference-group parameterisation removed the previous
## near-perfect log_b / betaFEED_b ridge?
## ============================================================

cor_a_b <- cor(
  s_REF10_FR[, "log_a10_fed"],
  s_REF10_FR[, "log_b_fed"]
)


cor_b_feed <- cor(
  s_REF10_FR[, "log_b_fed"],
  s_REF10_FR[, "betaFEED_b"]
)


cor_a_feed <- cor(
  s_REF10_FR[, "log_a10_fed"],
  s_REF10_FR[, "betaFEED_a"]
)


cat("\n")
cat("============================================================\n")
cat("KEY POSTERIOR CORRELATIONS\n")
cat("============================================================\n")

cat(
  "log_a10_fed vs log_b_fed:",
  round(cor_a_b, 3),
  "\n"
)

cat(
  "log_b_fed vs betaFEED_b:",
  round(cor_b_feed, 3),
  "\n"
)

cat(
  "log_a10_fed vs betaFEED_a:",
  round(cor_a_feed, 3),
  "\n"
)


## ============================================================
## 38. b-SIDE CORRELATION MATRIX
## ============================================================

cor_b_matrix <- cor(
  s_REF10_FR[
    ,
    c(
      "log_b_fed",
      "betaSEX[2]",
      "betaDENS_b",
      "betaFEED_b"
    )
  ]
)


cat("\n")
cat("============================================================\n")
cat("b-SIDE CORRELATION MATRIX\n")
cat("============================================================\n")


print(
  round(
    cor_b_matrix,
    3
  )
)


## ============================================================
## 39. CONVERGENCE SUMMARY
## ============================================================

conv_REF10_FR <- sm_REF10_FR %>%
  as.data.frame() %>%
  rownames_to_column(
    "parameter"
  ) %>%
  as_tibble() %>%
  arrange(
    desc(Rhat)
  )


cat("\n")
cat("============================================================\n")
cat("CONVERGENCE - WORST RHAT FIRST\n")
cat("============================================================\n")


print(
  conv_REF10_FR,
  n = Inf
)


cat(
  "\nMaximum Rhat:",
  max(
    conv_REF10_FR$Rhat,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Minimum ESS:",
  min(
    conv_REF10_FR$n.eff,
    na.rm = TRUE
  ),
  "\n"
)


## ============================================================
## 40. TEMPORAL EFFECTS
## ============================================================

sm_time_REF10_FR <- MCMCsummary(
  run_REF10_FR$samples,
  params = "u_time"
)


cat("\n")
cat("============================================================\n")
cat("TEMPORAL COHORT EFFECTS\n")
cat("============================================================\n")


print(
  sm_time_REF10_FR
)


## ============================================================
## 41. ISLAND EFFECTS
## ============================================================

sm_island_REF10_FR <- MCMCsummary(
  run_REF10_FR$samples,
  params = "island_a"
)


cat("\n")
cat("============================================================\n")
cat("RESIDUAL ISLAND EFFECTS ON a10\n")
cat("============================================================\n")


print(
  sm_island_REF10_FR
)


## ============================================================
## 42. SAVE COMPLETE RESULTS
## ============================================================

write_csv(
  effect_REF10_FR,
  "Model_PaceOfLife_REFAGE10_FEEDREF_100k_effects.csv"
)


saveRDS(
  list(
    fit =
      run_REF10_FR,
    
    summary =
      sm_REF10_FR,
    
    effects =
      effect_REF10_FR,
    
    convergence =
      conv_REF10_FR,
    
    temporal =
      sm_time_REF10_FR,
    
    island =
      sm_island_REF10_FR,
    
    correlation_a10_b =
      cor_a_b,
    
    correlation_b_feed =
      cor_b_feed,
    
    correlation_a_feed =
      cor_a_feed,
    
    b_correlation_matrix =
      cor_b_matrix,
    
    age_ref =
      age_ref,
    
    density_mean =
      dens_mu,
    
    density_sd =
      dens_sd,
    
    birth_change_audit =
      birth_change_audit,
    
    weak_birth_IDs =
      growth_info %>%
      filter(
        birth_weak_or_adjusted
      ) %>%
      select(
        ID,
        ID_raw,
        island,
        birth_source
      )
  ),
  
  "Model_PaceOfLife_REFAGE10_FEEDREF_100k_COMPLETE.rds"
)


cat("\n")
cat("============================================================\n")
cat("REFAGE10 FED-REFERENCE 100k TEST COMPLETE\n")
cat("============================================================\n")

print(sm_REF10_FR)

print(
  effect_REF10_FR,
  n = Inf
)

cat(
  "\nlog_a10_fed vs log_b_fed:",
  round(cor_a_b, 3),
  "\n"
)

cat(
  "log_b_fed vs betaFEED_b:",
  round(cor_b_feed, 3),
  "\n"
)

cat(
  "log_a10_fed vs betaFEED_a:",
  round(cor_a_feed, 3),
  "\n"
)

print(
  conv_REF10_FR,
  n = Inf
)