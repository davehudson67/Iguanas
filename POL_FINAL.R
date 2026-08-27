## ============================================================
## MODEL_PaceOfLife_REFAGE10_FEEDREF_TIMEPOOL_100k.R
##
## COMPLETE CORRECTED SCRIPT
##
## Gompertz parameterisation:
##
##   h(t) = a10 * exp[b * (t - 10)]
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
## Hierarchical heterogeneity on a10:
##   residual island effect
##   temporal cohort effect
##
## TEMPORAL STRUCTURE:
##
##   1946-1965 pooled
##   1966-1970
##   1971-1975
##   1976-1980
##   1981-1985
##   1986-1990
##   1991-1995
##   1996-2000
##   2001-2005
##   2006-2010
##   2011-2020 pooled
##
## TEST RUN:
##   niter  = 100000
##   burnin = 30000
##   chains = 3
##   thin   = 10
##
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
## DataPrep.R contains the older blanket birth/death correction.
## We rebuild the intended birth intervals below.
## ============================================================

cintB_adj_original <- cintB_adj

stopifnot(
  nrow(growth_info) == nind,
  nrow(cintB_adj_original) == nind,
  nrow(cintD_adj) == nind
)


## ============================================================
## 3. REBUILD CLEAN CALENDAR-SCALE BIRTH WINDOWS
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
    
    ref_yr <- if(has_svl) {
      growth_info$first_svl_year[i]
    } else {
      growth_info$tF_cal[i]
    }
    
    raw_low <- floor(
      ref_yr -
        growth_info$pred_max_age[i]
    )
    
    raw_high <- floor(
      ref_yr -
        growth_info$pred_min_age[i]
    )
    
    
    if(has_svl) {
      
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
      
      birth_source[i] <-
        "Unknown birth - no SVL"
    }
    
    
    ## --------------------------------------------------------
    ## FIRST CAPTURE IS HARD EVIDENCE THE ANIMAL WAS ALIVE
    ## --------------------------------------------------------
    
    final_high <- min(
      raw_high,
      growth_info$tF_cal[i]
    )
    
    final_low <- raw_low
    
    
    ## If SVL-derived interval lies entirely after first capture,
    ## use the one-year interval ending at first capture.
    
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
    "Rebuilt birth windows alter Origin."
  )
}


## ============================================================
## 5. SHIFT CLEAN BIRTH WINDOWS TO MODEL TIME
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
      "Invalid final clean birth interval for ",
      growth_info$ID_raw[i]
    )
  }
}


## ============================================================
## 6. AUDIT CHANGE FROM ORIGINAL DataPrep.R WINDOWS
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
  sum(
    growth_info$birth_weak_or_adjusted
  ),
  "\n"
)


## ============================================================
## 8. BIRTH VALIDITY
## ============================================================

stopifnot(
  all(
    cintB_adj[, 1] <
      cintB_adj[, 2]
  ),
  
  all(
    cintB_adj[, 1] >= 1
  ),
  
  all(
    growth_info$birth_upper_cal <=
      growth_info$tF_cal
  )
)


## ============================================================
## 9. TRUE INTERVAL-DEATH IMPOSSIBILITY CHECK
##
## Do NOT modify right-censored birth intervals.
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
## 11. SEX
##
## female = 0
## male   = 1
## ============================================================

stopifnot(
  length(sex) == nind,
  !anyNA(sex),
  all(sex %in% c(0, 1))
)


sex_mean <- mean(
  sex
)


sex_c <-
  sex -
  sex_mean


cat(
  "Sex mean:",
  sex_mean,
  "\n"
)


## ============================================================
## 12. FEEDING REFERENCE CODING
##
## fed   =  0
## unfed = -1
##
## betaFEED = fed - unfed on log scale.
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


cat("\nFeeding coding:\n")


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
      (
        usable_m2 /
          10000
      )
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


if(
  anyNA(
    density_complete$density_filled_raw
  )
) {
  
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
    "Invalid density SD."
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


for(r in seq_len(
  nrow(
    density_complete
  )
)) {
  
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
  !anyNA(
    density_matrix
  )
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
## 14. REVISED TEMPORAL COHORT STRUCTURE
##
## RAW 5-YEAR BLOCKS:
##
##  1 = 1946-1950
##  2 = 1951-1955
##  3 = 1956-1960
##  4 = 1961-1965
##  5 = 1966-1970
##  6 = 1971-1975
##  7 = 1976-1980
##  8 = 1981-1985
##  9 = 1986-1990
## 10 = 1991-1995
## 11 = 1996-2000
## 12 = 2001-2005
## 13 = 2006-2010
## 14 = 2011-2015
## 15 = 2016-2020
##
## MODEL BLOCKS:
##
## raw 1-4   -> model 1
## raw 5     -> model 2
## raw 6     -> model 3
## ...
## raw 13    -> model 10
## raw 14-15 -> model 11
##
## Algebraic mapping:
##
##   model_block =
##       min(11, max(1, raw_block - 3))
##
## This avoids dynamic indexing of a constant vector in NIMBLE.
## ============================================================

block_width <- 5


n_time_blocks_raw <- ceiling(
  tMax /
    block_width
)


if(n_time_blocks_raw != 15) {
  
  stop(
    "Temporal structure assumes 15 raw blocks; found ",
    n_time_blocks_raw
  )
}


n_time_blocks <- 11


## ============================================================
## 15. DESCRIPTIVE TEMPORAL TABLE
## ============================================================

raw_temporal_table <- tibble(
  raw_block =
    seq_len(
      n_time_blocks_raw
    ),
  
  raw_model_start =
    (
      seq_len(
        n_time_blocks_raw
      ) -
        1
    ) *
    block_width +
    1,
  
  raw_model_end =
    pmin(
      seq_len(
        n_time_blocks_raw
      ) *
        block_width,
      tMax
    )
) %>%
  mutate(
    model_block =
      pmin(
        n_time_blocks,
        pmax(
          1,
          raw_block - 3
        )
      ),
    
    calendar_start =
      raw_model_start +
      Origin -
      1,
    
    calendar_end =
      raw_model_end +
      Origin -
      1
  )


cat("\n")
cat("============================================================\n")
cat("TEMPORAL RAW -> MODEL BLOCK MAPPING\n")
cat("============================================================\n")


print(
  raw_temporal_table,
  n = Inf
)


model_temporal_table <- raw_temporal_table %>%
  group_by(model_block) %>%
  summarise(
    calendar_start =
      min(
        calendar_start
      ),
    
    calendar_end =
      max(
        calendar_end
      ),
    
    raw_blocks =
      paste(
        raw_block,
        collapse = ","
      ),
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("FINAL MODEL TEMPORAL PERIODS\n")
cat("============================================================\n")


print(
  model_temporal_table,
  n = Inf
)


cat(
  "\nTemporal structure:",
  n_time_blocks_raw,
  "raw blocks ->",
  n_time_blocks,
  "model blocks\n"
)


## ============================================================
## 16. DESCRIPTIVE SUPPORT FOR NEW TEMPORAL BLOCKS
##
## Midpoint used ONLY for this descriptive check.
## Birth remains latent in model.
## ============================================================

birth_mid_idx <- rowMeans(
  cintB_adj
)


birth_mid_raw_block <- pmax(
  1,
  pmin(
    n_time_blocks_raw,
    floor(
      (
        birth_mid_idx -
          1
      ) /
        block_width
    ) +
      1
  )
)


birth_mid_model_block <- pmin(
  n_time_blocks,
  pmax(
    1,
    birth_mid_raw_block - 3
  )
)


block_support <- tibble(
  model_block =
    birth_mid_model_block
) %>%
  count(
    model_block,
    name = "n_midpoint"
  ) %>%
  right_join(
    model_temporal_table,
    by = "model_block"
  ) %>%
  mutate(
    n_midpoint =
      replace_na(
        n_midpoint,
        0L
      )
  ) %>%
  arrange(
    model_block
  )


cat("\n")
cat("============================================================\n")
cat("MIDPOINT SUPPORT FOR NEW TEMPORAL BLOCKS\n")
cat("============================================================\n")


print(
  block_support,
  n = Inf
)


## ============================================================
## 17. FINAL SETUP CHECKS
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
  
  n_time_blocks == 11,
  
  nrow(density_matrix) == n_islands,
  ncol(density_matrix) == tMax,
  
  nrow(cum_effort_for_model) == n_islands,
  
  ncol(cum_effort_for_model) >=
    tMax + 1,
  
  all(
    cintB_adj[, 1] <
      cintB_adj[, 2]
  ),
  
  all(
    cintB_adj[, 1] >= 1
  ),
  
  all(
    island_idx >= 1
  ),
  
  all(
    island_idx <= n_islands
  )
)


cat("\n")
cat("============================================================\n")
cat("SETUP COMPLETE\n")
cat("============================================================\n")
cat("Individuals:", nind, "\n")
cat("Islands:", n_islands, "\n")
cat("Reference age:", age_ref, "\n")
cat("Temporal model blocks:", n_time_blocks, "\n")
cat("============================================================\n")


## ============================================================
## 18. NIMBLE MODEL
## ============================================================

code_REF10_FR_TP <- nimbleCode({
  
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
    
    
    ## --------------------------------------------------------
    ## RAW 5-YEAR COHORT BLOCK
    ## --------------------------------------------------------
    
    time_block_raw[i] <- max(
      1,
      min(
        n_time_blocks_raw,
        floor(
          (
            cohort_idx[i] -
              1
          ) /
            block_width
        ) +
          1
      )
    )
    
    
    ## --------------------------------------------------------
    ## REVISED TEMPORAL BLOCK
    ##
    ## raw 1-4   -> 1
    ## raw 5     -> 2
    ## ...
    ## raw 13    -> 10
    ## raw 14-15 -> 11
    ##
    ## No dynamic indexing of a constant is required.
    ## --------------------------------------------------------
    
    time_block[i] <- min(
      n_time_blocks,
      max(
        1,
        time_block_raw[i] - 3
      )
    )
    
    
    ## --------------------------------------------------------
    ## LATENT COHORT DENSITY
    ## --------------------------------------------------------
    
    dens_cohort[i] <-
      density_matrix[
        island_idx[i],
        cohort_idx[i]
      ]
    
    
    ## ========================================================
    ## MORTALITY LEVEL AT AGE 10
    ## ========================================================
    
    log_a10_ind[i] <-
      log_a10_fed +
      betaSEX[1] *
      sex_c[i] +
      betaDENS_a *
      dens_cohort[i] +
      betaFEED_a *
      feed_ref[i] +
      u_time[
        time_block[i]
      ] +
      island_a[
        island_idx[i]
      ]
    
    
    a10_ind[i] <- exp(
      log_a10_ind[i]
    )
    
    
    ## ========================================================
    ## ACTUARIAL SENESCENCE
    ## ========================================================
    
    log_b_ind[i] <-
      log_b_fed +
      betaSEX[2] *
      sex_c[i] +
      betaDENS_b *
      dens_cohort[i] +
      betaFEED_b *
      feed_ref[i]
    
    
    bmult[i] <- exp(
      log_b_ind[i]
    )
    
    
    ## ========================================================
    ## CONVERT a10 TO CONVENTIONAL AGE-ZERO GOMPERTZ a
    ##
    ## h(t) = a10 * exp[b(t-age_ref)]
    ##
    ## a0 = a10 * exp(-b*age_ref)
    ## ========================================================
    
    log_amult[i] <-
      log_a10_ind[i] -
      bmult[i] *
      age_ref
    
    
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
    ## DETECTION COMPONENT
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
  ## 19. FOCAL BIOLOGICAL EFFECTS
  ## ==========================================================
  
  for(k in 1:2) {
    
    betaSEX[k] ~ dnorm(
      0,
      sd = 1
    )
  }
  
  
  betaDENS_a ~ dnorm(
    0,
    sd = 1
  )
  
  
  betaDENS_b ~ dnorm(
    0,
    sd = 1
  )
  
  
  betaFEED_a ~ dnorm(
    0,
    sd = 1
  )
  
  
  betaFEED_b ~ dnorm(
    0,
    sd = 1
  )
  
  
  ## ==========================================================
  ## 20. RESIDUAL ISLAND HETEROGENEITY ON a10
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
      island_a_raw[
        1:n_islands
      ]
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
  ## 21. TEMPORAL COHORT HETEROGENEITY ON a10
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
      u_time_raw[
        1:n_time_blocks
      ]
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
  ## 22. FED REFERENCE PARAMETERS
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
  ## 23. DETECTION
  ## ==========================================================
  
  mean.p ~ dunif(
    0,
    1
  )
})


## ============================================================
## 24. CONSTANTS
##
## NOTE:
## There is deliberately NO time_block_map constant.
## ============================================================

consts_REF10_FR_TP <- list(
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
## 25. DATA
## ============================================================

data_REF10_FR_TP <- list(
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
## 26. INITIAL VALUES
## ============================================================

make_inits_REF10_FR_TP <- function(
    chain_id = 1
) {
  
  set.seed(
    10100 +
      chain_id
  )
  
  
  ## ----------------------------------------------------------
  ## LATENT BIRTH
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
  ## LATENT DEATH
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
          "No valid interval-death initial value for ",
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
    ## RIGHT-CENSORED
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


inits_REF10_FR_TP <- lapply(
  1:3,
  make_inits_REF10_FR_TP
)


## ============================================================
## 27. BUILD MODEL
## ============================================================

model_REF10_FR_TP <- nimbleModel(
  code = code_REF10_FR_TP,
  constants = consts_REF10_FR_TP,
  data = data_REF10_FR_TP,
  inits = inits_REF10_FR_TP[[1]],
  calculate = FALSE
)


lp_REF10_FR_TP <-
  model_REF10_FR_TP$calculate()


cat(
  "\nInitial log probability:",
  lp_REF10_FR_TP,
  "\n"
)


if(
  !is.finite(
    lp_REF10_FR_TP
  )
) {
  
  stop(
    "Initial model log probability is not finite."
  )
}


## ============================================================
## 28. CONFIGURE MCMC
## ============================================================

conf_REF10_FR_TP <- configureMCMC(
  model_REF10_FR_TP,
  
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
## 29. JOINT a10-SIDE SAMPLER
## ============================================================

a_pars <- c(
  "log_a10_fed",
  "betaSEX[1]",
  "betaDENS_a",
  "betaFEED_a"
)


conf_REF10_FR_TP$removeSamplers(
  a_pars
)


conf_REF10_FR_TP$addSampler(
  target = a_pars,
  type = "AF_slice"
)


## ============================================================
## 30. JOINT b-SIDE SAMPLER
## ============================================================

b_pars <- c(
  "log_b_fed",
  "betaSEX[2]",
  "betaDENS_b",
  "betaFEED_b"
)


conf_REF10_FR_TP$removeSamplers(
  b_pars
)


conf_REF10_FR_TP$addSampler(
  target = b_pars,
  type = "AF_slice"
)


## ============================================================
## 31. SCALE PARAMETERS
## ============================================================

for(v in c(
  "sigma_island_a",
  "sigma_time"
)) {
  
  conf_REF10_FR_TP$removeSamplers(
    v
  )
  
  
  conf_REF10_FR_TP$addSampler(
    target = v,
    type = "slice"
  )
}


## ============================================================
## 32. LATENT BIRTH + LIFESPAN JOINT SAMPLERS
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
  
  
  conf_REF10_FR_TP$removeSamplers(
    life_nodes
  )
  
  
  conf_REF10_FR_TP$addSampler(
    target = life_nodes,
    type = "AF_slice"
  )
}


cat("\nCustom samplers assigned:\n")


cat(
  "a10 block:",
  paste(
    a_pars,
    collapse = ", "
  ),
  "\n"
)


cat(
  "b block:",
  paste(
    b_pars,
    collapse = ", "
  ),
  "\n"
)


cat(
  "latent life-history blocks:",
  nind,
  "\n"
)


## ============================================================
## 33. BUILD + COMPILE
## ============================================================

mcmc_REF10_FR_TP <- buildMCMC(
  conf_REF10_FR_TP
)


cModel_REF10_FR_TP <- compileNimble(
  model_REF10_FR_TP,
  resetFunctions = TRUE
)


cMCMC_REF10_FR_TP <- compileNimble(
  mcmc_REF10_FR_TP,
  project = model_REF10_FR_TP
)


## ============================================================
## 34. RUN 100k TEST
## ============================================================

cat("\n")
cat("============================================================\n")
cat("RUNNING TIME-POOLED REFAGE10 MODEL\n")
cat("============================================================\n")
cat("Iterations: 100000\n")
cat("Burn-in:    30000\n")
cat("Chains:     3\n")
cat("Thin:       10\n")
cat("============================================================\n")


run_REF10_FR_TP <- runMCMC(
  cMCMC_REF10_FR_TP,
  
  niter = 100000,
  nburnin = 30000,
  nchains = 3,
  thin = 10,
  
  inits =
    inits_REF10_FR_TP,
  
  setSeed = TRUE,
  
  samplesAsCodaMCMC = TRUE,
  
  summary = TRUE,
  
  WAIC = FALSE,
  
  progressBar = TRUE
)


## ============================================================
## 35. SAVE RAW FIT IMMEDIATELY
## ============================================================

saveRDS(
  run_REF10_FR_TP,
  "Model_PaceOfLife_REFAGE10_FEEDREF_TIMEPOOL_100k.rds"
)


cat(
  "\nSaved: Model_PaceOfLife_REFAGE10_FEEDREF_TIMEPOOL_100k.rds\n"
)


## ============================================================
## 36. PRIMARY SUMMARY
## ============================================================

sm_REF10_FR_TP <- MCMCsummary(
  run_REF10_FR_TP$samples,
  
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
  sm_REF10_FR_TP
)


## ============================================================
## 37. POSTERIOR MATRIX
## ============================================================

s_REF10_FR_TP <- as.matrix(
  run_REF10_FR_TP$samples
)


## ============================================================
## 38. BIOLOGICAL EFFECT SUMMARY
##
## exp(beta):
##
## sex:
##   male / female
##
## density:
##   +1 SD cohort density
##
## feeding:
##   fed / unfed
## ============================================================

effect_REF10_FR_TP <- tibble(
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
        s_REF10_FR_TP[
          ,
          parameter
        ]
      ),
    
    low_beta =
      unname(
        quantile(
          s_REF10_FR_TP[
            ,
            parameter
          ],
          0.025
        )
      ),
    
    high_beta =
      unname(
        quantile(
          s_REF10_FR_TP[
            ,
            parameter
          ],
          0.975
        )
      ),
    
    median_ratio =
      median(
        exp(
          s_REF10_FR_TP[
            ,
            parameter
          ]
        )
      ),
    
    low_ratio =
      unname(
        quantile(
          exp(
            s_REF10_FR_TP[
              ,
              parameter
            ]
          ),
          0.025
        )
      ),
    
    high_ratio =
      unname(
        quantile(
          exp(
            s_REF10_FR_TP[
              ,
              parameter
            ]
          ),
          0.975
        )
      ),
    
    Pr_positive =
      mean(
        s_REF10_FR_TP[
          ,
          parameter
        ] > 0
      ),
    
    Pr_negative =
      mean(
        s_REF10_FR_TP[
          ,
          parameter
        ] < 0
      )
  ) %>%
  ungroup()


cat("\n")
cat("============================================================\n")
cat("BIOLOGICAL EFFECT SUMMARY\n")
cat("============================================================\n")


print(
  effect_REF10_FR_TP,
  n = Inf
)


## ============================================================
## 39. KEY POSTERIOR CORRELATIONS
## ============================================================

cor_a_b_TP <- cor(
  s_REF10_FR_TP[
    ,
    "log_a10_fed"
  ],
  
  s_REF10_FR_TP[
    ,
    "log_b_fed"
  ]
)


cor_b_feed_TP <- cor(
  s_REF10_FR_TP[
    ,
    "log_b_fed"
  ],
  
  s_REF10_FR_TP[
    ,
    "betaFEED_b"
  ]
)


cor_a_feed_TP <- cor(
  s_REF10_FR_TP[
    ,
    "log_a10_fed"
  ],
  
  s_REF10_FR_TP[
    ,
    "betaFEED_a"
  ]
)


cat("\n")
cat("============================================================\n")
cat("KEY POSTERIOR CORRELATIONS\n")
cat("============================================================\n")


cat(
  "log_a10_fed vs log_b_fed:",
  round(
    cor_a_b_TP,
    3
  ),
  "\n"
)


cat(
  "log_b_fed vs betaFEED_b:",
  round(
    cor_b_feed_TP,
    3
  ),
  "\n"
)


cat(
  "log_a10_fed vs betaFEED_a:",
  round(
    cor_a_feed_TP,
    3
  ),
  "\n"
)


## ============================================================
## 40. b-SIDE CORRELATION MATRIX
## ============================================================

cor_b_matrix_TP <- cor(
  s_REF10_FR_TP[
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
    cor_b_matrix_TP,
    3
  )
)


## ============================================================
## 41. CONVERGENCE SUMMARY
## ============================================================

conv_REF10_FR_TP <- sm_REF10_FR_TP %>%
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
  conv_REF10_FR_TP,
  n = Inf
)


cat(
  "\nMaximum Rhat:",
  max(
    conv_REF10_FR_TP$Rhat,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Minimum ESS:",
  min(
    conv_REF10_FR_TP$n.eff,
    na.rm = TRUE
  ),
  "\n"
)


## ============================================================
## 42. TEMPORAL EFFECTS
## ============================================================

sm_time_REF10_FR_TP <- MCMCsummary(
  run_REF10_FR_TP$samples,
  params = "u_time"
)


cat("\n")
cat("============================================================\n")
cat("TEMPORAL COHORT EFFECTS\n")
cat("============================================================\n")


print(
  sm_time_REF10_FR_TP
)


## ============================================================
## 43. ISLAND EFFECTS
## ============================================================

sm_island_REF10_FR_TP <- MCMCsummary(
  run_REF10_FR_TP$samples,
  params = "island_a"
)


cat("\n")
cat("============================================================\n")
cat("RESIDUAL ISLAND EFFECTS ON MORTALITY AT AGE 10\n")
cat("============================================================\n")


print(
  sm_island_REF10_FR_TP
)


## ============================================================
## 44. SAVE COMPLETE RESULTS
## ============================================================

write_csv(
  effect_REF10_FR_TP,
  "Model_PaceOfLife_REFAGE10_FEEDREF_TIMEPOOL_100k_effects.csv"
)


saveRDS(
  list(
    fit =
      run_REF10_FR_TP,
    
    summary =
      sm_REF10_FR_TP,
    
    effects =
      effect_REF10_FR_TP,
    
    convergence =
      conv_REF10_FR_TP,
    
    temporal =
      sm_time_REF10_FR_TP,
    
    temporal_mapping =
      model_temporal_table,
    
    temporal_support =
      block_support,
    
    island =
      sm_island_REF10_FR_TP,
    
    correlation_a10_b =
      cor_a_b_TP,
    
    correlation_b_feed =
      cor_b_feed_TP,
    
    correlation_a_feed =
      cor_a_feed_TP,
    
    b_correlation_matrix =
      cor_b_matrix_TP,
    
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
  
  "Model_PaceOfLife_REFAGE10_FEEDREF_TIMEPOOL_100k_COMPLETE.rds"
)


cat("\n")
cat("============================================================\n")
cat("TIME-POOLED REFAGE10 100k TEST COMPLETE\n")
cat("============================================================\n")

## ============================================================
## TIME-POOLED MODEL: FINAL 100k DIAGNOSTICS
## ============================================================

cat("\n")
cat("============================================================\n")
cat("PRIMARY POSTERIOR SUMMARY\n")
cat("============================================================\n")

print(
  sm_REF10_FR_TP
)


cat("\n")
cat("============================================================\n")
cat("BIOLOGICAL EFFECT SUMMARY\n")
cat("============================================================\n")

print(
  effect_REF10_FR_TP,
  n = Inf
)


cat("\n")
cat("============================================================\n")
cat("KEY POSTERIOR CORRELATIONS\n")
cat("============================================================\n")

cat(
  "log_a10_fed vs log_b_fed:",
  round(cor_a_b_TP, 3),
  "\n"
)

cat(
  "log_b_fed vs betaFEED_b:",
  round(cor_b_feed_TP, 3),
  "\n"
)

cat(
  "log_a10_fed vs betaFEED_a:",
  round(cor_a_feed_TP, 3),
  "\n"
)


cat("\n")
cat("============================================================\n")
cat("b-SIDE CORRELATION MATRIX\n")
cat("============================================================\n")

print(
  round(
    cor_b_matrix_TP,
    3
  )
)


cat("\n")
cat("============================================================\n")
cat("CONVERGENCE - WORST RHAT FIRST\n")
cat("============================================================\n")

print(
  conv_REF10_FR_TP,
  n = Inf
)


cat("\n")
cat("============================================================\n")
cat("TEMPORAL COHORT EFFECTS\n")
cat("============================================================\n")

print(
  sm_time_REF10_FR_TP
)


cat("\n")
cat("============================================================\n")
cat("TEMPORAL BLOCK SUPPORT\n")
cat("============================================================\n")

print(
  block_support,
  n = Inf
)


## ============================================================
## FINAL 500k RUN
## REFAGE10 + FED REFERENCE + POOLED TEMPORAL BLOCKS
## ============================================================

inits_REF10_FR_TP_FINAL <- lapply(
  1:3,
  function(chain_id) {
    make_inits_REF10_FR_TP(
      chain_id = chain_id + 100
    )
  }
)

cat("\n")
cat("============================================================\n")
cat("RUNNING FINAL TIME-POOLED REFAGE10 MODEL\n")
cat("============================================================\n")
cat("Iterations: 500000\n")
cat("Burn-in:    100000\n")
cat("Chains:     3\n")
cat("Thin:       10\n")
cat("============================================================\n")


run_REF10_FR_TP_FINAL <- runMCMC(
  cMCMC_REF10_FR_TP,
  
  niter = 500000,
  nburnin = 100000,
  nchains = 3,
  thin = 10,
  
  inits =
    inits_REF10_FR_TP_FINAL,
  
  setSeed = TRUE,
  
  samplesAsCodaMCMC = TRUE,
  
  summary = TRUE,
  
  WAIC = FALSE,
  
  progressBar = TRUE
)


## ============================================================
## SAVE RAW FIT IMMEDIATELY
## ============================================================

saveRDS(
  run_REF10_FR_TP_FINAL,
  "Model_PaceOfLife_REFAGE10_FEEDREF_TIMEPOOL_FINAL_500k.rds"
)


cat(
  "\nSaved: Model_PaceOfLife_REFAGE10_FEEDREF_TIMEPOOL_FINAL_500k.rds\n"
)


## ============================================================
## PRIMARY SUMMARY
## ============================================================

sm_REF10_FR_TP_FINAL <- MCMCsummary(
  run_REF10_FR_TP_FINAL$samples,
  
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
cat("FINAL PRIMARY POSTERIOR SUMMARY\n")
cat("============================================================\n")

print(
  sm_REF10_FR_TP_FINAL
)


## ============================================================
## POSTERIOR MATRIX
## ============================================================

s_REF10_FR_TP_FINAL <- as.matrix(
  run_REF10_FR_TP_FINAL$samples
)


## ============================================================
## BIOLOGICAL EFFECT SUMMARY
## ============================================================

effect_REF10_FR_TP_FINAL <- tibble(
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
        s_REF10_FR_TP_FINAL[, parameter]
      ),
    
    low_beta =
      unname(
        quantile(
          s_REF10_FR_TP_FINAL[, parameter],
          0.025
        )
      ),
    
    high_beta =
      unname(
        quantile(
          s_REF10_FR_TP_FINAL[, parameter],
          0.975
        )
      ),
    
    median_ratio =
      median(
        exp(
          s_REF10_FR_TP_FINAL[, parameter]
        )
      ),
    
    low_ratio =
      unname(
        quantile(
          exp(
            s_REF10_FR_TP_FINAL[, parameter]
          ),
          0.025
        )
      ),
    
    high_ratio =
      unname(
        quantile(
          exp(
            s_REF10_FR_TP_FINAL[, parameter]
          ),
          0.975
        )
      ),
    
    Pr_positive =
      mean(
        s_REF10_FR_TP_FINAL[, parameter] > 0
      ),
    
    Pr_negative =
      mean(
        s_REF10_FR_TP_FINAL[, parameter] < 0
      )
  ) %>%
  ungroup()


cat("\n")
cat("============================================================\n")
cat("FINAL BIOLOGICAL EFFECT SUMMARY\n")
cat("============================================================\n")

print(
  effect_REF10_FR_TP_FINAL,
  n = Inf
)


## ============================================================
## TEMPORAL EFFECTS
## ============================================================

sm_time_REF10_FR_TP_FINAL <- MCMCsummary(
  run_REF10_FR_TP_FINAL$samples,
  params = "u_time"
)


cat("\n")
cat("============================================================\n")
cat("FINAL TEMPORAL EFFECTS\n")
cat("============================================================\n")

print(
  sm_time_REF10_FR_TP_FINAL
)


## ============================================================
## ISLAND EFFECTS
## ============================================================

sm_island_REF10_FR_TP_FINAL <- MCMCsummary(
  run_REF10_FR_TP_FINAL$samples,
  params = "island_a"
)


cat("\n")
cat("============================================================\n")
cat("FINAL ISLAND EFFECTS\n")
cat("============================================================\n")

print(
  sm_island_REF10_FR_TP_FINAL
)


## ============================================================
## KEY CORRELATIONS
## ============================================================

cor_a_b_FINAL <- cor(
  s_REF10_FR_TP_FINAL[, "log_a10_fed"],
  s_REF10_FR_TP_FINAL[, "log_b_fed"]
)

cor_b_feed_FINAL <- cor(
  s_REF10_FR_TP_FINAL[, "log_b_fed"],
  s_REF10_FR_TP_FINAL[, "betaFEED_b"]
)

cor_a_feed_FINAL <- cor(
  s_REF10_FR_TP_FINAL[, "log_a10_fed"],
  s_REF10_FR_TP_FINAL[, "betaFEED_a"]
)


cat("\n")
cat("============================================================\n")
cat("FINAL KEY CORRELATIONS\n")
cat("============================================================\n")

cat(
  "log_a10_fed vs log_b_fed:",
  round(cor_a_b_FINAL, 3),
  "\n"
)

cat(
  "log_b_fed vs betaFEED_b:",
  round(cor_b_feed_FINAL, 3),
  "\n"
)

cat(
  "log_a10_fed vs betaFEED_a:",
  round(cor_a_feed_FINAL, 3),
  "\n"
)


## ============================================================
## b-SIDE CORRELATION MATRIX
## ============================================================

cor_b_matrix_FINAL <- cor(
  s_REF10_FR_TP_FINAL[
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
cat("FINAL b-SIDE CORRELATION MATRIX\n")
cat("============================================================\n")

print(
  round(
    cor_b_matrix_FINAL,
    3
  )
)


## ============================================================
## CONVERGENCE
## ============================================================

conv_REF10_FR_TP_FINAL <- sm_REF10_FR_TP_FINAL %>%
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
cat("FINAL CONVERGENCE - WORST RHAT FIRST\n")
cat("============================================================\n")

print(
  conv_REF10_FR_TP_FINAL,
  n = Inf
)


cat(
  "\nMaximum Rhat:",
  max(
    conv_REF10_FR_TP_FINAL$Rhat,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Minimum ESS:",
  min(
    conv_REF10_FR_TP_FINAL$n.eff,
    na.rm = TRUE
  ),
  "\n"
)


## ============================================================
## SAVE COMPLETE FINAL ANALYSIS
## ============================================================

write_csv(
  effect_REF10_FR_TP_FINAL,
  "Model_PaceOfLife_REFAGE10_FEEDREF_TIMEPOOL_FINAL_500k_effects.csv"
)


saveRDS(
  list(
    fit =
      run_REF10_FR_TP_FINAL,
    
    summary =
      sm_REF10_FR_TP_FINAL,
    
    effects =
      effect_REF10_FR_TP_FINAL,
    
    convergence =
      conv_REF10_FR_TP_FINAL,
    
    temporal =
      sm_time_REF10_FR_TP_FINAL,
    
    temporal_mapping =
      model_temporal_table,
    
    temporal_support =
      block_support,
    
    island =
      sm_island_REF10_FR_TP_FINAL,
    
    correlation_a10_b =
      cor_a_b_FINAL,
    
    correlation_b_feed =
      cor_b_feed_FINAL,
    
    correlation_a_feed =
      cor_a_feed_FINAL,
    
    b_correlation_matrix =
      cor_b_matrix_FINAL,
    
    age_ref =
      age_ref,
    
    density_mean =
      dens_mu,
    
    density_sd =
      dens_sd,
    
    birth_change_audit =
      birth_change_audit
  ),
  
  "Model_PaceOfLife_REFAGE10_FEEDREF_TIMEPOOL_FINAL_500k_COMPLETE.rds"
)


cat("\n")
cat("============================================================\n")
cat("FINAL 500k MODEL COMPLETE\n")
cat("============================================================\n")
