## ============================================================
## DataPrep_REFAGE10_Setup.R
##
## CANONICAL SETUP FOR REFERENCE-AGE PACE-OF-LIFE MODEL
##
## Starts from the existing canonical DataPrep.R, then:
##
## 1. Rebuilds birth intervals transparently from:
##      - known age, where available
##      - SVL-age 5th/95th quantile predictions otherwise
##      - broad 0-40 year range for individuals lacking SVL
##
## 2. Enforces the biologically certain rule:
##      birth <= first capture
##
## 3. DOES NOT apply the old blanket rule:
##      birth upper >= death/censor upper -> subtract 1
##
##    This was not appropriate for right-censored animals.
##
## 4. Explicitly flags:
##      - SVL/first-capture conflicts
##      - no-SVL birth histories
##
## 5. Builds:
##      - centred sex
##      - island-level centred feeding
##      - cohort density matrix
##      - temporal cohort blocks
##      - age_ref = 10
##
## ============================================================

source("DataPrep.R")

library(tidyverse)
library(nimble)

cat("\n")
cat("============================================================\n")
cat("REFERENCE-AGE MODEL SETUP\n")
cat("============================================================\n")

## ============================================================
## 1. SAVE ORIGINAL BIRTH WINDOWS FOR AUDIT
## ============================================================

cintB_adj_original <- cintB_adj

stopifnot(
  nrow(growth_info) == nind,
  nrow(cintB_adj_original) == nind,
  nrow(cintD_adj) == nind
)

## ============================================================
## 2. REBUILD CALENDAR-SCALE BIRTH WINDOWS
##
## This reproduces the intended biological construction without
## the later blanket birth/death correction.
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
  ## A. KNOWN AGE
  ## ----------------------------------------------------------
  
  if(!is.na(growth_info$known_birth_cal[i])) {
    
    cintB_cal_clean[i,] <-
      growth_info$known_birth_cal[i] +
      c(-1, 0)
    
    birth_source[i] <- "Known birth"
    
  } else {
    
    ## --------------------------------------------------------
    ## B. UNKNOWN AGE
    ##
    ## If SVL exists:
    ##   use SVL-derived age interval.
    ##
    ## If SVL does not exist:
    ##   pred_min_age = 0
    ##   pred_max_age = 40
    ##   and first capture becomes reference year.
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
    
    ## --------------------------------------------------------
    ## Does the SVL-derived birth window imply birth after an
    ## earlier recorded capture?
    ## --------------------------------------------------------
    
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
    ## First capture is direct evidence that the animal existed.
    ##
    ## Therefore birth cannot occur after first capture.
    ## --------------------------------------------------------
    
    final_high <- min(
      raw_high,
      growth_info$tF_cal[i]
    )
    
    final_low <- raw_low
    
    ## If constraining upper birth year invalidates the lower
    ## bound, use the one-year interval ending at first capture.
    if(final_low >= final_high) {
      final_low <- final_high - 1
    }
    
    cintB_cal_clean[i,] <- c(
      final_low,
      final_high
    )
  }
  
  if(any(!is.finite(cintB_cal_clean[i,]))) {
    stop(
      "Missing/invalid rebuilt birth interval for ID ",
      growth_info$ID_raw[i]
    )
  }
  
  if(
    cintB_cal_clean[i,1] >=
    cintB_cal_clean[i,2]
  ) {
    stop(
      "Non-positive rebuilt birth interval for ID ",
      growth_info$ID_raw[i]
    )
  }
}

## ============================================================
## 3. CHECK TIME ORIGIN
##
## Origin was originally created BEFORE the questionable
## birth/death adjustment. Rebuilding the intended intervals
## should therefore retain the same origin.
## ============================================================

Origin_rebuilt <- min(
  cintB_cal_clean[,1],
  na.rm = TRUE
)

cat(
  "\nExisting Origin:",
  Origin,
  "\n"
)

cat(
  "Rebuilt-birth Origin:",
  Origin_rebuilt,
  "\n"
)

if(Origin_rebuilt != Origin) {
  stop(
    "Rebuilt birth windows change Origin. ",
    "Do not proceed without reviewing this."
  )
}

## ============================================================
## 4. SHIFT CLEAN BIRTH WINDOWS TO MODEL TIME
## ============================================================

cintB_adj_clean <-
  cintB_cal_clean -
  Origin +
  1

## Enforce lower model boundary and first-capture constraint.

for(i in seq_len(nind)) {
  
  tF_adj_i <-
    growth_info$tF_cal[i] -
    Origin +
    1
  
  cintB_adj_clean[i,2] <- min(
    cintB_adj_clean[i,2],
    tF_adj_i
  )
  
  cintB_adj_clean[i,1] <- max(
    1,
    cintB_adj_clean[i,1]
  )
  
  if(
    cintB_adj_clean[i,1] >=
    cintB_adj_clean[i,2]
  ) {
    
    cintB_adj_clean[i,1] <- max(
      1,
      cintB_adj_clean[i,2] - 1
    )
  }
  
  if(
    cintB_adj_clean[i,1] >=
    cintB_adj_clean[i,2]
  ) {
    stop(
      "Invalid final clean birth interval for ID ",
      growth_info$ID_raw[i]
    )
  }
}

## ============================================================
## 5. COMPARE OLD AND CLEAN BIRTH WINDOWS
##
## This tells us exactly how many animals were affected by the
## old blanket birth/death correction.
## ============================================================

birth_changed <- which(
  rowSums(
    abs(
      cintB_adj_clean -
        cintB_adj_original
    ) > 1e-8
  ) > 0
)

cat(
  "\nIndividuals whose birth interval changes after removing",
  "the old birth/death correction:",
  length(birth_changed),
  "\n"
)

if(length(birth_changed) > 0) {
  
  birth_change_audit <- tibble(
    i = birth_changed,
    ID = growth_info$ID_raw[birth_changed],
    island = growth_info$island[birth_changed],
    censor = censoredD_vec[birth_changed],
    
    old_low =
      cintB_adj_original[birth_changed,1] +
      Origin - 1,
    
    old_high =
      cintB_adj_original[birth_changed,2] +
      Origin - 1,
    
    clean_low =
      cintB_adj_clean[birth_changed,1] +
      Origin - 1,
    
    clean_high =
      cintB_adj_clean[birth_changed,2] +
      Origin - 1
  )
  
  print(
    birth_change_audit,
    n = Inf
  )
  
  cat("\nChanges by censoring category:\n")
  
  print(
    birth_change_audit %>%
      count(censor)
  )
}

## ============================================================
## 6. REPLACE BIRTH WINDOWS WITH CLEAN VERSION
## ============================================================

cintB_adj <- cintB_adj_clean

## ============================================================
## 7. ADD EXPLICIT BIRTH-HISTORY FLAGS
## ============================================================

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
      cintB_adj[,1] +
      Origin - 1,
    
    birth_upper_cal =
      cintB_adj[,2] +
      Origin - 1,
    
    birth_window_width =
      cintB_adj[,2] -
      cintB_adj[,1]
  )

cat("\nBirth-source summary:\n")

print(
  growth_info %>%
    count(birth_source)
)

cat(
  "\nKnown birth:",
  sum(
    !is.na(
      growth_info$known_birth_cal
    )
  ),
  "\n"
)

cat(
  "SVL/first-capture conflicts:",
  sum(
    growth_info$birth_svl_capture_conflict
  ),
  "\n"
)

cat(
  "Unknown birth with no SVL:",
  sum(
    growth_info$birth_no_svl
  ),
  "\n"
)

cat(
  "Weak/adjusted birth histories for future sensitivity:",
  sum(
    growth_info$birth_weak_or_adjusted
  ),
  "\n"
)

## ============================================================
## 8. BIRTH-WINDOW VALIDITY
## ============================================================

stopifnot(
  all(
    cintB_adj[,1] <
      cintB_adj[,2]
  ),
  
  all(
    cintB_adj[,1] >= 1
  ),
  
  all(
    growth_info$birth_upper_cal <=
      growth_info$tF_cal
  )
)

## ============================================================
## 9. BIRTH VS DEATH/CENSORING AUDIT
##
## IMPORTANT:
##
## We DO NOT modify birth time merely because its upper bound
## overlaps a right-censoring point.
##
## For an interval-censored death, the truly impossible case is
## where the ENTIRE birth interval begins at/after the death
## interval's upper bound.
## ============================================================

impossible_interval_birth <- which(
  censoredD_vec == 1 &
    cintB_adj[,1] >=
    cintD_adj[,2]
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
    "At least one interval-censored death has an impossible ",
    "birth interval."
  )
}

birth_censor_overlap <- which(
  cintB_adj[,2] >=
    cintD_adj[,2]
)

cat(
  "\nBirth upper >= death/censor upper:",
  length(birth_censor_overlap),
  "\n"
)

if(length(birth_censor_overlap) > 0) {
  
  cat("By censoring category:\n")
  
  print(
    tibble(
      censor =
        censoredD_vec[
          birth_censor_overlap
        ]
    ) %>%
      count(censor)
  )
}

## ============================================================
## 10. REBUILD GENERIC INITIAL LIFE-HISTORY VALUES
##
## Model-specific chains will still generate their own random
## initial values. These are retained for canonical prep.
## ============================================================

tB_init <- rowMeans(
  cintB_adj
)

tD_init <- numeric(
  nind
)

for(i in seq_len(nind)) {
  
  if(censoredD_vec[i] == 1) {
    
    tD_init[i] <- mean(
      cintD_adj[i,]
    )
    
  } else if(censoredD_vec[i] == 2) {
    
    tD_init[i] <-
      cintD_adj[i,2] +
      1
  }
  
  if(tD_init[i] <= tB_init[i]) {
    
    if(censoredD_vec[i] == 1) {
      
      stop(
        "Initial birth >= interval death for ID ",
        growth_info$ID_raw[i]
      )
      
    } else {
      
      tD_init[i] <-
        tB_init[i] +
        1
    }
  }
}

tstar_init <- tD_init - tB_init

stopifnot(
  all(
    tB_init >=
      cintB_adj[,1]
  ),
  
  all(
    tB_init <=
      cintB_adj[,2]
  ),
  
  all(
    tstar_init > 0
  )
)

## ============================================================
## 11. REFERENCE AGE
##
## Chosen a priori from the observed age distribution:
##
## known-birth median observed midpoint ~9.5 years
## all-individual median first observed age ~10.5 years
##
## Therefore use a clear fixed reference age of 10 years.
## ============================================================

age_ref <- 10

cat(
  "\nReference age:",
  age_ref,
  "years\n"
)

## ============================================================
## 12. CENTRED SEX
## ============================================================

stopifnot(
  length(sex) == nind,
  !anyNA(sex)
)

sex_mean <- mean(
  sex
)

sex_c <- sex - sex_mean

cat(
  "Sex mean:",
  sex_mean,
  "\n"
)

## ============================================================
## 13. ISLAND-LEVEL FEEDING
##
## Treatment is an island/block-level variable.
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
          feeding != "none"
        )
      ),
    .groups = "drop"
  )

feed_island <- feeding_by_island$feed[
  match(
    ordered_islands,
    feeding_by_island$island
  )
]

stopifnot(
  !anyNA(
    feed_island
  )
)

feed_island_mean <- mean(
  feed_island
)

feed_c_island <-
  feed_island -
  feed_island_mean

feed_c <-
  feed_c_island[
    island_idx
  ]

cat("\nFeeding structure:\n")

print(
  tibble(
    island = ordered_islands,
    
    feeding = ifelse(
      feed_island == 1,
      "fed",
      "unfed"
    ),
    
    feed_c =
      feed_c_island
  )
)

## ============================================================
## 14. COHORT DENSITY
##
## MNA per usable hectare.
##
## Observed survey years use annual MNA.
## Missing historical cohort years use that island's mean
## observed density.
##
## log1p transform + global standardisation.
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
  unnest(
    year_idx
  ) %>%
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
    "Missing density remains after filling."
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
  
  density_matrix[jj,tt] <-
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
## 15. TEMPORAL COHORT STRUCTURE
##
## 5-year blocks.
## Final two raw blocks pooled as before.
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
  "modelled blocks\n"
)

## ============================================================
## 16. FINAL SETUP CHECKS
## ============================================================

stopifnot(
  nrow(growth_info) == nind,
  
  nrow(cintB_adj) == nind,
  ncol(cintB_adj) == 2,
  
  nrow(cintD_adj) == nind,
  ncol(cintD_adj) == 2,
  
  length(censoredD_vec) == nind,
  length(sex_c) == nind,
  length(feed_c) == nind,
  
  length(island_idx) == nind,
  length(island_start_time) == nind,
  length(island_start_idx) == nind,
  
  nrow(density_matrix) == n_islands,
  ncol(density_matrix) == tMax,
  
  nrow(cum_effort_for_model) == n_islands,
  ncol(cum_effort_for_model) >= tMax + 1,
  
  all(cintB_adj[,1] < cintB_adj[,2]),
  all(cintB_adj[,1] >= 1),
  
  all(island_idx >= 1),
  all(island_idx <= n_islands)
)

cat("\n")
cat("============================================================\n")
cat("REFERENCE-AGE SETUP COMPLETE\n")
cat("============================================================\n")
cat("Individuals:", nind, "\n")
cat("Islands:", n_islands, "\n")
cat("Reference age:", age_ref, "\n")
cat(
  "Weak/adjusted birth histories:",
  sum(
    growth_info$birth_weak_or_adjusted
  ),
  "\n"
)
cat("============================================================\n")



###############################################################################
###############################################################################
###############################################################################
###############################################################################
## ============================================================
## MODEL_PaceOfLife_REFAGE10_100k.R
##
## BIOLOGICALLY SPECIFIED PACE-OF-LIFE MODEL
##
## Gompertz parameterisation:
##
##   h(t) = a10 * exp[b * (t - 10)]
##
## a10:
##   mortality hazard at age 10
##
## b:
##   rate of actuarial senescence
##
## FOCAL BIOLOGICAL EFFECTS:
##
##   SEX
##   COHORT DENSITY
##   FEEDING
##
## Each focal biological predictor affects BOTH:
##
##   a10 = mortality level at age 10
##   b   = actuarial senescence
##
## HIERARCHICAL HETEROGENEITY:
##
##   residual island effect -> a10
##   temporal cohort effect -> a10
##
## TEST RUN:
##
##   100,000 iterations
##   30,000 burn-in
##   3 chains
##   thin = 10
##
## ============================================================

source("DataPrep_REFAGE10_Setup.R")

source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

library(nimble)
library(MCMCvis)
library(tidyverse)
library(coda)

q
