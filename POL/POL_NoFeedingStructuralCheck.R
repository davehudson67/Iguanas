## ============================================================
## IGUANA PACE-OF-LIFE MODEL
##
## NO-FEEDING STRUCTURAL CHECK
##
## - Fresh-session setup
## - Corrected birth intervals
## - Reference-age Gompertz: age 10
## - Sex: female = 0, male = 1
## - Cohort density:
##      island x cohort-year specific
##      first 5 survey years treated as start-up
##      established annual MNA used from survey year 6 onward
##      established island mean used otherwise
## - Shared pooled temporal cohort effects
## - Residual island effect on mortality at age 10
## - NO feeding term
##
## DEVELOPMENT RUN:
##   100000 iterations
##   30000 burn-in
##   3 chains
##   thin = 10
##
## ============================================================


## ============================================================
## 1. CANONICAL DATA PREPARATION
## ============================================================

source("DataPrep.R")


## ============================================================
## 2. LOAD CUSTOM DISTRIBUTIONS
##
## Handles either:
##   Dist_GompertzLB.R
##
## or:
##   Distributions/Dist_GompertzLB.R
## ============================================================

source_first_existing <- function(paths) {
  
  hit <- paths[
    file.exists(
      paths
    )
  ]
  
  if(
    length(
      hit
    ) == 0
  ) {
    
    stop(
      "Could not find any of: ",
      paste(
        paths,
        collapse = ", "
      )
    )
  }
  
  source(
    hit[1]
  )
}


source_first_existing(
  c(
    "Dist_GompertzLB.R",
    "Distributions/Dist_GompertzLB.R"
  )
)


source_first_existing(
  c(
    "Dist_Gompertz.R",
    "Distributions/Dist_Gompertz.R"
  )
)


source_first_existing(
  c(
    "Dist_GompertzNim.R",
    "Distributions/Dist_GompertzNim.R"
  )
)


## ============================================================
## 3. LIBRARIES
## ============================================================

library(tidyverse)
library(nimble)
library(MCMCvis)
library(coda)


cat("\n")
cat("============================================================\n")
cat("NO-FEEDING REFAGE10 MODEL\n")
cat("5-YEAR CORRECTED COHORT DENSITY\n")
cat("============================================================\n")


## ============================================================
## 4. BASIC CHECKS
## ============================================================

stopifnot(
  exists("growth_info"),
  exists("cintB_adj"),
  exists("cintD_adj"),
  exists("censoredD_vec"),
  exists("actual_surveys"),
  exists("ordered_islands"),
  exists("island_idx"),
  exists("island_start_time"),
  exists("island_start_idx"),
  exists("cum_effort_for_model"),
  exists("y"),
  exists("sex"),
  exists("Origin"),
  exists("tMax"),
  exists("nind"),
  exists("n_islands")
)


stopifnot(
  nrow(
    growth_info
  ) == nind
)


## ============================================================
## 5. REBUILD CLEAN BIRTH WINDOWS
##
## This removes the old blanket adjustment that incorrectly
## shortened birth intervals for some right-censored animals.
## ============================================================

cintB_adj_original <- cintB_adj


cintB_cal_clean <- matrix(
  NA_real_,
  nrow = nind,
  ncol = 2
)


birth_source <- character(
  nind
)


birth_svl_capture_conflict <- logical(
  nind
)


for(i in seq_len(nind)) {
  
  ## ----------------------------------------------------------
  ## KNOWN BIRTH / KNOWN AGE
  ## ----------------------------------------------------------
  
  if(
    !is.na(
      growth_info$known_birth_cal[i]
    )
  ) {
    
    cintB_cal_clean[i, ] <-
      growth_info$known_birth_cal[i] +
      c(
        -1,
        0
      )
    
    
    birth_source[i] <-
      "Known birth"
    
  } else {
    
    ## --------------------------------------------------------
    ## UNKNOWN BIRTH
    ## --------------------------------------------------------
    
    has_svl <-
      is.finite(
        growth_info$first_svl[i]
      ) &&
      is.finite(
        growth_info$first_svl_year[i]
      )
    
    
    if(
      has_svl
    ) {
      
      ref_yr <-
        growth_info$first_svl_year[i]
      
    } else {
      
      ref_yr <-
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
    
    
    if(
      has_svl
    ) {
      
      birth_svl_capture_conflict[i] <-
        raw_high >
        growth_info$tF_cal[i]
      
      
      if(
        birth_svl_capture_conflict[i]
      ) {
        
        birth_source[i] <-
          "SVL informed - constrained by first capture"
        
      } else {
        
        birth_source[i] <-
          "SVL informed"
      }
      
    } else {
      
      birth_source[i] <-
        "Unknown birth - no SVL"
    }
    
    
    ## --------------------------------------------------------
    ## FIRST CAPTURE IS HARD EVIDENCE OF EXISTENCE
    ## --------------------------------------------------------
    
    final_high <- min(
      raw_high,
      growth_info$tF_cal[i]
    )
    
    
    final_low <- raw_low
    
    
    if(
      final_low >=
      final_high
    ) {
      
      final_low <-
        final_high - 1
    }
    
    
    cintB_cal_clean[i, ] <- c(
      final_low,
      final_high
    )
  }
  
  
  if(
    any(
      !is.finite(
        cintB_cal_clean[i, ]
      )
    )
  ) {
    
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
## 6. CHECK ORIGIN
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


if(
  Origin_rebuilt !=
  Origin
) {
  
  stop(
    "Rebuilt birth windows alter Origin."
  )
}


## ============================================================
## 7. SHIFT BIRTH WINDOWS TO MODEL TIME
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
## 8. BIRTH-WINDOW AUDIT
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


if(
  length(
    birth_changed
  ) > 0
) {
  
  birth_change_audit <- tibble(
    i =
      birth_changed,
    
    ID =
      growth_info$ID_raw[
        birth_changed
      ],
    
    island =
      growth_info$island[
        birth_changed
      ],
    
    censor =
      censoredD_vec[
        birth_changed
      ],
    
    old_low =
      cintB_adj_original[
        birth_changed,
        1
      ] +
      Origin -
      1,
    
    old_high =
      cintB_adj_original[
        birth_changed,
        2
      ] +
      Origin -
      1,
    
    clean_low =
      cintB_adj_clean[
        birth_changed,
        1
      ] +
      Origin -
      1,
    
    clean_high =
      cintB_adj_clean[
        birth_changed,
        2
      ] +
      Origin -
      1
  )
}


cat(
  "\nBirth windows changed:",
  length(
    birth_changed
  ),
  "\n"
)


if(
  nrow(
    birth_change_audit
  ) > 0
) {
  
  birth_change_audit %>%
    count(
      censor
    ) %>%
    print()
}


## ============================================================
## 9. USE CLEAN BIRTH WINDOWS
## ============================================================

cintB_adj <-
  cintB_adj_clean


growth_info <- growth_info %>%
  mutate(
    birth_source =
      birth_source,
    
    birth_svl_capture_conflict =
      birth_svl_capture_conflict,
    
    birth_no_svl =
      is.na(
        known_birth_cal
      ) &
      !is.finite(
        first_svl
      ),
    
    birth_weak_or_adjusted =
      birth_no_svl |
      birth_svl_capture_conflict,
    
    birth_lower_cal =
      cintB_adj[, 1] +
      Origin -
      1,
    
    birth_upper_cal =
      cintB_adj[, 2] +
      Origin -
      1,
    
    birth_window_width =
      cintB_adj[, 2] -
      cintB_adj[, 1]
  )


cat("\nBirth-source summary:\n")


growth_info %>%
  count(
    birth_source
  ) %>%
  print(
    n = Inf
  )


cat(
  "\nWeak/adjusted birth histories:",
  sum(
    growth_info$birth_weak_or_adjusted
  ),
  "\n"
)


## ============================================================
## 10. VALIDATE BIRTH / DEATH HISTORIES
## ============================================================

stopifnot(
  all(
    cintB_adj[, 1] <
      cintB_adj[, 2]
  ),
  
  all(
    cintB_adj[, 1] >=
      1
  ),
  
  all(
    growth_info$birth_upper_cal <=
      growth_info$tF_cal
  )
)


impossible_interval_birth <- which(
  censoredD_vec == 1 &
    cintB_adj[, 1] >=
    cintD_adj[, 2]
)


if(
  length(
    impossible_interval_birth
  ) > 0
) {
  
  stop(
    "Impossible interval-censored birth/death history."
  )
}


## ============================================================
## 11. REFERENCE AGE
## ============================================================

age_ref <- 10


cat(
  "\nReference age:",
  age_ref,
  "years\n"
)


## ============================================================
## 12. SEX
##
## Female = 0
## Male   = 1
##
## No centering.
## ============================================================

sex01 <- as.numeric(
  sex
)


stopifnot(
  length(
    sex01
  ) == nind,
  
  !anyNA(
    sex01
  ),
  
  all(
    sex01 %in%
      c(
        0,
        1
      )
  )
)


cat("\nSex distribution:\n")


print(
  table(
    sex01
  )
)


## ============================================================
## 13. HABITAT AREA
## ============================================================

habitat <- read_csv(
  "Data/Iguana_islands_habitat.csv",
  show_col_types = FALSE
) %>%
  transmute(
    island =
      islands,
    
    usable_m2 =
      Useable,
    
    usable_ha =
      usable_m2 /
      10000
  )


stopifnot(
  all(
    ordered_islands %in%
      habitat$island
  )
)


## ============================================================
## 14. INDIVIDUAL KNOWN LIFE-HISTORY YEARS
## ============================================================

life_years <- growth_info %>%
  transmute(
    ID,
    ID_raw,
    island,
    
    first_year =
      tF_cal,
    
    last_year =
      tL_cal
  )


stopifnot(
  !anyNA(
    life_years$first_year
  ),
  
  !anyNA(
    life_years$last_year
  ),
  
  all(
    life_years$last_year >=
      life_years$first_year
  )
)


## ============================================================
## 15. RECONSTRUCT MINIMUM NUMBER ALIVE
##
## Each individual contributes from first known year
## through last known year inclusive.
## ============================================================

mna_calendar <- life_years %>%
  rowwise() %>%
  mutate(
    calendar_year = list(
      seq(
        first_year,
        last_year
      )
    )
  ) %>%
  unnest(
    calendar_year
  ) %>%
  ungroup() %>%
  count(
    island,
    calendar_year,
    name = "MNA"
  )


## ============================================================
## 16. ACTUAL SURVEY YEARS
## ============================================================

survey_years <- actual_surveys %>%
  transmute(
    island =
      Island,
    
    year_idx =
      yr_adj,
    
    calendar_year =
      yr_adj +
      Origin -
      1
  ) %>%
  distinct(
    island,
    year_idx,
    calendar_year
  )


## ============================================================
## 17. MNA DENSITY FOR SURVEY YEARS
## ============================================================

mna_survey <- survey_years %>%
  left_join(
    mna_calendar,
    by = c(
      "island",
      "calendar_year"
    )
  ) %>%
  left_join(
    habitat,
    by = "island"
  ) %>%
  mutate(
    MNA =
      replace_na(
        MNA,
        0L
      ),
    
    density_raw =
      MNA /
      usable_ha,
    
    log_density_raw =
      log1p(
        density_raw
      )
  ) %>%
  group_by(
    island
  ) %>%
  arrange(
    calendar_year,
    .by_group = TRUE
  ) %>%
  mutate(
    survey_number =
      row_number()
  ) %>%
  ungroup()


## ============================================================
## 18. DEFINE START-UP PERIOD
##
## First five ACTUAL survey years on each island are treated
## as monitoring / ascertainment start-up.
##
## Annual MNA density is trusted from survey year 6 onward.
## ============================================================

startup_n_surveys <- 5


mna_survey <- mna_survey %>%
  mutate(
    monitoring_phase =
      if_else(
        survey_number <=
          startup_n_surveys,
        
        "Start-up",
        
        "Established"
      )
  )


## ============================================================
## 19. ESTABLISHED ANNUAL DENSITY DATA
## ============================================================

density_established_observed <- mna_survey %>%
  filter(
    survey_number >
      startup_n_surveys
  )


established_support <- density_established_observed %>%
  group_by(
    island
  ) %>%
  summarise(
    n_established_years =
      n(),
    
    first_established_year =
      min(
        calendar_year
      ),
    
    last_established_year =
      max(
        calendar_year
      ),
    
    mean_MNA =
      mean(
        MNA
      ),
    
    median_MNA =
      median(
        MNA
      ),
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("ESTABLISHED MONITORING SUPPORT\n")
cat("============================================================\n")


print(
  established_support,
  n = Inf
)


if(
  nrow(
    established_support
  ) !=
  n_islands
) {
  
  stop(
    "At least one island has no established monitoring years."
  )
}


if(
  any(
    established_support$n_established_years <
    2
  )
) {
  
  stop(
    "At least one island has fewer than two established survey years."
  )
}


## ============================================================
## 20. ESTABLISHED ISLAND DENSITY BASELINE
##
## Calculated on log1p(density) scale.
## ============================================================

stable_island_density <- density_established_observed %>%
  group_by(
    island
  ) %>%
  summarise(
    stable_log_density =
      mean(
        log_density_raw,
        na.rm = TRUE
      ),
    
    stable_raw_density =
      expm1(
        stable_log_density
      ),
    
    n_established =
      n(),
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("ESTABLISHED ISLAND DENSITY BASELINES\n")
cat("============================================================\n")


print(
  stable_island_density,
  n = Inf
)


stopifnot(
  nrow(
    stable_island_density
  ) == n_islands,
  
  !anyNA(
    stable_island_density$stable_log_density
  )
)


## ============================================================
## 21. BUILD CORRECTED ISLAND x YEAR DENSITY TABLE
##
## Established survey year:
##   observed annual MNA density
##
## Otherwise:
##   established island mean
##
## Includes:
##   pre-monitoring years
##   first five survey years
##   unsurveyed gaps
## ============================================================

density_complete <- expand_grid(
  island =
    ordered_islands,
  
  year_idx =
    seq_len(
      tMax
    )
) %>%
  mutate(
    calendar_year =
      year_idx +
      Origin -
      1
  ) %>%
  left_join(
    mna_survey %>%
      select(
        island,
        year_idx,
        calendar_year,
        survey_number,
        MNA,
        density_raw,
        log_density_raw,
        monitoring_phase
      ),
    by = c(
      "island",
      "year_idx",
      "calendar_year"
    )
  ) %>%
  left_join(
    stable_island_density,
    by = "island"
  ) %>%
  mutate(
    trusted_annual_density =
      !is.na(
        survey_number
      ) &
      survey_number >
      startup_n_surveys,
    
    density_source =
      if_else(
        trusted_annual_density,
        "Established annual MNA",
        "Established island mean"
      ),
    
    log_density_used =
      if_else(
        trusted_annual_density,
        log_density_raw,
        stable_log_density
      ),
    
    density_used_raw =
      expm1(
        log_density_used
      )
  )


stopifnot(
  !anyNA(
    density_complete$log_density_used
  )
)


## ============================================================
## 22. STANDARDISE DENSITY
##
## Mean and SD come ONLY from trusted established annual
## observations.
## ============================================================

dens_mu <- mean(
  density_established_observed$
    log_density_raw,
  na.rm = TRUE
)


dens_sd <- sd(
  density_established_observed$
    log_density_raw,
  na.rm = TRUE
)


if(
  !is.finite(
    dens_sd
  ) ||
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
        log_density_used -
          dens_mu
      ) /
      dens_sd
  )


cat("\n")
cat("============================================================\n")
cat("CORRECTED DENSITY SCALE\n")
cat("============================================================\n")


cat(
  "Mean established log-density:",
  dens_mu,
  "\n"
)


cat(
  "SD established log-density:",
  dens_sd,
  "\n"
)


print(
  summary(
    density_complete$density_z
  )
)


## ============================================================
## 23. DENSITY SOURCE AUDIT
## ============================================================

density_source_summary <- density_complete %>%
  count(
    island,
    density_source
  ) %>%
  group_by(
    island
  ) %>%
  mutate(
    proportion =
      n /
      sum(n)
  ) %>%
  ungroup()


cat("\n")
cat("============================================================\n")
cat("DENSITY SOURCE SUMMARY\n")
cat("============================================================\n")


print(
  density_source_summary,
  n = Inf
)


## ============================================================
## 24. START-UP VERSUS ESTABLISHED DENSITY
## ============================================================

startup_density_comparison <- mna_survey %>%
  left_join(
    stable_island_density,
    by = "island"
  ) %>%
  filter(
    survey_number <=
      startup_n_surveys
  ) %>%
  group_by(
    island
  ) %>%
  summarise(
    first5_mean_raw_density =
      mean(
        density_raw
      ),
    
    first5_median_raw_density =
      median(
        density_raw
      ),
    
    established_baseline =
      first(
        stable_raw_density
      ),
    
    first5_to_established_ratio =
      first5_mean_raw_density /
      established_baseline,
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("START-UP VS ESTABLISHED DENSITY\n")
cat("============================================================\n")


print(
  startup_density_comparison,
  n = Inf
)


## ============================================================
## 25. BUILD DENSITY MATRIX
## ============================================================

density_matrix <- matrix(
  NA_real_,
  nrow = n_islands,
  ncol = tMax,
  dimnames = list(
    ordered_islands,
    seq_len(
      tMax
    )
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


## ============================================================
## 26. SAVE DENSITY AUDIT TABLES
## ============================================================

write_csv(
  density_complete,
  "Corrected_Cohort_Density_Island_Year_NOFEED.csv"
)


write_csv(
  stable_island_density,
  "Corrected_Cohort_Density_Island_Baselines_NOFEED.csv"
)


write_csv(
  startup_density_comparison,
  "Corrected_Cohort_Density_Startup_Comparison_NOFEED.csv"
)


## ============================================================
## 27. TEMPORAL COHORT STRUCTURE
##
## Raw 5-year blocks:
##   15 blocks from 1946 onward
##
## Model blocks:
##   1  = 1946-1965
##   2  = 1966-1970
##   3  = 1971-1975
##   4  = 1976-1980
##   5  = 1981-1985
##   6  = 1986-1990
##   7  = 1991-1995
##   8  = 1996-2000
##   9  = 2001-2005
##   10 = 2006-2010
##   11 = 2011-2019
## ============================================================

block_width <- 5


n_time_blocks_raw <- ceiling(
  tMax /
    block_width
)


if(
  n_time_blocks_raw != 15
) {
  
  stop(
    "Expected 15 raw temporal blocks but found ",
    n_time_blocks_raw
  )
}


n_time_blocks <- 11


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
          raw_block -
            3
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


model_temporal_table <- raw_temporal_table %>%
  group_by(
    model_block
  ) %>%
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
cat("TEMPORAL BLOCKS\n")
cat("============================================================\n")


print(
  model_temporal_table,
  n = Inf
)


## ============================================================
## 28. FINAL SETUP CHECKS
## ============================================================

stopifnot(
  length(
    sex01
  ) == nind,
  
  nrow(
    density_matrix
  ) == n_islands,
  
  ncol(
    density_matrix
  ) == tMax,
  
  !anyNA(
    density_matrix
  ),
  
  n_time_blocks == 11,
  
  nrow(
    cum_effort_for_model
  ) == n_islands,
  
  ncol(
    cum_effort_for_model
  ) >=
    tMax + 1
)


cat("\n")
cat("============================================================\n")
cat("SETUP COMPLETE\n")
cat("============================================================\n")
cat("Individuals:", nind, "\n")
cat("Islands:", n_islands, "\n")
cat("Reference age:", age_ref, "\n")
cat("Density start-up surveys:", startup_n_surveys, "\n")
cat("Temporal blocks:", n_time_blocks, "\n")
cat("============================================================\n")


## ============================================================
## 29. NIMBLE MODEL
## ============================================================

code_REF10_NOFEED <- nimbleCode({
  
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
    ## LATENT COHORT YEAR
    ## --------------------------------------------------------
    
    cohort_idx[i] <- max(
      1,
      min(
        tMax,
        floor(
          tB[i]
        ) +
          1
      )
    )
    
    
    ## --------------------------------------------------------
    ## TEMPORAL BLOCK
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
    
    
    time_block[i] <- min(
      n_time_blocks,
      max(
        1,
        time_block_raw[i] -
          3
      )
    )
    
    
    ## --------------------------------------------------------
    ## CORRECTED ISLAND x COHORT-YEAR DENSITY
    ## --------------------------------------------------------
    
    dens_cohort[i] <-
      density_matrix[
        island_idx[i],
        cohort_idx[i]
      ]
    
    
    ## ========================================================
    ## MORTALITY AT AGE 10
    ##
    ## log(a10)
    ## ========================================================
    
    log_a10_ind[i] <-
      log_a10 +
      betaSEX[1] *
      sex01[i] +
      betaDENS_a *
      dens_cohort[i] +
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
    ##
    ## log(b)
    ## ========================================================
    
    log_b_ind[i] <-
      log_b +
      betaSEX[2] *
      sex01[i] +
      betaDENS_b *
      dens_cohort[i]
    
    
    bmult[i] <- exp(
      log_b_ind[i]
    )
    
    
    ## --------------------------------------------------------
    ## CONVERT AGE-10 PARAMETER TO STANDARD GOMPERTZ a0
    ##
    ## a0 = a10 * exp(-b * 10)
    ## --------------------------------------------------------
    
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
        floor(
          tB[i]
        ) +
          1
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
        floor(
          tD[i]
        ) +
          1
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
  ## 30. FOCAL BIOLOGICAL EFFECTS
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
  
  
  ## ==========================================================
  ## 31. RESIDUAL ISLAND HETEROGENEITY ON a10
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
  ## 32. SHARED TEMPORAL COHORT HETEROGENEITY ON a10
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
  ## 33. REFERENCE GOMPERTZ PARAMETERS
  ##
  ## Reference:
  ## female
  ## average corrected cohort density
  ## zero island deviation
  ## zero temporal deviation
  ## ==========================================================
  
  log_a10 ~ dnorm(
    log(0.015),
    sd = 2
  )
  
  
  log_b ~ dnorm(
    log(0.20),
    sd = 1
  )
  
  
  a10 <- exp(
    log_a10
  )
  
  
  b <- exp(
    log_b
  )
  
  
  ## ==========================================================
  ## 34. DETECTION
  ## ==========================================================
  
  mean.p ~ dunif(
    0,
    1
  )
})


## ============================================================
## 35. CONSTANTS
## ============================================================

consts_REF10_NOFEED <- list(
  nind =
    nind,
  
  tMax =
    tMax,
  
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
  
  sex01 =
    as.numeric(
      sex01
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
## 36. DATA
## ============================================================

data_REF10_NOFEED <- list(
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
## 37. INITIAL VALUE FUNCTION
## ============================================================

make_inits_REF10_NOFEED <- function(
    chain_id = 1
) {
  
  set.seed(
    16000 +
      chain_id
  )
  
  
  tB_start <- runif(
    nind,
    min = cintB_adj[, 1],
    max = cintB_adj[, 2]
  )
  
  
  tD_start <- numeric(
    nind
  )
  
  
  for(i in seq_len(nind)) {
    
    ## --------------------------------------------------------
    ## INTERVAL-CENSORED DEATH
    ## --------------------------------------------------------
    
    if(
      censoredD_vec[i] == 1
    ) {
      
      loD <-
        cintD_adj[i, 1] +
        0.001
      
      
      hiD <-
        cintD_adj[i, 2] -
        0.001
      
      
      if(
        hiD <=
        loD
      ) {
        
        stop(
          "Invalid death interval for ",
          growth_info$ID_raw[i]
        )
      }
      
      
      if(
        tB_start[i] >=
        hiD
      ) {
        
        possible_birth_high <- min(
          cintB_adj[i, 2],
          hiD -
            0.01
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
        tB_start[i] +
          0.01
      )
      
      
      if(
        loD_use >=
        hiD
      ) {
        
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
    
    if(
      censoredD_vec[i] == 2
    ) {
      
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
    
    
    if(
      !censoredD_vec[i] %in%
      c(
        1,
        2
      )
    ) {
      
      stop(
        "Unexpected censor code for ",
        growth_info$ID_raw[i],
        ": ",
        censoredD_vec[i]
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
    
    log_a10 =
      rnorm(
        1,
        log(0.005),
        0.20
      ),
    
    log_b =
      rnorm(
        1,
        log(0.18),
        0.08
      ),
    
    betaSEX =
      c(
        rnorm(
          1,
          0.35,
          0.08
        ),
        
        rnorm(
          1,
          0.17,
          0.05
        )
      ),
    
    betaDENS_a =
      rnorm(
        1,
        -0.8,
        0.10
      ),
    
    betaDENS_b =
      rnorm(
        1,
        -0.10,
        0.05
      ),
    
    sigma_island_a =
      runif(
        1,
        0.7,
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
        0.23,
        0.28
      )
  )
}


inits_REF10_NOFEED <- lapply(
  1:3,
  make_inits_REF10_NOFEED
)


## ============================================================
## 38. BUILD MODEL
## ============================================================

model_REF10_NOFEED <- nimbleModel(
  code =
    code_REF10_NOFEED,
  
  constants =
    consts_REF10_NOFEED,
  
  data =
    data_REF10_NOFEED,
  
  inits =
    inits_REF10_NOFEED[[1]],
  
  calculate = FALSE
)


lp_REF10_NOFEED <-
  model_REF10_NOFEED$calculate()


cat(
  "\nInitial log probability:",
  lp_REF10_NOFEED,
  "\n"
)


if(
  !is.finite(
    lp_REF10_NOFEED
  )
) {
  
  stop(
    "Initial model log probability is not finite."
  )
}


## ============================================================
## 39. CONFIGURE MCMC
## ============================================================

conf_REF10_NOFEED <- configureMCMC(
  model_REF10_NOFEED,
  
  monitors = c(
    "a10",
    "b",
    
    "log_a10",
    "log_b",
    
    "betaSEX",
    
    "betaDENS_a",
    "betaDENS_b",
    
    "sigma_island_a",
    "sigma_time",
    
    "island_a",
    "u_time",
    
    "mean.p"
  ),
  
  enableWAIC = FALSE
)


## ============================================================
## 40. JOINT a10-SIDE SAMPLER
## ============================================================

a_pars_NOFEED <- c(
  "log_a10",
  "betaSEX[1]",
  "betaDENS_a"
)


conf_REF10_NOFEED$removeSamplers(
  a_pars_NOFEED
)


conf_REF10_NOFEED$addSampler(
  target =
    a_pars_NOFEED,
  
  type =
    "AF_slice"
)


## ============================================================
## 41. JOINT b-SIDE SAMPLER
## ============================================================

b_pars_NOFEED <- c(
  "log_b",
  "betaSEX[2]",
  "betaDENS_b"
)


conf_REF10_NOFEED$removeSamplers(
  b_pars_NOFEED
)


conf_REF10_NOFEED$addSampler(
  target =
    b_pars_NOFEED,
  
  type =
    "AF_slice"
)


## ============================================================
## 42. VARIANCE PARAMETERS
## ============================================================

for(v in c(
  "sigma_island_a",
  "sigma_time"
)) {
  
  conf_REF10_NOFEED$removeSamplers(
    v
  )
  
  
  conf_REF10_NOFEED$addSampler(
    target =
      v,
    
    type =
      "slice"
  )
}


## ============================================================
## 43. JOINT LATENT BIRTH + LIFESPAN SAMPLERS
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
  
  
  conf_REF10_NOFEED$removeSamplers(
    life_nodes
  )
  
  
  conf_REF10_NOFEED$addSampler(
    target =
      life_nodes,
    
    type =
      "AF_slice"
  )
}


## ============================================================
## 44. BUILD + COMPILE
## ============================================================

mcmc_REF10_NOFEED <- buildMCMC(
  conf_REF10_NOFEED
)


cModel_REF10_NOFEED <- compileNimble(
  model_REF10_NOFEED,
  resetFunctions = TRUE
)


cMCMC_REF10_NOFEED <- compileNimble(
  mcmc_REF10_NOFEED,
  project = model_REF10_NOFEED
)


## ============================================================
## 45. RUN 100k NO-FEEDING MODEL
## ============================================================

cat("\n")
cat("============================================================\n")
cat("RUNNING NO-FEEDING 100k MODEL\n")
cat("============================================================\n")
cat("Reference age: 10 years\n")
cat("Density: corrected 5-survey-year start-up\n")
cat("Iterations: 100000\n")
cat("Burn-in: 30000\n")
cat("Chains: 3\n")
cat("Thin: 10\n")
cat("============================================================\n")


run_REF10_NOFEED <- runMCMC(
  cMCMC_REF10_NOFEED,
  
  niter = 100000,
  nburnin = 30000,
  nchains = 3,
  thin = 10,
  
  inits =
    inits_REF10_NOFEED,
  
  setSeed = TRUE,
  
  samplesAsCodaMCMC = TRUE,
  
  summary = TRUE,
  
  WAIC = FALSE,
  
  progressBar = TRUE
)


## ============================================================
## 46. SAVE RAW FIT IMMEDIATELY
## ============================================================

saveRDS(
  run_REF10_NOFEED,
  "Model_PaceOfLife_REFAGE10_DENSITY_CORRECTED_NOFEED_100k.rds"
)


cat(
  "\nSaved raw fit.\n"
)


## ============================================================
## 47. PRIMARY SUMMARY
## ============================================================

sm_REF10_NOFEED <- MCMCsummary(
  run_REF10_NOFEED$samples,
  
  params = c(
    "a10",
    "b",
    
    "log_a10",
    "log_b",
    
    "betaSEX",
    
    "betaDENS_a",
    "betaDENS_b",
    
    "sigma_island_a",
    "sigma_time",
    
    "mean.p"
  )
)


cat("\n")
cat("============================================================\n")
cat("NO-FEEDING POSTERIOR SUMMARY\n")
cat("============================================================\n")


print(
  sm_REF10_NOFEED
)


## ============================================================
## 48. POSTERIOR MATRIX
## ============================================================

s_REF10_NOFEED <- as.matrix(
  run_REF10_NOFEED$samples
)


## ============================================================
## 49. BIOLOGICAL EFFECT SUMMARY
## ============================================================

effect_REF10_NOFEED <- tibble(
  effect = c(
    "Sex -> mortality at age 10",
    "Sex -> senescence",
    "Density -> mortality at age 10",
    "Density -> senescence"
  ),
  
  parameter = c(
    "betaSEX[1]",
    "betaSEX[2]",
    "betaDENS_a",
    "betaDENS_b"
  )
) %>%
  rowwise() %>%
  mutate(
    median_beta =
      median(
        s_REF10_NOFEED[
          ,
          parameter
        ]
      ),
    
    low_beta =
      unname(
        quantile(
          s_REF10_NOFEED[
            ,
            parameter
          ],
          0.025
        )
      ),
    
    high_beta =
      unname(
        quantile(
          s_REF10_NOFEED[
            ,
            parameter
          ],
          0.975
        )
      ),
    
    median_ratio =
      median(
        exp(
          s_REF10_NOFEED[
            ,
            parameter
          ]
        )
      ),
    
    low_ratio =
      unname(
        quantile(
          exp(
            s_REF10_NOFEED[
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
            s_REF10_NOFEED[
              ,
              parameter
            ]
          ),
          0.975
        )
      ),
    
    Pr_positive =
      mean(
        s_REF10_NOFEED[
          ,
          parameter
        ] >
          0
      ),
    
    Pr_negative =
      mean(
        s_REF10_NOFEED[
          ,
          parameter
        ] <
          0
      )
  ) %>%
  ungroup()


cat("\n")
cat("============================================================\n")
cat("NO-FEEDING BIOLOGICAL EFFECTS\n")
cat("============================================================\n")


print(
  effect_REF10_NOFEED,
  n = Inf
)


## ============================================================
## 50. POSTERIOR CORRELATIONS
## ============================================================

cor_a_b_NOFEED <- cor(
  s_REF10_NOFEED[
    ,
    "log_a10"
  ],
  
  s_REF10_NOFEED[
    ,
    "log_b"
  ]
)


cor_b_matrix_NOFEED <- cor(
  s_REF10_NOFEED[
    ,
    c(
      "log_b",
      "betaSEX[2]",
      "betaDENS_b"
    )
  ]
)


cat("\n")
cat("============================================================\n")
cat("NO-FEEDING POSTERIOR CORRELATIONS\n")
cat("============================================================\n")


cat(
  "log_a10 vs log_b:",
  round(
    cor_a_b_NOFEED,
    3
  ),
  "\n"
)


cat("\n")


print(
  round(
    cor_b_matrix_NOFEED,
    3
  )
)


## ============================================================
## 51. CONVERGENCE
## ============================================================

conv_REF10_NOFEED <- sm_REF10_NOFEED %>%
  as.data.frame() %>%
  rownames_to_column(
    "parameter"
  ) %>%
  as_tibble() %>%
  arrange(
    desc(
      Rhat
    )
  )


cat("\n")
cat("============================================================\n")
cat("NO-FEEDING CONVERGENCE\n")
cat("============================================================\n")


print(
  conv_REF10_NOFEED,
  n = Inf
)


cat(
  "\nMaximum Rhat:",
  max(
    conv_REF10_NOFEED$Rhat,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Minimum ESS:",
  min(
    conv_REF10_NOFEED$n.eff,
    na.rm = TRUE
  ),
  "\n"
)


## ============================================================
## 52. TEMPORAL EFFECTS
## ============================================================

sm_time_REF10_NOFEED <- MCMCsummary(
  run_REF10_NOFEED$samples,
  params = "u_time"
)


cat("\n")
cat("============================================================\n")
cat("NO-FEEDING TEMPORAL EFFECTS\n")
cat("============================================================\n")


print(
  sm_time_REF10_NOFEED
)


## ============================================================
## 53. ISLAND EFFECTS
## ============================================================

sm_island_REF10_NOFEED <- MCMCsummary(
  run_REF10_NOFEED$samples,
  params = "island_a"
)


cat("\n")
cat("============================================================\n")
cat("NO-FEEDING ISLAND EFFECTS\n")
cat("============================================================\n")


print(
  sm_island_REF10_NOFEED
)


## ============================================================
## 54. OPTIONAL:
## COMPARE AGAINST EXISTING 5-YEAR FEEDING MODEL
##
## Only runs if the previous raw-fit RDS exists.
## ============================================================

feeding_fit_file <-
  "Model_PaceOfLife_REFAGE10_DENSITY_CORRECTED_100k.rds"


if(
  file.exists(
    feeding_fit_file
  )
) {
  
  cat("\n")
  cat("============================================================\n")
  cat("COMPARING FEEDING VS NO-FEEDING MODELS\n")
  cat("============================================================\n")
  
  
  run_REF10_FEED_COMPARE <- readRDS(
    feeding_fit_file
  )
  
  
  s_REF10_FEED_COMPARE <- as.matrix(
    run_REF10_FEED_COMPARE$samples
  )
  
  
  compare_model_parameter <- function(
    samples,
    parameter,
    model_name
  ) {
    
    x <-
      samples[
        ,
        parameter
      ]
    
    
    tibble(
      model =
        model_name,
      
      parameter =
        parameter,
      
      median =
        median(
          x
        ),
      
      low =
        unname(
          quantile(
            x,
            0.025
          )
        ),
      
      high =
        unname(
          quantile(
            x,
            0.975
          )
        ),
      
      Pr_positive =
        mean(
          x >
            0
        ),
      
      Pr_negative =
        mean(
          x <
            0
        )
    )
  }
  
  
  comparison_FEED_NOFEED <- bind_rows(
    compare_model_parameter(
      s_REF10_FEED_COMPARE,
      "betaSEX[1]",
      "Feeding included"
    ),
    
    compare_model_parameter(
      s_REF10_NOFEED,
      "betaSEX[1]",
      "No feeding"
    ),
    
    compare_model_parameter(
      s_REF10_FEED_COMPARE,
      "betaSEX[2]",
      "Feeding included"
    ),
    
    compare_model_parameter(
      s_REF10_NOFEED,
      "betaSEX[2]",
      "No feeding"
    ),
    
    compare_model_parameter(
      s_REF10_FEED_COMPARE,
      "betaDENS_a",
      "Feeding included"
    ),
    
    compare_model_parameter(
      s_REF10_NOFEED,
      "betaDENS_a",
      "No feeding"
    ),
    
    compare_model_parameter(
      s_REF10_FEED_COMPARE,
      "betaDENS_b",
      "Feeding included"
    ),
    
    compare_model_parameter(
      s_REF10_NOFEED,
      "betaDENS_b",
      "No feeding"
    )
  )
  
  
  print(
    comparison_FEED_NOFEED,
    n = Inf
  )
  
  
  write_csv(
    comparison_FEED_NOFEED,
    "Sensitivity_Feeding_vs_NoFeeding.csv"
  )
  
} else {
  
  comparison_FEED_NOFEED <-
    NULL
  
  
  cat(
    "\nPrevious feeding-model RDS not found; ",
    "skipping direct comparison.\n"
  )
}


## ============================================================
## 55. SAVE COMPLETE ANALYSIS
## ============================================================

write_csv(
  effect_REF10_NOFEED,
  "Model_PaceOfLife_REFAGE10_DENSITY_CORRECTED_NOFEED_100k_effects.csv"
)


saveRDS(
  list(
    fit =
      run_REF10_NOFEED,
    
    summary =
      sm_REF10_NOFEED,
    
    effects =
      effect_REF10_NOFEED,
    
    convergence =
      conv_REF10_NOFEED,
    
    temporal =
      sm_time_REF10_NOFEED,
    
    island =
      sm_island_REF10_NOFEED,
    
    density_table =
      density_complete,
    
    density_baselines =
      stable_island_density,
    
    density_support =
      established_support,
    
    density_startup_comparison =
      startup_density_comparison,
    
    density_source_summary =
      density_source_summary,
    
    density_mean =
      dens_mu,
    
    density_sd =
      dens_sd,
    
    startup_n_surveys =
      startup_n_surveys,
    
    age_ref =
      age_ref,
    
    birth_change_audit =
      birth_change_audit,
    
    weak_birth_n =
      sum(
        growth_info$birth_weak_or_adjusted
      ),
    
    correlation_a10_b =
      cor_a_b_NOFEED,
    
    b_correlation_matrix =
      cor_b_matrix_NOFEED,
    
    feeding_comparison =
      comparison_FEED_NOFEED,
    
    temporal_mapping =
      model_temporal_table
  ),
  
  "Model_PaceOfLife_REFAGE10_DENSITY_CORRECTED_NOFEED_100k_COMPLETE.rds"
)


cat("\n")
cat("============================================================\n")
cat("NO-FEEDING 100k MODEL COMPLETE\n")
cat("============================================================\n")