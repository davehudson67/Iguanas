## ============================================================
## IGUANA PACE-OF-LIFE MODEL
##
## WITHIN-ISLAND COHORT DENSITY
## + EXPLICIT TEMPORAL TREND
## + RESIDUAL TEMPORAL PERIOD VARIATION
##
## - Fresh-session setup
## - Corrected latent birth intervals
## - Reference-age Gompertz: age 10
## - Sex:
##      female = 0
##      male   = 1
## - NO feeding term
## - Density:
##      within-island deviation from established mean
##      first 5 survey years excluded from annual signal
##      unobserved/start-up years assigned deviation = 0
## - Explicit birth-cohort trend on a10:
##      effect per decade
## - Residual pooled temporal effects on a10
## - Residual island effect on a10
##
## DEVELOPMENT / STRUCTURAL CHECK:
##   niter   = 100000
##   burn-in = 30000
##   chains  = 3
##   thin    = 10
##
## ============================================================


## ============================================================
## 1. CANONICAL DATA PREPARATION
## ============================================================

source("DataPrep.R")


## ============================================================
## 2. LOAD CUSTOM DISTRIBUTIONS
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
cat("WITHIN-ISLAND DENSITY + TEMPORAL TREND MODEL\n")
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
## Removes the old blanket birth/death correction that
## incorrectly shortened some right-censored histories.
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
    ## FIRST CAPTURE IS HARD EVIDENCE THAT THE ANIMAL
    ## WAS ALREADY ALIVE
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
        final_high -
        1
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
## 7. SHIFT BIRTH WINDOWS INTO MODEL TIME
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
      cintB_adj_clean[i, 2] -
        1
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


cat("\n")
cat("Birth-source summary:\n")


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
##
## Age 10 remains fixed across all primary/sensitivity models.
##
## h(t) = a10 * exp[b * (t - 10)]
##
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


cat("\n")
cat("Sex coding:\n")


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
## 15. RECONSTRUCT MNA
##
## Individual contributes from first through last known year.
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
## 18. MONITORING START-UP RULE
##
## First five ACTUAL survey years on each island are treated as
## ascertainment/start-up.
##
## Trusted annual density begins at survey year 6.
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
## 19. TRUSTED ESTABLISHED DENSITY DATA
## ============================================================

density_established <- mna_survey %>%
  filter(
    survey_number >
      startup_n_surveys
  )


established_support <- density_established %>%
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


## ============================================================
## 20. ESTABLISHED ISLAND MEAN LOG-DENSITY
##
## This is the island-specific baseline.
## ============================================================

island_density_baseline <- density_established %>%
  group_by(
    island
  ) %>%
  summarise(
    mean_log_density =
      mean(
        log_density_raw,
        na.rm = TRUE
      ),
    
    mean_raw_density =
      expm1(
        mean_log_density
      ),
    
    n_established =
      n(),
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("ISLAND ESTABLISHED DENSITY BASELINES\n")
cat("============================================================\n")


print(
  island_density_baseline,
  n = Inf
)


stopifnot(
  nrow(
    island_density_baseline
  ) == n_islands,
  
  !anyNA(
    island_density_baseline$mean_log_density
  )
)


## ============================================================
## 21. WITHIN-ISLAND DENSITY DEVIATION
##
## For TRUSTED established survey years:
##
## D_within(j,t) =
##
## log(1 + density_jt)
## -
## mean_j[log(1 + density)]
##
##
## For:
##
## - pre-monitoring years
## - first five survey years
## - unsurveyed gaps
##
## deviation = 0.
##
## This means:
##
## "we have no evidence that this year's density differs
##  from the island's established typical density"
##
## rather than inventing annual information.
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
    island_density_baseline,
    by = "island"
  ) %>%
  mutate(
    trusted_annual_density =
      !is.na(
        survey_number
      ) &
      survey_number >
      startup_n_surveys,
    
    density_within_raw =
      if_else(
        trusted_annual_density,
        
        log_density_raw -
          mean_log_density,
        
        0
      ),
    
    density_source =
      if_else(
        trusted_annual_density,
        
        "Observed within-island deviation",
        
        "No annual deviation information"
      )
  )


stopifnot(
  !anyNA(
    density_complete$density_within_raw
  )
)


## ============================================================
## 22. STANDARDISE WITHIN-ISLAND DENSITY
##
## Scale uses only the actual observed established deviations.
##
## Because deviations are calculated around each island mean,
## their pooled mean is effectively zero.
## ============================================================

density_within_sd <- sd(
  density_complete$density_within_raw[
    density_complete$trusted_annual_density
  ],
  na.rm = TRUE
)


if(
  !is.finite(
    density_within_sd
  ) ||
  density_within_sd <= 0
) {
  
  stop(
    "Invalid within-island density SD."
  )
}


density_complete <- density_complete %>%
  mutate(
    density_within_z =
      density_within_raw /
      density_within_sd
  )


cat("\n")
cat("============================================================\n")
cat("WITHIN-ISLAND DENSITY SCALE\n")
cat("============================================================\n")


cat(
  "SD of trusted within-island log-density deviations:",
  density_within_sd,
  "\n"
)


print(
  summary(
    density_complete$density_within_z[
      density_complete$trusted_annual_density
    ]
  )
)


## ============================================================
## 23. CONFIRM WITHIN-ISLAND MEANS ARE APPROXIMATELY ZERO
## ============================================================

within_density_check <- density_complete %>%
  filter(
    trusted_annual_density
  ) %>%
  group_by(
    island
  ) %>%
  summarise(
    mean_deviation =
      mean(
        density_within_raw
      ),
    
    sd_deviation =
      sd(
        density_within_raw
      ),
    
    n =
      n(),
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("WITHIN-ISLAND DENSITY CHECK\n")
cat("============================================================\n")


print(
  within_density_check,
  n = Inf
)


## ============================================================
## 24. DENSITY-TIME AUDIT
##
## Trusted observations only.
## ============================================================

density_time_audit <- density_complete %>%
  filter(
    trusted_annual_density
  ) %>%
  group_by(
    island
  ) %>%
  summarise(
    n =
      n(),
    
    cor_density_year =
      if(
        n() >= 3
      ) {
        
        cor(
          calendar_year,
          density_within_raw
        )
        
      } else {
        
        NA_real_
      },
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("WITHIN-ISLAND DENSITY / TIME CORRELATION\n")
cat("============================================================\n")


print(
  density_time_audit,
  n = Inf
)


## ============================================================
## 25. BUILD WITHIN-ISLAND DENSITY MATRIX
## ============================================================

density_within_matrix <- matrix(
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
  
  
  density_within_matrix[jj, tt] <-
    density_complete$density_within_z[r]
}


stopifnot(
  !anyNA(
    density_within_matrix
  )
)


## ============================================================
## 26. TEMPORAL BLOCK STRUCTURE
##
## Residual temporal variation still uses:
##
## 1  = 1946-1965
## 2  = 1966-1970
## 3  = 1971-1975
## 4  = 1976-1980
## 5  = 1981-1985
## 6  = 1986-1990
## 7  = 1991-1995
## 8  = 1996-2000
## 9  = 2001-2005
## 10 = 2006-2010
## 11 = 2011-2019
##
## But these now represent deviations around the explicit
## temporal trend rather than the entire temporal signal.
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
    "Expected 15 raw temporal blocks; found ",
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
  ) %>%
  mutate(
    midpoint_year =
      (
        calendar_start +
          calendar_end
      ) /
      2
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
## 27. EXPLICIT TEMPORAL TREND
##
## Time is expressed in DECADES.
##
## Reference year = 1990.
##
## betaTIME therefore means:
##
## change in log mortality hazard at age 10 for animals
## born one decade later.
##
## exp(betaTIME) =
## hazard ratio per decade.
## ============================================================

time_ref_year <- 1990


time_scale_years <- 10


cat("\n")
cat("Temporal reference year:", time_ref_year, "\n")
cat("Temporal effect scale: per decade\n")


## ============================================================
## 28. FINAL SETUP CHECKS
## ============================================================

stopifnot(
  length(
    sex01
  ) == nind,
  
  nrow(
    density_within_matrix
  ) == n_islands,
  
  ncol(
    density_within_matrix
  ) == tMax,
  
  !anyNA(
    density_within_matrix
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
cat("Temporal reference year:", time_ref_year, "\n")
cat("Density start-up surveys:", startup_n_surveys, "\n")
cat("Within-density SD:", density_within_sd, "\n")
cat("Residual temporal blocks:", n_time_blocks, "\n")
cat("============================================================\n")


## ============================================================
## 29. NIMBLE MODEL
## ============================================================

code_REF10_WITHIN_TREND <- nimbleCode({
  
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
    ## CALENDAR BIRTH YEAR
    ##
    ## tB is on shifted model-time scale.
    ## --------------------------------------------------------
    
    birth_year[i] <-
      tB[i] +
      Origin -
      1
    
    
    ## --------------------------------------------------------
    ## CONTINUOUS TEMPORAL TREND
    ##
    ## Units = decades from 1990
    ## --------------------------------------------------------
    
    birth_decade_c[i] <-
      (
        birth_year[i] -
          time_ref_year
      ) /
      time_scale_years
    
    
    ## --------------------------------------------------------
    ## RESIDUAL TEMPORAL BLOCK
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
    ## WITHIN-ISLAND COHORT DENSITY
    ## --------------------------------------------------------
    
    dens_within[i] <-
      density_within_matrix[
        island_idx[i],
        cohort_idx[i]
      ]
    
    
    ## ========================================================
    ## MORTALITY AT AGE 10
    ##
    ## Effects:
    ##
    ## sex
    ## within-island cohort density
    ## continuous cohort trend
    ## residual temporal period effect
    ## residual island effect
    ## ========================================================
    
    log_a10_ind[i] <-
      log_a10 +
      betaSEX[1] *
      sex01[i] +
      betaDENS_a *
      dens_within[i] +
      betaTIME *
      birth_decade_c[i] +
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
    ## No temporal trend is added to b in this sensitivity.
    ##
    ## We are testing the temporal pattern that was observed
    ## specifically in mortality level.
    ## ========================================================
    
    log_b_ind[i] <-
      log_b +
      betaSEX[2] *
      sex01[i] +
      betaDENS_b *
      dens_within[i]
    
    
    bmult[i] <- exp(
      log_b_ind[i]
    )
    
    
    ## --------------------------------------------------------
    ## CONVERT REFERENCE-AGE GOMPERTZ TO STANDARD a0
    ##
    ## a0 = a10 * exp(-b * age_ref)
    ## --------------------------------------------------------
    
    log_amult[i] <-
      log_a10_ind[i] -
      bmult[i] *
      age_ref
    
    
    amult[i] <- exp(
      log_amult[i]
    )
    
    
    ## --------------------------------------------------------
    ## LATENT LIFESPAN / DEATH
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
  ## 30. BIOLOGICAL EFFECTS
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
  ## 31. LINEAR BIRTH-COHORT TREND
  ##
  ## betaTIME:
  ## log hazard ratio at age 10 per decade later birth.
  ## ==========================================================
  
  betaTIME ~ dnorm(
    0,
    sd = 1
  )
  
  
  HR_TIME_DECADE <- exp(
    betaTIME
  )
  
  
  ## ==========================================================
  ## 32. RESIDUAL ISLAND EFFECT ON a10
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
  ## 33. RESIDUAL TEMPORAL-PERIOD VARIATION
  ##
  ## These effects now represent departures from the explicit
  ## secular trend.
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
  ## 34. REFERENCE GOMPERTZ PARAMETERS
  ##
  ## log_a10 corresponds approximately to:
  ##
  ## female
  ## within-island density deviation = 0
  ## birth year 1990
  ## zero residual island effect
  ## zero residual temporal-period effect
  ##
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
  ## 35. DETECTION
  ## ==========================================================
  
  mean.p ~ dunif(
    0,
    1
  )
})


## ============================================================
## 36. CONSTANTS
## ============================================================

consts_REF10_WITHIN_TREND <- list(
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
  
  Origin =
    Origin,
  
  time_ref_year =
    time_ref_year,
  
  time_scale_years =
    time_scale_years,
  
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
## 37. DATA
## ============================================================

data_REF10_WITHIN_TREND <- list(
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
  
  density_within_matrix =
    as.matrix(
      density_within_matrix
    )
)


## ============================================================
## 38. INITIAL VALUES
## ============================================================

make_inits_REF10_WITHIN_TREND <- function(
    chain_id = 1
) {
  
  set.seed(
    21000 +
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
        log(0.004),
        0.20
      ),
    
    log_b =
      rnorm(
        1,
        log(0.17),
        0.08
      ),
    
    betaSEX =
      c(
        rnorm(
          1,
          0.36,
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
        -0.2,
        0.10
      ),
    
    betaDENS_b =
      rnorm(
        1,
        -0.05,
        0.05
      ),
    
    ## Approximate audit suggested a strong negative trend,
    ## but deliberately initialise chains differently.
    
    betaTIME =
      rnorm(
        1,
        -0.7,
        0.15
      ),
    
    sigma_island_a =
      runif(
        1,
        0.6,
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
        0.3,
        1.5
      ),
    
    u_time_raw =
      rnorm(
        n_time_blocks,
        0,
        0.3
      ),
    
    mean.p =
      runif(
        1,
        0.23,
        0.28
      )
  )
}


inits_REF10_WITHIN_TREND <- lapply(
  1:3,
  make_inits_REF10_WITHIN_TREND
)


## ============================================================
## 39. BUILD MODEL
## ============================================================

model_REF10_WITHIN_TREND <- nimbleModel(
  code =
    code_REF10_WITHIN_TREND,
  
  constants =
    consts_REF10_WITHIN_TREND,
  
  data =
    data_REF10_WITHIN_TREND,
  
  inits =
    inits_REF10_WITHIN_TREND[[1]],
  
  calculate = FALSE
)


lp_REF10_WITHIN_TREND <-
  model_REF10_WITHIN_TREND$calculate()


cat(
  "\nInitial log probability:",
  lp_REF10_WITHIN_TREND,
  "\n"
)


if(
  !is.finite(
    lp_REF10_WITHIN_TREND
  )
) {
  
  stop(
    "Initial model log probability is not finite."
  )
}


## ============================================================
## 40. CONFIGURE MCMC
## ============================================================

conf_REF10_WITHIN_TREND <- configureMCMC(
  model_REF10_WITHIN_TREND,
  
  monitors = c(
    "a10",
    "b",
    
    "log_a10",
    "log_b",
    
    "betaSEX",
    
    "betaDENS_a",
    "betaDENS_b",
    
    "betaTIME",
    "HR_TIME_DECADE",
    
    "sigma_island_a",
    "sigma_time",
    
    "island_a",
    "u_time",
    
    "mean.p"
  ),
  
  enableWAIC = FALSE
)


## ============================================================
## 41. JOINT a10-SIDE SAMPLER
##
## Includes temporal trend because it contributes directly
## to mortality level.
## ============================================================

a_pars_TREND <- c(
  "log_a10",
  "betaSEX[1]",
  "betaDENS_a",
  "betaTIME"
)


conf_REF10_WITHIN_TREND$removeSamplers(
  a_pars_TREND
)


conf_REF10_WITHIN_TREND$addSampler(
  target =
    a_pars_TREND,
  
  type =
    "AF_slice"
)


## ============================================================
## 42. JOINT b-SIDE SAMPLER
## ============================================================

b_pars_TREND <- c(
  "log_b",
  "betaSEX[2]",
  "betaDENS_b"
)


conf_REF10_WITHIN_TREND$removeSamplers(
  b_pars_TREND
)


conf_REF10_WITHIN_TREND$addSampler(
  target =
    b_pars_TREND,
  
  type =
    "AF_slice"
)


## ============================================================
## 43. VARIANCE PARAMETERS
## ============================================================

for(v in c(
  "sigma_island_a",
  "sigma_time"
)) {
  
  conf_REF10_WITHIN_TREND$removeSamplers(
    v
  )
  
  
  conf_REF10_WITHIN_TREND$addSampler(
    target =
      v,
    
    type =
      "slice"
  )
}


## ============================================================
## 44. LATENT BIRTH + LIFESPAN SAMPLERS
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
  
  
  conf_REF10_WITHIN_TREND$removeSamplers(
    life_nodes
  )
  
  
  conf_REF10_WITHIN_TREND$addSampler(
    target =
      life_nodes,
    
    type =
      "AF_slice"
  )
}


## ============================================================
## 45. BUILD + COMPILE
## ============================================================

mcmc_REF10_WITHIN_TREND <- buildMCMC(
  conf_REF10_WITHIN_TREND
)


cModel_REF10_WITHIN_TREND <- compileNimble(
  model_REF10_WITHIN_TREND,
  resetFunctions = TRUE
)


cMCMC_REF10_WITHIN_TREND <- compileNimble(
  mcmc_REF10_WITHIN_TREND,
  project = model_REF10_WITHIN_TREND
)


## ============================================================
## 46. RUN 100k STRUCTURAL CHECK
## ============================================================

cat("\n")
cat("============================================================\n")
cat("RUNNING WITHIN-DENSITY + TEMPORAL-TREND MODEL\n")
cat("============================================================\n")
cat("Iterations: 100000\n")
cat("Burn-in:    30000\n")
cat("Chains:     3\n")
cat("Thin:       10\n")
cat("============================================================\n")


run_REF10_WITHIN_TREND <- runMCMC(
  cMCMC_REF10_WITHIN_TREND,
  
  niter = 100000,
  nburnin = 30000,
  nchains = 3,
  thin = 10,
  
  inits =
    inits_REF10_WITHIN_TREND,
  
  setSeed = TRUE,
  
  samplesAsCodaMCMC = TRUE,
  
  summary = TRUE,
  
  WAIC = FALSE,
  
  progressBar = TRUE
)


## ============================================================
## 47. SAVE RAW FIT IMMEDIATELY
## ============================================================

saveRDS(
  run_REF10_WITHIN_TREND,
  "Model_PaceOfLife_REF10_WITHIN_DENSITY_TEMPORAL_TREND_100k.rds"
)


## ============================================================
## 48. PRIMARY SUMMARY
## ============================================================

sm_REF10_WITHIN_TREND <- MCMCsummary(
  run_REF10_WITHIN_TREND$samples,
  
  params = c(
    "a10",
    "b",
    
    "log_a10",
    "log_b",
    
    "betaSEX",
    
    "betaDENS_a",
    "betaDENS_b",
    
    "betaTIME",
    "HR_TIME_DECADE",
    
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
  sm_REF10_WITHIN_TREND
)


## ============================================================
## 49. POSTERIOR MATRIX
## ============================================================

s_REF10_WITHIN_TREND <- as.matrix(
  run_REF10_WITHIN_TREND$samples
)


## ============================================================
## 50. BIOLOGICAL EFFECT SUMMARY
##
## Sex ratios:
##   male / female
##
## Density ratios:
##   +1 SD within-island density deviation
##
## Time ratio:
##   one decade later birth cohort
## ============================================================

effect_REF10_WITHIN_TREND <- tibble(
  effect = c(
    "Sex -> mortality at age 10",
    "Sex -> senescence",
    
    "Within-island density -> mortality at age 10",
    "Within-island density -> senescence",
    
    "Birth cohort -> mortality at age 10 per decade"
  ),
  
  parameter = c(
    "betaSEX[1]",
    "betaSEX[2]",
    
    "betaDENS_a",
    "betaDENS_b",
    
    "betaTIME"
  )
) %>%
  rowwise() %>%
  mutate(
    median_beta =
      median(
        s_REF10_WITHIN_TREND[
          ,
          parameter
        ]
      ),
    
    low_beta =
      unname(
        quantile(
          s_REF10_WITHIN_TREND[
            ,
            parameter
          ],
          0.025
        )
      ),
    
    high_beta =
      unname(
        quantile(
          s_REF10_WITHIN_TREND[
            ,
            parameter
          ],
          0.975
        )
      ),
    
    median_ratio =
      median(
        exp(
          s_REF10_WITHIN_TREND[
            ,
            parameter
          ]
        )
      ),
    
    low_ratio =
      unname(
        quantile(
          exp(
            s_REF10_WITHIN_TREND[
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
            s_REF10_WITHIN_TREND[
              ,
              parameter
            ]
          ),
          0.975
        )
      ),
    
    Pr_positive =
      mean(
        s_REF10_WITHIN_TREND[
          ,
          parameter
        ] >
          0
      ),
    
    Pr_negative =
      mean(
        s_REF10_WITHIN_TREND[
          ,
          parameter
        ] <
          0
      )
  ) %>%
  ungroup()


cat("\n")
cat("============================================================\n")
cat("BIOLOGICAL EFFECT SUMMARY\n")
cat("============================================================\n")


print(
  effect_REF10_WITHIN_TREND,
  n = Inf
)


## ============================================================
## 51. EXTRA TIME INTERPRETATION
## ============================================================

betaTIME_draw <-
  s_REF10_WITHIN_TREND[
    ,
    "betaTIME"
  ]


time_effect_interpretation <- tibble(
  quantity = c(
    "Hazard ratio per decade later birth",
    "Percent change in age-10 hazard per decade"
  ),
  
  median = c(
    median(
      exp(
        betaTIME_draw
      )
    ),
    
    median(
      (
        exp(
          betaTIME_draw
        ) -
          1
      ) *
        100
    )
  ),
  
  low = c(
    unname(
      quantile(
        exp(
          betaTIME_draw
        ),
        0.025
      )
    ),
    
    unname(
      quantile(
        (
          exp(
            betaTIME_draw
          ) -
            1
        ) *
          100,
        0.025
      )
    )
  ),
  
  high = c(
    unname(
      quantile(
        exp(
          betaTIME_draw
        ),
        0.975
      )
    ),
    
    unname(
      quantile(
        (
          exp(
            betaTIME_draw
          ) -
            1
        ) *
          100,
        0.975
      )
    )
  )
)


cat("\n")
cat("============================================================\n")
cat("TEMPORAL TREND INTERPRETATION\n")
cat("============================================================\n")


print(
  time_effect_interpretation,
  n = Inf
)


## ============================================================
## 52. PARAMETER CORRELATIONS
## ============================================================

cor_a_b_TREND <- cor(
  s_REF10_WITHIN_TREND[
    ,
    "log_a10"
  ],
  
  s_REF10_WITHIN_TREND[
    ,
    "log_b"
  ]
)


a_side_cor_TREND <- cor(
  s_REF10_WITHIN_TREND[
    ,
    c(
      "log_a10",
      "betaSEX[1]",
      "betaDENS_a",
      "betaTIME"
    )
  ]
)


b_side_cor_TREND <- cor(
  s_REF10_WITHIN_TREND[
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
cat("POSTERIOR CORRELATIONS\n")
cat("============================================================\n")


cat(
  "log_a10 vs log_b:",
  round(
    cor_a_b_TREND,
    3
  ),
  "\n\n"
)


cat("a10-side correlations:\n")


print(
  round(
    a_side_cor_TREND,
    3
  )
)


cat("\n")


cat("b-side correlations:\n")


print(
  round(
    b_side_cor_TREND,
    3
  )
)


## ============================================================
## 53. CONVERGENCE
## ============================================================

conv_REF10_WITHIN_TREND <- sm_REF10_WITHIN_TREND %>%
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
cat("CONVERGENCE - WORST RHAT FIRST\n")
cat("============================================================\n")


print(
  conv_REF10_WITHIN_TREND,
  n = Inf
)


cat(
  "\nMaximum Rhat:",
  max(
    conv_REF10_WITHIN_TREND$Rhat,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Minimum ESS:",
  min(
    conv_REF10_WITHIN_TREND$n.eff,
    na.rm = TRUE
  ),
  "\n"
)


## ============================================================
## 54. RESIDUAL TEMPORAL EFFECTS
##
## Crucial question:
##
## Has sigma_time become smaller after including the trend?
##
## Are the residual periods now centred around zero rather than
## recreating the entire secular trajectory?
## ============================================================

sm_time_REF10_WITHIN_TREND <- MCMCsummary(
  run_REF10_WITHIN_TREND$samples,
  params = "u_time"
) %>%
  as.data.frame() %>%
  rownames_to_column(
    "parameter"
  ) %>%
  as_tibble() %>%
  mutate(
    model_block =
      readr::parse_number(
        parameter
      )
  ) %>%
  left_join(
    model_temporal_table,
    by = "model_block"
  )


cat("\n")
cat("============================================================\n")
cat("RESIDUAL TEMPORAL EFFECTS AFTER LINEAR TREND\n")
cat("============================================================\n")


print(
  sm_time_REF10_WITHIN_TREND,
  n = Inf
)


## ============================================================
## 55. ISLAND EFFECTS
## ============================================================

sm_island_REF10_WITHIN_TREND <- MCMCsummary(
  run_REF10_WITHIN_TREND$samples,
  params = "island_a"
)


cat("\n")
cat("============================================================\n")
cat("RESIDUAL ISLAND EFFECTS\n")
cat("============================================================\n")


print(
  sm_island_REF10_WITHIN_TREND
)


## ============================================================
## 56. OPTIONAL COMPARISON WITH PREVIOUS NO-FEEDING MODEL
##
## Sex is directly comparable.
##
## Density is NOT exactly the same estimand:
##
## previous = total between + within island density
## current  = within-island deviation only
##
## Therefore density comparison is descriptive, not a direct
## coefficient-stability test.
## ============================================================

old_fit_file <-
  "Model_PaceOfLife_REFAGE10_DENSITY_CORRECTED_NOFEED_100k_COMPLETE.rds"


comparison_OLD_NEW <-
  NULL


if(
  file.exists(
    old_fit_file
  )
) {
  
  old_obj <- readRDS(
    old_fit_file
  )
  
  
  s_old <- as.matrix(
    old_obj$fit$samples
  )
  
  
  compare_one <- function(
    samples,
    parameter,
    model
  ) {
    
    x <-
      samples[
        ,
        parameter
      ]
    
    
    tibble(
      model =
        model,
      
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
          x > 0
        ),
      
      Pr_negative =
        mean(
          x < 0
        )
    )
  }
  
  
  comparison_OLD_NEW <- bind_rows(
    compare_one(
      s_old,
      "betaSEX[1]",
      "Previous overall-density model"
    ),
    
    compare_one(
      s_REF10_WITHIN_TREND,
      "betaSEX[1]",
      "Within-density + trend model"
    ),
    
    compare_one(
      s_old,
      "betaSEX[2]",
      "Previous overall-density model"
    ),
    
    compare_one(
      s_REF10_WITHIN_TREND,
      "betaSEX[2]",
      "Within-density + trend model"
    ),
    
    compare_one(
      s_old,
      "betaDENS_a",
      "Previous overall-density model"
    ),
    
    compare_one(
      s_REF10_WITHIN_TREND,
      "betaDENS_a",
      "Within-density + trend model"
    ),
    
    compare_one(
      s_old,
      "betaDENS_b",
      "Previous overall-density model"
    ),
    
    compare_one(
      s_REF10_WITHIN_TREND,
      "betaDENS_b",
      "Within-density + trend model"
    )
  )
  
  
  cat("\n")
  cat("============================================================\n")
  cat("PREVIOUS VS WITHIN-DENSITY + TREND MODEL\n")
  cat("============================================================\n")
  
  
  print(
    comparison_OLD_NEW,
    n = Inf
  )
}


## ============================================================
## 57. COMMON RAW LOG-DENSITY SCALE
##
## This still does NOT make the estimands identical.
##
## But it helps distinguish a scaling change from a substantive
## change in the density coefficient.
## ============================================================

density_raw_scale_summary <- tibble(
  effect = c(
    "Within density -> a10",
    "Within density -> b"
  ),
  
  standardized_beta = c(
    median(
      s_REF10_WITHIN_TREND[
        ,
        "betaDENS_a"
      ]
    ),
    
    median(
      s_REF10_WITHIN_TREND[
        ,
        "betaDENS_b"
      ]
    )
  ),
  
  within_density_sd =
    density_within_sd
) %>%
  mutate(
    beta_per_unit_log_density =
      standardized_beta /
      within_density_sd
  )


cat("\n")
cat("============================================================\n")
cat("WITHIN-DENSITY EFFECTS ON RAW LOG-DENSITY SCALE\n")
cat("============================================================\n")


print(
  density_raw_scale_summary,
  n = Inf
)


## ============================================================
## 58. PLOT RESIDUAL TEMPORAL EFFECTS
## ============================================================

temporal_residual_plot <- ggplot(
  sm_time_REF10_WITHIN_TREND,
  aes(
    x =
      midpoint_year,
    
    y =
      mean
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_errorbar(
    aes(
      ymin =
        `2.5%`,
      
      ymax =
        `97.5%`
    ),
    width = 0
  ) +
  geom_point(
    size = 2
  ) +
  geom_line() +
  labs(
    x =
      "Birth-cohort period midpoint",
    
    y =
      "Residual temporal effect on log age-10 mortality",
    
    title =
      "Residual cohort-period effects after accounting for temporal trend"
  ) +
  theme_bw()


print(
  temporal_residual_plot
)


ggsave(
  "Temporal_Residuals_After_Trend.png",
  temporal_residual_plot,
  width = 9,
  height = 6,
  dpi = 300
)


## ============================================================
## 59. PLOT WITHIN-ISLAND DENSITY THROUGH TIME
## ============================================================

within_density_plot <- density_complete %>%
  filter(
    trusted_annual_density
  ) %>%
  ggplot(
    aes(
      x =
        calendar_year,
      
      y =
        density_within_z
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_line() +
  geom_point() +
  facet_wrap(
    ~ island
  ) +
  labs(
    x =
      "Calendar year",
    
    y =
      "Within-island cohort density deviation (SD)",
    
    title =
      "Established within-island density variation"
  ) +
  theme_bw()


print(
  within_density_plot
)


ggsave(
  "Within_Island_Density_Deviation_Model.png",
  within_density_plot,
  width = 12,
  height = 8,
  dpi = 300
)


## ============================================================
## 60. SAVE TABLES
## ============================================================

write_csv(
  density_complete,
  "Within_Island_Density_Model_Table.csv"
)


write_csv(
  within_density_check,
  "Within_Island_Density_Check.csv"
)


write_csv(
  density_time_audit,
  "Within_Island_Density_Time_Audit.csv"
)


write_csv(
  effect_REF10_WITHIN_TREND,
  "Model_REF10_WITHIN_DENSITY_TEMPORAL_TREND_effects.csv"
)


write_csv(
  time_effect_interpretation,
  "Model_REF10_WITHIN_DENSITY_TEMPORAL_TREND_time_effect.csv"
)


write_csv(
  sm_time_REF10_WITHIN_TREND,
  "Model_REF10_WITHIN_DENSITY_TEMPORAL_TREND_residual_time.csv"
)


write_csv(
  density_raw_scale_summary,
  "Model_REF10_WITHIN_DENSITY_raw_scale_effects.csv"
)


if(
  !is.null(
    comparison_OLD_NEW
  )
) {
  
  write_csv(
    comparison_OLD_NEW,
    "Comparison_Previous_vs_WithinDensityTrend.csv"
  )
}


## ============================================================
## 61. SAVE COMPLETE ANALYSIS
## ============================================================

saveRDS(
  list(
    fit =
      run_REF10_WITHIN_TREND,
    
    summary =
      sm_REF10_WITHIN_TREND,
    
    effects =
      effect_REF10_WITHIN_TREND,
    
    convergence =
      conv_REF10_WITHIN_TREND,
    
    temporal_residual =
      sm_time_REF10_WITHIN_TREND,
    
    island =
      sm_island_REF10_WITHIN_TREND,
    
    time_interpretation =
      time_effect_interpretation,
    
    density_table =
      density_complete,
    
    density_baselines =
      island_density_baseline,
    
    density_support =
      established_support,
    
    density_within_check =
      within_density_check,
    
    density_time_audit =
      density_time_audit,
    
    density_within_sd =
      density_within_sd,
    
    density_raw_scale_summary =
      density_raw_scale_summary,
    
    comparison_previous =
      comparison_OLD_NEW,
    
    temporal_mapping =
      model_temporal_table,
    
    startup_n_surveys =
      startup_n_surveys,
    
    age_ref =
      age_ref,
    
    time_ref_year =
      time_ref_year,
    
    time_scale_years =
      time_scale_years,
    
    birth_change_audit =
      birth_change_audit,
    
    weak_birth_n =
      sum(
        growth_info$birth_weak_or_adjusted
      ),
    
    correlation_a10_b =
      cor_a_b_TREND,
    
    a_side_correlation_matrix =
      a_side_cor_TREND,
    
    b_side_correlation_matrix =
      b_side_cor_TREND
  ),
  
  "Model_PaceOfLife_REF10_WITHIN_DENSITY_TEMPORAL_TREND_100k_COMPLETE.rds"
)


cat("\n")
cat("============================================================\n")
cat("WITHIN-DENSITY + TEMPORAL-TREND MODEL COMPLETE\n")
cat("============================================================\n")