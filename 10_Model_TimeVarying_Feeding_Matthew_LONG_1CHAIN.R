## ============================================================
## 10_Model_TimeVarying_Feeding_Matthew_LONG_1CHAIN.R
##
## NINE-SITE TIME-VARYING SURVIVAL-RECAPTURE MODEL
##
## Survival:
##   - continuous-age Gompertz hazard parameterised at age 10
##   - left truncation at each individual's first live capture
##   - latent / partially observed birth time and sex
##   - within-site birth-cohort density
##   - linear birth-cohort effect on age-10 mortality
##   - residual centred site effect on age-10 mortality
##   - calendar-varying high/regular provisioning
##   - FFRC post-Matthew calendar indicator
##   - centred site-specific calendar-time nuisance slopes
##
## Detection:
##   - reconstructed field visits separated by >14-day gaps
##   - first-capture event excluded; no pre-entry zeros
##   - sex and age (saturating at age 15)
##   - centred site and field-visit effects
##
## NOTE:
##   The final publishable individual-detection model is generated
##   from this source by Build_IndivDetect_Model_Files.R, which adds
##   persistent individual heterogeneity in detection without changing
##   the survival model.
## ============================================================
library(tidyverse)
library(lubridate)
library(nimble)
library(MCMCvis)
library(coda)


## ============================================================
## 1. LOAD FROZEN DATA PREP
## ============================================================

source(
  "DataPrep_Canonical_LatentSex_V5.R"
)


covariate_file <-
  "Data/Iguana_ModelCovariates_Through2019.rds"


if(!file.exists(covariate_file)) {
  stop(
    "Cannot find ",
    covariate_file,
    ". Run 05_Build_Density_and_Exposure_Histories.R first."
  )
}


cov_obj <- readRDS(
  covariate_file
)


stopifnot(
  identical(
    cov_obj$ordered_islands,
    ordered_islands
  ),
  cov_obj$Origin ==
    Origin,
  cov_obj$tMax ==
    tMax
)


density_matrix <-
  cov_obj$density_matrix


birth_density_support <-
  cov_obj$birth_density_support


## ============================================================
## 2. LOAD CALENDAR-VARYING LEFT-TRUNCATED GOMPERTZ
## ============================================================

gompertz_candidates <- c(
  "Distributions/Dist_GompertzCalendarLT_V2.R",
  "Dist_GompertzCalendarLT_V2.R"
)


gompertz_file <- gompertz_candidates[
  file.exists(
    gompertz_candidates
  )
][1]


if(is.na(gompertz_file)) {
  stop(
    "Cannot find Dist_GompertzCalendarLT.R. Put it either in the ",
    "project root or Distributions/."
  )
}


source(
  gompertz_file
)


## ============================================================
## 2B. LOAD FINAL CALENDAR EXPOSURES
##
## Primary feeding contrast:
##   high / regular provisioning = 1
##   none or low / intermittent = 0
##
## Hurricane Matthew primary contrast:
##   FFRC 2017 onward = 1
##
## Large-male removal is deliberately omitted.
## ============================================================

exposure_file <-
  "Data/Iguana_SurveyOccasion_and_CalendarExposure_Audit.rds"


if(!file.exists(exposure_file)) {
  stop(
    "Cannot find ",
    exposure_file,
    ". Run 07_Audit_SurveyOccasions_and_CalendarExposures_V2.R first."
  )
}


exposure_obj <- readRDS(
  exposure_file
)


feeding_high_matrix <-
  exposure_obj$feeding_high_matrix


matthew_post_matrix <-
  exposure_obj$matthew_post_matrix


stopifnot(
  all(
    dim(
      feeding_high_matrix
    ) ==
      c(
        n_islands,
        tMax
      )
  ),
  all(
    dim(
      matthew_post_matrix
    ) ==
      c(
        n_islands,
        tMax
      )
  )
)


## Duplicate the small site-year matrices to individual rows so the
## custom NIMBLE distribution receives an ordinary data vector for
## each individual.
feed_individual_matrix <-
  feeding_high_matrix[
    island_idx,
    ,
    drop =
      FALSE
  ]


matthew_individual_matrix <-
  matthew_post_matrix[
    island_idx,
    ,
    drop =
      FALSE
  ]


period_ref_year <- 2010


## ============================================================
## 3. DENSITY ELIGIBILITY
##
## Conservative primary rule:
##   density contributes only when the complete possible birth
##   interval is within supported density years.
##
## This prevents missing historical density from feeding back into
## latent birth-year inference.
## ============================================================

density_eligible <- birth_density_support %>%
  arrange(
    ID
  ) %>%
  pull(
    prop_density_supported
  ) ==
  1


density_eligible <- as.integer(
  density_eligible
)


stopifnot(
  length(
    density_eligible
  ) ==
    nind
)


cat("\n")
cat("============================================================\n")
cat("DENSITY ELIGIBILITY FOR PRIMARY BASELINE\n")
cat("============================================================\n")
cat(
  "Fully supported individuals contributing to density effect: ",
  sum(
    density_eligible
  ),
  " / ",
  nind,
  "\n",
  sep =
    ""
)
cat(
  "Other individuals retained with density contribution = 0: ",
  sum(
    density_eligible ==
      0
  ),
  "\n",
  sep =
    ""
)
cat("============================================================\n")


## ============================================================
## 4. RECONSTRUCT FIELD-VISIT SURVEY OCCASIONS
##
## The audit showed that a site-year is often too coarse:
## 30/138 site-years have capture activity separated by large gaps.
##
## The candidate occasion count stabilised at a gap threshold of
## about 14-21 days. We therefore use >14 days as the primary rule:
##
##   new field visit = next capture day is >14 days after the
##                     previous capture day at that site/year.
##
## This gives a conservative visit definition while preserving
## clearly separated spring/summer/autumn trips as distinct events.
## ============================================================

visit_gap_days <- 14


capture_days <- AllData %>%
  filter(
    alive_encounter,
    !is.na(
      CapDate
    )
  ) %>%
  distinct(
    Island,
    CalYear,
    CapDate
  ) %>%
  arrange(
    Island,
    CalYear,
    CapDate
  ) %>%
  group_by(
    Island,
    CalYear
  ) %>%
  mutate(
    gap_days =
      as.integer(
        CapDate -
        lag(
          CapDate
        )
      ),

    new_visit =
      if_else(
        row_number() ==
          1L,
        1L,
        as.integer(
          gap_days >
            visit_gap_days
        )
      ),

    visit_within_year =
      cumsum(
        new_visit
      )
  ) %>%
  ungroup()


live_event_rows <- AllData %>%
  filter(
    alive_encounter,
    !is.na(
      CapDate
    )
  ) %>%
  inner_join(
    capture_days %>%
      select(
        Island,
        CalYear,
        CapDate,
        visit_within_year
      ),
    by =
      c(
        "Island",
        "CalYear",
        "CapDate"
      )
  )


survey_events_all <- live_event_rows %>%
  group_by(
    Island,
    CalYear,
    visit_within_year
  ) %>%
  summarise(
    event_start_date =
      min(
        CapDate
      ),

    event_end_date =
      max(
        CapDate
      ),

    event_start_cal =
      decimal_date(
        event_start_date
      ),

    event_end_cal =
      decimal_date(
        event_end_date
      ),

    event_time_cal =
      mean(
        decimal_date(
          CapDate
        )
      ),

    event_span_days =
      as.integer(
        event_end_date -
        event_start_date
      ),

    n_capture_days =
      n_distinct(
        CapDate
      ),

    n_live_captures =
      n_distinct(
        ID
      ),

    .groups =
      "drop"
  ) %>%
  arrange(
    CalYear,
    event_start_date,
    Island
  ) %>%
  mutate(
    event_id_all =
      row_number(),

    island_idx =
      match(
        Island,
        ordered_islands
      )
  )


## Attach reconstructed event IDs back to live capture rows.
live_event_rows <- live_event_rows %>%
  left_join(
    survey_events_all %>%
      select(
        Island,
        CalYear,
        visit_within_year,
        event_id_all,
        event_start_cal,
        event_end_cal,
        event_time_cal
      ),
    by =
      c(
        "Island",
        "CalYear",
        "visit_within_year"
      )
  )


## ============================================================
## 5. FIRST LIVE-CAPTURE EVENT BY INDIVIDUAL
##
## First capture itself is conditioned upon and therefore excluded
## from the detection likelihood. Earlier visits have no zeros.
## Later visits in the same calendar year ARE allowed.
## ============================================================

first_live <- live_event_rows %>%
  arrange(
    ID,
    CapDate
  ) %>%
  group_by(
    ID
  ) %>%
  slice(
    1L
  ) %>%
  ungroup() %>%
  transmute(
    ID =
      ID,

    first_event_id =
      event_id_all,

    first_live_time_cal =
      decimal_date(
        CapDate
      )
  )


stopifnot(
  nrow(
    first_live
  ) ==
    nind
)


## ============================================================
## 6. OBSERVED LIVE DETECTIONS BY INDIVIDUAL x VISIT
## ============================================================

live_detection_by_event <- live_event_rows %>%
  group_by(
    ID,
    Island,
    event_id_all
  ) %>%
  summarise(
    detected =
      1L,

    capture_time_cal =
      min(
        decimal_date(
          CapDate
        )
      ),

    .groups =
      "drop"
  )


## ============================================================
## 7. POST-ENTRY ENCOUNTER PAIRS
##
## Rules:
##   - no pre-first-capture zeros
##   - first-capture visit excluded
##   - all subsequent reconstructed visits are opportunities
##   - known permanent removals contribute no opportunities after
##     removal
##   - known observed deaths contribute no opportunities after the
##     observed death time
##
## For detected pairs use actual capture time.
## For non-detections use the mean capture time of that field visit.
## ============================================================

removal_times <- AllData %>%
  filter(
    permanent_removal
  ) %>%
  group_by(
    ID
  ) %>%
  summarise(
    removal_time_cal =
      min(
        decimal_date(
          CapDate
        )
      ),

    .groups =
      "drop"
  )


individual_lookup <- growth_info %>%
  transmute(
    ID =
      ID,

    Island =
      island,

    observed_death_time_cal =
      tD_cal
  ) %>%
  left_join(
    first_live,
    by =
      "ID"
  ) %>%
  left_join(
    removal_times,
    by =
      "ID"
  )


encounter_pairs <- individual_lookup %>%
  inner_join(
    survey_events_all,
    by =
      "Island",
    relationship =
      "many-to-many"
  ) %>%
  filter(
    event_start_cal >
      first_live_time_cal
  ) %>%
  filter(
    event_id_all !=
      first_event_id
  ) %>%
  filter(
    is.na(
      removal_time_cal
    ) |
      event_start_cal <=
        removal_time_cal
  ) %>%
  filter(
    is.na(
      observed_death_time_cal
    ) |
      event_start_cal <=
        observed_death_time_cal
  ) %>%
  left_join(
    live_detection_by_event,
    by =
      c(
        "ID",
        "Island",
        "event_id_all"
      )
  ) %>%
  mutate(
    detected =
      replace_na(
        detected,
        0L
      ),

    event_time_cal_pair =
      if_else(
        detected ==
          1L,
        capture_time_cal,
        event_time_cal
      )
  )


used_event_ids <- sort(
  unique(
    encounter_pairs$event_id_all
  )
)


used_events <- survey_events_all %>%
  filter(
    event_id_all %in%
      used_event_ids
  ) %>%
  arrange(
    event_id_all
  ) %>%
  mutate(
    event_id =
      row_number()
  )


event_id_map <- used_events %>%
  select(
    event_id_all,
    event_id
  )


encounter_pairs <- encounter_pairs %>%
  inner_join(
    event_id_map,
    by =
      "event_id_all"
  ) %>%
  arrange(
    ID,
    event_time_cal_pair
  )


nObs <- nrow(
  encounter_pairs
)


n_events <- nrow(
  used_events
)


cat("\n")
cat("============================================================\n")
cat("RECONSTRUCTED FIELD-VISIT DETECTION HISTORY\n")
cat("============================================================\n")
cat(
  "Visit-gap threshold: >",
  visit_gap_days,
  " days\n",
  sep =
    ""
)
cat(
  "Recognised field-visit events: ",
  nrow(
    survey_events_all
  ),
  "\n",
  sep =
    ""
)
cat(
  "Field-visit events used after first capture: ",
  n_events,
  "\n",
  sep =
    ""
)
cat(
  "Post-entry individual-event pairs: ",
  nObs,
  "\n",
  sep =
    ""
)
cat(
  "Post-entry live detections: ",
  sum(
    encounter_pairs$detected
  ),
  "\n",
  sep =
    ""
)
cat(
  "Raw post-entry encounter rate: ",
  round(
    mean(
      encounter_pairs$detected
    ),
    4
  ),
  "\n",
  sep =
    ""
)
cat(
  "Permanent-removal individuals: ",
  nrow(
    removal_times
  ),
  "\n",
  sep =
    ""
)
cat("============================================================\n")


## ============================================================
## 8. MODEL INPUT VECTORS
## ============================================================

obs_ind <- as.integer(
  encounter_pairs$ID
)


obs_event <- as.integer(
  encounter_pairs$event_id
)


obs_time <-
  encounter_pairs$event_time_cal_pair -
  Origin +
  1


obs_detected <- as.integer(
  encounter_pairs$detected
)


age_ref <- 10


detection_adult_age <- 15


time_ref_year <- 1990


time_scale_years <- 10


## First live capture on adjusted model-time scale.
entry_time <-
  growth_info$tF_cal -
  Origin +
  1


## ============================================================
## 9. MODEL
## ============================================================

code_BASE <- nimbleCode({

  ## ----------------------------------------------------------
  ## POPULATION PARAMETERS
  ## ----------------------------------------------------------

  log_a10 ~ dnorm(
    log(
      0.03
    ),
    sd =
      0.75
  )

  log_b ~ dnorm(
    log(
      0.08
    ),
    sd =
      0.50
  )

  a10_ref <- exp(
    log_a10
  )

  b_ref <- exp(
    log_b
  )


  ## ----------------------------------------------------------
  ## SEX / DENSITY / COHORT EFFECTS
  ## ----------------------------------------------------------

  betaSEX[1] ~ dnorm(
    0,
    sd =
      0.75
  )

  betaSEX[2] ~ dnorm(
    0,
    sd =
      0.50
  )

  betaDENS_a ~ dnorm(
    0,
    sd =
      0.75
  )

  betaDENS_b ~ dnorm(
    0,
    sd =
      0.50
  )

  betaCOHORT ~ dnorm(
    0,
    sd =
      0.50
  )


  ## ----------------------------------------------------------
  ## CALENDAR-VARYING SURVIVAL EXPOSURES
  ##
  ## betaFEED:
  ##   proportional-hazard effect while a site-year is classified
  ##   as high / regular provisioning.
  ##
  ## betaMATTHEW:
  ##   proportional-hazard effect for FFRC from 2017 onward.
  ## ----------------------------------------------------------

  betaFEED ~ dnorm(
    0,
    sd =
      0.75
  )

  betaMATTHEW ~ dnorm(
    0,
    sd =
      0.75
  )


  ## ----------------------------------------------------------
  ## CENTERED SITE-SPECIFIC CALENDAR SLOPES
  ##
  ## No unrestricted global period trend is fitted because age,
  ## birth cohort and period are structurally related.
  ##
  ## Instead, sites can deviate around a zero-average calendar
  ## slope. This helps prevent feeding onset at Bitter / LeafSouth
  ## from absorbing every unrelated temporal change at those sites.
  ## ----------------------------------------------------------

  sigma_site_period ~ dunif(
    0,
    0.75
  )

  for(j in 1:n_islands) {

    site_period_raw[j] ~ dnorm(
      0,
      sd =
        1
    )
  }

  mean_site_period_raw <-
    sum(
      site_period_raw[
        1:n_islands
      ]
    ) /
    n_islands

  for(j in 1:n_islands) {

    site_period[j] <-
      sigma_site_period *
      (
        site_period_raw[j] -
        mean_site_period_raw
      )
  }


  ## ----------------------------------------------------------
  ## RESIDUAL SITE EFFECT ON AGE-10 MORTALITY
  ## ----------------------------------------------------------

  sigma_island_a ~ dunif(
    0,
    2
  )

  for(j in 1:n_islands) {

    island_a_raw[j] ~ dnorm(
      0,
      sd =
        1
    )
  }

  mean_island_a_raw <-
    sum(
      island_a_raw[
        1:n_islands
      ]
    ) /
    n_islands

  for(j in 1:n_islands) {

    island_a[j] <-
      sigma_island_a *
      (
        island_a_raw[j] -
        mean_island_a_raw
      )
  }


  ## ----------------------------------------------------------
  ## DETECTION PARAMETERS
  ## ----------------------------------------------------------

  alpha_p ~ dnorm(
    qlogis(
      0.20
    ),
    sd =
      1.5
  )

  betaP_SEX ~ dnorm(
    0,
    sd =
      1
  )

  betaP_AGE ~ dnorm(
    0,
    sd =
      1
  )


  sigma_island_p ~ dunif(
    0,
    3
  )

  for(j in 1:n_islands) {

    island_p_raw[j] ~ dnorm(
      0,
      sd =
        1
    )
  }

  mean_island_p_raw <-
    sum(
      island_p_raw[
        1:n_islands
      ]
    ) /
    n_islands

  for(j in 1:n_islands) {

    island_p[j] <-
      sigma_island_p *
      (
        island_p_raw[j] -
        mean_island_p_raw
      )
  }


  sigma_survey_p ~ dunif(
    0,
    3
  )

  for(e in 1:n_events) {

    survey_p_raw[e] ~ dnorm(
      0,
      sd =
        1
    )
  }

  mean_survey_p_raw <-
    sum(
      survey_p_raw[
        1:n_events
      ]
    ) /
    n_events

  for(e in 1:n_events) {

    survey_p[e] <-
      sigma_survey_p *
      (
        survey_p_raw[e] -
        mean_survey_p_raw
      )
  }


  ## ----------------------------------------------------------
  ## INDIVIDUAL LIFE HISTORIES
  ## ----------------------------------------------------------

  for(i in 1:nind) {

    ## Partially observed sex:
    ## observed 0/1 values are data;
    ## NA values are sampled.
    sex[i] ~ dbern(
      0.5
    )


    ## Latent birth time.
    tB[i] ~ dunif(
      cintB[
        i,
        1
      ],
      cintB[
        i,
        2
      ]
    )


    ## Correct delayed entry:
    ## condition on survival to THIS INDIVIDUAL'S first live capture.
    entry_age[i] <-
      max(
        0,
        entry_time[i] -
        tB[i]
      )


    ## Birth cohort index.
    cohort_idx[i] <-
      max(
        1,
        min(
          tMax,
          floor(
            tB[i]
          )
        )
      )


    ## Calendar birth year and centred decade.
    birth_year_cal[i] <-
      tB[i] +
      Origin -
      1

    birth_decade_c[i] <-
      (
        birth_year_cal[i] -
        time_ref_year
      ) /
      time_scale_years


    ## Only fully supported birth windows inform density.
    dens_cohort[i] <-
      density_eligible[i] *
      density_matrix[
        island_idx[i],
        cohort_idx[i]
      ]


    ## Age-10 mortality level.
    log_a10_i[i] <-
      log_a10 +
      betaSEX[1] *
      sex[i] +
      betaDENS_a *
      dens_cohort[i] +
      betaCOHORT *
      birth_decade_c[i] +
      island_a[
        island_idx[i]
      ]


    ## Senescence rate.
    log_b_i[i] <-
      log_b +
      betaSEX[2] *
      sex[i] +
      betaDENS_b *
      dens_cohort[i]

    b_i[i] <-
      exp(
        log_b_i[i]
      )


    ## Age-10 hazard on the natural scale.
    a10_i[i] <-
      exp(
        log_a10_i[i]
      )


    ## Age at death, conditional on surviving to first capture.
    ##
    ## The underlying Gompertz ageing curve remains age-based,
    ## while feeding, Matthew and the centered site-period slope
    ## modify instantaneous mortality in the calendar years when
    ## those exposures apply.
    tstar[i] ~ dGompertzCalendarLT(
      a10_i[i],
      b_i[i],
      tB[i],
      lowerBound =
        entry_age[i],
      feed =
        feed_individual_matrix[
          i,
          1:tMax
        ],
      matthew =
        matthew_individual_matrix[
          i,
          1:tMax
        ],
      betaFeed =
        betaFEED,
      betaMatthew =
        betaMATTHEW,
      sitePeriodSlope =
        site_period[
          island_idx[i]
        ],
      Origin =
        Origin,
      periodRefYear =
        period_ref_year,
      tMax =
        tMax
    )


    ## Calendar/model time of death.
    tD[i] <-
      tB[i] +
      tstar[i]


    ## Observed interval/right censoring.
    censoredD[i] ~ dinterval(
      tD[i],
      cintD[
        i,
        1:2
      ]
    )
  }


  ## ----------------------------------------------------------
  ## POST-ENTRY DETECTION LIKELIHOOD
  ## ----------------------------------------------------------

  for(k in 1:nObs) {

    age_obs[k] <-
      max(
        0,
        obs_time[k] -
        tB[
          obs_ind[k]
        ]
      )

    age_score[k] <-
      min(
        1,
        age_obs[k] /
        detection_adult_age
      )


    logit(
      p_det[k]
    ) <-
      alpha_p +
      betaP_SEX *
      sex[
        obs_ind[k]
      ] +
      betaP_AGE *
      age_score[k] +
      island_p[
        island_idx[
          obs_ind[k]
        ]
      ] +
      survey_p[
        obs_event[k]
      ]


    ## If death precedes a later survey, detection probability is zero.
    alive_at_event[k] <-
      step(
        tD[
          obs_ind[k]
        ] -
        obs_time[k]
      )


    p_obs[k] <-
      alive_at_event[k] *
      p_det[k]


    detected[k] ~ dbern(
      p_obs[k]
    )
  }
})


## ============================================================
## 10. CONSTANTS / DATA
## ============================================================

consts_BASE <- list(
  nind =
    nind,

  n_islands =
    n_islands,

  nObs =
    nObs,

  n_events =
    n_events,

  tMax =
    tMax,

  Origin =
    Origin,

  detection_adult_age =
    detection_adult_age,

  time_ref_year =
    time_ref_year,

  time_scale_years =
    time_scale_years,

  period_ref_year =
    period_ref_year,

  cintB =
    cintB_adj,

  cintD =
    cintD_adj,

  island_idx =
    island_idx,

  density_eligible =
    density_eligible,

  entry_time =
    entry_time,

  obs_ind =
    obs_ind,

  obs_event =
    obs_event,

  obs_time =
    obs_time
)


data_BASE <- list(
  sex =
    sex_obs,

  censoredD =
    censoredD_vec,

  detected =
    obs_detected,

  ## Must be supplied as DATA rather than a constant because
  ## cohort_idx[i] is stochastic/dynamic in:
  ##   density_matrix[island_idx[i], cohort_idx[i]]
  density_matrix =
    density_matrix,

  feed_individual_matrix =
    feed_individual_matrix,

  matthew_individual_matrix =
    matthew_individual_matrix
)


## ============================================================
## 11. INITIAL VALUES
## ============================================================

make_inits <- function(
  chain_id
) {

  set.seed(
    6000 +
      chain_id
  )


  tB_start <- runif(
    nind,
    min =
      cintB_adj[
        ,
        1
      ] +
      0.0001,
    max =
      cintB_adj[
        ,
        2
      ] -
      0.0001
  )


  tD_start <- numeric(
    nind
  )


  for(i in seq_len(nind)) {

    if(
      censoredD_vec[i] ==
        1L
    ) {

      ## Use the exact midpoint of the observed death interval.
      ## This is robust even for extremely short intervals such as
      ## animals that died ~1 hour after a live capture.
      tD_start[i] <-
        mean(
          cintD_adj[
            i,
          ]
        )

    } else {

      tD_start[i] <-
        cintD_adj[
          i,
          2
        ] +
        runif(
          1,
          0.5,
          3
        )
    }


    if(
      !is.finite(
        tD_start[i]
      )
    ) {
      stop(
        "Non-finite initial death time for individual ",
        i,
        " (",
        growth_info$ID_raw[i],
        ")."
      )
    }

    if(
      !is.finite(
        entry_time[i]
      )
    ) {
      stop(
        "Non-finite entry time for individual ",
        i,
        " (",
        growth_info$ID_raw[i],
        ")."
      )
    }

    if(
      tD_start[i] <=
        entry_time[i]
    ) {

      ## For a known death interval we must remain INSIDE that
      ## interval. A midpoint should already satisfy this; if it
      ## does not, stop rather than invent an impossible value.
      if(
        censoredD_vec[i] ==
          1L
      ) {
        stop(
          "Interval-censored death is not after first live capture for ",
          growth_info$ID_raw[i],
          ". Inspect its capture/death timing."
        )
      }

      tD_start[i] <-
        entry_time[i] +
        runif(
          1,
          0.5,
          2
        )
    }
  }


  sex_start <- rep(
    NA_real_,
    nind
  )


  sex_start[
    sex_is_latent
  ] <- rbinom(
    sum(
      sex_is_latent
    ),
    1,
    0.5
  )


  list(
    sex =
      sex_start,

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
        log(
          0.03
        ),
        0.08
      ),

    log_b =
      rnorm(
        1,
        log(
          0.08
        ),
        0.08
      ),

    betaSEX =
      rnorm(
        2,
        0,
        0.05
      ),

    betaDENS_a =
      rnorm(
        1,
        0,
        0.05
      ),

    betaDENS_b =
      rnorm(
        1,
        0,
        0.05
      ),

    betaCOHORT =
      rnorm(
        1,
        0,
        0.05
      ),

    betaFEED =
      rnorm(
        1,
        0,
        0.05
      ),

    betaMATTHEW =
      rnorm(
        1,
        0,
        0.05
      ),

    sigma_site_period =
      runif(
        1,
        0.02,
        0.15
      ),

    site_period_raw =
      rnorm(
        n_islands,
        0,
        0.10
      ),

    sigma_island_a =
      runif(
        1,
        0.2,
        0.8
      ),

    island_a_raw =
      rnorm(
        n_islands,
        0,
        0.3
      ),

    alpha_p =
      rnorm(
        1,
        qlogis(
          0.20
        ),
        0.10
      ),

    betaP_SEX =
      rnorm(
        1,
        0,
        0.05
      ),

    betaP_AGE =
      rnorm(
        1,
        0,
        0.05
      ),

    sigma_island_p =
      runif(
        1,
        0.3,
        1.2
      ),

    island_p_raw =
      rnorm(
        n_islands,
        0,
        0.3
      ),

    sigma_survey_p =
      runif(
        1,
        0.3,
        1.2
      ),

    survey_p_raw =
      rnorm(
        n_events,
        0,
        0.3
      )
  )
}


## ============================================================
## ONE CHAIN PER TERMINAL
##
## Set in the shell, e.g.
##   CHAIN_ID=1 Rscript 10_Model_TimeVarying_Feeding_Matthew_LONG_1CHAIN.R
##
## Each terminal/process gets independent initial values and an
## independent MCMC random seed.
## ============================================================

CHAIN_ID <- as.integer(
  Sys.getenv(
    "CHAIN_ID",
    "1"
  )
)


if(
  is.na(
    CHAIN_ID
  ) |
  CHAIN_ID <
    1
) {
  stop(
    "CHAIN_ID must be a positive integer."
  )
}


nchains <- 1L


inits_BASE <- list(
  make_inits(
    CHAIN_ID
  )
)


## ============================================================
## 12. BUILD MODEL + INITIAL LOG PROBABILITY
## ============================================================

cat("\nDefining baseline model...\n")


model_BASE <- nimbleModel(
  code =
    code_BASE,

  constants =
    consts_BASE,

  data =
    data_BASE,

  inits =
    inits_BASE[[1]],

  calculate =
    FALSE
)


lp_BASE <- model_BASE$calculate()


cat(
  "Initial log probability: ",
  lp_BASE,
  "\n",
  sep =
    ""
)


if(!is.finite(lp_BASE)) {

  cat("\nINITIALIZATION INFO\n")

  print(
    model_BASE$initializeInfo()
  )

  stop(
    "Initial log probability is not finite."
  )
}


## ============================================================
## 13. MCMC CONFIGURATION
## ============================================================

latent_sex_nodes <- paste0(
  "sex[",
  which(
    sex_is_latent
  ),
  "]"
)


monitors_BASE <- c(
  "a10_ref",
  "b_ref",
  "log_a10",
  "log_b",
  "betaSEX",
  "betaDENS_a",
  "betaDENS_b",
  "betaCOHORT",
  "betaFEED",
  "betaMATTHEW",
  "sigma_site_period",
  "site_period",
  "sigma_island_a",
  "island_a",
  "alpha_p",
  "betaP_SEX",
  "betaP_AGE",
  "sigma_island_p",
  "island_p",
  "sigma_survey_p",
  latent_sex_nodes
)


conf_BASE <- configureMCMC(
  model_BASE,
  monitors =
    monitors_BASE,
  enableWAIC =
    FALSE
)


## ------------------------------------------------------------
## Joint a-side sampler
## ------------------------------------------------------------

a_side <- c(
  "log_a10",
  "betaSEX[1]",
  "betaDENS_a",
  "betaCOHORT",
  "betaFEED",
  "betaMATTHEW"
)


conf_BASE$removeSamplers(
  a_side
)


conf_BASE$addSampler(
  target =
    a_side,
  type =
    "AF_slice"
)


## ------------------------------------------------------------
## Joint b-side sampler
## ------------------------------------------------------------

b_side <- c(
  "log_b",
  "betaSEX[2]",
  "betaDENS_b"
)


conf_BASE$removeSamplers(
  b_side
)


conf_BASE$addSampler(
  target =
    b_side,
  type =
    "AF_slice"
)


## ------------------------------------------------------------
## Scale parameters
## ------------------------------------------------------------

for(v in c(
  "sigma_site_period",
  "sigma_island_a",
  "sigma_island_p",
  "sigma_survey_p"
)) {

  conf_BASE$removeSamplers(
    v
  )

  conf_BASE$addSampler(
    target =
      v,
    type =
      "slice"
  )
}


## ------------------------------------------------------------
## Joint latent birth + age-at-death sampler
## ------------------------------------------------------------

cat(
  "Replacing tB/tstar samplers with joint AF_slice samplers...\n"
)


for(i in seq_len(nind)) {

  nodes_i <- c(
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

  conf_BASE$removeSamplers(
    nodes_i
  )

  conf_BASE$addSampler(
    target =
      nodes_i,
    type =
      "AF_slice"
  )
}


## ============================================================
## 14. COMPILE
## ============================================================

mcmc_BASE <- buildMCMC(
  conf_BASE
)


cat(
  "Compiling model...\n"
)


cModel_BASE <- compileNimble(
  model_BASE,
  resetFunctions =
    TRUE
)


cat(
  "Compiling MCMC...\n"
)


cMCMC_BASE <- compileNimble(
  mcmc_BASE,
  project =
    model_BASE
)


## ============================================================
## 15. LONG-RUN SETTINGS: ONE CHAIN PER TERMINAL
## ============================================================

niter <- as.integer(
  Sys.getenv(
    "NITER",
    "150000"
  )
)


nburnin <- as.integer(
  Sys.getenv(
    "NBURN",
    "50000"
  )
)


thin <- as.integer(
  Sys.getenv(
    "THIN",
    "10"
  )
)


if(
  is.na(
    niter
  ) |
  is.na(
    nburnin
  ) |
  is.na(
    thin
  ) |
  niter <=
    nburnin |
  nburnin <
    0 |
  thin <
    1
) {
  stop(
    "Invalid NITER/NBURN/THIN settings."
  )
}


chain_seed <-
  90000L +
  CHAIN_ID


cat("\n")
cat("============================================================\n")
cat("TIME-VARYING FEEDING + MATTHEW: LONG SINGLE-CHAIN RUN\n")
cat("============================================================\n")
cat(
  "Chain ID: ",
  CHAIN_ID,
  "\n",
  sep =
    ""
)
cat(
  "MCMC seed: ",
  chain_seed,
  "\n",
  sep =
    ""
)
cat(
  "Iterations: ",
  niter,
  "\n",
  sep =
    ""
)
cat(
  "Burn-in: ",
  nburnin,
  "\n",
  sep =
    ""
)
cat(
  "Thin: ",
  thin,
  "\n",
  sep =
    ""
)
cat(
  "Retained draws: ",
  floor(
    (
      niter -
      nburnin
    ) /
      thin
  ),
  "\n",
  sep =
    ""
)
cat(
  "Individuals: ",
  nind,
  "\n",
  sep =
    ""
)
cat(
  "Latent sex individuals: ",
  sum(
    sex_is_latent
  ),
  "\n",
  sep =
    ""
)
cat(
  "Density-eligible individuals: ",
  sum(
    density_eligible
  ),
  "\n",
  sep =
    ""
)
cat(
  "Survey events: ",
  n_events,
  "\n",
  sep =
    ""
)
cat(
  "Encounter pairs: ",
  nObs,
  "\n",
  sep =
    ""
)
cat("============================================================\n")


## ============================================================
## 16. RUN THIS CHAIN
## ============================================================

run_BASE <- runMCMC(
  cMCMC_BASE,
  niter =
    niter,
  nburnin =
    nburnin,
  nchains =
    1,
  thin =
    thin,
  inits =
    inits_BASE,
  setSeed =
    chain_seed,
  samplesAsCodaMCMC =
    TRUE,
  summary =
    TRUE,
  WAIC =
    FALSE,
  progressBar =
    TRUE
)


## ============================================================
## 17. SAVE CHAIN-SPECIFIC RESULT
##
## Do NOT calculate Rhat here: convergence is a between-chain
## diagnostic and will be calculated only after all chains are
## combined.
## ============================================================

out_dir <-
  "ModelRuns/TimeVarying_Feeding_Matthew_LONG"


dir.create(
  out_dir,
  recursive =
    TRUE,
  showWarnings =
    FALSE
)


latent_lookup <- tibble(
  ID =
    which(
      sex_is_latent
    ),

  ID_raw =
    growth_info$ID_raw[
      sex_is_latent
    ],

  island =
    growth_info$island[
      sex_is_latent
    ],

  sex_state =
    growth_info$sex_state[
      sex_is_latent
    ]
)


chain_file <- file.path(
  out_dir,
  paste0(
    "Model_TimeVarying_Feeding_Matthew_LONG_CHAIN",
    CHAIN_ID,
    ".rds"
  )
)


saveRDS(
  list(
    chain_id =
      CHAIN_ID,

    chain_seed =
      chain_seed,

    niter =
      niter,

    nburnin =
      nburnin,

    thin =
      thin,

    fit =
      run_BASE,

    latent_lookup =
      latent_lookup,

    density_eligible =
      density_eligible,

    birth_density_support =
      birth_density_support,

    feeding_high_matrix =
      feeding_high_matrix,

    matthew_post_matrix =
      matthew_post_matrix,

    period_ref_year =
      period_ref_year,

    visit_gap_days =
      visit_gap_days,

    survey_events_used =
      used_events,

    survey_events_all =
      survey_events_all,

    detection_adult_age =
      detection_adult_age,

    Origin =
      Origin,

    tMax =
      tMax,

    ordered_islands =
      ordered_islands
  ),
  chain_file
)


cat("\n")
cat("============================================================\n")
cat("CHAIN COMPLETE\n")
cat("============================================================\n")
cat(
  "Chain ",
  CHAIN_ID,
  " saved to:\n  ",
  chain_file,
  "\n",
  sep =
    ""
)
cat("Do not assess Rhat until all independent chains are combined.\n")
cat("============================================================\n")
