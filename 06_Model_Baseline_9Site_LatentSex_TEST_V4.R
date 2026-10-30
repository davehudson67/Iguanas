## ============================================================
## 06_Model_Baseline_9Site_LatentSex_TEST_V4.R
##
## CORRECTED NINE-SITE BASELINE SURVIVAL MODEL
##
## Purpose:
##   Establish a clean biological baseline BEFORE adding calendar-
##   varying provisioning, Hurricane Matthew, and male removals.
##
## Survival:
##   - continuous-age Gompertz hazard
##   - correct left truncation at FIRST LIVE CAPTURE
##   - age-10 reference parameterisation
##   - partially observed / latent sex
##   - within-site birth-cohort density
##   - linear birth-cohort effect
##   - residual site effect on age-10 mortality
##
## Detection:
##   - first capture excluded
##   - subsequent actual site-year survey events only
##   - sex
##   - age, saturating at 15 y
##   - residual site detection effect
##   - residual survey-event detection effect
##
## Density:
##   Primary baseline is deliberately conservative:
##   only individuals whose ENTIRE possible birth interval lies
##   within supported density years contribute to betaDENS.
##
##   Other individuals remain in the survival/detection model,
##   but their density contribution is fixed to zero (= site mean).
##
## IMPORTANT:
##   Provisioning, hurricane and removal terms are NOT added here.
##   Those vary by CALENDAR YEAR and require the next time-varying
##   hazard model. Adding them as static individual covariates here
##   would be biologically wrong.
##
## Default RUN_MODE = TEST.
##
## To run longer later:
##   Sys.setenv(RUN_MODE = "LONG")
##   source("06_Model_Baseline_9Site_LatentSex_TEST.R")
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
## 2. LOAD CORRECT LEFT-TRUNCATED GOMPERTZ DISTRIBUTION
## ============================================================

gompertz_candidates <- c(
  "Distributions/Dist_GompertzLT.R",
  "Dist_GompertzLT.R"
)


gompertz_file <- gompertz_candidates[
  file.exists(
    gompertz_candidates
  )
][1]


if(is.na(gompertz_file)) {
  stop(
    "Cannot find Dist_GompertzLT.R. Put the corrected distribution ",
    "either in the project root or Distributions/."
  )
}


source(
  gompertz_file
)


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
## 4. BUILD ACTUAL SURVEY-EVENT TABLE
##
## One event = one monitored site x calendar year.
##
## The event time is the mean decimal capture date for live
## encounters in that site-year.
## ============================================================

survey_events_all <- AllData %>%
  filter(
    alive_encounter
  ) %>%
  group_by(
    Island,
    CalYear
  ) %>%
  summarise(
    event_time_cal =
      mean(
        decimal_date(
          CapDate
        )
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


## ============================================================
## 5. FIRST LIVE CAPTURE YEAR BY INDIVIDUAL
## ============================================================

first_live <- AllData %>%
  filter(
    alive_encounter
  ) %>%
  group_by(
    ID
  ) %>%
  summarise(
    first_live_year =
      min(
        CalYear
      ),

    first_live_time_cal =
      min(
        decimal_date(
          CapDate
        )
      ),

    .groups =
      "drop"
  )


stopifnot(
  nrow(
    first_live
  ) ==
    nind
)


## ============================================================
## 6. OBSERVED LIVE DETECTIONS BY INDIVIDUAL x SITE-YEAR
## ============================================================

live_detection_by_year <- AllData %>%
  filter(
    alive_encounter
  ) %>%
  group_by(
    ID,
    Island,
    CalYear
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
## We exclude the entire first-capture site-year from detection
## likelihood. There are no pre-first-capture zeros.
##
## For a detected pair we use that individual's actual capture date.
## For a non-detection we use the site-year event time.
##
## This avoids making a capture-then-death animal appear dead before
## the survey at which it was actually caught alive.
## ============================================================

individual_lookup <- growth_info %>%
  transmute(
    ID =
      ID,

    Island =
      island
  ) %>%
  left_join(
    first_live,
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
    CalYear >
      first_live_year
  ) %>%
  left_join(
    live_detection_by_year,
    by =
      c(
        "ID",
        "Island",
        "CalYear"
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
cat("DELAYED-ENTRY SURVEY DETECTION HISTORY\n")
cat("============================================================\n")
cat(
  "Recognised site-year survey events: ",
  nrow(
    survey_events_all
  ),
  "\n",
  sep =
    ""
)
cat(
  "Survey events used after first capture: ",
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


    ## dGompertzLT is parameterised as h(age) = a * exp(b*age).
    ## Convert from age-10 hazard to birth-scale a.
    a_birth_i[i] <-
      exp(
        log_a10_i[i] -
        b_i[i] *
        age_ref
      )


    ## Age at death, conditional on surviving to first capture.
    tstar[i] ~ dGompertzLT(
      a_birth_i[i],
      b_i[i],
      lowerBound =
        entry_age[i]
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

  age_ref =
    age_ref,

  detection_adult_age =
    detection_adult_age,

  time_ref_year =
    time_ref_year,

  time_scale_years =
    time_scale_years,

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
    density_matrix
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


nchains <- 3


inits_BASE <- lapply(
  seq_len(
    nchains
  ),
  make_inits
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
  "betaCOHORT"
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
## 15. RUN SETTINGS
## ============================================================

RUN_MODE <- toupper(
  Sys.getenv(
    "RUN_MODE",
    "TEST"
  )
)


if(
  RUN_MODE ==
    "LONG"
) {

  niter <- 150000
  nburnin <- 50000
  thin <- 10

} else {

  RUN_MODE <- "TEST"

  niter <- 3000
  nburnin <- 1000
  thin <- 2
}


cat("\n")
cat("============================================================\n")
cat("RUNNING CORRECTED NINE-SITE BASELINE MODEL\n")
cat("============================================================\n")
cat(
  "Mode: ",
  RUN_MODE,
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
  "Chains: ",
  nchains,
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


run_BASE <- runMCMC(
  cMCMC_BASE,
  niter =
    niter,
  nburnin =
    nburnin,
  nchains =
    nchains,
  thin =
    thin,
  inits =
    inits_BASE,
  setSeed =
    TRUE,
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
## 16. PRIMARY SUMMARY
## ============================================================

primary_parameters <- c(
  "a10_ref",
  "b_ref",
  "log_a10",
  "log_b",
  "betaSEX",
  "betaDENS_a",
  "betaDENS_b",
  "betaCOHORT",
  "sigma_island_a",
  "alpha_p",
  "betaP_SEX",
  "betaP_AGE",
  "sigma_island_p",
  "sigma_survey_p"
)


sm_BASE <- MCMCsummary(
  run_BASE$samples,
  params =
    primary_parameters
) %>%
  as.data.frame() %>%
  rownames_to_column(
    "parameter"
  ) %>%
  as_tibble()


cat("\n")
cat("============================================================\n")
cat("PRIMARY POSTERIOR SUMMARY\n")
cat("============================================================\n")

print(
  sm_BASE,
  n = Inf
)


## ============================================================
## 17. LATENT-SEX POSTERIOR
## ============================================================

sample_matrix <- as.matrix(
  run_BASE$samples
)


latent_sex_summary <- tibble(
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


latent_sex_summary$Pr_male <- vapply(
  latent_sex_summary$ID,
  function(i) {

    node_name <- paste0(
      "sex[",
      i,
      "]"
    )

    mean(
      sample_matrix[
        ,
        node_name
      ]
    )
  },
  numeric(
    1
  )
)


latent_sex_summary <- latent_sex_summary %>%
  mutate(
    posterior_class =
      case_when(
        Pr_male >=
          0.90 ~
          "strong_male",

        Pr_male <=
          0.10 ~
          "strong_female",

        TRUE ~
          "uncertain"
      )
  ) %>%
  arrange(
    desc(
      abs(
        Pr_male -
          0.5
      )
    )
  )


cat("\n")
cat("============================================================\n")
cat("LATENT SEX POSTERIOR\n")
cat("============================================================\n")
cat(
  "Strong male P>=0.90: ",
  sum(
    latent_sex_summary$Pr_male >=
      0.90
  ),
  "\n",
  sep =
    ""
)
cat(
  "Strong female P<=0.10: ",
  sum(
    latent_sex_summary$Pr_male <=
      0.10
  ),
  "\n",
  sep =
    ""
)
cat(
  "Remaining uncertain: ",
  sum(
    latent_sex_summary$Pr_male >
      0.10 &
      latent_sex_summary$Pr_male <
      0.90
  ),
  "\n",
  sep =
    ""
)


## ============================================================
## 18. SAVE
## ============================================================

out_dir <-
  "DataAudits/Baseline9Site"


dir.create(
  out_dir,
  recursive =
    TRUE,
  showWarnings =
    FALSE
)


write_csv(
  sm_BASE,
  file.path(
    out_dir,
    paste0(
      "Baseline9Site_",
      RUN_MODE,
      "_PrimarySummary.csv"
    )
  )
)


write_csv(
  latent_sex_summary,
  file.path(
    out_dir,
    paste0(
      "Baseline9Site_",
      RUN_MODE,
      "_LatentSexPosterior.csv"
    )
  )
)


write_csv(
  used_events,
  file.path(
    out_dir,
    "SurveyDetection_UsedEvents.csv"
  )
)


write_csv(
  encounter_pairs,
  file.path(
    out_dir,
    "SurveyDetection_EncounterPairs.csv"
  )
)


saveRDS(
  list(
    fit =
      run_BASE,

    summary =
      sm_BASE,

    latent_sex =
      latent_sex_summary,

    density_eligible =
      density_eligible,

    density_support =
      birth_density_support,

    survey_events_used =
      used_events,

    survey_events_all =
      survey_events_all,

    encounter_pairs =
      encounter_pairs,

    detection_adult_age =
      detection_adult_age,

    age_ref =
      age_ref,

    time_ref_year =
      time_ref_year,

    run_mode =
      RUN_MODE
  ),
  paste0(
    "Model_Baseline9Site_LatentSex_",
    RUN_MODE,
    "_COMPLETE.rds"
  )
)


cat("\n")
cat("============================================================\n")
cat("BASELINE MODEL COMPLETE\n")
cat("============================================================\n")
cat(
  "Saved model: Model_Baseline9Site_LatentSex_",
  RUN_MODE,
  "_COMPLETE.rds\n",
  sep =
    ""
)
cat(
  "Audit output: ",
  out_dir,
  "\n",
  sep =
    ""
)
cat("============================================================\n")
