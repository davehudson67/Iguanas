## ============================================================
## 14_PosteriorPredictiveChecks_Primary_V2.R
##
## Posterior predictive checks for the converged PRIMARY model.
##
## This script performs a joint survival + detection PPC on the
## subset with reliable known-age information and resolved sex.
##
## Why this subset?
## The long primary run did not monitor all 4,152 latent birth/death
## times. For the known-age subset, birth windows are already tight,
## so we can propagate that small remaining uncertainty directly
## while generating replicated capture histories.
##
## Replicated data preserve:
##   - each individual's observed first-capture time
##   - their island
##   - the actual reconstructed 14-day field-visit schedule
##   - birth-density treatment used in the primary model
##   - high/regular provisioning history
##   - FFRC post-2016 Matthew indicator
##   - posterior site mortality/detection effects
##   - new survey-event random effects drawn from the fitted
##     survey-effect distribution
##
## Main checks:
##   - total post-entry detections
##   - proportion never recaptured
##   - detections per individual
##   - long recapture spans (5, 10, 20 years)
##   - event-level capture-count variation
##   - site-specific detection totals
##
## Usage:
##   N_PPC=1000 PPC_CORES=4 Rscript 14_PosteriorPredictiveChecks_Primary.R
##
## ============================================================


library(tidyverse)
library(lubridate)
library(coda)
library(parallel)


## ============================================================
## 1. LOAD PRIMARY POSTERIOR + FROZEN DATA
## ============================================================

## IMPORTANT:
## DataPrep_Canonical_LatentSex_V5.R starts by clearing the R
## workspace. Therefore it MUST be sourced before loading the
## posterior object or defining objects needed below.
source(
  "DataPrep_Canonical_LatentSex_V5.R"
)


fit_file <-
  "Model_TimeVarying_Feeding_Matthew_LONG_3CHAINS_COMPLETE.rds"


if(
  !file.exists(
    fit_file
  )
) {
  stop(
    "Cannot find ",
    fit_file
  )
}


fit_obj <- readRDS(
  fit_file
)


cov_obj <- readRDS(
  "Data/Iguana_ModelCovariates_Through2019.rds"
)


density_matrix <-
  cov_obj$density_matrix


density_eligible <-
  as.integer(
    fit_obj$density_eligible
  )


feeding_matrix <-
  fit_obj$feeding_high_matrix


matthew_matrix <-
  fit_obj$matthew_post_matrix


ordered_islands <-
  fit_obj$ordered_islands


n_islands <-
  length(
    ordered_islands
  )


visit_gap_days <-
  fit_obj$visit_gap_days


period_ref_year <-
  fit_obj$period_ref_year


sample_matrix <- as.matrix(
  fit_obj$samples
)


## ============================================================
## 2. SETTINGS
## ============================================================

N_PPC <- as.integer(
  Sys.getenv(
    "N_PPC",
    "1000"
  )
)


PPC_CORES <- as.integer(
  Sys.getenv(
    "PPC_CORES",
    "4"
  )
)


if(
  is.na(
    N_PPC
  ) |
  N_PPC <
    50
) {
  stop(
    "N_PPC must be at least 50."
  )
}


if(
  is.na(
    PPC_CORES
  ) |
  PPC_CORES <
    1
) {
  PPC_CORES <-
    1L
}


PPC_CORES <- min(
  PPC_CORES,
  parallel::detectCores(
    logical =
      FALSE
  )
)


set.seed(
  73021
)


draw_ids <- sample(
  seq_len(
    nrow(
      sample_matrix
    )
  ),
  size =
    N_PPC,
  replace =
    N_PPC >
      nrow(
        sample_matrix
      )
)


## ============================================================
## 3. RECONSTRUCT THE SAME 14-DAY FIELD VISITS
## ============================================================

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


encounter_pairs <- encounter_pairs %>%
  inner_join(
    used_events %>%
      select(
        event_id_all,
        event_id
      ),
    by =
      "event_id_all"
  )


if(
  nrow(
    used_events
  ) !=
    nrow(
      fit_obj$survey_events_used
    )
) {
  stop(
    "Reconstructed field-visit count does not match fitted model."
  )
}


## ============================================================
## 4. PPC SUBSET: KNOWN AGE + RESOLVED SEX
## ============================================================

ppc_individuals <- growth_info %>%
  filter(
    known_age_consistent,
    !is.na(
      known_birth_low
    ),
    !is.na(
      known_birth_high
    ),
    sex_state %in%
      c(
        "observed_F",
        "observed_M"
      )
  ) %>%
  transmute(
    ID =
      ID,

    ID_raw =
      ID_raw,

    Island =
      island,

    island_idx =
      match(
        island,
        ordered_islands
      ),

    sex =
      if_else(
        sex_state ==
          "observed_M",
        1L,
        0L
      ),

    entry_time =
      tF_cal -
      Origin +
      1,

    birth_low =
      cintB_adj[
        ID,
        1
      ],

    birth_high =
      cintB_adj[
        ID,
        2
      ]
  ) %>%
  arrange(
    ID
  ) %>%
  mutate(
    ppc_ind =
      row_number()
  )


if(
  nrow(
    ppc_individuals
  ) <
    500
) {
  stop(
    "Unexpectedly small known-age/resolved-sex PPC subset."
  )
}


pairs_ppc <- encounter_pairs %>%
  semi_join(
    ppc_individuals,
    by =
      "ID"
  ) %>%
  left_join(
    ppc_individuals %>%
      select(
        ID,
        ppc_ind
      ),
    by =
      "ID"
  ) %>%
  mutate(
    obs_time =
      event_time_cal_pair -
      Origin +
      1
  ) %>%
  arrange(
    ppc_ind,
    event_id
  )


## ============================================================
## 5. HELPERS
## ============================================================

get_indexed <- function(
    draw,
    stem,
    n
) {

  draw[
    paste0(
      stem,
      "[",
      seq_len(
        n
      ),
      "]"
    )
  ]
}


simulate_death_age <- function(
    a10,
    b,
    tB,
    lower_bound,
    feed,
    matthew,
    beta_feed,
    beta_matthew,
    site_period_slope,
    Origin,
    period_ref_year,
    tMax
) {

  remaining_H <- rexp(
    1,
    rate =
      1
  )


  current_cal <-
    tB +
    lower_bound


  study_end_boundary <-
    tMax +
    1


  while(
    current_cal <
      study_end_boundary
  ) {

    year_idx <- floor(
      current_cal
    )


    year_idx <- max(
      1,
      min(
        tMax,
        year_idx
      )
    )


    next_boundary <-
      min(
        floor(
          current_cal
        ) +
          1,
        study_end_boundary
      )


    age_start <-
      current_cal -
      tB


    age_end <-
      next_boundary -
      tB


    calendar_year <-
      Origin +
      year_idx -
      1


    period_decade <-
      (
        calendar_year -
        period_ref_year
      ) /
      10


    eta <-
      beta_feed *
      feed[
        year_idx
      ] +
      beta_matthew *
      matthew[
        year_idx
      ] +
      site_period_slope *
      period_decade


    multiplier <- exp(
      eta
    )


    H_seg <-
      multiplier *
      (
        a10 /
        b
      ) *
      (
        exp(
          b *
          (
            age_end -
            10
          )
        ) -
        exp(
          b *
          (
            age_start -
            10
          )
        )
      )


    if(
      remaining_H <=
        H_seg
    ) {

      base_start <-
        exp(
          b *
          (
            age_start -
            10
          )
        )


      inside <-
        base_start +
        (
          remaining_H *
          b
        ) /
        (
          multiplier *
          a10
        )


      return(
        10 +
          log(
            inside
          ) /
          b
      )
    }


    remaining_H <-
      remaining_H -
      H_seg


    current_cal <-
      next_boundary
  }


  ## Proper tail beyond final observed year: final calendar
  ## multiplier held fixed, matching the fitted distribution.
  year_idx <-
    tMax


  age_start <-
    current_cal -
    tB


  calendar_year <-
    Origin +
    tMax -
    1


  period_decade <-
    (
      calendar_year -
      period_ref_year
    ) /
    10


  eta <-
    beta_feed *
    feed[
      year_idx
    ] +
    beta_matthew *
    matthew[
      year_idx
    ] +
    site_period_slope *
    period_decade


  multiplier <- exp(
    eta
  )


  base_start <-
    exp(
      b *
      (
        age_start -
        10
      )
    )


  inside <-
    base_start +
    (
      remaining_H *
      b
    ) /
    (
      multiplier *
      a10
    )


  10 +
    log(
      inside
    ) /
    b
}


capture_history_stats <- function(
    y,
    pair_df,
    indiv_df
) {

  det_counts <- tabulate(
    pair_df$ppc_ind[
      y ==
        1L
    ],
    nbins =
      nrow(
        indiv_df
      )
  )


  first_time <- indiv_df$entry_time


  last_rep_time <- rep(
    NA_real_,
    nrow(
      indiv_df
    )
  )


  detected_rows <- which(
    y ==
      1L
  )


  if(
    length(
      detected_rows
    ) >
      0
  ) {

    split_times <- split(
      pair_df$obs_time[
        detected_rows
      ],
      pair_df$ppc_ind[
        detected_rows
      ]
    )


    idx <- as.integer(
      names(
        split_times
      )
    )


    last_rep_time[
      idx
    ] <- vapply(
      split_times,
      max,
      numeric(
        1
      )
    )
  }


  recapture_span <- ifelse(
    is.na(
      last_rep_time
    ),
    0,
    pmax(
      0,
      last_rep_time -
        first_time
    )
  )


  event_counts <- tabulate(
    pair_df$event_id[
      y ==
        1L
    ],
    nbins =
      nrow(
        used_events
      )
  )


  c(
    total_detections =
      sum(
        y
      ),

    prop_never_recaptured =
      mean(
        det_counts ==
          0
      ),

    median_detections_per_individual =
      median(
        det_counts
      ),

    q90_detections_per_individual =
      unname(
        quantile(
          det_counts,
          0.90
        )
      ),

    n_recap_span_ge_5y =
      sum(
        recapture_span >=
          5
      ),

    n_recap_span_ge_10y =
      sum(
        recapture_span >=
          10
      ),

    n_recap_span_ge_20y =
      sum(
        recapture_span >=
          20
      ),

    max_recapture_span_years =
      max(
        recapture_span
      ),

    sd_event_detection_count =
      sd(
        event_counts
      ),

    max_event_detection_count =
      max(
        event_counts
      )
  )
}


count_categories <- function(
    y,
    pair_df,
    indiv_df
) {

  det_counts <- tabulate(
    pair_df$ppc_ind[
      y ==
        1L
    ],
    nbins =
      nrow(
        indiv_df
      )
  )


  c(
    `0` =
      mean(
        det_counts ==
          0
      ),

    `1` =
      mean(
        det_counts ==
          1
      ),

    `2` =
      mean(
        det_counts ==
          2
      ),

    `3` =
      mean(
        det_counts ==
          3
      ),

    `4` =
      mean(
        det_counts ==
          4
      ),

    `5+` =
      mean(
        det_counts >=
          5
      )
  )
}


site_counts <- function(
    y,
    pair_df
) {

  out <- tapply(
    y,
    factor(
      pair_df$Island,
      levels =
        ordered_islands
    ),
    sum
  )


  as.numeric(
    out
  )
}


## ============================================================
## 6. OBSERVED PPC STATISTICS
## ============================================================

y_obs <- as.integer(
  pairs_ppc$detected
)


observed_stats <- capture_history_stats(
  y_obs,
  pairs_ppc,
  ppc_individuals
)


observed_count_categories <- count_categories(
  y_obs,
  pairs_ppc,
  ppc_individuals
)


observed_site_counts <- site_counts(
  y_obs,
  pairs_ppc
)


## ============================================================
## 7. ONE POSTERIOR-PREDICTIVE REPLICATE
## ============================================================

run_one_ppc <- function(
    k
) {

  set.seed(
    810000 +
      k
  )


  draw <- sample_matrix[
    draw_ids[
      k
    ],
    ,
    drop =
      TRUE
  ]


  tB <- runif(
    nrow(
      ppc_individuals
    ),
    min =
      ppc_individuals$birth_low,
    max =
      ppc_individuals$birth_high
  )


  cohort_idx <- pmax(
    1,
    pmin(
      tMax,
      floor(
        tB
      ) +
        1
    )
  )


  dens <- numeric(
    nrow(
      ppc_individuals
    )
  )


  for(ii in seq_len(nrow(ppc_individuals))) {

    original_id <-
      ppc_individuals$ID[
        ii
      ]


    dens[
      ii
    ] <-
      density_eligible[
        original_id
      ] *
      density_matrix[
        ppc_individuals$island_idx[
          ii
        ],
        cohort_idx[
          ii
        ]
      ]
  }


  birth_year_cal <-
    tB +
    Origin -
    1


  birth_decade_c <-
    (
      birth_year_cal -
      1990
    ) /
    10


  island_a <- get_indexed(
    draw,
    "island_a",
    n_islands
  )


  island_p <- get_indexed(
    draw,
    "island_p",
    n_islands
  )


  site_period <- get_indexed(
    draw,
    "site_period",
    n_islands
  )


  sex <- ppc_individuals$sex


  log_a10_i <-
    draw[
      "log_a10"
    ] +
    draw[
      "betaSEX[1]"
    ] *
    sex +
    draw[
      "betaDENS_a"
    ] *
    dens +
    draw[
      "betaCOHORT"
    ] *
    birth_decade_c +
    island_a[
      ppc_individuals$island_idx
    ]


  log_b_i <-
    draw[
      "log_b"
    ] +
    draw[
      "betaSEX[2]"
    ] *
    sex +
    draw[
      "betaDENS_b"
    ] *
    dens


  a10_i <- exp(
    log_a10_i
  )


  b_i <- exp(
    log_b_i
  )


  entry_age <- pmax(
    0,
    ppc_individuals$entry_time -
      tB
  )


  death_age <- numeric(
    nrow(
      ppc_individuals
    )
  )


  for(ii in seq_len(nrow(ppc_individuals))) {

    j <-
      ppc_individuals$island_idx[
        ii
      ]


    death_age[
      ii
    ] <- simulate_death_age(
      a10 =
        a10_i[
          ii
        ],
      b =
        b_i[
          ii
        ],
      tB =
        tB[
          ii
        ],
      lower_bound =
        entry_age[
          ii
        ],
      feed =
        feeding_matrix[
          j,
          ,
          drop =
            TRUE
        ],
      matthew =
        matthew_matrix[
          j,
          ,
          drop =
            TRUE
        ],
      beta_feed =
        draw[
          "betaFEED"
        ],
      beta_matthew =
        draw[
          "betaMATTHEW"
        ],
      site_period_slope =
        site_period[
          j
        ],
      Origin =
        Origin,
      period_ref_year =
        period_ref_year,
      tMax =
        tMax
    )
  }


  death_time <-
    tB +
    death_age


  event_raw <- rnorm(
    nrow(
      used_events
    )
  )


  event_raw <-
    event_raw -
    mean(
      event_raw
    )


  event_re <-
    draw[
      "sigma_survey_p"
    ] *
    event_raw


  pair_ind <-
    pairs_ppc$ppc_ind


  pair_age <- pmax(
    0,
    pairs_ppc$obs_time -
      tB[
        pair_ind
      ]
  )


  age_score <- pmin(
    1,
    pair_age /
      15
  )


  eta_p <-
    draw[
      "alpha_p"
    ] +
    draw[
      "betaP_SEX"
    ] *
    sex[
      pair_ind
    ] +
    draw[
      "betaP_AGE"
    ] *
    age_score +
    island_p[
      pairs_ppc$island_idx
    ] +
    event_re[
      pairs_ppc$event_id
    ]


  alive <-
    death_time[
      pair_ind
    ] >=
    pairs_ppc$obs_time


  p_det <-
    as.numeric(
      alive
    ) *
    plogis(
      eta_p
    )


  y_rep <- rbinom(
    length(
      p_det
    ),
    size =
      1,
    prob =
      p_det
  )


  list(
    stats =
      capture_history_stats(
        y_rep,
        pairs_ppc,
        ppc_individuals
      ),

    count_categories =
      count_categories(
        y_rep,
        pairs_ppc,
        ppc_individuals
      ),

    site_counts =
      site_counts(
        y_rep,
        pairs_ppc
      )
  )
}


## ============================================================
## 8. RUN PPC REPLICATES
## ============================================================

cat("\n")
cat("============================================================\n")
cat("PRIMARY MODEL POSTERIOR PREDICTIVE CHECKS\n")
cat("============================================================\n")
cat(
  "Posterior draws available: ",
  nrow(
    sample_matrix
  ),
  "\n",
  sep =
    ""
)
cat(
  "PPC replicates: ",
  N_PPC,
  "\n",
  sep =
    ""
)
cat(
  "Parallel cores: ",
  PPC_CORES,
  "\n",
  sep =
    ""
)
cat(
  "Known-age/resolved-sex individuals: ",
  nrow(
    ppc_individuals
  ),
  "\n",
  sep =
    ""
)
cat(
  "Post-entry encounter pairs in PPC subset: ",
  nrow(
    pairs_ppc
  ),
  "\n",
  sep =
    ""
)
cat("============================================================\n")


if(
  .Platform$OS.type ==
    "unix" &
  PPC_CORES >
    1
) {

  ppc_results <- parallel::mclapply(
    seq_len(
      N_PPC
    ),
    run_one_ppc,
    mc.cores =
      PPC_CORES,
    mc.preschedule =
      TRUE
  )

} else {

  ppc_results <- lapply(
    seq_len(
      N_PPC
    ),
    run_one_ppc
  )
}


## ============================================================
## 9. SUMMARISE PPC
## ============================================================

rep_stats <- bind_rows(
  lapply(
    seq_along(
      ppc_results
    ),
    function(k) {

      as_tibble_row(
        ppc_results[
          [
            k
          ]
        ]$stats
      ) %>%
        mutate(
          replicate =
            k,
          .before =
            1
        )
    }
  )
)


global_summary <- tibble(
  statistic =
    names(
      observed_stats
    ),

  observed =
    as.numeric(
      observed_stats
    )
) %>%
  mutate(
    pred_median =
      map_dbl(
        statistic,
        function(s) {
          median(
            rep_stats[
              [
                s
              ]
            ]
          )
        }
      ),

    pred_low =
      map_dbl(
        statistic,
        function(s) {
          unname(
            quantile(
              rep_stats[
                [
                  s
                ]
              ],
              0.025
            )
          )
        }
      ),

    pred_high =
      map_dbl(
        statistic,
        function(s) {
          unname(
            quantile(
              rep_stats[
                [
                  s
                ]
              ],
              0.975
            )
          )
        }
      ),

    Pr_rep_ge_obs =
      map_dbl(
        statistic,
        function(s) {
          mean(
            rep_stats[
              [
                s
              ]
            ] >=
              observed[
                statistic ==
                  s
              ]
          )
        }
      ),

    two_sided_ppp =
      2 *
      pmin(
        Pr_rep_ge_obs,
        1 -
          Pr_rep_ge_obs
      )
  )


site_rep_matrix <- do.call(
  rbind,
  lapply(
    ppc_results,
    function(x) {
      x$site_counts
    }
  )
)


colnames(
  site_rep_matrix
) <- ordered_islands


site_summary <- tibble(
  Island =
    ordered_islands,

  observed =
    observed_site_counts,

  pred_median =
    apply(
      site_rep_matrix,
      2,
      median
    ),

  pred_low =
    apply(
      site_rep_matrix,
      2,
      quantile,
      probs =
        0.025
    ),

  pred_high =
    apply(
      site_rep_matrix,
      2,
      quantile,
      probs =
        0.975
    ),

  Pr_rep_ge_obs =
    vapply(
      seq_along(
        ordered_islands
      ),
      function(j) {
        mean(
          site_rep_matrix[
            ,
            j
          ] >=
            observed_site_counts[
              j
            ]
        )
      },
      numeric(
        1
      )
    )
) %>%
  mutate(
    two_sided_ppp =
      2 *
      pmin(
        Pr_rep_ge_obs,
        1 -
          Pr_rep_ge_obs
      )
  )


count_rep_matrix <- do.call(
  rbind,
  lapply(
    ppc_results,
    function(x) {
      x$count_categories
    }
  )
)


count_category_summary <- tibble(
  category =
    colnames(
      count_rep_matrix
    ),

  observed =
    as.numeric(
      observed_count_categories
    ),

  pred_median =
    apply(
      count_rep_matrix,
      2,
      median
    ),

  pred_low =
    apply(
      count_rep_matrix,
      2,
      quantile,
      probs =
        0.025
    ),

  pred_high =
    apply(
      count_rep_matrix,
      2,
      quantile,
      probs =
        0.975
    )
)


## ============================================================
## 10. OUTPUTS + PLOTS
## ============================================================

out_dir <-
  "DataAudits/PosteriorPredictive_Primary"


dir.create(
  out_dir,
  recursive =
    TRUE,
  showWarnings =
    FALSE
)


write_csv(
  global_summary,
  file.path(
    out_dir,
    "PPC_GlobalStatistics.csv"
  )
)


write_csv(
  site_summary,
  file.path(
    out_dir,
    "PPC_SiteDetectionTotals.csv"
  )
)


write_csv(
  count_category_summary,
  file.path(
    out_dir,
    "PPC_DetectionsPerIndividual.csv"
  )
)


write_csv(
  rep_stats,
  file.path(
    out_dir,
    "PPC_ReplicateStatistics.csv"
  )
)


p_global <- global_summary %>%
  mutate(
    statistic =
      forcats::fct_reorder(
        statistic,
        observed
      )
  ) %>%
  ggplot(
    aes(
      x =
        statistic
    )
  ) +
  geom_linerange(
    aes(
      ymin =
        pred_low,
      ymax =
        pred_high
    ),
    linewidth =
      1
  ) +
  geom_point(
    aes(
      y =
        pred_median
    ),
    size =
      2
  ) +
  geom_point(
    aes(
      y =
        observed
    ),
    shape =
      4,
    size =
      3,
    stroke =
      1
  ) +
  coord_flip() +
  labs(
    x =
      NULL,
    y =
      "Statistic",
    title =
      "Posterior predictive checks: global capture-history statistics",
    subtitle =
      "Line = 95% posterior predictive interval; dot = predictive median; cross = observed"
  ) +
  theme_bw()


ggsave(
  file.path(
    out_dir,
    "PPC_GlobalStatistics.png"
  ),
  p_global,
  width =
    10,
  height =
    7,
  dpi =
    160
)


p_site <- site_summary %>%
  ggplot(
    aes(
      x =
        Island
    )
  ) +
  geom_linerange(
    aes(
      ymin =
        pred_low,
      ymax =
        pred_high
    ),
    linewidth =
      1
  ) +
  geom_point(
    aes(
      y =
        pred_median
    ),
    size =
      2
  ) +
  geom_point(
    aes(
      y =
        observed
    ),
    shape =
      4,
    size =
      3,
    stroke =
      1
  ) +
  coord_flip() +
  labs(
    x =
      NULL,
    y =
      "Post-entry detections",
    title =
      "Posterior predictive checks by site"
  ) +
  theme_bw()


ggsave(
  file.path(
    out_dir,
    "PPC_SiteDetectionTotals.png"
  ),
  p_site,
  width =
    8,
  height =
    6,
  dpi =
    160
)


p_counts <- count_category_summary %>%
  ggplot(
    aes(
      x =
        category
    )
  ) +
  geom_linerange(
    aes(
      ymin =
        pred_low,
      ymax =
        pred_high
    ),
    linewidth =
      1
  ) +
  geom_point(
    aes(
      y =
        pred_median
    ),
    size =
      2
  ) +
  geom_point(
    aes(
      y =
        observed
    ),
    shape =
      4,
    size =
      3,
    stroke =
      1
  ) +
  labs(
    x =
      "Number of post-entry detections",
    y =
      "Proportion of individuals",
    title =
      "Posterior predictive capture-frequency distribution"
  ) +
  theme_bw()


ggsave(
  file.path(
    out_dir,
    "PPC_DetectionsPerIndividual.png"
  ),
  p_counts,
  width =
    8,
  height =
    5,
  dpi =
    160
)


saveRDS(
  list(
    global_summary =
      global_summary,

    site_summary =
      site_summary,

    count_category_summary =
      count_category_summary,

    replicate_statistics =
      rep_stats,

    ppc_individuals =
      ppc_individuals,

    encounter_pairs =
      pairs_ppc,

    N_PPC =
      N_PPC,

    PPC_CORES =
      PPC_CORES
  ),
  file.path(
    out_dir,
    "PosteriorPredictive_Primary_COMPLETE.rds"
  )
)


cat("\n")
cat("GLOBAL PPC SUMMARY\n")
cat("------------------------------------------------------------\n")

print(
  global_summary,
  n = Inf
)


cat("\nSITE PPC SUMMARY\n")
cat("------------------------------------------------------------\n")

print(
  site_summary,
  n = Inf
)


cat("\n")
cat("============================================================\n")
cat("POSTERIOR PREDICTIVE CHECKS COMPLETE\n")
cat("============================================================\n")
cat(
  "Outputs: ",
  out_dir,
  "\n",
  sep =
    ""
)
cat("============================================================\n")
