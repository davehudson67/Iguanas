## ============================================================
## DataPrep_Canonical_LatentSex.R
##
## Revised data-prep pipeline for the cleaned/canonical iguana data.
##
## KEY CHANGES FROM THE OLD DataPrep.R
## ----------------------------------
## 1. Starts from:
##      Data/Iguana_Canonical_PreDataPrep_Through2019.rds
##
## 2. Uses 9 biological sites:
##      Bitter, FFRC, Gaulin, LeafNorth, LeafSouth,
##      NAdder, Noddy, UCay, WhiteB
##
## 3. Keeps Leaf North and Leaf South separate.
##
## 4. Sex handling:
##      * unknown + F only  -> F
##      * unknown + M only  -> M
##      * all unknown       -> latent
##      * both F and M      -> latent
##
##    sex_obs is therefore:
##      0 = observed female
##      1 = observed male
##     NA = latent sex in NIMBLE
##
## 5. SVL->age prediction respects latent sex:
##      * known sex: use sex-specific prediction
##      * latent sex: use UNION of F and M prediction intervals
##
## 6. Uses only consistent "known" northern ages as age anchors.
##    Individuals flagged known_age_consistent == FALSE are retained,
##    but their birth year is inferred from SVL like an unknown-age animal.
##
## 7. Permanent removals to Pasture Cay remain LIVE observations and
##    are right-censored after removal, never treated as deaths.
##
## 8. Does NOT create the old static individual-level feeding factor.
##    Time-varying site-year provisioning/tourism will be constructed
##    separately for the final annual survival model.
##
## 9. Does NOT apply the old blanket birth/death interval shrink.
##
## 10. Produces sex_obs and sex_init for partially observed latent sex.
## ============================================================

library(tidyverse)
library(lubridate)
library(quantreg)
library(nimble)

rm(list = ls())

## ============================================================
## 1. LOAD CANONICAL DATA
## ============================================================

canonical_file <- "Data/Iguana_Canonical_PreDataPrep_Through2019.rds"

if(!file.exists(canonical_file)) {
  stop(
    "Cannot find ",
    canonical_file,
    "\nRun 01_Combine_Clean_Raw_Iguana_Data_V2.R and ",
    "03_Resolve_Identity_and_Prepare_ModelData.R first."
  )
}

AllData_raw <- readRDS(canonical_file)

required_columns <- c(
  "Site",
  "Species",
  "individual_id_clean",
  "sex",
  "capture_date",
  "calendar_year",
  "svl_cm",
  "age_known_for_model",
  "known_age_consistent",
  "known_birth_low",
  "known_birth_high",
  "known_birth_mid",
  "comment_dead",
  "permanent_removal"
)

missing_columns <- setdiff(
  required_columns,
  names(AllData_raw)
)

if(length(missing_columns) > 0) {
  stop(
    "Canonical input is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

## ============================================================
## 2. STANDARDISE NAMES USED BY THE EXISTING MODEL CODE
## ============================================================

AllData <- AllData_raw %>%
  transmute(
    Island = Site,
    Species = Species,
    ID_raw = individual_id_clean,
    animal_id = animal_id,
    CapDate = as.Date(capture_date),
    CalYear = as.integer(calendar_year),
    svl = as.numeric(svl_cm),
    body_mass = as.numeric(body_mass_g),
    sex_record = case_when(
      sex == "F" ~ "F",
      sex == "M" ~ "M",
      TRUE ~ NA_character_
    ),
    age_known = as.numeric(age_known_for_model),
    known_age_consistent = coalesce(
      as.logical(known_age_consistent),
      FALSE
    ),
    known_birth_low = as.numeric(known_birth_low),
    known_birth_high = as.numeric(known_birth_high),
    known_birth_mid = as.numeric(known_birth_mid),
    comment_dead_raw = coalesce(
      as.logical(comment_dead),
      FALSE
    ),
    permanent_removal = coalesce(
      as.logical(permanent_removal),
      FALSE
    ),
    feeding_raw = if(
      "feeding_raw" %in% names(AllData_raw)
    ) {
      AllData_raw$feeding_raw
    } else {
      NA_character_
    },
    comments = if(
      "comments" %in% names(AllData_raw)
    ) {
      AllData_raw$comments
    } else {
      NA_character_
    },
    source_sheet = if(
      "source_sheet" %in% names(AllData_raw)
    ) {
      AllData_raw$source_sheet
    } else {
      NA_character_
    },
    source_row = if(
      "source_row" %in% names(AllData_raw)
    ) {
      AllData_raw$source_row
    } else {
      NA_integer_
    }
  ) %>%
  filter(
    !is.na(ID_raw),
    !is.na(CapDate),
    !is.na(CalYear),
    CalYear <= 2019
  )

## ============================================================
## 3. TRUE ANIMAL-DEATH CLASSIFICATION
##
## IMPORTANT:
## The original broad "dead" keyword flag also caught comments such as:
##   "pit dead"
##   "dead pit, new pit"
##   "old tag dead"
##
## Those refer to failed PIT tags, NOT dead iguanas.
##
## We therefore re-classify death from the original comments here.
## Strong animal-death wording overrides tag-failure wording.
## Historical/referential notes such as "[DEAD in 2004]" on an
## earlier capture row are not treated as a death on that row.
## ============================================================

death_text <- str_to_lower(
  str_squish(
    replace_na(
      AllData$comments,
      ""
    )
  )
)

has_dead_word <- str_detect(
  death_text,
  "\\bdead\\b"
)

strong_animal_death <- str_detect(
  death_text,
  paste0(
    "\\bcarcass\\b|",
    "\\bdeceased\\b|",
    "\\bdied\\b|",
    "\\bfound\\s+(?:newly\\s+|freshly\\s+)?dead\\b|",
    "\\bfresh(?:ly)?\\s+dead\\b|",
    "\\bdead\\s+(?:on|at|in|by|just|after|before)\\b|",
    "\\(dead\\)|",
    "\\bdead\\b.*\\bpit\\s+(?:ok|confirmed|recovered)\\b|",
    "\\bpit\\s+(?:ok|confirmed|recovered)\\b.*\\bdead\\b"
  )
)

pit_or_tag_failure <- str_detect(
  death_text,
  paste0(
    "\\bpit\\s+(?:tag\\s+)?(?:was\\s+)?dead\\b|",
    "\\bold\\s+(?:pit|tag)\\s+dead\\b|",
    "\\btag\\s+(?:was\\s+)?dead\\b|",
    "\\btransponder\\s+(?:was\\s+)?dead\\b|",
    "\\bdead\\s+pit(?:\\s*\\d+)?\\b.*",
    "(?:new\\s+pit|re[- ]?pit|palp|feel|broken|not\\s+palpable)"
  )
)

other_nonanimal_dead <- str_detect(
  death_text,
  "\\bdead\\s+tissue\\b|\\bdead\\s+skin\\b|\\bdead\\s+transmitter\\b"
)

## Notes that refer to a death in another year / another animal,
## rather than reporting death on this capture row.
historical_dead_reference <- str_detect(
  death_text,
  paste0(
    "\\[\\s*dead\\s+in\\s+[12][0-9]{3}\\s*\\]|",
    "\\bthis\\s+or\\s+next\\s+one\\s+dead\\s+in\\s+[12][0-9]{3}\\b|",
    "\\bfrom\\s+female\\b.*\\bdead\\s+in\\s+[12][0-9]{3}\\b"
  )
)

## If a comment explicitly says "dead in YYYY", only let that phrase
## classify the current row as death when YYYY matches the row year.
dead_year_in_text <- suppressWarnings(
  as.integer(
    str_match(
      death_text,
      "\\bdead\\s+in\\s+([12][0-9]{3})\\b"
    )[, 2]
  )
)

historical_dead_reference <- historical_dead_reference |
  (
    !is.na(
      dead_year_in_text
    ) &
      dead_year_in_text !=
        AllData$CalYear
  )

AllData <- AllData %>%
  mutate(
    animal_dead =
      (
        strong_animal_death |
          (
            has_dead_word &
              !pit_or_tag_failure &
              !other_nonanimal_dead
          )
      ) &
      !historical_dead_reference &
      !permanent_removal,

    death_classification =
      case_when(
        animal_dead ~
          "animal_death",

        historical_dead_reference ~
          "historical_or_referential_death_note",

        pit_or_tag_failure ~
          "pit_or_tag_failure_not_animal_death",

        other_nonanimal_dead ~
          "other_nonanimal_dead_text",

        comment_dead_raw ~
          "broad_old_dead_flag_not_confirmed",

        TRUE ~
          "no_death_text"
      ),

    Status =
      if_else(
        animal_dead,
        2L,
        1L
      )
  )


## ============================================================
## 4. DEATH AUDIT + LATER-RECAPTURE OVERRIDE
##
## A true death cannot be followed by a later live recapture of the
## same cleaned individual. Some historic comments use "dead" for
## something other than the animal itself.
##
## Example now verified in the raw workbook:
##   LeafNorth_100C, 2002-06-26:
##       "840/157 (dead) NOOO"
##   Same animal/PIT (77D91E) was recaptured repeatedly through 2013.
##
## Therefore:
##   1. classify candidate deaths from comments;
##   2. look for a later live record for that SAME cleaned ID;
##   3. if one exists, override the earlier candidate death to LIVE;
##   4. retain all such overrides in an audit file.
##
## This is safer than deleting later captures or treating an
## impossible early "death" as biological truth.
## ============================================================

prelim_death_rows <- AllData %>%
  filter(
    animal_dead
  ) %>%
  select(
    ID_raw,
    candidate_death_date =
      CapDate,
    candidate_death_year =
      CalYear,
    candidate_comments =
      comments
  )


later_recapture_override <- prelim_death_rows %>%
  inner_join(
    AllData %>%
      filter(
        !animal_dead
      ) %>%
      select(
        ID_raw,
        later_capture_date =
          CapDate,
        later_capture_year =
          CalYear,
        later_comments =
          comments
      ),
    by =
      "ID_raw"
  ) %>%
  filter(
    later_capture_date >
      candidate_death_date
  ) %>%
  group_by(
    ID_raw,
    candidate_death_date,
    candidate_death_year,
    candidate_comments
  ) %>%
  summarise(
    first_later_capture_date =
      min(
        later_capture_date
      ),

    last_later_capture_date =
      max(
        later_capture_date
      ),

    n_later_live_records =
      n(),

    .groups =
      "drop"
  )


if(nrow(later_recapture_override) > 0) {

  message(
    "Overriding ",
    nrow(
      later_recapture_override
    ),
    " candidate death row(s) because the same cleaned individual ",
    "was captured alive later."
  )

  AllData <- AllData %>%
    left_join(
      later_recapture_override %>%
        transmute(
          ID_raw,
          CapDate =
            candidate_death_date,
          death_overridden_by_later_recapture =
            TRUE
        ),
      by =
        c(
          "ID_raw",
          "CapDate"
        )
    ) %>%
    mutate(
      death_overridden_by_later_recapture =
        coalesce(
          death_overridden_by_later_recapture,
          FALSE
        ),

      animal_dead =
        animal_dead &
        !death_overridden_by_later_recapture,

      Status =
        if_else(
          animal_dead,
          2L,
          1L
        ),

      death_classification =
        if_else(
          death_overridden_by_later_recapture,
          "death_text_overridden_by_later_live_recapture",
          death_classification
        )
    )

} else {

  AllData <- AllData %>%
    mutate(
      death_overridden_by_later_recapture =
        FALSE
    )
}


## ------------------------------------------------------------
## Some genuine death records describe an animal that was captured
## alive and then died during/shortly after handling.
##
## These rows contribute BOTH:
##   * a live encounter for capture/detection, and
##   * a death observation shortly afterwards.
##
## Clear examples in the source data include:
##   "found dead 2 days later"
##   "dead 1 hr later"
##   "died while on boat"
##   "died in bag overnight"
## ------------------------------------------------------------

AllData <- AllData %>%
  mutate(
    capture_then_death =
      animal_dead &
      str_detect(
        death_text,
        paste0(
          "\\bdied\\b|",
          "\\bdeceased\\s+in\\s+(?:the\\s+)?bag\\b|",
          "\\bdead\\s+in\\s+(?:the\\s+)?bag\\b|",
          "\\bfound\\s+dead\\s+[0-9]+\\s*(?:hr|hrs|hour|hours|day|days)\\s+later\\b|",
          "\\bdead\\s+[0-9]+\\s*(?:hr|hrs|hour|hours|day|days)\\s+later\\b"
        )
      ),

    ## For annual-scale survival, all we need is a small positive
    ## separation between the live encounter and subsequent death.
    ## Where "2 days later" is explicitly stated, preserve that;
    ## otherwise use one day as a conservative within-year offset.
    death_delay_years =
      case_when(
        capture_then_death &
          str_detect(
            death_text,
            "\\b2\\s*days?\\s+later\\b"
          ) ~
          2 / 365.25,

        capture_then_death &
          str_detect(
            death_text,
            "\\b1\\s*(?:hr|hrs|hour|hours)\\s+later\\b"
          ) ~
          1 / (24 * 365.25),

        capture_then_death ~
          1 / 365.25,

        TRUE ~
          0
      ),

    ## A capture-then-death row is a genuine live encounter even
    ## though Status remains 2 for the death observation.
    alive_encounter =
      Status ==
        1L |
      capture_then_death,

    ## Continuous calendar-scale death time used for censoring.
    death_time_cal =
      if_else(
        animal_dead,
        decimal_date(
          CapDate
        ) +
          death_delay_years,
        NA_real_
      )
  )


death_audit <- AllData %>%
  filter(
    animal_dead |
      comment_dead_raw |
      death_overridden_by_later_recapture |
      str_detect(
        death_text,
        "\\bdead\\b|\\bcarcass\\b|\\bdeceased\\b|\\bdied\\b"
      )
  ) %>%
  select(
    ID_raw,
    Island,
    animal_id,
    CapDate,
    CalYear,
    sex_record,
    Status,
    animal_dead,
    capture_then_death,
    alive_encounter,
    death_delay_years,
    death_overridden_by_later_recapture,
    comment_dead_raw,
    death_classification,
    comments,
    source_sheet,
    source_row
  ) %>%
  arrange(
    Island,
    ID_raw,
    CapDate
  )


death_dates <- AllData %>%
  filter(
    Status ==
      2L
  ) %>%
  group_by(
    ID_raw
  ) %>%
  summarise(
    first_dead_date =
      min(
        CapDate,
        na.rm = TRUE
      ),

    first_dead_year =
      year(
        first_dead_date
      ),

    .groups =
      "drop"
  )


post_death_records <- AllData %>%
  left_join(
    death_dates,
    by =
      "ID_raw"
  ) %>%
  filter(
    !is.na(
      first_dead_date
    ),
    CapDate >
      first_dead_date,
    Status ==
      1L
  ) %>%
  arrange(
    ID_raw,
    CapDate
  )


dir.create(
  "DataAudits/DataPrep",
  recursive = TRUE,
  showWarnings = FALSE
)


write_csv(
  death_audit,
  "DataAudits/DataPrep/DeathClassificationAudit.csv"
)


write_csv(
  later_recapture_override,
  "DataAudits/DataPrep/DeathOverrides_LaterRecapture.csv"
)


write_csv(
  post_death_records,
  "DataAudits/DataPrep/PostDeathLiveRecords.csv"
)


message(
  "Confirmed animal-death rows after override: ",
  sum(
    AllData$Status ==
      2L
  )
)

message(
  "Death-text rows overridden by later live recapture: ",
  sum(
    AllData$death_overridden_by_later_recapture
  )
)
message(
  "Capture-then-death rows counted as live encounters: ",
  sum(
    AllData$capture_then_death
  )
)

message(
  "PIT/tag-death text correctly retained as live: ",
  sum(
    AllData$death_classification ==
      "pit_or_tag_failure_not_animal_death"
  )
)

message(
  "Historical/referential death notes retained as live: ",
  sum(
    AllData$death_classification ==
      "historical_or_referential_death_note"
  )
)


if(nrow(post_death_records) > 0) {

  print(
    post_death_records %>%
      select(
        ID_raw,
        Island,
        animal_id,
        CapDate,
        CalYear,
        first_dead_date,
        sex_record,
        svl,
        comments
      ),
    n = Inf
  )

  stop(
    nrow(
      post_death_records
    ),
    " live records still occur after a classified animal death. ",
    "Nothing has been deleted. Inspect ",
    "DataAudits/DataPrep/DeathClassificationAudit.csv and ",
    "PostDeathLiveRecords.csv."
  )
}

## ============================================================
## 5. INDIVIDUAL SEX STATE
##
## Agreed rule:
##   unknown + female only -> female
##   unknown + male only   -> male
##   both female and male  -> latent
##   all unknown           -> latent
## ============================================================

sex_by_individual <- AllData %>%
  group_by(ID_raw) %>%
  summarise(
    n_F = sum(
      sex_record == "F",
      na.rm = TRUE
    ),
    n_M = sum(
      sex_record == "M",
      na.rm = TRUE
    ),
    n_unknown = sum(
      is.na(sex_record)
    ),
    sex_state = case_when(
      n_F > 0 & n_M == 0 ~ "observed_F",
      n_M > 0 & n_F == 0 ~ "observed_M",
      n_F > 0 & n_M > 0 ~ "conflict_latent",
      TRUE ~ "all_unknown_latent"
    ),
    sex_clean = case_when(
      sex_state == "observed_F" ~ "F",
      sex_state == "observed_M" ~ "M",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  )

AllData <- AllData %>%
  left_join(
    sex_by_individual,
    by = "ID_raw"
  )

## ============================================================
## 6. NUMERIC INDIVIDUAL ID
## ============================================================

id_levels <- sort(
  unique(AllData$ID_raw)
)

AllData <- AllData %>%
  mutate(
    ID = as.integer(
      factor(
        ID_raw,
        levels = id_levels
      )
    )
  )

nind <- length(id_levels)

## ============================================================
## 7. GROWTH CALIBRATION: SVL -> AGE
##
## Use only:
##   * genuinely known/near-known age records,
##   * age-consistent individuals,
##   * known sex.
## ============================================================

calibration_data <- AllData %>%
  filter(
    known_age_consistent,
    !is.na(age_known),
    !is.na(svl),
    age_known > 0,
    !is.na(sex_clean)
  ) %>%
  add_count(
    ID,
    name = "n_calibration_records_id"
  ) %>%
  mutate(
    ## Retain all longitudinal known-age measurements, but give each
    ## biological individual equal TOTAL influence in the quantile
    ## calibration regardless of how often it was recaptured.
    calibration_weight =
      1 / n_calibration_records_id,

    ## Constant rescaling does not alter weighted quantile-regression
    ## coefficients and keeps weights numerically well behaved.
    calibration_weight =
      calibration_weight /
      mean(calibration_weight),

    age_for_calibration = pmin(
      age_known,
      40
    ),
    log_age = log(
      age_for_calibration
    ),
    sex_f = factor(
      sex_clean,
      levels = c(
        "F",
        "M"
      )
    )
  )

if(nrow(calibration_data) == 0) {
  stop(
    "No usable known-age/known-sex records for SVL-age calibration."
  )
}

if(length(unique(calibration_data$sex_f)) < 2) {
  stop(
    "SVL-age calibration requires both female and male known-sex records."
  )
}

fit_min_log <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.05,
  data = calibration_data,
  weights = calibration_weight,
  method = "fn"
)

fit_max_log <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.95,
  data = calibration_data,
  weights = calibration_weight,
  method = "fn"
)

## ============================================================
## 8. INDIVIDUAL LIFE-HISTORY SUMMARY
## ============================================================

growth_info <- AllData %>%
  arrange(
    ID,
    CapDate
  ) %>%
  group_by(ID) %>%
  summarise(
    ID_raw = ID_raw[1],
    island = Island[1],
    species = Species[1],
    sex_state = sex_state[1],
    sex_clean = sex_clean[1],
    n_F = n_F[1],
    n_M = n_M[1],
    n_unknown_sex = n_unknown[1],

    first_svl = {
      x <- svl[
        !is.na(svl)
      ]

      if(length(x) == 0) {
        NA_real_
      } else {
        x[1]
      }
    },

    first_svl_year = {
      x <- CalYear[
        !is.na(svl)
      ]

      if(length(x) == 0) {
        NA_real_
      } else {
        x[1]
      }
    },

    first_capture_year =
      min(
        CalYear,
        na.rm = TRUE
      ),

    total_growth = {
      x <- svl[
        !is.na(svl)
      ]

      if(length(x) < 2) {
        0
      } else {
        max(x) - min(x)
      }
    },

    known_age_consistent = any(
      known_age_consistent
    ),

    known_birth_low = {
      x <- known_birth_low[
        !is.na(known_birth_low) &
          known_age_consistent
      ]

      if(length(x) == 0) {
        NA_real_
      } else {
        x[1]
      }
    },

    known_birth_high = {
      x <- known_birth_high[
        !is.na(known_birth_high) &
          known_age_consistent
      ]

      if(length(x) == 0) {
        NA_real_
      } else {
        x[1]
      }
    },

    known_birth_mid = {
      x <- known_birth_mid[
        !is.na(known_birth_mid) &
          known_age_consistent
      ]

      if(length(x) == 0) {
        NA_real_
      } else {
        x[1]
      }
    },

    ## Continuous calendar time is used for survival/censoring.
    ## This matters when an animal is captured alive and then dies
    ## later in the SAME calendar year.
    tF_cal = {
      x <- decimal_date(
        CapDate[
          alive_encounter
        ]
      )

      if(length(x) == 0) {
        NA_real_
      } else {
        min(
          x,
          na.rm = TRUE
        )
      }
    },

    tL_cal = {
      x <- decimal_date(
        CapDate[
          alive_encounter
        ]
      )

      if(length(x) == 0) {
        NA_real_
      } else {
        max(
          x,
          na.rm = TRUE
        )
      }
    },

    tD_cal = {
      x <- death_time_cal[
        animal_dead
      ]

      if(length(x) == 0) {
        NA_real_
      } else {
        min(
          x,
          na.rm = TRUE
        )
      }
    },

    permanent_removal = any(
      permanent_removal
    ),

    removal_year = {
      x <- CalYear[
        permanent_removal
      ]

      if(length(x) == 0) {
        NA_real_
      } else {
        min(x)
      }
    },

    .groups = "drop"
  ) %>%
  arrange(ID)

## ============================================================
## 9. SAFETY CHECKS
## ============================================================

if(any(is.na(growth_info$tF_cal))) {
  stop(
    "Some individuals have missing first-capture year."
  )
}

if(any(is.na(growth_info$tL_cal))) {
  bad <- growth_info %>%
    filter(
      is.na(tL_cal)
    )

  print(
    bad,
    n = Inf
  )

  stop(
    "Some individuals have no live observation."
  )
}

same_year_live_death <- growth_info %>%
  filter(
    !is.na(
      tD_cal
    ),
    floor(
      tD_cal
    ) ==
      floor(
        tL_cal
      ),
    tD_cal >
      tL_cal
  )

message(
  "Same-calendar-year live-then-death histories resolved with exact dates: ",
  nrow(
    same_year_live_death
  )
)

if(
  any(
    !is.na(growth_info$tD_cal) &
      growth_info$tD_cal <= growth_info$tL_cal
  )
) {
  bad <- growth_info %>%
    filter(
      !is.na(tD_cal),
      tD_cal <= tL_cal
    )

  print(
    bad %>%
      select(
        ID,
        ID_raw,
        island,
        sex_state,
        tF_cal,
        tL_cal,
        tD_cal,
        permanent_removal,
        removal_year
      ),
    n = Inf
  )

  stop(
    "Some individuals still have tD_cal <= tL_cal even after using exact dates."
  )
}

if(
  any(
    growth_info$permanent_removal &
      !is.na(growth_info$tD_cal)
  )
) {
  bad <- growth_info %>%
    filter(
      permanent_removal,
      !is.na(tD_cal)
    )

  print(
    bad,
    n = Inf
  )

  stop(
    "An individual is both permanently removed and later recorded dead."
  )
}

## ============================================================
## 10. SEX-SPECIFIC SVL -> AGE BOUNDS
##
## Latent-sex individuals use the union of the F and M bounds.
## ============================================================

predict_age_bounds <- function(
  svl_value,
  sex_value
) {
  if(is.na(svl_value)) {
    return(
      c(
        min_age = 0,
        max_age = 40
      )
    )
  }

  nd <- data.frame(
    svl = svl_value,
    sex_f = factor(
      sex_value,
      levels = c(
        "F",
        "M"
      )
    )
  )

  p_min <- as.numeric(
    exp(
      predict(
        fit_min_log,
        newdata = nd
      )
    )
  )

  p_max <- as.numeric(
    exp(
      predict(
        fit_max_log,
        newdata = nd
      )
    )
  )

  p_max <- min(
    p_max,
    40
  )

  p_min <- min(
    p_min,
    p_max - 1
  )

  if(svl_value < 25) {
    p_max <- min(
      p_max,
      4
    )
  }

  min_age <- max(
    0,
    p_min
  )

  max_age <- max(
    min_age + 1,
    p_max
  )

  c(
    min_age = min_age,
    max_age = max_age
  )
}

female_bounds <- t(
  vapply(
    growth_info$first_svl,
    predict_age_bounds,
    sex_value = "F",
    FUN.VALUE = c(
      min_age = 0,
      max_age = 0
    )
  )
)

male_bounds <- t(
  vapply(
    growth_info$first_svl,
    predict_age_bounds,
    sex_value = "M",
    FUN.VALUE = c(
      min_age = 0,
      max_age = 0
    )
  )
)

growth_info <- growth_info %>%
  mutate(
    pred_min_age_F = female_bounds[, "min_age"],
    pred_max_age_F = female_bounds[, "max_age"],
    pred_min_age_M = male_bounds[, "min_age"],
    pred_max_age_M = male_bounds[, "max_age"],

    pred_min_age = case_when(
      sex_clean == "F" ~ pred_min_age_F,
      sex_clean == "M" ~ pred_min_age_M,
      TRUE ~ pmin(
        pred_min_age_F,
        pred_min_age_M
      )
    ),

    pred_max_age = case_when(
      sex_clean == "F" ~ pred_max_age_F,
      sex_clean == "M" ~ pred_max_age_M,
      TRUE ~ pmax(
        pred_max_age_F,
        pred_max_age_M
      )
    ),

    pred_max_age = if_else(
      total_growth > 5,
      pmin(
        pred_max_age,
        15
      ),
      pred_max_age
    ),

    pred_min_age = pmax(
      0,
      pred_min_age
    ),

    pred_max_age = pmax(
      pred_min_age + 1,
      pred_max_age
    )
  )

## ============================================================
## 11. BUILD CALENDAR-SCALE BIRTH WINDOWS
## ============================================================

cintB_cal <- matrix(
  NA_real_,
  nrow = nind,
  ncol = 2
)

for(i in seq_len(nind)) {
  use_known_birth <-
    growth_info$known_age_consistent[i] &&
    !is.na(
      growth_info$known_birth_low[i]
    ) &&
    !is.na(
      growth_info$known_birth_high[i]
    )

  if(use_known_birth) {
    lower <- growth_info$known_birth_low[i]
    upper <- growth_info$known_birth_high[i]

    if(lower >= upper) {
      mid <- growth_info$known_birth_mid[i]
      lower <- mid - 0.5
      upper <- mid + 0.5
    }

    cintB_cal[i, ] <- c(
      lower,
      upper
    )
  } else {
    ref_yr <- if(
      !is.na(
        growth_info$first_svl_year[i]
      )
    ) {
      growth_info$first_svl_year[i]
    } else {
      growth_info$first_capture_year[i]
    }

    cintB_cal[i, 1] <- floor(
      ref_yr -
        growth_info$pred_max_age[i]
    )

    cintB_cal[i, 2] <- floor(
      ref_yr -
        growth_info$pred_min_age[i]
    )
  }

  if(any(is.na(cintB_cal[i, ]))) {
    stop(
      "Missing birth interval for individual ",
      growth_info$ID_raw[i]
    )
  }

  if(cintB_cal[i, 1] >= cintB_cal[i, 2]) {
    cintB_cal[i, 1] <- cintB_cal[i, 2] - 1
  }

  cintB_cal[i, 2] <- min(
    cintB_cal[i, 2],
    growth_info$tF_cal[i]
  )

  if(cintB_cal[i, 1] >= cintB_cal[i, 2]) {
    cintB_cal[i, 1] <- cintB_cal[i, 2] - 1
  }
}

## ============================================================
## 12. DYNAMIC TIME ORIGIN
## ============================================================

Origin <- floor(
  min(
    cintB_cal[, 1],
    na.rm = TRUE
  )
)

cintB_adj <- cintB_cal -
  Origin +
  1

if(any(cintB_adj[, 1] <= 0)) {
  stop(
    "Adjusted birth lower bound is non-positive."
  )
}

if(any(cintB_adj[, 1] >= cintB_adj[, 2])) {
  stop(
    "Invalid adjusted birth interval."
  )
}

## ============================================================
## 13. DEATH / RIGHT-CENSORING INTERVALS
## ============================================================

cintD_cal <- matrix(
  NA_real_,
  nrow = nind,
  ncol = 2
)

censoredD_vec <- integer(
  nind
)

for(i in seq_len(nind)) {
  if(is.na(growth_info$tD_cal[i])) {
    cintD_cal[i, 1] <- 0
    cintD_cal[i, 2] <- growth_info$tL_cal[i]
    censoredD_vec[i] <- 2L
  } else {
    cintD_cal[i, 1] <- growth_info$tL_cal[i]
    cintD_cal[i, 2] <- growth_info$tD_cal[i]
    censoredD_vec[i] <- 1L
  }
}

cintD_adj <- cintD_cal

cintD_adj[, 1] <- ifelse(
  cintD_cal[, 1] == 0,
  0,
  cintD_cal[, 1] - Origin + 1
)

cintD_adj[, 2] <- cintD_cal[, 2] - Origin + 1

if(any(is.na(cintD_adj))) {
  stop(
    "NA values found in death/censoring intervals."
  )
}

if(
  any(
    censoredD_vec == 1L &
      cintD_adj[, 2] <= cintD_adj[, 1]
  )
) {
  bad <- which(
    censoredD_vec == 1L &
      cintD_adj[, 2] <= cintD_adj[, 1]
  )

  print(
    tibble(
      i = bad,
      ID_raw = growth_info$ID_raw[bad],
      island = growth_info$island[bad],
      last_alive = growth_info$tL_cal[bad],
      found_dead = growth_info$tD_cal[bad]
    ),
    n = Inf
  )

  stop(
    "Invalid interval-censored death window."
  )
}

bad_birth_death <- which(
  censoredD_vec == 1L &
    cintB_adj[, 1] >= cintD_adj[, 2]
)

if(length(bad_birth_death) > 0) {
  print(
    tibble(
      i = bad_birth_death,
      ID_raw = growth_info$ID_raw[bad_birth_death],
      island = growth_info$island[bad_birth_death],
      birth_low_cal =
        cintB_adj[bad_birth_death, 1] +
        Origin -
        1,
      birth_high_cal =
        cintB_adj[bad_birth_death, 2] +
        Origin -
        1,
      death_upper_cal =
        cintD_adj[bad_birth_death, 2] +
        Origin -
        1
    ),
    n = Inf
  )

  stop(
    "Impossible birth/death support found for interval-censored deaths."
  )
}

## ============================================================
## 14. SEX DATA FOR NIMBLE
##
## sex_obs:
##   0  observed female
##   1  observed male
##   NA latent
##
## In NIMBLE:
##   sex[i] ~ dbern(0.5)
##
## Pass:
##   data = list(sex = sex_obs, ...)
## ============================================================

sex_obs <- case_when(
  growth_info$sex_clean == "F" ~ 0,
  growth_info$sex_clean == "M" ~ 1,
  TRUE ~ NA_real_
)

sex_is_latent <- is.na(
  sex_obs
)

sex_init <- ifelse(
  sex_is_latent,
  0,
  sex_obs
)

sex_state <- growth_info$sex_state

## ============================================================
## 15. SPECIES AND ISLAND INDICES
## ============================================================

ordered_islands <- c(
  "Bitter",
  "FFRC",
  "Gaulin",
  "LeafNorth",
  "LeafSouth",
  "NAdder",
  "Noddy",
  "UCay",
  "WhiteB"
)

if(
  !setequal(
    unique(growth_info$island),
    ordered_islands
  )
) {
  stop(
    "Unexpected site set. Found: ",
    paste(
      sort(
        unique(growth_info$island)
      ),
      collapse = ", "
    )
  )
}

n_islands <- length(
  ordered_islands
)

island_idx <- as.integer(
  factor(
    growth_info$island,
    levels = ordered_islands
  )
)

ordered_species <- c(
  "figginsi",
  "inornata"
)

species_idx <- as.integer(
  factor(
    growth_info$species,
    levels = ordered_species
  )
)

## ============================================================
## 16. CAPTURE HISTORY MATRIX
## ============================================================

tMax <- max(
  AllData$CalYear,
  na.rm = TRUE
) -
  Origin +
  1

CH <- matrix(
  0L,
  nrow = nind,
  ncol = tMax
)

for(r in seq_len(nrow(AllData))) {
  i <- AllData$ID[r]
  t <- AllData$CalYear[r] - Origin + 1

  if(AllData$alive_encounter[r]) {
    CH[i, t] <- 1L
  }
}

y <- as.numeric(
  rowSums(CH)
)

## ============================================================
## 17. SITE-YEAR SURVEY EVENTS / EFFORT
## ============================================================

survey_events <- AllData %>%
  filter(
    alive_encounter
  ) %>%
  distinct(
    Island,
    CalYear
  ) %>%
  arrange(
    CalYear,
    Island
  ) %>%
  mutate(
    survey_event_id = row_number(),
    island_idx = match(
      Island,
      ordered_islands
    ),
    year_adj = CalYear -
      Origin +
      1
  )

n_survey_events <- nrow(
  survey_events
)

effort_matrix <- matrix(
  0L,
  nrow = n_islands,
  ncol = tMax
)

for(r in seq_len(nrow(survey_events))) {
  effort_matrix[
    survey_events$island_idx[r],
    survey_events$year_adj[r]
  ] <- 1L
}

cum_effort_adj <- t(
  apply(
    effort_matrix,
    1,
    cumsum
  )
)

cum_effort_for_model <- cbind(
  0,
  cum_effort_adj
)

island_starts_cal <- survey_events %>%
  group_by(
    Island
  ) %>%
  summarise(
    S = min(
      CalYear
    ),
    .groups = "drop"
  )

island_start <- island_starts_cal$S[
  match(
    growth_info$island,
    island_starts_cal$Island
  )
] -
  Origin +
  1

island_start_time <- island_start

island_start_idx <- pmax(
  1,
  pmin(
    tMax,
    floor(
      island_start_time
    ) +
      1
  )
)

## ============================================================
## 18. FIRST / LAST / DEATH TIMES
## ============================================================

tF <- growth_info$tF_cal -
  Origin +
  1

tL <- growth_info$tL_cal -
  Origin +
  1

tD <- ifelse(
  is.na(growth_info$tD_cal),
  NA,
  growth_info$tD_cal -
    Origin +
    1
)

## ============================================================
## 19. INITIAL VALUES FOR LATENT BIRTH AND DEATH
## ============================================================

tB_init <- rowMeans(
  cintB_adj
)

tD_init <- numeric(
  nind
)

for(i in seq_len(nind)) {
  if(censoredD_vec[i] == 1L) {
    tD_init[i] <- mean(
      cintD_adj[i, ]
    )
  } else if(censoredD_vec[i] == 2L) {
    tD_init[i] <- cintD_adj[i, 2] + 1
  } else {
    stop(
      "Unexpected censoring category for individual ",
      i
    )
  }

  if(tD_init[i] <= tB_init[i]) {
    tD_init[i] <- tB_init[i] + 1
  }
}

tstar_init <- pmax(
  1,
  tD_init -
    tB_init
)

## ============================================================
## 20. LATENT-SEX AUDIT
## ============================================================

sex_audit <- growth_info %>%
  transmute(
    ID = ID,
    ID_raw = ID_raw,
    island = island,
    species = species,
    sex_state = sex_state,
    n_F = n_F,
    n_M = n_M,
    n_unknown = n_unknown_sex,
    sex_obs = sex_obs,
    sex_is_latent = sex_is_latent,
    first_svl = first_svl,
    pred_min_age_F = pred_min_age_F,
    pred_max_age_F = pred_max_age_F,
    pred_min_age_M = pred_min_age_M,
    pred_max_age_M = pred_max_age_M,
    birth_low_cal =
      cintB_adj[, 1] +
      Origin -
      1,
    birth_high_cal =
      cintB_adj[, 2] +
      Origin -
      1
  )

dir.create(
  "DataAudits/DataPrep",
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  sex_audit,
  "DataAudits/DataPrep/LatentSexAudit.csv"
)

write_csv(
  growth_info,
  "DataAudits/DataPrep/GrowthInfo.csv"
)

write_csv(
  survey_events,
  "DataAudits/DataPrep/SurveyEvents.csv"
)

## ============================================================
## 21. FINAL SANITY CHECKS
## ============================================================

stopifnot(
  nrow(growth_info) == nind,
  nrow(cintB_adj) == nind,
  ncol(cintB_adj) == 2,
  nrow(cintD_adj) == nind,
  ncol(cintD_adj) == 2,
  length(censoredD_vec) == nind,
  length(sex_obs) == nind,
  length(sex_init) == nind,
  length(sex_is_latent) == nind,
  length(species_idx) == nind,
  length(island_idx) == nind,
  length(y) == nind,
  length(tB_init) == nind,
  length(tstar_init) == nind,
  all(cintB_adj[, 1] < cintB_adj[, 2]),
  all(cintD_adj[, 2] > 0),
  all(
    censoredD_vec %in%
      c(
        1L,
        2L
      )
  ),
  all(
    island_idx >= 1
  ),
  all(
    island_idx <= n_islands
  ),
  all(
    species_idx %in%
      c(
        1L,
        2L
      )
  ),
  all(
    tB_init >= cintB_adj[, 1]
  ),
  all(
    tB_init <= cintB_adj[, 2]
  ),
  all(
    tstar_init > 0
  )
)

## ============================================================
## 22. CONSOLE SUMMARY
## ============================================================

message(
  "------------------------------------------------------------"
)
message(
  "Revised canonical DataPrep verified."
)
message(
  "Origin Year / Model Year 1: ",
  Origin
)
message(
  "Individuals: ",
  nind
)
message(
  "Model duration: ",
  tMax,
  " years"
)
message(
  "Sites: ",
  n_islands
)
message(
  "Survey events (site-years): ",
  n_survey_events
)
message(
  "Confirmed animal-death rows after override: ",
  sum(
    AllData$Status == 2L
  )
)
message(
  "Same-calendar-year live-then-death histories: ",
  nrow(
    same_year_live_death
  )
)
message(
  "Death-text rows overridden by later live recapture: ",
  sum(
    AllData$death_overridden_by_later_recapture
  )
)
message(
  "PIT/tag death wording retained as live: ",
  sum(
    AllData$death_classification ==
      "pit_or_tag_failure_not_animal_death"
  )
)
message(
  "Observed female: ",
  sum(
    sex_obs == 0,
    na.rm = TRUE
  )
)
message(
  "Observed male: ",
  sum(
    sex_obs == 1,
    na.rm = TRUE
  )
)
message(
  "Latent sex total: ",
  sum(
    sex_is_latent
  )
)
message(
  "  conflicting F/M: ",
  sum(
    sex_state == "conflict_latent"
  )
)
message(
  "  all sex records unknown: ",
  sum(
    sex_state == "all_unknown_latent"
  )
)
message(
  "Average birth-window width: ",
  round(
    mean(
      cintB_adj[, 2] -
        cintB_adj[, 1]
    ),
    2
  ),
  " years"
)
message(
  "Death censoring categories:"
)
print(
  table(
    censoredD_vec
  )
)
message(
  "Birth interval summary:"
)
print(
  summary(
    cintB_adj
  )
)
message(
  "Initial lifespan summary:"
)
print(
  summary(
    tstar_init
  )
)
message(
  "------------------------------------------------------------"
)

## ============================================================
## 23. HOW TO USE sex_obs IN NIMBLE
##
## In the NIMBLE model code:
##
##   for(i in 1:nind) {
##     sex[i] ~ dbern(0.5)
##     ...
##   }
##
## Pass:
##
##   data = list(
##     sex = sex_obs,
##     ...
##   )
##
## and initialise:
##
##   inits = list(
##     sex = sex_init,
##     ...
##   )
##
## Observed 0/1 values stay fixed; NA values are sampled.
## ============================================================
