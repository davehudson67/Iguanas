## ============================================================
## 07_Audit_SurveyOccasions_and_CalendarExposures_V2.R
##
## PURPOSE
## -------
## Before the long baseline / time-varying survival model:
##
## 1. Audit whether a "survey occasion" should remain a site-year
##    or whether some site-years clearly contain multiple distinct
##    field visits.
##
## 2. Quantify how sensitive reconstructed field-visit counts are
##    to plausible gaps between capture dates.
##
## 3. Build the final site-year provisioning history to be used in
##    the calendar-varying survival model.
##
## 4. Build Hurricane Matthew covariates for FFRC.
##
## IMPORTANT
## ---------
## - Provisioning is a SITE x YEAR exposure, not an individual
##   feeding treatment.
##
## - Raw feeding labels are audited, not used to assign individual
##   causal treatment.
##
## - The reported 2019-2023 decline in visitation on Leaf North /
##   U Cay is NOT coded as a 2019 treatment change because our
##   permitted data stop in 2019 and there is effectively no
##   post-change follow-up.
##
## - Targeted large-male removal is deliberately omitted from all
##   model covariates. We do not have strong enough evidence to
##   parameterise it defensibly.
##
## ============================================================


library(tidyverse)
library(lubridate)


## ============================================================
## 1. LOAD FROZEN CLEAN DATA
## ============================================================

source(
  "DataPrep_Canonical_LatentSex_V5.R"
)


out_dir <- "DataAudits/SurveyOccasions_Feeding"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 2. LIVE CAPTURE DAYS
## ============================================================

live_rows <- AllData %>%
  filter(
    alive_encounter,
    !is.na(
      CapDate
    ),
    !is.na(
      Island
    ),
    !is.na(
      CalYear
    )
  )


capture_days <- live_rows %>%
  distinct(
    Island,
    CalYear,
    CapDate
  ) %>%
  arrange(
    Island,
    CalYear,
    CapDate
  )


## ============================================================
## 3. SITE-YEAR DATE-SPAN AUDIT
## ============================================================

site_year_date_audit <- capture_days %>%
  group_by(
    Island,
    CalYear
  ) %>%
  summarise(
    n_capture_days =
      n(),

    first_capture_date =
      min(
        CapDate
      ),

    last_capture_date =
      max(
        CapDate
      ),

    span_days =
      as.integer(
        last_capture_date -
        first_capture_date
      ),

    max_gap_days =
      {
        dd <- sort(
          unique(
            CapDate
          )
        )

        if(
          length(
            dd
          ) <=
            1
        ) {
          0
        } else {
          max(
            as.integer(
              diff(
                dd
              )
            )
          )
        }
      },

    median_gap_days =
      {
        dd <- sort(
          unique(
            CapDate
          )
        )

        if(
          length(
            dd
          ) <=
            1
        ) {
          0
        } else {
          median(
            as.integer(
              diff(
                dd
              )
            )
          )
        }
      },

    .groups =
      "drop"
  ) %>%
  arrange(
    desc(
      max_gap_days
    ),
    desc(
      span_days
    )
  )


## Every actual gap between successive live-capture dates.
date_gaps <- capture_days %>%
  group_by(
    Island,
    CalYear
  ) %>%
  arrange(
    CapDate,
    .by_group =
      TRUE
  ) %>%
  mutate(
    previous_date =
      lag(
        CapDate
      ),

    gap_days =
      as.integer(
        CapDate -
        previous_date
      )
  ) %>%
  ungroup() %>%
  filter(
    !is.na(
      gap_days
    )
  ) %>%
  arrange(
    desc(
      gap_days
    )
  )


## ============================================================
## 4. RECONSTRUCT CANDIDATE FIELD VISITS
##
## We do not pick a threshold in advance.
##
## For each candidate threshold:
##   a new visit begins whenever the next capture day is MORE
##   than threshold_days after the previous capture day.
##
## We compare:
##   3, 5, 7, 10, 14, 21, 30 days
##
## If results stabilise over a sensible range, that gives us an
## empirical basis for defining field visits.
## ============================================================

candidate_thresholds <- c(
  3,
  5,
  7,
  10,
  14,
  21,
  30
)


build_candidate_visits <- function(
    threshold_days
) {

  day_assignments <- capture_days %>%
    group_by(
      Island,
      CalYear
    ) %>%
    arrange(
      CapDate,
      .by_group =
        TRUE
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
              threshold_days
          )
        ),

      visit_within_year =
        cumsum(
          new_visit
        )
    ) %>%
    ungroup()


  visit_rows <- live_rows %>%
    inner_join(
      day_assignments %>%
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


  visits <- visit_rows %>%
    group_by(
      Island,
      CalYear,
      visit_within_year
    ) %>%
    summarise(
      visit_start =
        min(
          CapDate
        ),

      visit_end =
        max(
          CapDate
        ),

      visit_span_days =
        as.integer(
          visit_end -
          visit_start
        ),

      n_capture_days =
        n_distinct(
          CapDate
        ),

      n_live_rows =
        n(),

      n_individuals =
        n_distinct(
          ID
        ),

      .groups =
        "drop"
    ) %>%
    arrange(
      Island,
      CalYear,
      visit_start
    ) %>%
    mutate(
      threshold_days =
        threshold_days,

      candidate_visit_id =
        row_number()
    )


  site_year_visit_counts <- visits %>%
    count(
      Island,
      CalYear,
      name =
        "n_visits"
    )


  summary_row <- tibble(
    threshold_days =
      threshold_days,

    n_candidate_visits =
      nrow(
        visits
      ),

    n_site_years =
      n_distinct(
        paste(
          visits$Island,
          visits$CalYear
        )
      ),

    n_site_years_split_into_multiple_visits =
      sum(
        site_year_visit_counts$n_visits >
          1
      ),

    median_visit_span_days =
      median(
        visits$visit_span_days
      ),

    q95_visit_span_days =
      unname(
        quantile(
          visits$visit_span_days,
          0.95
        )
      ),

    max_visit_span_days =
      max(
        visits$visit_span_days
      ),

    median_capture_days_per_visit =
      median(
        visits$n_capture_days
      ),

    median_individuals_per_visit =
      median(
        visits$n_individuals
      )
  )


  list(
    visits =
      visits,

    day_assignments =
      day_assignments,

    summary =
      summary_row
  )
}


candidate_results <- map(
  candidate_thresholds,
  build_candidate_visits
)


threshold_comparison <- map_dfr(
  candidate_results,
  "summary"
)


names(
  candidate_results
) <- paste0(
  "gap_",
  candidate_thresholds,
  "d"
)


## Save detailed candidate visit tables for the most plausible
## middle thresholds so we can inspect them without rerunning.
for(th in c(
  7,
  10,
  14,
  21
)) {

  obj <- candidate_results[[
    paste0(
      "gap_",
      th,
      "d"
    )
  ]]

  write_csv(
    obj$visits,
    file.path(
      out_dir,
      paste0(
        "CandidateVisits_GapGT",
        th,
        "Days.csv"
      )
    )
  )
}


## ============================================================
## 5. SITE-YEAR VERSUS FIELD-VISIT SUPPORT
##
## Flag the site-years where a site-year detection occasion is
## most questionable because capture activity spans a long period
## or contains a large internal gap.
## ============================================================

flagged_site_years <- site_year_date_audit %>%
  filter(
    span_days >
      21 |
      max_gap_days >
        14
  ) %>%
  arrange(
    desc(
      max_gap_days
    ),
    desc(
      span_days
    )
  )


## ============================================================
## 6. FINAL COAUTHOR-BASED PROVISIONING HISTORY THROUGH 2019
##
## Primary state:
##   none
##   low_intermit
##   high_regular
##
## Primary model contrast:
##   high_regular vs other
##
## Sensitivity:
##   any provisioning vs none
##
## Note:
## LeafNorth and UCay stay high_regular through 2019.
## We do not force the later 2019-2023 visitation decline into
## the final permitted year.
## ============================================================

feeding_history <- tribble(
  ~Site,        ~start_year, ~end_year, ~feeding_state,  ~evidence_note,

  "Bitter",          1980,       2005,  "none",
  "No organised boat provisioning before the reported 2006 onset.",

  "Bitter",          2006,       2019,  "high_regular",
  "Organised tourist boats began around 2006; intense/daily feeding thereafter.",

  "FFRC",            1980,       2019,  "none",
  "No provisioning history reported.",

  "Gaulin",          1980,       2019,  "low_intermit",
  "Longstanding lower-intensity visitation/provisioning; reportedly waned as Bitter became more popular.",

  "LeafNorth",       1980,       2019,  "high_regular",
  "Tourism/provisioning already established when long-term monitoring began; no defensible within-data cessation before 2019.",

  "LeafSouth",       1980,       2011,  "none",
  "Provisioning not documented before approximately 2012.",

  "LeafSouth",       2012,       2019,  "high_regular",
  "Intense/daily feeding reported from approximately 2012 onward.",

  "NAdder",          1980,       2019,  "none",
  "Reported not visited/provisioned.",

  "Noddy",           1980,       2019,  "none",
  "Reported not visited/provisioned.",

  "UCay",            1980,       2019,  "high_regular",
  "Regular provisioning/tourism throughout the permitted monitoring period.",

  "WhiteB",          1980,       1997,  "none",
  "Provisioning not documented before the late 1990s.",

  "WhiteB",          1998,       2019,  "low_intermit",
  "Modest/localised provisioning from the late 1990s; later waning retained qualitatively rather than assigned a numerical decline."
)


calendar_years <- seq.int(
  Origin,
  2019
)


expand_feeding_history <- function(
    site_name
) {

  map_dfr(
    calendar_years,
    function(yr) {

      hit <- feeding_history %>%
        filter(
          Site ==
            site_name,
          start_year <=
            yr,
          end_year >=
            yr
        )


      if(
        nrow(
          hit
        ) >
          1
      ) {
        stop(
          "Overlapping feeding history for ",
          site_name,
          " in ",
          yr
        )
      }


      if(
        nrow(
          hit
        ) ==
          0
      ) {

        return(
          tibble(
            Site =
              site_name,

            calendar_year =
              yr,

            feeding_state =
              NA_character_,

            evidence_note =
              NA_character_
          )
        )
      }


      tibble(
        Site =
          site_name,

        calendar_year =
          yr,

        feeding_state =
          hit$feeding_state[1],

        evidence_note =
          hit$evidence_note[1]
      )
    }
  )
}


feeding_site_year <- map_dfr(
  ordered_islands,
  expand_feeding_history
) %>%
  mutate(
    feeding_high =
      as.integer(
        feeding_state ==
          "high_regular"
      ),

    feeding_any =
      as.integer(
        feeding_state %in%
          c(
            "low_intermit",
            "high_regular"
          )
      ),

    year_idx =
      calendar_year -
      Origin +
      1
  )


## ============================================================
## 7. FEEDING SUPPORT IN THE OBSERVED DATA
## ============================================================

observed_site_years <- live_rows %>%
  distinct(
    Island,
    CalYear
  ) %>%
  rename(
    Site =
      Island,
    calendar_year =
      CalYear
  )


feeding_survey_support <- observed_site_years %>%
  left_join(
    feeding_site_year,
    by =
      c(
        "Site",
        "calendar_year"
      )
  ) %>%
  group_by(
    Site,
    feeding_state,
    feeding_high,
    feeding_any
  ) %>%
  summarise(
    n_survey_years =
      n(),

    first_year =
      min(
        calendar_year
      ),

    last_year =
      max(
        calendar_year
      ),

    .groups =
      "drop"
  ) %>%
  arrange(
    Site,
    first_year
  )


feeding_capture_support <- live_rows %>%
  transmute(
    Site =
      Island,

    calendar_year =
      CalYear,

    ID
  ) %>%
  left_join(
    feeding_site_year,
    by =
      c(
        "Site",
        "calendar_year"
      )
  ) %>%
  group_by(
    Site,
    feeding_state,
    feeding_high,
    feeding_any
  ) %>%
  summarise(
    n_live_rows =
      n(),

    n_individuals =
      n_distinct(
        ID
      ),

    .groups =
      "drop"
  ) %>%
  arrange(
    Site,
    feeding_state
  )


## Which sites actually provide within-site feeding change?
feeding_transition_support <- observed_site_years %>%
  left_join(
    feeding_site_year,
    by =
      c(
        "Site",
        "calendar_year"
      )
  ) %>%
  group_by(
    Site
  ) %>%
  summarise(
    n_primary_states_observed =
      n_distinct(
        feeding_high,
        na.rm =
          TRUE
      ),

    n_detailed_states_observed =
      n_distinct(
        feeding_state,
        na.rm =
          TRUE
      ),

    first_observed_year =
      min(
        calendar_year
      ),

    last_observed_year =
      max(
        calendar_year
      ),

    first_high_year =
      if(
        any(
          feeding_high ==
            1,
          na.rm =
            TRUE
        )
      ) {
        min(
          calendar_year[
            feeding_high ==
              1
          ],
          na.rm =
            TRUE
        )
      } else {
        NA_integer_
      },

    last_nonhigh_year =
      if(
        any(
          feeding_high ==
            0,
          na.rm =
            TRUE
        )
      ) {
        max(
          calendar_year[
            feeding_high ==
              0
          ],
          na.rm =
            TRUE
        )
      } else {
        NA_integer_
      },

    .groups =
      "drop"
  ) %>%
  arrange(
    desc(
      n_primary_states_observed
    ),
    Site
  )


## ============================================================
## 8. RAW FEEDING-LABEL AUDIT
##
## Retain these as observational metadata only.
## They are NOT used to assign individual treatment.
## ============================================================

if(
  "feeding_raw" %in%
    names(
      live_rows
    )
) {

  raw_feeding_audit <- live_rows %>%
    mutate(
      feeding_raw_clean =
        str_to_lower(
          str_squish(
            as.character(
              feeding_raw
            )
          )
        ),

      feeding_raw_clean =
        na_if(
          feeding_raw_clean,
          ""
        )
    ) %>%
    count(
      Island,
      CalYear,
      feeding_raw_clean,
      name =
        "n_live_rows"
    ) %>%
    arrange(
      Island,
      CalYear,
      desc(
        n_live_rows
      )
    )

} else {

  raw_feeding_audit <- tibble(
    Island =
      character(),
    CalYear =
      integer(),
    feeding_raw_clean =
      character(),
    n_live_rows =
      integer()
  )

  warning(
    "feeding_raw is not present in AllData."
  )
}


## ============================================================
## 9. HURRICANE MATTHEW SUPPORT
##
## Matthew occurred in 2016.
##
## Primary survival intervention later:
##   FFRC post-2016 = 2017 onward
##
## Sensitivity:
##   acute 2016-2018
## ============================================================

matthew_site_year <- expand_grid(
  Site =
    ordered_islands,

  calendar_year =
    calendar_years
) %>%
  mutate(
    matthew_post2016 =
      as.integer(
        Site ==
          "FFRC" &
          calendar_year >=
            2017
      ),

    matthew_acute_2016_2018 =
      as.integer(
        Site ==
          "FFRC" &
          calendar_year >=
            2016 &
          calendar_year <=
            2018
      ),

    year_idx =
      calendar_year -
      Origin +
      1
  )


matthew_support <- live_rows %>%
  filter(
    Island ==
      "FFRC"
  ) %>%
  mutate(
    matthew_period =
      case_when(
        CalYear <=
          2015 ~
          "pre_2016",

        CalYear ==
          2016 ~
          "2016",

        CalYear >=
          2017 ~
          "post_2016"
      )
  ) %>%
  group_by(
    matthew_period
  ) %>%
  summarise(
    n_live_rows =
      n(),

    n_individuals =
      n_distinct(
        ID
      ),

    n_survey_years =
      n_distinct(
        CalYear
      ),

    first_year =
      min(
        CalYear
      ),

    last_year =
      max(
        CalYear
      ),

    .groups =
      "drop"
  )


## ============================================================
## 10. BUILD FINAL CALENDAR-EXPOSURE MATRICES
##
## These will feed the next annual time-varying survival model.
## ============================================================

feeding_high_matrix <- matrix(
  0L,
  nrow =
    n_islands,
  ncol =
    tMax,
  dimnames =
    list(
      ordered_islands,
      seq_len(
        tMax
      )
    )
)


feeding_any_matrix <- feeding_high_matrix

matthew_post_matrix <- feeding_high_matrix

matthew_acute_matrix <- feeding_high_matrix


for(r in seq_len(nrow(feeding_site_year))) {

  jj <- match(
    feeding_site_year$Site[r],
    ordered_islands
  )

  tt <- feeding_site_year$year_idx[r]

  if(
    tt >=
      1 &
      tt <=
        tMax
  ) {

    feeding_high_matrix[
      jj,
      tt
    ] <-
      replace_na(
        feeding_site_year$feeding_high[r],
        0L
      )

    feeding_any_matrix[
      jj,
      tt
    ] <-
      replace_na(
        feeding_site_year$feeding_any[r],
        0L
      )
  }
}


for(r in seq_len(nrow(matthew_site_year))) {

  jj <- match(
    matthew_site_year$Site[r],
    ordered_islands
  )

  tt <- matthew_site_year$year_idx[r]

  if(
    tt >=
      1 &
      tt <=
        tMax
  ) {

    matthew_post_matrix[
      jj,
      tt
    ] <-
      matthew_site_year$matthew_post2016[r]

    matthew_acute_matrix[
      jj,
      tt
    ] <-
      matthew_site_year$matthew_acute_2016_2018[r]
  }
}


## ============================================================
## 11. SAVE AUDITS + FINAL CALENDAR EXPOSURES
## ============================================================

write_csv(
  site_year_date_audit,
  file.path(
    out_dir,
    "01_SiteYear_CaptureDateSpanAudit.csv"
  )
)


write_csv(
  date_gaps,
  file.path(
    out_dir,
    "02_All_WithinSiteYear_DateGaps.csv"
  )
)


write_csv(
  threshold_comparison,
  file.path(
    out_dir,
    "03_CandidateVisitThreshold_Comparison.csv"
  )
)


write_csv(
  flagged_site_years,
  file.path(
    out_dir,
    "04_Flagged_SiteYears_LongSpan_or_LargeGap.csv"
  )
)


write_csv(
  feeding_history,
  file.path(
    out_dir,
    "05_Final_Coauthor_FeedingHistory.csv"
  )
)


write_csv(
  feeding_site_year,
  file.path(
    out_dir,
    "06_Final_Feeding_SiteYear.csv"
  )
)


write_csv(
  feeding_survey_support,
  file.path(
    out_dir,
    "07_Feeding_SurveyYearSupport.csv"
  )
)


write_csv(
  feeding_capture_support,
  file.path(
    out_dir,
    "08_Feeding_CaptureSupport.csv"
  )
)


write_csv(
  feeding_transition_support,
  file.path(
    out_dir,
    "09_Feeding_WithinSiteTransitionSupport.csv"
  )
)


write_csv(
  raw_feeding_audit,
  file.path(
    out_dir,
    "10_Raw_FeedingLabel_Audit.csv"
  )
)


write_csv(
  matthew_support,
  file.path(
    out_dir,
    "11_FFRC_Matthew_Support.csv"
  )
)


saveRDS(
  list(
    candidate_thresholds =
      candidate_thresholds,

    threshold_comparison =
      threshold_comparison,

    site_year_date_audit =
      site_year_date_audit,

    flagged_site_years =
      flagged_site_years,

    candidate_visits =
      candidate_results,

    feeding_history =
      feeding_history,

    feeding_site_year =
      feeding_site_year,

    feeding_high_matrix =
      feeding_high_matrix,

    feeding_any_matrix =
      feeding_any_matrix,

    matthew_site_year =
      matthew_site_year,

    matthew_post_matrix =
      matthew_post_matrix,

    matthew_acute_matrix =
      matthew_acute_matrix
  ),
  "Data/Iguana_SurveyOccasion_and_CalendarExposure_Audit.rds"
)


## ============================================================
## 12. CONSOLE SUMMARY
## ============================================================

cat("\n")
cat("============================================================\n")
cat("SURVEY-OCCASION + CALENDAR-EXPOSURE AUDIT\n")
cat("============================================================\n")


cat("\nSITE-YEAR CAPTURE-DATE SPANS\n")
cat("------------------------------------------------------------\n")

print(
  summary(
    site_year_date_audit %>%
      select(
        n_capture_days,
        span_days,
        max_gap_days
      )
  )
)


cat("\nLARGEST WITHIN-SITE-YEAR GAPS\n")
cat("------------------------------------------------------------\n")

print(
  date_gaps %>%
    select(
      Island,
      CalYear,
      previous_date,
      CapDate,
      gap_days
    ) %>%
    slice_head(
      n =
        25
    ),
  n = 25
)


cat("\nCANDIDATE FIELD-VISIT THRESHOLDS\n")
cat("------------------------------------------------------------\n")

print(
  threshold_comparison,
  n = Inf
)


cat("\nFLAGGED SITE-YEARS: SPAN >21 DAYS OR MAX GAP >14 DAYS\n")
cat("------------------------------------------------------------\n")

print(
  flagged_site_years,
  n = Inf
)


cat("\nFINAL FEEDING HISTORY THROUGH 2019\n")
cat("------------------------------------------------------------\n")

print(
  feeding_history,
  n = Inf
)


cat("\nFEEDING SUPPORT IN OBSERVED SURVEY YEARS\n")
cat("------------------------------------------------------------\n")

print(
  feeding_survey_support,
  n = Inf
)


cat("\nWITHIN-SITE FEEDING TRANSITION SUPPORT\n")
cat("------------------------------------------------------------\n")

print(
  feeding_transition_support,
  n = Inf
)


cat("\nFFRC / HURRICANE MATTHEW SUPPORT\n")
cat("------------------------------------------------------------\n")

print(
  matthew_support,
  n = Inf
)


cat("\nMODELLING DECISIONS CARRIED FORWARD\n")
cat("------------------------------------------------------------\n")
cat("1. Feeding is site x calendar-year, never individual treatment.\n")
cat("2. Primary feeding contrast: high/regular vs all other states.\n")
cat("3. Sensitivity: any provisioning vs none.\n")
cat("4. LeafNorth and UCay remain high/regular through 2019.\n")
cat("5. Matthew primary: FFRC 2017 onward; sensitivity 2016-2018.\n")
cat("6. Large-male removal is omitted entirely.\n")
cat("7. Survey-occasion definition is NOT chosen automatically here.\n")
cat("   Choose it from the date-gap/threshold audit before long MCMC.\n")


cat("\nSaved:\n")
cat("  Data/Iguana_SurveyOccasion_and_CalendarExposure_Audit.rds\n")
cat("  DataAudits/SurveyOccasions_Feeding/*.csv\n")
cat("============================================================\n")
