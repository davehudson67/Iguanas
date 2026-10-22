## ============================================================
## 02_Audit_Canonical_Cleaning_Issues.R
##
## PURPOSE
## -------
## Focused audit of the canonical combined dataset before DataPrep.R.
##
## This script DOES NOT automatically alter biological records.
## It classifies:
##   1. duplicate animal-date records
##   2. conflicting sex labels
##   3. multiple PITs within one reconciled individual
##   4. the same PIT appearing under multiple individuals
##   5. permanent-removal histories
##   6. known-age consistency
##   7. raw feeding labels by site/year
##
## INPUT
## -----
## Data/Iguana_Canonical_CaptureData_Through2019.rds
##
## OUTPUT
## ------
## DataAudits/CleaningStage2/
## ============================================================


library(tidyverse)
library(lubridate)


## ============================================================
## 1. LOAD
## ============================================================

infile <- "Data/Iguana_Canonical_CaptureData_Through2019.rds"

if(!file.exists(infile)) {
  stop(
    "Could not find: ",
    infile,
    "\nRun 01_Combine_Clean_Raw_Iguana_Data_V2.R first."
  )
}


dat <- readRDS(
  infile
)


out_dir <- "DataAudits/CleaningStage2"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 2. DUPLICATE ANIMAL-DATE RECORDS
##
## Classify whether same-day duplicates are exact/near-exact or
## contain genuinely different measurements/comments.
## ============================================================

dup_keys <- dat %>%
  filter(
    !is.na(
      capture_date
    )
  ) %>%
  count(
    Site,
    individual_id,
    capture_date,
    name = "n_rows"
  ) %>%
  filter(
    n_rows > 1
  )


duplicate_detail <- dat %>%
  semi_join(
    dup_keys,
    by = c(
      "Site",
      "individual_id",
      "capture_date"
    )
  ) %>%
  arrange(
    Site,
    individual_id,
    capture_date,
    source_sheet,
    source_row
  )


duplicate_summary <- duplicate_detail %>%
  group_by(
    Site,
    individual_id,
    capture_date
  ) %>%
  summarise(
    n_rows = n(),

    n_distinct_svl =
      n_distinct(
        svl_cm,
        na.rm = TRUE
      ),

    n_distinct_mass =
      n_distinct(
        body_mass_g,
        na.rm = TRUE
      ),

    n_distinct_sex =
      n_distinct(
        sex,
        na.rm = TRUE
      ),

    n_distinct_pit =
      n_distinct(
        pit_tag,
        na.rm = TRUE
      ),

    same_source_sheet =
      n_distinct(
        source_sheet
      ) ==
        1,

    source_rows =
      paste(
        source_row,
        collapse = ", "
      ),

    comments =
      paste(
        unique(
          na.omit(
            comments
          )
        ),
        collapse = " || "
      ),

    .groups = "drop"
  ) %>%
  mutate(
    likely_exact_duplicate =
      n_distinct_svl <= 1 &
      n_distinct_mass <= 1 &
      n_distinct_sex <= 1 &
      n_distinct_pit <= 1
  )


## ============================================================
## 3. SEX CONFLICTS
##
## Do not simply majority-vote automatically.
##
## Summarise:
##   * counts of F/M records
##   * first and last recorded sex
##   * whether conflict is a single discordant record
##   * whether one sex dominates >=80% of sexed captures
## ============================================================

sex_detail <- dat %>%
  filter(
    !is.na(
      sex
    )
  ) %>%
  group_by(
    Site,
    individual_id
  ) %>%
  filter(
    n_distinct(
      sex
    ) >
      1
  ) %>%
  ungroup() %>%
  arrange(
    Site,
    individual_id,
    capture_date,
    calendar_year
  )


sex_conflict_summary <- sex_detail %>%
  group_by(
    Site,
    individual_id
  ) %>%
  summarise(
    n_sexed_records =
      n(),

    n_F =
      sum(
        sex ==
          "F"
      ),

    n_M =
      sum(
        sex ==
          "M"
      ),

    first_year =
      min(
        calendar_year
      ),

    last_year =
      max(
        calendar_year
      ),

    first_sex =
      sex[
        order(
          coalesce(
            capture_date,
            as.Date(
              paste0(
                calendar_year,
                "-06-15"
              )
            )
          )
        )
      ][1],

    last_sex =
      sex[
        order(
          coalesce(
            capture_date,
            as.Date(
              paste0(
                calendar_year,
                "-06-15"
              )
            )
          )
        )
      ][
        n()
      ],

    F_years =
      paste(
        sort(
          unique(
            calendar_year[
              sex ==
                "F"
            ]
          )
        ),
        collapse = ","
      ),

    M_years =
      paste(
        sort(
          unique(
            calendar_year[
              sex ==
                "M"
            ]
          )
        ),
        collapse = ","
      ),

    .groups = "drop"
  ) %>%
  mutate(
    dominant_sex =
      case_when(
        n_F >
          n_M ~
          "F",

        n_M >
          n_F ~
          "M",

        TRUE ~
          "tie"
      ),

    dominant_fraction =
      pmax(
        n_F,
        n_M
      ) /
      n_sexed_records,

    singleton_discordance =
      pmin(
        n_F,
        n_M
      ) ==
        1,

    strong_majority =
      dominant_fraction >=
        0.8
  ) %>%
  arrange(
    Site,
    individual_id
  )


## ============================================================
## 4. MULTIPLE PIT TAGS WITHIN INDIVIDUAL
##
## Expected in some animals because collaborators reconciled IDs.
## ============================================================

multi_pit_detail <- dat %>%
  filter(
    !is.na(
      pit_tag
    )
  ) %>%
  group_by(
    Site,
    individual_id
  ) %>%
  filter(
    n_distinct(
      pit_tag
    ) >
      1
  ) %>%
  ungroup() %>%
  arrange(
    Site,
    individual_id,
    capture_date,
    calendar_year
  )


multi_pit_summary <- multi_pit_detail %>%
  group_by(
    Site,
    individual_id
  ) %>%
  summarise(
    n_records =
      n(),

    n_pit_tags =
      n_distinct(
        pit_tag
      ),

    pit_tags =
      paste(
        sort(
          unique(
            pit_tag
          )
        ),
        collapse = " | "
      ),

    first_year =
      min(
        calendar_year
      ),

    last_year =
      max(
        calendar_year
      ),

    comments =
      paste(
        unique(
          na.omit(
            comments
          )
        ),
        collapse = " || "
      ),

    .groups = "drop"
  ) %>%
  mutate(
    comment_mentions_tag_issue =
      str_detect(
        str_to_lower(
          comments
        ),
        "pit|tag|lost|replace|re[- ]?tag|fail"
      )
  )


## ============================================================
## 5. REVERSE PIT CHECK
##
## More important than multiple PITs within one reconciled animal:
## does the SAME PIT appear under >1 stable individual?
## ============================================================

pit_to_multiple_ids <- dat %>%
  filter(
    !is.na(
      pit_tag
    )
  ) %>%
  distinct(
    Site,
    individual_id,
    pit_tag,
    calendar_year
  ) %>%
  group_by(
    pit_tag
  ) %>%
  summarise(
    n_individuals =
      n_distinct(
        individual_id
      ),

    n_sites =
      n_distinct(
        Site
      ),

    sites =
      paste(
        sort(
          unique(
            Site
          )
        ),
        collapse = " | "
      ),

    individual_ids =
      paste(
        sort(
          unique(
            individual_id
          )
        ),
        collapse = " | "
      ),

    first_year =
      min(
        calendar_year
      ),

    last_year =
      max(
        calendar_year
      ),

    .groups = "drop"
  ) %>%
  filter(
    n_individuals >
      1
  ) %>%
  arrange(
    desc(
      n_sites
    ),
    pit_tag
  )


pit_to_multiple_ids_detail <- dat %>%
  semi_join(
    pit_to_multiple_ids,
    by = "pit_tag"
  ) %>%
  arrange(
    pit_tag,
    Site,
    individual_id,
    capture_date,
    calendar_year
  )


## ============================================================
## 6. PERMANENT REMOVAL CHECK
##
## Confirm:
##   * exactly four removals
##   * no later capture on LeafSouth after the removal date
## ============================================================

removals <- dat %>%
  filter(
    permanent_removal
  ) %>%
  select(
    Site,
    individual_id,
    animal_id,
    pit_tag,
    capture_date,
    calendar_year,
    removal_destination,
    comments,
    source_sheet,
    source_row
  )


removal_followup <- removals %>%
  select(
    individual_id,
    removal_date =
      capture_date
  ) %>%
  left_join(
    dat %>%
      select(
        individual_id,
        later_date =
          capture_date,
        later_site =
          Site,
        later_comments =
          comments
      ),
    by =
      "individual_id"
  ) %>%
  filter(
    !is.na(
      later_date
    ),
    later_date >
      removal_date
  ) %>%
  arrange(
    individual_id,
    later_date
  )


## ============================================================
## 7. KNOWN-AGE CONSISTENCY
##
## Susannah noted ages may be close rather than perfectly exact.
## We therefore report implied birth-year spread rather than
## enforcing exact equality.
## ============================================================

known_age_detail <- dat %>%
  filter(
    !is.na(
      age_known_for_model
    )
  ) %>%
  mutate(
    implied_birth_year =
      calendar_year -
      age_known_for_model
  ) %>%
  arrange(
    Site,
    individual_id,
    calendar_year,
    capture_date
  )


known_age_summary <- known_age_detail %>%
  group_by(
    Site,
    individual_id
  ) %>%
  summarise(
    n_known_records =
      n(),

    birth_year_min =
      min(
        implied_birth_year
      ),

    birth_year_max =
      max(
        implied_birth_year
      ),

    birth_year_range =
      birth_year_max -
      birth_year_min,

    median_implied_birth_year =
      median(
        implied_birth_year
      ),

    .groups = "drop"
  ) %>%
  arrange(
    desc(
      birth_year_range
    )
  )


known_age_site_summary <- known_age_summary %>%
  group_by(
    Site
  ) %>%
  summarise(
    n_known_age_individuals =
      n(),

    n_range_0 =
      sum(
        birth_year_range ==
          0
      ),

    n_range_le_1 =
      sum(
        birth_year_range <=
          1
      ),

    n_range_le_2 =
      sum(
        birth_year_range <=
          2
      ),

    n_range_gt_2 =
      sum(
        birth_year_range >
          2
      ),

    median_birth_year_range =
      median(
        birth_year_range
      ),

    max_birth_year_range =
      max(
        birth_year_range
      ),

    .groups = "drop"
  )


## ============================================================
## 8. RAW FEEDING LABELS BY SITE/YEAR
##
## This remains an audit only.
## ============================================================

feeding_site_year <- dat %>%
  count(
    Site,
    calendar_year,
    feeding_raw,
    name =
      "n_records"
  ) %>%
  arrange(
    Site,
    calendar_year,
    feeding_raw
  )


feeding_transition_summary <- feeding_site_year %>%
  group_by(
    Site
  ) %>%
  summarise(
    first_year =
      min(
        calendar_year
      ),

    last_year =
      max(
        calendar_year
      ),

    raw_states =
      paste(
        unique(
          na.omit(
            feeding_raw
          )
        ),
        collapse = " -> "
      ),

    n_distinct_raw_states =
      n_distinct(
        feeding_raw[
          !is.na(
            feeding_raw
          )
        ]
      ),

    .groups = "drop"
  )


## ============================================================
## 9. WRITE AUDITS
## ============================================================

write_csv(
  duplicate_summary,
  file.path(
    out_dir,
    "01_DuplicateCaptureSummary.csv"
  )
)

write_csv(
  duplicate_detail,
  file.path(
    out_dir,
    "01_DuplicateCaptureDetail.csv"
  )
)

write_csv(
  sex_conflict_summary,
  file.path(
    out_dir,
    "02_SexConflictSummary.csv"
  )
)

write_csv(
  sex_detail,
  file.path(
    out_dir,
    "02_SexConflictDetail.csv"
  )
)

write_csv(
  multi_pit_summary,
  file.path(
    out_dir,
    "03_MultiplePITsWithinIndividual.csv"
  )
)

write_csv(
  multi_pit_detail,
  file.path(
    out_dir,
    "03_MultiplePITsWithinIndividual_Detail.csv"
  )
)

write_csv(
  pit_to_multiple_ids,
  file.path(
    out_dir,
    "04_SamePIT_MultipleIndividuals.csv"
  )
)

write_csv(
  pit_to_multiple_ids_detail,
  file.path(
    out_dir,
    "04_SamePIT_MultipleIndividuals_Detail.csv"
  )
)

write_csv(
  removals,
  file.path(
    out_dir,
    "05_PermanentRemovals.csv"
  )
)

write_csv(
  removal_followup,
  file.path(
    out_dir,
    "05_PostRemovalCaptures.csv"
  )
)

write_csv(
  known_age_summary,
  file.path(
    out_dir,
    "06_KnownAgeConsistency_Individuals.csv"
  )
)

write_csv(
  known_age_site_summary,
  file.path(
    out_dir,
    "06_KnownAgeConsistency_Sites.csv"
  )
)

write_csv(
  feeding_site_year,
  file.path(
    out_dir,
    "07_RawFeedingBySiteYear.csv"
  )
)

write_csv(
  feeding_transition_summary,
  file.path(
    out_dir,
    "07_RawFeedingTransitionSummary.csv"
  )
)


## ============================================================
## 10. CONSOLE SUMMARY
## ============================================================

cat("\n")
cat("============================================================\n")
cat("STAGE-2 CLEANING AUDIT\n")
cat("============================================================\n")

cat("\nDUPLICATE ANIMAL-DATE GROUPS\n")
cat("------------------------------------------------------------\n")
print(
  duplicate_summary,
  n = Inf
)

cat("\nSEX CONFLICT SUMMARY\n")
cat("------------------------------------------------------------\n")
print(
  sex_conflict_summary,
  n = Inf
)

cat("\nSEX CONFLICT COUNTS\n")
cat("------------------------------------------------------------\n")
cat(
  "Total conflicting individuals:",
  nrow(
    sex_conflict_summary
  ),
  "\n"
)
cat(
  "Singleton discordance:",
  sum(
    sex_conflict_summary$singleton_discordance
  ),
  "\n"
)
cat(
  "Strong >=80% majority:",
  sum(
    sex_conflict_summary$strong_majority
  ),
  "\n"
)

cat("\nMULTIPLE PITs WITHIN RECONCILED INDIVIDUAL\n")
cat("------------------------------------------------------------\n")
print(
  multi_pit_summary,
  n = Inf
)

cat("\nSAME PIT ASSIGNED TO MULTIPLE INDIVIDUALS\n")
cat("------------------------------------------------------------\n")
print(
  pit_to_multiple_ids,
  n = Inf
)

cat("\nPERMANENT REMOVAL CHECK\n")
cat("------------------------------------------------------------\n")
print(
  removals,
  n = Inf
)
cat(
  "Captures after permanent removal:",
  nrow(
    removal_followup
  ),
  "\n"
)

cat("\nKNOWN-AGE CONSISTENCY BY SITE\n")
cat("------------------------------------------------------------\n")
print(
  known_age_site_summary,
  n = Inf
)

cat("\nRAW FEEDING LABEL TRANSITIONS\n")
cat("------------------------------------------------------------\n")
print(
  feeding_transition_summary,
  n = Inf
)

cat("\n============================================================\n")
cat("AUDIT COMPLETE\n")
cat("============================================================\n")
cat(
  "Outputs:",
  out_dir,
  "\n"
)
cat("============================================================\n")
