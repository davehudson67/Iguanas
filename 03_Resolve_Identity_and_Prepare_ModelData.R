## ============================================================
## 03_Resolve_Identity_and_Prepare_ModelData.R
##
## Conservative identity-cleaning pass for the canonical dataset.
##
## IMPORTANT:
##   * Northern groomed ID is trusted.
##   * Southern Ig# is trusted EXCEPT where the data prove it was reused.
##   * PIT is NOT used to merge otherwise distinct named animals.
##   * Placeholder IDs are split by valid PIT where possible.
##   * The four Leaf South animals translocated to Pasture Cay are
##     kept as their original histories; later reuse of Ig# 19-22
##     is split into new individuals.
##   * Any remaining sex conflict is left unresolved (sex_clean = NA).
##   * Known ages are converted to a birth-year support interval;
##     individuals whose known-age records disagree by >2 yr are
##     flagged and NOT used as known-age anchors yet.
## ============================================================


library(tidyverse)
library(lubridate)


## ============================================================
## 1. LOAD
## ============================================================

infile <- "Data/Iguana_Canonical_CaptureData_Through2019.rds"

if(!file.exists(infile)) {
  stop(
    "Missing ",
    infile,
    ". Run the canonical importer first."
  )
}


dat <- readRDS(
  infile
)


out_dir <- "DataAudits/IdentityResolution"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 2. HELPERS
## ============================================================

clean_pit <- function(x) {

  z <- as.character(x) %>%
    str_squish()

  z[
    is.na(z) |
      z %in% c(
        "",
        "NA",
        "N/A"
      )
  ] <- NA_character_

  z
}


is_valid_pit <- function(x) {

  z <- clean_pit(x)

  !is.na(z) &
    str_detect(
      z,
      "^[0-9A-Fa-f]{6,}$"
    )
}


is_placeholder_id <- function(x) {

  z <- str_to_lower(
    str_squish(
      replace_na(
        as.character(x),
        ""
      )
    )
  )

  z %in% c(
    "",
    "?",
    "dead",
    "notes",
    "note",
    "see notes",
    "see note",
    "no #",
    "no#",
    "no number",
    "unknown",
    "unk"
  ) |
    str_detect(
      z,
      "^either .*check"
    )
}


## ============================================================
## 3. IDENTIFY THE PROVEN LEAF SOUTH Ig# REUSE
##
## Collaborators confirmed that Ig# 19-22 were permanently removed
## to Pasture Cay in 2006. Therefore any later row carrying the same
## Ig# cannot be the same biological individual.
## ============================================================

removal_dates <- dat %>%
  filter(
    permanent_removal
  ) %>%
  group_by(
    individual_id
  ) %>%
  summarise(
    removal_date =
      min(
        capture_date,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


dat2 <- dat %>%
  left_join(
    removal_dates,
    by =
      "individual_id"
  ) %>%
  mutate(
    pit_tag =
      clean_pit(
        pit_tag
      ),

    valid_pit =
      is_valid_pit(
        pit_tag
      ),

    placeholder_id =
      is_placeholder_id(
        animal_id
      ),

    proven_post_removal_id_reuse =
      Site ==
        "LeafSouth" &
      !is.na(
        removal_date
      ) &
      !is.na(
        capture_date
      ) &
      capture_date >
        removal_date,

    identity_resolution =
      case_when(
        proven_post_removal_id_reuse &
          valid_pit ~
          "split_post_translocation_Ig_reuse_by_PIT",

        proven_post_removal_id_reuse &
          !valid_pit ~
          "split_post_translocation_Ig_reuse_no_valid_PIT",

        placeholder_id &
          valid_pit ~
          "placeholder_ID_resolved_by_PIT",

        placeholder_id &
          !valid_pit ~
          "placeholder_ID_unresolved",

        TRUE ~
          "trusted_source_ID"
      ),

    individual_id_clean =
      case_when(
        proven_post_removal_id_reuse &
          valid_pit ~
          paste0(
            Site,
            "_PIT_",
            pit_tag
          ),

        proven_post_removal_id_reuse &
          !valid_pit ~
          paste0(
            Site,
            "_UNRESOLVED_ROW_",
            source_sheet,
            "_",
            source_row
          ),

        placeholder_id &
          valid_pit ~
          paste0(
            Site,
            "_PIT_",
            pit_tag
          ),

        placeholder_id &
          !valid_pit ~
          NA_character_,

        TRUE ~
          individual_id
      )
  )


## ============================================================
## 4. RECORDS WITHOUT A USABLE INDIVIDUAL ID
##
## Example: unmarked carcasses with ID = DEAD and no PIT.
## Preserve them in an audit, but do not use them as CMR histories.
## ============================================================

unresolved_identity_rows <- dat2 %>%
  filter(
    is.na(
      individual_id_clean
    )
  )


model_dat <- dat2 %>%
  filter(
    !is.na(
      individual_id_clean
    )
  )


## ============================================================
## 5. COLLAPSE SAME-INDIVIDUAL / SAME-DATE DUPLICATES
##
## For CMR there is only one encounter on a given date.
## Keep the most complete row and retain a record of all source rows.
## ============================================================

model_dat <- model_dat %>%
  mutate(
    completeness_score =
      rowSums(
        cbind(
          !is.na(
            pit_tag
          ),
          !is.na(
            sex
          ),
          !is.na(
            svl_cm
          ),
          !is.na(
            body_mass_g
          ),
          !is.na(
            age_recorded
          ),
          !is.na(
            comments
          )
        )
      )
  )


duplicate_groups_before <- model_dat %>%
  filter(
    !is.na(
      capture_date
    )
  ) %>%
  count(
    Site,
    individual_id_clean,
    capture_date,
    name =
      "n_rows"
  ) %>%
  filter(
    n_rows >
      1
  )


duplicate_detail_before <- model_dat %>%
  semi_join(
    duplicate_groups_before,
    by = c(
      "Site",
      "individual_id_clean",
      "capture_date"
    )
  ) %>%
  arrange(
    Site,
    individual_id_clean,
    capture_date,
    desc(
      completeness_score
    ),
    source_row
  )


same_day_provenance <- duplicate_detail_before %>%
  group_by(
    Site,
    individual_id_clean,
    capture_date
  ) %>%
  summarise(
    duplicate_source_rows =
      paste(
        paste0(
          source_sheet,
          ":",
          source_row
        ),
        collapse = " | "
      ),

    duplicate_comments =
      paste(
        unique(
          na.omit(
            comments
          )
        ),
        collapse = " || "
      ),

    .groups =
      "drop"
  )


model_dat <- model_dat %>%
  arrange(
    Site,
    individual_id_clean,
    capture_date,
    desc(
      completeness_score
    ),
    source_row
  ) %>%
  group_by(
    Site,
    individual_id_clean,
    capture_date
  ) %>%
  slice_head(
    n =
      1
  ) %>%
  ungroup() %>%
  left_join(
    same_day_provenance,
    by = c(
      "Site",
      "individual_id_clean",
      "capture_date"
    )
  )


## ============================================================
## 6. RESOLVE SEX CONSERVATIVELY
##
## If an individual is consistently sexed, use that sex.
## If both F and M occur, leave sex_clean = NA rather than guessing.
## A majority candidate is retained for later sensitivity checks.
## ============================================================

sex_resolution <- model_dat %>%
  filter(
    !is.na(
      sex
    )
  ) %>%
  group_by(
    Site,
    individual_id_clean
  ) %>%
  summarise(
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

    n_sexed =
      n_F +
      n_M,

    n_sex_values =
      n_distinct(
        sex
      ),

    .groups =
      "drop"
  ) %>%
  mutate(
    sex_conflict =
      n_sex_values >
        1,

    sex_clean =
      case_when(
        n_F >
          0 &
          n_M ==
            0 ~
          "F",

        n_M >
          0 &
          n_F ==
            0 ~
          "M",

        TRUE ~
          NA_character_
      ),

    sex_majority_candidate =
      case_when(
        n_F >
          n_M ~
          "F",

        n_M >
          n_F ~
          "M",

        TRUE ~
          NA_character_
      ),

    sex_majority_fraction =
      pmax(
        n_F,
        n_M
      ) /
      n_sexed
  )


model_dat <- model_dat %>%
  left_join(
    sex_resolution %>%
      select(
        individual_id_clean,
        sex_conflict,
        sex_clean,
        sex_majority_candidate,
        sex_majority_fraction
      ),
    by =
      "individual_id_clean"
  )


## ============================================================
## 7. KNOWN-AGE -> BIRTH-YEAR SUPPORT
##
## The field team noted that these known ages can differ by about
## +/-1 year and at most about 2 years.
##
## Rather than treating every recorded age as exact, derive the
## individual's implied birth-year support from ALL known-age records.
##
## If the spread exceeds 2 years, do not use as a known-age anchor
## until inspected.
## ============================================================

known_birth_support <- model_dat %>%
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
  group_by(
    Site,
    individual_id_clean
  ) %>%
  summarise(
    n_known_age_records =
      n(),

    implied_birth_min =
      min(
        implied_birth_year,
        na.rm = TRUE
      ),

    implied_birth_max =
      max(
        implied_birth_year,
        na.rm = TRUE
      ),

    implied_birth_median =
      median(
        implied_birth_year,
        na.rm = TRUE
      ),

    implied_birth_range =
      implied_birth_max -
      implied_birth_min,

    .groups =
      "drop"
  ) %>%
  mutate(
    known_age_consistent =
      implied_birth_range <=
        2,

    known_birth_low =
      if_else(
        known_age_consistent,
        floor(
          implied_birth_min
        ),
        NA_real_
      ),

    known_birth_high =
      if_else(
        known_age_consistent,
        ceiling(
          implied_birth_max
        ),
        NA_real_
      ),

    known_birth_mid =
      if_else(
        known_age_consistent,
        implied_birth_median,
        NA_real_
      )
  )


model_dat <- model_dat %>%
  left_join(
    known_birth_support %>%
      select(
        individual_id_clean,
        known_age_consistent,
        known_birth_low,
        known_birth_high,
        known_birth_mid,
        implied_birth_range
      ),
    by =
      "individual_id_clean"
  )


## ============================================================
## 8. PERMANENT REMOVAL / RIGHT-CENSOR FLAG
## ============================================================

model_dat <- model_dat %>%
  mutate(
    permanent_removal_right_censor =
      permanent_removal
  )


## ============================================================
## 9. CHECK SAME PIT ACROSS CLEAN INDIVIDUALS
##
## This is an audit only. We DO NOT merge named animals merely
## because PIT values collide.
## ============================================================

same_pit_clean_ids <- model_dat %>%
  filter(
    valid_pit
  ) %>%
  distinct(
    Site,
    individual_id_clean,
    pit_tag,
    calendar_year
  ) %>%
  group_by(
    pit_tag
  ) %>%
  summarise(
    n_clean_individuals =
      n_distinct(
        individual_id_clean
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

    clean_individuals =
      paste(
        sort(
          unique(
            individual_id_clean
          )
        ),
        collapse = " | "
      ),

    .groups =
      "drop"
  ) %>%
  filter(
    n_clean_individuals >
      1
  ) %>%
  arrange(
    desc(
      n_sites
    ),
    pit_tag
  )


## ============================================================
## 10. AUDITS AFTER IDENTITY RESOLUTION
## ============================================================

identity_changes <- dat2 %>%
  filter(
    identity_resolution !=
      "trusted_source_ID"
  ) %>%
  select(
    Site,
    animal_id,
    individual_id,
    individual_id_clean,
    pit_tag,
    capture_date,
    calendar_year,
    identity_resolution,
    permanent_removal,
    removal_date,
    comments,
    source_sheet,
    source_row
  ) %>%
  arrange(
    Site,
    animal_id,
    capture_date
  )


sex_conflicts_after <- sex_resolution %>%
  filter(
    sex_conflict
  ) %>%
  arrange(
    Site,
    individual_id_clean
  )


known_age_conflicts <- known_birth_support %>%
  filter(
    !known_age_consistent
  ) %>%
  arrange(
    desc(
      implied_birth_range
    )
  )


known_age_site_summary <- known_birth_support %>%
  group_by(
    Site
  ) %>%
  summarise(
    n_known_age_individuals =
      n(),

    n_consistent_le_2yr =
      sum(
        known_age_consistent
      ),

    n_conflict_gt_2yr =
      sum(
        !known_age_consistent
      ),

    .groups =
      "drop"
  )


## Confirm there are no later records under the SAME CLEAN ID
## after a permanent-removal event.
clean_removal_dates <- model_dat %>%
  filter(
    permanent_removal
  ) %>%
  select(
    individual_id_clean,
    removal_date_clean =
      capture_date
  )


post_removal_same_clean_id <- clean_removal_dates %>%
  left_join(
    model_dat %>%
      select(
        individual_id_clean,
        later_date =
          capture_date,
        later_site =
          Site,
        later_pit =
          pit_tag
      ),
    by =
      "individual_id_clean"
  ) %>%
  filter(
    !is.na(
      later_date
    ),
    later_date >
      removal_date_clean
  )


## ============================================================
## 11. SAVE MODEL-PREPARATION DATASET
## ============================================================

model_dat <- model_dat %>%
  select(
    -completeness_score
  ) %>%
  arrange(
    Site,
    individual_id_clean,
    capture_date,
    calendar_year
  )


write_csv(
  model_dat,
  "Data/Iguana_Canonical_PreDataPrep_Through2019.csv",
  na =
    ""
)

saveRDS(
  model_dat,
  "Data/Iguana_Canonical_PreDataPrep_Through2019.rds"
)


write_csv(
  identity_changes,
  file.path(
    out_dir,
    "01_IdentityChanges.csv"
  )
)

write_csv(
  unresolved_identity_rows,
  file.path(
    out_dir,
    "02_UnresolvedIdentityRows_ExcludedFromCMR.csv"
  )
)

write_csv(
  duplicate_detail_before,
  file.path(
    out_dir,
    "03_SameDayDuplicates_BeforeCollapse.csv"
  )
)

write_csv(
  sex_conflicts_after,
  file.path(
    out_dir,
    "04_SexConflicts_AfterIdentityResolution.csv"
  )
)

write_csv(
  known_birth_support,
  file.path(
    out_dir,
    "05_KnownBirthSupport_All.csv"
  )
)

write_csv(
  known_age_conflicts,
  file.path(
    out_dir,
    "05_KnownAgeConflicts_GT2yr.csv"
  )
)

write_csv(
  known_age_site_summary,
  file.path(
    out_dir,
    "05_KnownAgeSupport_BySite.csv"
  )
)

write_csv(
  same_pit_clean_ids,
  file.path(
    out_dir,
    "06_SamePIT_MultipleCleanIndividuals.csv"
  )
)

write_csv(
  post_removal_same_clean_id,
  file.path(
    out_dir,
    "07_PostRemovalCaptures_SameCleanID.csv"
  )
)


## ============================================================
## 12. CONSOLE SUMMARY
## ============================================================

cat("\n")
cat("============================================================\n")
cat("IDENTITY RESOLUTION COMPLETE\n")
cat("============================================================\n")

cat(
  "Input rows:",
  nrow(
    dat
  ),
  "\n"
)

cat(
  "Rows excluded from CMR because identity unresolved:",
  nrow(
    unresolved_identity_rows
  ),
  "\n"
)

cat(
  "Same-day duplicate groups collapsed:",
  nrow(
    duplicate_groups_before
  ),
  "\n"
)

cat(
  "Output capture rows:",
  nrow(
    model_dat
  ),
  "\n"
)

cat(
  "Clean individuals:",
  n_distinct(
    model_dat$individual_id_clean
  ),
  "\n"
)

cat(
  "Identity rows changed/split:",
  nrow(
    identity_changes
  ),
  "\n"
)

cat(
  "Remaining individuals with sex conflicts:",
  nrow(
    sex_conflicts_after
  ),
  "\n"
)

cat(
  "Known-age individuals with <=2 yr implied birth spread:",
  sum(
    known_birth_support$known_age_consistent
  ),
  "\n"
)

cat(
  "Known-age individuals with >2 yr conflict:",
  sum(
    !known_birth_support$known_age_consistent
  ),
  "\n"
)

cat(
  "Same PIT still attached to >1 clean individual:",
  nrow(
    same_pit_clean_ids
  ),
  "\n"
)

cat(
  "Captures after permanent removal under same CLEAN identity:",
  nrow(
    post_removal_same_clean_id
  ),
  "\n"
)

cat("\nKNOWN-AGE SUPPORT BY SITE\n")
cat("------------------------------------------------------------\n")

print(
  known_age_site_summary,
  n = Inf
)

cat("\nREMAINING SEX CONFLICTS\n")
cat("------------------------------------------------------------\n")

print(
  sex_conflicts_after,
  n = Inf
)

cat("\nKNOWN-AGE CONFLICTS >2 YEARS\n")
cat("------------------------------------------------------------\n")

print(
  known_age_conflicts,
  n = Inf
)

cat("\n============================================================\n")
cat("Saved:\n")
cat("  Data/Iguana_Canonical_PreDataPrep_Through2019.csv\n")
cat("  Data/Iguana_Canonical_PreDataPrep_Through2019.rds\n")
cat("  DataAudits/IdentityResolution/*.csv\n")
cat("============================================================\n")
