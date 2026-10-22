## ============================================================
## 01_Combine_Clean_Raw_Iguana_Data.R
##
## PURPOSE
## -------
## Combine the two agreed source workbooks into one traceable,
## minimally interpreted capture-level dataset through 2019.
##
## This is deliberately a DATA-CLEANING step only.
## It does NOT yet create survival histories, latent ages,
## densities, survey-event likelihoods, or final provisioning
## covariates.
##
## AGREED SOURCES
## ---------------
## Northern populations:
##   AllensIguana.PitTag.Combined Data 082026.xlsx
##
## Central/southern populations:
##   C figginsi morphometric data_MASTER vers 12.12.2019_
##   Chuck_ID with feeding.xlsx
##
## AGREED SITE MAPPING
## -------------------
## Northern:
##   FlatRock2019-2024  -> FFRC       -> inornata
##   LeafCay2019-2026   -> LeafNorth  -> inornata
##   Ucay2019-2026      -> UCay       -> inornata
##
## Southern:
##   Bitter Guana Cay   -> Bitter     -> figginsi
##   Gaulin Cay         -> Gaulin     -> figginsi
##   Noddy Cay          -> Noddy      -> figginsi
##   White Bay          -> WhiteB     -> figginsi
##   North Adderly      -> NAdder     -> figginsi
##   Leaf Cay           -> LeafSouth  -> figginsi
##
## EXCLUDED
## --------
##   Alligator Cay      -> explicitly excluded by collaborators
##   Allen Cay          -> not part of the current analysis
##   Old northern sheets in Chuck's workbook are NOT used because
##   they are superseded by the groomed northern workbook.
##
## KEY RULES
## ---------
## * Restrict all analysis data to <= 2019.
## * Northern ages:
##       k = known/near-known -> can later anchor birth year
##       e = estimated        -> retained, but NOT treated as known
## * Southern Ig# is the reconciled stable individual identifier.
## * Northern ID field is the groomed stable individual identifier.
## * Four Leaf South animals moved to Pasture Cay are FLAGGED as
##   permanent removals. They are NOT converted to deaths here.
## * Feeding is retained exactly as recorded at capture.
##   No final tourism/provisioning exposure is imposed in this file.
## ============================================================


library(tidyverse)
library(readxl)
library(lubridate)


## ============================================================
## 1. FILE LOCATIONS
## ============================================================

find_file <- function(filename) {

  candidates <- c(
    filename,
    file.path("Data", filename)
  )

  hit <- candidates[file.exists(candidates)]

  if(length(hit) == 0) {
    stop(
      "Could not find: ",
      filename,
      "\nPut the workbook in either the repository root or Data/."
    )
  }

  hit[1]
}


north_file <- find_file(
  "AllensIguana.PitTag.Combined Data 082026.xlsx"
)

south_file <- find_file(
  "C figginsi morphometric data_MASTER vers 12.12.2019_Chuck_ID with feeding.xlsx"
)


out_dir <- "DataAudits/CanonicalImport"

dir.create(
  "Data",
  showWarnings = FALSE
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 2. GENERAL HELPERS
## ============================================================

clean_chr <- function(x) {

  x <- as.character(x)

  x <- str_replace_all(
    x,
    "\u00a0",
    " "
  )

  x <- str_squish(x)

  x[
    is.na(x) |
      x %in% c(
        "",
        "NA",
        "N/A",
        "na",
        "n/a",
        "NULL",
        "."
      )
  ] <- NA_character_

  x
}


clean_num <- function(x) {

  suppressWarnings(
    as.numeric(
      clean_chr(x)
    )
  )
}


normalise_sex <- function(x) {

  z <- str_to_lower(
    clean_chr(x)
  )

  case_when(
    str_detect(z, "^f") ~ "F",
    str_detect(z, "^m") ~ "M",
    TRUE ~ NA_character_
  )
}


pick_column <- function(
  d,
  aliases,
  required = FALSE
) {

  nm <- names(d)

  normalised_names <- nm %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "")

  normalised_aliases <- aliases %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "")

  idx <- match(
    normalised_aliases,
    normalised_names
  )

  idx <- idx[
    !is.na(idx)
  ]

  if(length(idx) == 0) {

    if(required) {
      stop(
        "Could not find required column. Tried: ",
        paste(
          aliases,
          collapse = ", "
        )
      )
    }

    return(
      rep(
        NA_character_,
        nrow(d)
      )
    )
  }

  d[[idx[1]]]
}


parse_capture_date <- function(x) {

  raw <- clean_chr(x)

  ## Some historic cells display a 2019 date as 2919.
  ## This may exist either as literal text ("2919") OR as an
  ## Excel serial representing a genuine year-2919 date.
  corrected <- raw %>%
    str_replace_all(
      "2919",
      "2019"
    )

  numeric_date <- suppressWarnings(
    as.numeric(
      corrected
    )
  )

  out <- as.Date(
    rep(
      NA_character_,
      length(corrected)
    )
  )

  ## Excel serial dates. Allow large serials temporarily so that
  ## accidental year-2919 dates can be recovered rather than lost.
  is_serial <- !is.na(numeric_date) &
    numeric_date >= 20000 &
    numeric_date <= 400000

  out[is_serial] <- as.Date(
    numeric_date[is_serial],
    origin = "1899-12-30"
  )

  ## Character dates.
  remaining <- is.na(out) &
    !is.na(corrected)

  parsed <- suppressWarnings(
    parse_date_time(
      corrected[remaining],
      orders = c(
        "ymd",
        "mdy",
        "dmy",
        "Ymd",
        "mdY",
        "dmY",
        "d-b-Y",
        "d-b-y",
        "m/d/Y",
        "d/m/Y",
        "m/d/y",
        "d/m/y"
      ),
      quiet = TRUE
    )
  )

  out[remaining] <- as.Date(
    parsed
  )

  ## Repair the known 2919 -> 2019 year-entry error while retaining
  ## the original month/day.
  bad_2919 <- !is.na(out) &
    year(out) ==
      2919

  if(any(bad_2919)) {

    out[bad_2919] <- make_date(
      year =
        2019,
      month =
        month(
          out[bad_2919]
        ),
      day =
        day(
          out[bad_2919]
        )
    )
  }

  out
}


year_from_date_text <- function(x) {

  raw <- clean_chr(x) %>%
    str_replace_all(
      "2919",
      "2019"
    )

  four_digit <- str_extract(
    raw,
    "(19|20)[0-9]{2}"
  )

  suppressWarnings(
    as.integer(
      four_digit
    )
  )
}


comment_flag <- function(
  x,
  pattern
) {

  str_detect(
    str_to_lower(
      replace_na(
        clean_chr(x),
        ""
      )
    ),
    pattern
  )
}


## ============================================================
## 3. IMPORT NORTHERN POPULATIONS
## ============================================================

read_north <- function(
  sheet,
  site
) {

  d <- read_excel(
    north_file,
    sheet = sheet,
    col_types = "text"
  ) %>%
    as_tibble() %>%
    mutate(
      source_row =
        row_number() +
        1L
    )

  animal_id <- clean_chr(
    pick_column(
      d,
      c(
        "ID",
        "Iguana ID",
        "Animal ID"
      ),
      required = TRUE
    )
  )

  pit <- clean_chr(
    pick_column(
      d,
      c(
        "Pit Tag",
        "Pit Tag Id",
        "PIT",
        "PIT Tag ID"
      )
    )
  )

  year_raw <- clean_num(
    pick_column(
      d,
      c(
        "Year",
        "CalYear",
        "Calendar Year"
      )
    )
  )

  date_raw <- clean_chr(
    pick_column(
      d,
      c(
        "Date",
        "Capture Date"
      )
    )
  )

  capture_date <- parse_capture_date(
    date_raw
  )

  date_year <- year(
    capture_date
  )

  text_year <- year_from_date_text(
    date_raw
  )

  ## Prefer the explicit Year column when it is plausible.
  calendar_year <- case_when(
    !is.na(year_raw) &
      year_raw >= 1980 &
      year_raw <= 2019 ~
      as.integer(
        year_raw
      ),

    !is.na(date_year) &
      date_year >= 1980 &
      date_year <= 2019 ~
      as.integer(
        date_year
      ),

    !is.na(text_year) &
      text_year >= 1980 &
      text_year <= 2019 ~
      text_year,

    TRUE ~
      NA_integer_
  )

  age_recorded <- clean_num(
    pick_column(
      d,
      c(
        "Age"
      )
    )
  )

  age_code <- str_to_lower(
    clean_chr(
      pick_column(
        d,
        c(
          "Age known/estimated (k/e)",
          "Age known estimated",
          "Age k/e"
        )
      )
    )
  )

  comments <- clean_chr(
    pick_column(
      d,
      c(
        "Comments",
        "Comment",
        "Notes"
      )
    )
  )

  feeding_raw <- clean_chr(
    pick_column(
      d,
      c(
        "Fed or not",
        "Fed",
        "feeding",
        "Feeding"
      )
    )
  )

  tibble(
    Site =
      site,

    Species =
      "inornata",

    animal_id =
      animal_id,

    individual_id =
      paste0(
        site,
        "_",
        animal_id
      ),

    pit_tag =
      pit,

    sex =
      normalise_sex(
        pick_column(
          d,
          c(
            "Sex"
          )
        )
      ),

    capture_date =
      capture_date,

    calendar_year =
      calendar_year,

    date_raw =
      date_raw,

    date_corrected_2919 =
      (
        !is.na(
          date_raw
        ) &
        str_detect(
          date_raw,
          "2919"
        )
      ) |
      (
        !is.na(
          suppressWarnings(
            as.numeric(
              date_raw
            )
          )
        ) &
        suppressWarnings(
          as.numeric(
            date_raw
          )
        ) >
          60000 &
        calendar_year ==
          2019
      ),

    svl_cm =
      clean_num(
        pick_column(
          d,
          c(
            "Svl",
            "SVL",
            "Svl (cm)",
            "SVL (cm)"
          )
        )
      ),

    tail_length_cm =
      clean_num(
        pick_column(
          d,
          c(
            "Tail-length",
            "Tail length",
            "Vtl",
            "Vtl (cm)"
          )
        )
      ),

    body_mass_g =
      clean_num(
        pick_column(
          d,
          c(
            "Body-mass",
            "Body mass",
            "BM",
            "BM (g)"
          )
        )
      ),

    age_recorded =
      age_recorded,

    age_source =
      case_when(
        age_code ==
          "k" ~
          "known",

        age_code ==
          "e" ~
          "estimated",

        !is.na(
          age_recorded
        ) ~
          "unspecified",

        TRUE ~
          NA_character_
      ),

    ## Only known northern ages can anchor later birth reconstruction.
    age_known_for_model =
      if_else(
        age_code ==
          "k",
        age_recorded,
        NA_real_
      ),

    feeding_raw =
      feeding_raw,

    comments =
      comments,

    permanent_removal =
      FALSE,

    removal_destination =
      NA_character_,

    source_group =
      "Northern_cleaned",

    source_file =
      basename(
        north_file
      ),

    source_sheet =
      sheet,

    source_row =
      d$source_row
  ) %>%
    filter(
      !is.na(
        animal_id
      ),
      !is.na(
        calendar_year
      ),
      calendar_year <=
        2019
    )
}


north_data <- bind_rows(
  read_north(
    "FlatRock2019-2024",
    "FFRC"
  ),

  read_north(
    "LeafCay2019-2026",
    "LeafNorth"
  ),

  read_north(
    "Ucay2019-2026",
    "UCay"
  )
)


## ============================================================
## 4. IMPORT CENTRAL / SOUTHERN POPULATIONS
## ============================================================

read_south <- function(
  sheet,
  site
) {

  d <- read_excel(
    south_file,
    sheet = sheet,
    col_types = "text"
  ) %>%
    as_tibble() %>%
    mutate(
      source_row =
        row_number() +
        1L
    )

  animal_id <- clean_chr(
    pick_column(
      d,
      c(
        "Ig #",
        "Ig#",
        "Iguana #",
        "Iguana ID",
        "ID"
      ),
      required = TRUE
    )
  )

  pit <- clean_chr(
    pick_column(
      d,
      c(
        "Pit #",
        "PIT #",
        "Pit Tag",
        "PIT Tag"
      )
    )
  )

  date_raw <- clean_chr(
    pick_column(
      d,
      c(
        "Capture Date",
        "Date"
      ),
      required = TRUE
    )
  )

  capture_date <- parse_capture_date(
    date_raw
  )

  calendar_year <- year(
    capture_date
  )

  ## Retain any legacy age field for provenance, but do NOT assume
  ## it is genuinely known age.
  age_recorded <- clean_num(
    pick_column(
      d,
      c(
        "Age",
        "age"
      )
    )
  )

  comments <- clean_chr(
    pick_column(
      d,
      c(
        "Comments",
        "Comment",
        "Notes"
      )
    )
  )

  translocated_pasture <- site ==
    "LeafSouth" &
    comment_flag(
      comments,
      "translocated\\s+to\\s+pasture\\s+cay"
    )

  tibble(
    Site =
      site,

    Species =
      "figginsi",

    animal_id =
      animal_id,

    individual_id =
      paste0(
        site,
        "_",
        animal_id
      ),

    pit_tag =
      pit,

    sex =
      normalise_sex(
        pick_column(
          d,
          c(
            "Sex"
          )
        )
      ),

    capture_date =
      capture_date,

    calendar_year =
      calendar_year,

    date_raw =
      date_raw,

    date_corrected_2919 =
      (
        !is.na(
          date_raw
        ) &
        str_detect(
          date_raw,
          "2919"
        )
      ) |
      (
        !is.na(
          suppressWarnings(
            as.numeric(
              date_raw
            )
          )
        ) &
        suppressWarnings(
          as.numeric(
            date_raw
          )
        ) >
          60000 &
        calendar_year ==
          2019
      ),

    svl_cm =
      clean_num(
        pick_column(
          d,
          c(
            "Svl (cm)",
            "SVL (cm)",
            "SVL",
            "Svl"
          )
        )
      ),

    tail_length_cm =
      clean_num(
        pick_column(
          d,
          c(
            "Vtl (cm)",
            "VTL (cm)",
            "Tail-length",
            "Tail length"
          )
        )
      ),

    body_mass_g =
      clean_num(
        pick_column(
          d,
          c(
            "BM (g)",
            "Body mass",
            "Body-mass",
            "Mass (g)"
          )
        )
      ),

    age_recorded =
      age_recorded,

    age_source =
      if_else(
        !is.na(
          age_recorded
        ),
        "legacy_unverified",
        NA_character_
      ),

    ## Do not use these as known-age anchors until independently
    ## verified.
    age_known_for_model =
      NA_real_,

    feeding_raw =
      clean_chr(
        pick_column(
          d,
          c(
            "feeding",
            "Feeding",
            "Fed or not",
            "Fed"
          )
        )
      ),

    comments =
      comments,

    permanent_removal =
      translocated_pasture,

    removal_destination =
      if_else(
        translocated_pasture,
        "Pasture Cay",
        NA_character_
      ),

    source_group =
      "Southern_original",

    source_file =
      basename(
        south_file
      ),

    source_sheet =
      sheet,

    source_row =
      d$source_row
  ) %>%
    filter(
      !is.na(
        animal_id
      ),
      !is.na(
        calendar_year
      ),
      calendar_year <=
        2019
    )
}


south_data <- bind_rows(
  read_south(
    "Bitter Guana Cay",
    "Bitter"
  ),

  read_south(
    "Gaulin Cay",
    "Gaulin"
  ),

  read_south(
    "Noddy Cay",
    "Noddy"
  ),

  read_south(
    "White Bay",
    "WhiteB"
  ),

  read_south(
    "North Adderly",
    "NAdder"
  ),

  read_south(
    "Leaf Cay",
    "LeafSouth"
  )
)


## ============================================================
## 5. COMBINE
## ============================================================

AllData_clean <- bind_rows(
  north_data,
  south_data
) %>%
  arrange(
    Site,
    animal_id,
    capture_date,
    calendar_year
  )


## ============================================================
## 6. DERIVED CLEANING FLAGS
##
## These flags do NOT delete records automatically.
## They tell us what needs inspection before DataPrep/model fitting.
## ============================================================

AllData_clean <- AllData_clean %>%
  mutate(
    comment_dead =
      comment_flag(
        comments,
        "\\bdead\\b|\\bdied\\b|\\bdeceased\\b|\\bcarcass\\b"
      ),

    comment_trap =
      comment_flag(
        comments,
        "\\btrap|trapped"
      ),

    comment_nesting =
      comment_flag(
        comments,
        "nest|nesting"
      ),

    comment_feeding_beach =
      comment_flag(
        comments,
        "feeding\\s+beach|tourist\\s+beach"
      )
  )


## ============================================================
## 7. AUDITS
## ============================================================

site_summary <- AllData_clean %>%
  group_by(
    Site,
    Species,
    source_group
  ) %>%
  summarise(
    n_records =
      n(),

    n_individuals =
      n_distinct(
        individual_id
      ),

    first_year =
      min(
        calendar_year,
        na.rm = TRUE
      ),

    last_year =
      max(
        calendar_year,
        na.rm = TRUE
      ),

    n_known_age_individuals =
      n_distinct(
        individual_id[
          !is.na(
            age_known_for_model
          )
        ]
      ),

    n_permanent_removals =
      sum(
        permanent_removal,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) %>%
  arrange(
    Species,
    Site
  )


## Same animal/date appearing more than once.
possible_duplicate_capture <- AllData_clean %>%
  filter(
    !is.na(
      capture_date
    )
  ) %>%
  group_by(
    Site,
    individual_id,
    capture_date
  ) %>%
  filter(
    n() >
      1
  ) %>%
  ungroup() %>%
  arrange(
    Site,
    individual_id,
    capture_date
  )


## Conflicting sex labels for the same stable individual.
sex_conflicts <- AllData_clean %>%
  filter(
    !is.na(
      sex
    )
  ) %>%
  distinct(
    individual_id,
    Site,
    sex
  ) %>%
  count(
    individual_id,
    Site,
    name =
      "n_sex_values"
  ) %>%
  filter(
    n_sex_values >
      1
  )


## Individuals associated with more than one PIT value.
## This is not automatically an error because tags can be replaced.
pit_history <- AllData_clean %>%
  filter(
    !is.na(
      pit_tag
    )
  ) %>%
  group_by(
    Site,
    individual_id
  ) %>%
  summarise(
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
        collapse =
          " | "
      ),

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
  filter(
    n_pit_tags >
      1
  ) %>%
  arrange(
    Site,
    individual_id
  )


## Northern known-age consistency.
known_age_consistency <- AllData_clean %>%
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
    individual_id
  ) %>%
  summarise(
    n_known_age_records =
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

    .groups =
      "drop"
  ) %>%
  arrange(
    desc(
      birth_year_range
    )
  )


## Raw feeding labels through time.
feeding_by_site_year <- AllData_clean %>%
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


## Permanent removals.
permanent_removals <- AllData_clean %>%
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


## Date corrections.
date_corrections <- AllData_clean %>%
  filter(
    date_corrected_2919
  ) %>%
  select(
    Site,
    individual_id,
    date_raw,
    capture_date,
    calendar_year,
    source_sheet,
    source_row
  )


## ============================================================
## 8. HARD CHECKS
## ============================================================

expected_sites <- c(
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


missing_sites <- setdiff(
  expected_sites,
  unique(
    AllData_clean$Site
  )
)


if(length(missing_sites) > 0) {
  stop(
    "Expected site(s) missing after import: ",
    paste(
      missing_sites,
      collapse = ", "
    )
  )
}


if(any(AllData_clean$calendar_year > 2019, na.rm = TRUE)) {
  stop(
    "Post-2019 records remain in the canonical dataset."
  )
}


if(any(is.na(AllData_clean$individual_id))) {
  stop(
    "Missing stable individual IDs remain."
  )
}


if(nrow(permanent_removals) != 4) {
  warning(
    "Expected 4 Leaf South -> Pasture Cay permanent removals; found ",
    nrow(
      permanent_removals
    ),
    ". Inspect PermanentRemovals.csv."
  )
}


## ============================================================
## 9. SAVE
## ============================================================

write_csv(
  AllData_clean,
  "Data/Iguana_Canonical_CaptureData_Through2019.csv",
  na = ""
)

saveRDS(
  AllData_clean,
  "Data/Iguana_Canonical_CaptureData_Through2019.rds"
)


write_csv(
  site_summary,
  file.path(
    out_dir,
    "01_SiteSummary.csv"
  )
)

write_csv(
  possible_duplicate_capture,
  file.path(
    out_dir,
    "02_PossibleDuplicateCaptures.csv"
  )
)

write_csv(
  sex_conflicts,
  file.path(
    out_dir,
    "03_SexConflicts.csv"
  )
)

write_csv(
  pit_history,
  file.path(
    out_dir,
    "04_MultiplePITsPerIndividual.csv"
  )
)

write_csv(
  known_age_consistency,
  file.path(
    out_dir,
    "05_KnownAgeConsistency.csv"
  )
)

write_csv(
  feeding_by_site_year,
  file.path(
    out_dir,
    "06_RawFeedingBySiteYear.csv"
  )
)

write_csv(
  permanent_removals,
  file.path(
    out_dir,
    "07_PermanentRemovals.csv"
  )
)

write_csv(
  date_corrections,
  file.path(
    out_dir,
    "08_DateCorrections.csv"
  )
)


## ============================================================
## 10. CONSOLE SUMMARY
## ============================================================

cat("\n")
cat("============================================================\n")
cat("CANONICAL CAPTURE DATASET CREATED\n")
cat("============================================================\n")
cat("Rows:", nrow(AllData_clean), "\n")
cat("Individuals:", n_distinct(AllData_clean$individual_id), "\n")
cat("Sites:", n_distinct(AllData_clean$Site), "\n")
cat("Years:", min(AllData_clean$calendar_year), "-", max(AllData_clean$calendar_year), "\n")
cat("\n")

print(
  site_summary,
  n = Inf
)

cat("\n")
cat("Possible duplicate animal-date records:",
    nrow(possible_duplicate_capture),
    "\n"
)

cat("Individuals with conflicting sex labels:",
    nrow(sex_conflicts),
    "\n"
)

cat("Individuals with >1 recorded PIT:",
    nrow(pit_history),
    "\n"
)

cat("Permanent removals flagged:",
    nrow(permanent_removals),
    "\n"
)

cat("2919 -> 2019 date corrections:",
    nrow(date_corrections),
    "\n"
)

cat("\nSaved:\n")
cat("  Data/Iguana_Canonical_CaptureData_Through2019.csv\n")
cat("  Data/Iguana_Canonical_CaptureData_Through2019.rds\n")
cat("  DataAudits/CanonicalImport/*.csv\n")
cat("============================================================\n")
