## ============================================================
## 05_Build_Density_and_Exposure_Histories.R
##
## Rebuilds:
##   1. within-site MNA-derived cohort density for 9 sites
##   2. site-year provisioning/tourism histories from coauthors
##   3. FFRC Hurricane Matthew intervention covariates
##   4. an empty exact-history template for targeted male removals
##
## Run after freezing:
##   DataPrep_Canonical_LatentSex_V5.R
## ============================================================

library(tidyverse)

source(
  "DataPrep_Canonical_LatentSex_V5.R"
)

out_dir <- "DataAudits/ModelCovariates"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 1. MINIMUM KNOWN ALIVE FROM CLEANED LIFE HISTORIES
##
## MNA is based on being known alive between first and last
## confirmed live observation. Direct density information is
## anchored only on actual survey years.
## ============================================================

first_last_year <- growth_info %>%
  transmute(
    ID =
      ID,

    ID_raw =
      ID_raw,

    island =
      island,

    first_year =
      floor(
        tF_cal
      ),

    last_year =
      floor(
        tL_cal
      )
  )

if(
  any(
    first_last_year$last_year <
      first_last_year$first_year
  )
) {
  stop(
    "At least one individual has last_year < first_year."
  )
}

mna_all_years <- first_last_year %>%
  rowwise() %>%
  mutate(
    calendar_year =
      list(
        seq.int(
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
    name =
      "MNA"
  )

survey_years <- survey_events %>%
  transmute(
    island =
      Island,

    calendar_year =
      CalYear
  ) %>%
  distinct() %>%
  arrange(
    island,
    calendar_year
  )

mna_survey <- survey_years %>%
  left_join(
    mna_all_years,
    by =
      c(
        "island",
        "calendar_year"
      )
  ) %>%
  mutate(
    MNA =
      replace_na(
        MNA,
        0L
      ),

    log_MNA =
      log1p(
        MNA
      )
  )


## ============================================================
## 2. WITHIN-SITE DENSITY
##
## We model relative cohort environment WITHIN each population,
## using deviation in log(1 + MNA) from that site's survey-year
## mean. This avoids carrying the old pooled Leaf Cay structure
## into the corrected analysis and does not require a between-site
## habitat-area comparison.
## ============================================================

site_density_baseline <- mna_survey %>%
  group_by(
    island
  ) %>%
  summarise(
    n_survey_years =
      n(),

    first_survey_year =
      min(
        calendar_year
      ),

    last_survey_year =
      max(
        calendar_year
      ),

    mean_log_MNA =
      mean(
        log_MNA
      ),

    sd_log_MNA =
      sd(
        log_MNA
      ),

    .groups =
      "drop"
  )

mna_survey <- mna_survey %>%
  left_join(
    site_density_baseline,
    by =
      "island"
  ) %>%
  mutate(
    density_within =
      log_MNA -
      mean_log_MNA
  )

density_within_sd <- sd(
  mna_survey$density_within
)

if(
  !is.finite(
    density_within_sd
  ) ||
    density_within_sd <=
      0
) {
  stop(
    "Invalid within-site density SD."
  )
}

mna_survey <- mna_survey %>%
  mutate(
    density_within_z =
      density_within /
      density_within_sd
  )


## ============================================================
## 3. COMPLETE ANNUAL DENSITY TRAJECTORY
##
## Survey year                 -> observed
## Between survey years        -> linear interpolation
## Outside observed support    -> 0 = site mean
##
## We explicitly avoid extrapolating trends before/after support.
## A support matrix is retained for sensitivity analyses.
## ============================================================

calendar_years <- seq.int(
  Origin,
  2019
)

interpolate_site_density <- function(
  site_name
) {

  observed <- mna_survey %>%
    filter(
      island ==
        site_name
    ) %>%
    arrange(
      calendar_year
    )

  if(nrow(observed) < 2) {
    stop(
      "Need >=2 observed survey years for ",
      site_name
    )
  }

  interp <- approx(
    x =
      observed$calendar_year,

    y =
      observed$density_within_z,

    xout =
      calendar_years,

    rule =
      1
  )$y

  tibble(
    island =
      site_name,

    calendar_year =
      calendar_years,

    density_interpolated =
      interp
  ) %>%
  left_join(
    observed %>%
      select(
        island,
        calendar_year,
        MNA,
        log_MNA,
        density_observed =
          density_within_z
      ),
    by =
      c(
        "island",
        "calendar_year"
      )
  ) %>%
  mutate(
    density_source =
      case_when(
        !is.na(
          density_observed
        ) ~
          "observed_survey_year",

        !is.na(
          density_interpolated
        ) ~
          "interpolated_between_surveys",

        TRUE ~
          "outside_observed_support_site_mean"
      ),

    density_supported =
      density_source !=
      "outside_observed_support_site_mean",

    density_within_z_model =
      replace_na(
        density_interpolated,
        0
      )
  )
}

density_complete <- map_dfr(
  ordered_islands,
  interpolate_site_density
) %>%
  mutate(
    year_idx =
      calendar_year -
      Origin +
      1
  ) %>%
  arrange(
    island,
    calendar_year
  )

density_matrix <- matrix(
  0,
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

density_support_matrix <- matrix(
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

for(r in seq_len(nrow(density_complete))) {

  jj <- match(
    density_complete$island[r],
    ordered_islands
  )

  tt <- density_complete$year_idx[r]

  density_matrix[
    jj,
    tt
  ] <-
    density_complete$density_within_z_model[r]

  density_support_matrix[
    jj,
    tt
  ] <-
    as.integer(
      density_complete$density_supported[r]
    )
}


## ============================================================
## 4. COAUTHOR / EXPERT SITE-YEAR HISTORY
##
## Categorical states only:
##   none
##   low_intermit
##   high_regular
##
## These are site x calendar-year exposures. They must NOT be
## interpreted as individual treatment assignments.
## ============================================================

exposure_history <- tribble(
  ~Site, ~start_year, ~end_year, ~provisioning_state, ~tourism_state, ~history_note,

  "Bitter", 1980, 2005, "none", "none",
  "No organised tourism/provisioning described before 2006.",

  "Bitter", 2006, 2019, "high_regular", "high_regular",
  "Organised tour boats began in 2006; intense/daily feeding thereafter.",

  "FFRC", 1980, 2019, "none", "none",
  "No tourism/provisioning reported.",

  "Gaulin", 1980, 2009, "low_intermit", "low_intermit",
  "Longstanding visitation/provisioning at relatively low intensity.",

  "Gaulin", 2010, 2019, "low_intermit", "low_intermit",
  "Activity reportedly declined after Bitter became more popular; kept categorical rather than inventing a numeric decline.",

  "LeafNorth", 1980, 2018, "high_regular", "high_regular",
  "Tourism/provisioning present at the start of monitoring and increased through the long-term study.",

  "LeafNorth", 2019, 2019, "low_intermit", "low_intermit",
  "Beginning of the reported marked reduction in visitation during 2019-2023.",

  "LeafSouth", 1980, 2011, "none", "none",
  "Provisioning not documented before the reported onset around 2012.",

  "LeafSouth", 2012, 2019, "high_regular", "high_regular",
  "Intense/daily feeding reported from approximately 2012 onward.",

  "NAdder", 1980, 2019, "none", "none",
  "Reported not visited/provisioned.",

  "Noddy", 1980, 2019, "none", "none",
  "Reported not visited/provisioned.",

  "UCay", 1980, 2018, "high_regular", "high_regular",
  "Regular provisioning/tourism throughout long monitoring; lower visitation than Leaf North but established and increasing.",

  "UCay", 2019, 2019, "low_intermit", "low_intermit",
  "Beginning of the reported marked reduction in visitation during 2019-2023.",

  "WhiteB", 1980, 1997, "none", "none",
  "Feeding not documented before the late 1990s.",

  "WhiteB", 1998, 2009, "low_intermit", "low_intermit",
  "Limited/local provisioning beginning in the late 1990s.",

  "WhiteB", 2010, 2019, "low_intermit", "low_intermit",
  "Provisioning/visitation reportedly waned after about 2010; kept categorical."
)


## ============================================================
## 5. EXPAND EXPOSURE HISTORY TO SITE x YEAR
## ============================================================

expand_history_for_site <- function(
  site_name
) {

  map_dfr(
    calendar_years,
    function(yr) {

      hit <- exposure_history %>%
        filter(
          Site ==
            site_name,
          start_year <=
            yr,
          end_year >=
            yr
        )

      if(nrow(hit) > 1) {
        stop(
          "Overlapping exposure histories for ",
          site_name,
          " in ",
          yr
        )
      }

      if(nrow(hit) == 0) {
        return(
          tibble(
            Site =
              site_name,
            calendar_year =
              yr,
            exposure_supported =
              FALSE,
            provisioning_state =
              NA_character_,
            tourism_state =
              NA_character_,
            history_note =
              NA_character_
          )
        )
      }

      tibble(
        Site =
          site_name,
        calendar_year =
          yr,
        exposure_supported =
          TRUE,
        provisioning_state =
          hit$provisioning_state[1],
        tourism_state =
          hit$tourism_state[1],
        history_note =
          hit$history_note[1]
      )
    }
  )
}

site_year_exposure <- map_dfr(
  ordered_islands,
  expand_history_for_site
) %>%
  mutate(
    provision_high =
      as.integer(
        provisioning_state ==
          "high_regular"
      ),

    provision_any =
      as.integer(
        provisioning_state %in%
          c(
            "low_intermit",
            "high_regular"
          )
      ),

    tourism_any =
      as.integer(
        tourism_state %in%
          c(
            "low_intermit",
            "high_regular"
          )
      ),

    provision_high_model =
      replace_na(
        provision_high,
        0L
      ),

    provision_any_model =
      replace_na(
        provision_any,
        0L
      ),

    tourism_any_model =
      replace_na(
        tourism_any,
        0L
      ),

    ## Hurricane Matthew was in 2016.
    ## Primary = regime after 2016.
    matthew_post2016 =
      as.integer(
        Site ==
          "FFRC" &
          calendar_year >=
            2017
      ),

    ## Sensitivity = acute disturbance window.
    matthew_acute_2017_2018 =
      as.integer(
        Site ==
          "FFRC" &
          calendar_year >=
            2017 &
          calendar_year <=
            2018
      ),

    year_idx =
      calendar_year -
      Origin +
      1
  )


## ============================================================
## 6. EXPOSURE MATRICES FOR NIMBLE
## ============================================================

new_zero_matrix <- function() {

  matrix(
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
}

provision_high_matrix <- new_zero_matrix()
provision_any_matrix <- new_zero_matrix()
tourism_any_matrix <- new_zero_matrix()
matthew_post_matrix <- new_zero_matrix()
matthew_acute_matrix <- new_zero_matrix()
exposure_support_matrix <- new_zero_matrix()

for(r in seq_len(nrow(site_year_exposure))) {

  jj <- match(
    site_year_exposure$Site[r],
    ordered_islands
  )

  tt <- site_year_exposure$year_idx[r]

  provision_high_matrix[
    jj,
    tt
  ] <-
    site_year_exposure$provision_high_model[r]

  provision_any_matrix[
    jj,
    tt
  ] <-
    site_year_exposure$provision_any_model[r]

  tourism_any_matrix[
    jj,
    tt
  ] <-
    site_year_exposure$tourism_any_model[r]

  matthew_post_matrix[
    jj,
    tt
  ] <-
    site_year_exposure$matthew_post2016[r]

  matthew_acute_matrix[
    jj,
    tt
  ] <-
    site_year_exposure$matthew_acute_2017_2018[r]

  exposure_support_matrix[
    jj,
    tt
  ] <-
    as.integer(
      site_year_exposure$exposure_supported[r]
    )
}


## ============================================================
## 7. TARGETED LARGE-MALE REMOVAL HISTORY
##
## Exact site/year history is not currently available in the
## retained project material. Do not infer deliberate removals
## from subsequent disappearance of large males.
##
## When exact records are recovered, populate this file. The
## intended model term is:
##
##   betaREMOVAL * sex[i] * removal_matrix[site, year]
##
## ============================================================

removal_file <- "Data/ExpertHistory_MaleRemoval.csv"

if(!file.exists(removal_file)) {

  write_csv(
    tibble(
      Site =
        character(),
      start_year =
        integer(),
      end_year =
        integer(),
      removal_state =
        character(),
      source =
        character(),
      notes =
        character()
    ),
    removal_file
  )

  message(
    "Created empty male-removal history template: ",
    removal_file
  )
}


## ============================================================
## 8. COMBINED SITE-YEAR TABLE
## ============================================================

model_site_year <- density_complete %>%
  rename(
    Site =
      island
  ) %>%
  left_join(
    site_year_exposure %>%
      select(
        Site,
        calendar_year,
        exposure_supported,
        provisioning_state,
        tourism_state,
        provision_high,
        provision_any,
        tourism_any,
        matthew_post2016,
        matthew_acute_2017_2018,
        history_note
      ),
    by =
      c(
        "Site",
        "calendar_year"
      )
  ) %>%
  arrange(
    Site,
    calendar_year
  )


## ============================================================
## 9. HOW MUCH OF EACH BIRTH WINDOW HAS DENSITY SUPPORT?
## ============================================================

birth_density_support <- map_dfr(
  seq_len(
    nind
  ),
  function(i) {

    low_year <- floor(
      cintB_adj[
        i,
        1
      ] +
        Origin -
        1
    )

    high_year <- floor(
      cintB_adj[
        i,
        2
      ] +
        Origin -
        1
    )

    possible_years <- seq.int(
      low_year,
      high_year
    )

    site_name <- growth_info$island[i]

    support <- density_complete %>%
      filter(
        island ==
          site_name,
        calendar_year %in%
          possible_years
      )

    tibble(
      ID =
        i,

      ID_raw =
        growth_info$ID_raw[i],

      island =
        site_name,

      birth_low_year =
        low_year,

      birth_high_year =
        high_year,

      n_possible_birth_years =
        length(
          possible_years
        ),

      n_density_supported_years =
        sum(
          support$density_supported
        ),

      prop_density_supported =
        n_density_supported_years /
        n_possible_birth_years
    )
  }
)


## ============================================================
## 10. SAVE
## ============================================================

write_csv(
  mna_survey,
  file.path(
    out_dir,
    "01_MNA_ObservedSurveyYears.csv"
  )
)

write_csv(
  site_density_baseline,
  file.path(
    out_dir,
    "01_Density_SiteBaselines.csv"
  )
)

write_csv(
  density_complete,
  file.path(
    out_dir,
    "01_Density_SiteYearModelTable.csv"
  )
)

write_csv(
  exposure_history,
  file.path(
    out_dir,
    "02_ExpertExposure_HistoryIntervals.csv"
  )
)

write_csv(
  site_year_exposure,
  file.path(
    out_dir,
    "02_ExpertExposure_SiteYear.csv"
  )
)

write_csv(
  model_site_year,
  file.path(
    out_dir,
    "03_Model_SiteYear_Covariates.csv"
  )
)

write_csv(
  birth_density_support,
  file.path(
    out_dir,
    "04_BirthWindow_DensitySupport.csv"
  )
)

saveRDS(
  list(
    Origin =
      Origin,

    tMax =
      tMax,

    ordered_islands =
      ordered_islands,

    density_matrix =
      density_matrix,

    density_support_matrix =
      density_support_matrix,

    density_within_sd =
      density_within_sd,

    density_table =
      density_complete,

    provision_high_matrix =
      provision_high_matrix,

    provision_any_matrix =
      provision_any_matrix,

    tourism_any_matrix =
      tourism_any_matrix,

    matthew_post_matrix =
      matthew_post_matrix,

    matthew_acute_matrix =
      matthew_acute_matrix,

    exposure_support_matrix =
      exposure_support_matrix,

    exposure_history =
      exposure_history,

    site_year =
      model_site_year,

    birth_density_support =
      birth_density_support
  ),
  "Data/Iguana_ModelCovariates_Through2019.rds"
)


## ============================================================
## 11. CONSOLE SUMMARY
## ============================================================

cat("\n")
cat("============================================================\n")
cat("DENSITY + EXPOSURE COVARIATES COMPLETE\n")
cat("============================================================\n")

cat("\nMNA / DENSITY SUPPORT BY SITE\n")
cat("------------------------------------------------------------\n")

print(
  mna_survey %>%
    group_by(
      island
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

      min_MNA =
        min(
          MNA
        ),

      median_MNA =
        median(
          MNA
        ),

      max_MNA =
        max(
          MNA
        ),

      .groups =
        "drop"
    ),
  n = Inf
)

cat("\nSITE-YEAR EXPOSURE INTERVALS\n")
cat("------------------------------------------------------------\n")

print(
  exposure_history,
  n = Inf
)

cat("\nBIRTH-WINDOW DENSITY SUPPORT\n")
cat("------------------------------------------------------------\n")

print(
  summary(
    birth_density_support$prop_density_supported
  )
)

cat(
  "\nIndividuals with zero supported density years in possible birth interval: ",
  sum(
    birth_density_support$prop_density_supported ==
      0
  ),
  "\n",
  sep =
    ""
)

cat(
  "Individuals with fully supported density birth interval: ",
  sum(
    birth_density_support$prop_density_supported ==
      1
  ),
  "\n",
  sep =
    ""
)

cat("\nMODEL CODING DECISIONS\n")
cat("------------------------------------------------------------\n")
cat("Primary provisioning: high/regular vs other states\n")
cat("Sensitivity provisioning: any provisioning vs none\n")
cat("FFRC primary disturbance: 2017 onward (post-2016)\n")
cat("FFRC sensitivity disturbance: 2016-2018\n")
cat("Targeted male removals: not yet coded without exact site/year history\n")

cat("\nSaved:\n")
cat("  Data/Iguana_ModelCovariates_Through2019.rds\n")
cat("  DataAudits/ModelCovariates/*.csv\n")
cat("  Data/ExpertHistory_MaleRemoval.csv\n")
cat("============================================================\n")
