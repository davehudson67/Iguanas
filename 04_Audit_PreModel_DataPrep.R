## ============================================================
## 04_Audit_PreModel_DataPrep.R
##
## Final audit before adapting/running the NIMBLE survival model.
##
## Audits:
##   1. confirmed deaths by site/year
##   2. multiple death records per individual
##   3. latent-sex composition by site
##   4. SVL-age calibration sample sizes
##   5. birth-window widths by information source
##   6. whether the rq() non-unique-solution warning materially
##      changes predicted age bounds
##
## Run after DataPrep_Canonical_LatentSex_V5.R, or this script will
## source it automatically if the required objects are absent.
## ============================================================


library(tidyverse)
library(quantreg)


if(
  !exists(
    "AllData"
  ) ||
    !exists(
      "growth_info"
    ) ||
    !exists(
      "fit_min_log"
    ) ||
    !exists(
      "fit_max_log"
    )
) {

  source(
    "DataPrep_Canonical_LatentSex_V5.R"
  )
}


out_dir <- "DataAudits/PreModel"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


## ============================================================
## 1. DEATH AUDIT
## ============================================================

death_detail <- AllData %>%
  filter(
    animal_dead
  ) %>%
  arrange(
    Island,
    ID_raw,
    CapDate
  ) %>%
  group_by(
    ID_raw
  ) %>%
  mutate(
    death_number_for_id =
      row_number()
  ) %>%
  ungroup() %>%
  select(
    Island,
    ID_raw,
    animal_id,
    CapDate,
    CalYear,
    capture_then_death,
    death_delay_years,
    sex_clean,
    comments,
    source_sheet,
    source_row,
    death_number_for_id
  )


death_by_site <- death_detail %>%
  group_by(
    Island
  ) %>%
  summarise(
    n_death_rows =
      n(),

    n_death_individuals =
      n_distinct(
        ID_raw
      ),

    first_death_year =
      min(
        CalYear
      ),

    last_death_year =
      max(
        CalYear
      ),

    n_capture_then_death =
      sum(
        capture_then_death
      ),

    .groups =
      "drop"
  )


multiple_deaths_per_id <- death_detail %>%
  count(
    Island,
    ID_raw,
    name =
      "n_death_rows"
  ) %>%
  filter(
    n_death_rows >
      1
  ) %>%
  arrange(
    desc(
      n_death_rows
    ),
    Island,
    ID_raw
  )


## ============================================================
## 2. SEX AUDIT
## ============================================================

sex_by_site <- growth_info %>%
  count(
    island,
    sex_state,
    name =
      "n_individuals"
  ) %>%
  group_by(
    island
  ) %>%
  mutate(
    proportion =
      n_individuals /
      sum(
        n_individuals
      )
  ) %>%
  ungroup() %>%
  arrange(
    island,
    sex_state
  )


latent_sex_individuals <- growth_info %>%
  filter(
    is.na(
      sex_clean
    )
  ) %>%
  select(
    ID,
    ID_raw,
    island,
    species,
    sex_state,
    n_F,
    n_M,
    n_unknown_sex,
    first_svl,
    pred_min_age_F,
    pred_max_age_F,
    pred_min_age_M,
    pred_max_age_M
  ) %>%
  arrange(
    island,
    sex_state,
    ID_raw
  )


## ============================================================
## 3. CALIBRATION DATA SUPPORT
## ============================================================

calibration_support <- calibration_data %>%
  count(
    Island,
    sex_clean,
    name =
      "n_records"
  ) %>%
  arrange(
    Island,
    sex_clean
  )


calibration_individual_support <- calibration_data %>%
  distinct(
    Island,
    ID_raw,
    sex_clean
  ) %>%
  count(
    Island,
    sex_clean,
    name =
      "n_individuals"
  ) %>%
  arrange(
    Island,
    sex_clean
  )


## ============================================================
## 4. BIRTH-WINDOW SOURCE AND WIDTH
## ============================================================

birth_window_audit <- growth_info %>%
  mutate(
    birth_basis =
      case_when(
        known_age_consistent &
          !is.na(
            known_birth_low
          ) &
          !is.na(
            known_birth_high
          ) ~
          "known_age",

        !is.na(
          first_svl
        ) &
          !is.na(
            sex_clean
          ) ~
          "SVL_known_sex",

        !is.na(
          first_svl
        ) &
          is.na(
            sex_clean
          ) ~
          "SVL_latent_sex_union",

        TRUE ~
          "no_SVL"
      ),

    birth_low_cal =
      cintB_adj[
        ,
        1
      ] +
      Origin -
      1,

    birth_high_cal =
      cintB_adj[
        ,
        2
      ] +
      Origin -
      1,

    birth_width =
      birth_high_cal -
      birth_low_cal
  ) %>%
  select(
    ID,
    ID_raw,
    island,
    species,
    sex_state,
    birth_basis,
    first_svl,
    known_age_consistent,
    birth_low_cal,
    birth_high_cal,
    birth_width
  )


birth_window_summary <- birth_window_audit %>%
  group_by(
    birth_basis
  ) %>%
  summarise(
    n =
      n(),

    median_width =
      median(
        birth_width
      ),

    mean_width =
      mean(
        birth_width
      ),

    q95_width =
      quantile(
        birth_width,
        0.95
      ),

    max_width =
      max(
        birth_width
      ),

    .groups =
      "drop"
  )


birth_window_by_site <- birth_window_audit %>%
  group_by(
    island,
    birth_basis
  ) %>%
  summarise(
    n =
      n(),

    median_width =
      median(
        birth_width
      ),

    mean_width =
      mean(
        birth_width
      ),

    max_width =
      max(
        birth_width
      ),

    .groups =
      "drop"
  ) %>%
  arrange(
    island,
    birth_basis
  )


## ============================================================
## 5. QUANTILE-REGRESSION SOLUTION STABILITY
##
## rq() warned that the Barrodale-Roberts solution may be nonunique.
## Fit the same models with the Frisch-Newton algorithm and compare
## PREDICTIONS rather than coefficients.
## ============================================================

fit_min_fn <- rq(
  log_age ~
    svl +
    I(
      svl^2
    ) +
    sex_f,
  tau =
    0.05,
  data =
    calibration_data,
  method =
    "fn"
)


fit_max_fn <- rq(
  log_age ~
    svl +
    I(
      svl^2
    ) +
    sex_f,
  tau =
    0.95,
  data =
    calibration_data,
  method =
    "fn"
)


svl_grid <- seq(
  min(
    calibration_data$svl,
    na.rm = TRUE
  ),
  max(
    calibration_data$svl,
    na.rm = TRUE
  ),
  length.out =
    200
)


prediction_grid <- expand_grid(
  svl =
    svl_grid,
  sex_f =
    factor(
      c(
        "F",
        "M"
      ),
      levels =
        c(
          "F",
          "M"
        )
    )
) %>%
  mutate(
    p05_br =
      exp(
        predict(
          fit_min_log,
          newdata =
            .
        )
      ),

    p95_br =
      exp(
        predict(
          fit_max_log,
          newdata =
            .
        )
      ),

    p05_fn =
      exp(
        predict(
          fit_min_fn,
          newdata =
            .
        )
      ),

    p95_fn =
      exp(
        predict(
          fit_max_fn,
          newdata =
            .
        )
      ),

    diff_p05 =
      p05_fn -
      p05_br,

    diff_p95 =
      p95_fn -
      p95_br,

    abs_diff_p05 =
      abs(
        diff_p05
      ),

    abs_diff_p95 =
      abs(
        diff_p95
      ),

    quantiles_cross_br =
      p05_br >
      p95_br,

    quantiles_cross_fn =
      p05_fn >
      p95_fn
  )


rq_stability_summary <- prediction_grid %>%
  summarise(
    median_abs_diff_p05 =
      median(
        abs_diff_p05
      ),

    max_abs_diff_p05 =
      max(
        abs_diff_p05
      ),

    median_abs_diff_p95 =
      median(
        abs_diff_p95
      ),

    max_abs_diff_p95 =
      max(
        abs_diff_p95
      ),

    n_crossing_br =
      sum(
        quantiles_cross_br
      ),

    n_crossing_fn =
      sum(
        quantiles_cross_fn
      )
  )


## ============================================================
## 6. WRITE OUTPUTS
## ============================================================

write_csv(
  death_detail,
  file.path(
    out_dir,
    "01_DeathDetail.csv"
  )
)


write_csv(
  death_by_site,
  file.path(
    out_dir,
    "01_DeathsBySite.csv"
  )
)


write_csv(
  multiple_deaths_per_id,
  file.path(
    out_dir,
    "01_MultipleDeathRowsPerIndividual.csv"
  )
)


write_csv(
  sex_by_site,
  file.path(
    out_dir,
    "02_SexStateBySite.csv"
  )
)


write_csv(
  latent_sex_individuals,
  file.path(
    out_dir,
    "02_LatentSexIndividuals.csv"
  )
)


write_csv(
  calibration_support,
  file.path(
    out_dir,
    "03_CalibrationRecordsBySiteSex.csv"
  )
)


write_csv(
  calibration_individual_support,
  file.path(
    out_dir,
    "03_CalibrationIndividualsBySiteSex.csv"
  )
)


write_csv(
  birth_window_audit,
  file.path(
    out_dir,
    "04_BirthWindowAudit.csv"
  )
)


write_csv(
  birth_window_summary,
  file.path(
    out_dir,
    "04_BirthWindowSummary.csv"
  )
)


write_csv(
  birth_window_by_site,
  file.path(
    out_dir,
    "04_BirthWindowBySite.csv"
  )
)


write_csv(
  prediction_grid,
  file.path(
    out_dir,
    "05_RQPredictionStabilityGrid.csv"
  )
)


write_csv(
  rq_stability_summary,
  file.path(
    out_dir,
    "05_RQPredictionStabilitySummary.csv"
  )
)


## ============================================================
## 7. CONSOLE SUMMARY
## ============================================================

cat("\n")
cat("============================================================\n")
cat("PRE-MODEL AUDIT\n")
cat("============================================================\n")


cat("\nDEATHS BY SITE\n")
cat("------------------------------------------------------------\n")

print(
  death_by_site,
  n = Inf
)


cat("\nMULTIPLE DEATH ROWS PER INDIVIDUAL\n")
cat("------------------------------------------------------------\n")

print(
  multiple_deaths_per_id,
  n = Inf
)


cat("\nSEX STATE BY SITE\n")
cat("------------------------------------------------------------\n")

print(
  sex_by_site,
  n = Inf
)


cat("\nCALIBRATION RECORDS BY SITE / SEX\n")
cat("------------------------------------------------------------\n")

print(
  calibration_support,
  n = Inf
)


cat("\nCALIBRATION INDIVIDUALS BY SITE / SEX\n")
cat("------------------------------------------------------------\n")

print(
  calibration_individual_support,
  n = Inf
)


cat("\nBIRTH-WINDOW SUMMARY\n")
cat("------------------------------------------------------------\n")

print(
  birth_window_summary,
  n = Inf
)


cat("\nRQ SOLUTION-STABILITY CHECK\n")
cat("------------------------------------------------------------\n")

print(
  rq_stability_summary,
  n = Inf
)


cat("\n============================================================\n")
cat("AUDIT COMPLETE\n")
cat("Outputs: DataAudits/PreModel/*.csv\n")
cat("============================================================\n")
