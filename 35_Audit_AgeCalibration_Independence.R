## ============================================================
## 35_Audit_AgeCalibration_Independence.R
##
## Sensitivity audit for SVL -> age calibration.
##
## The production DataPrep quantile regressions use all eligible
## known-age records. Because known-age individuals may be captured
## repeatedly, this audit asks whether repeated records materially
## influence the age bounds by refitting the same 5th/95th percentile
## models using ONE eligible record per biological individual.
##
## This is an audit only; it does not alter production DataPrep.
## ============================================================

library(tidyverse)
library(quantreg)

source("DataPrep_Canonical_LatentSex_V5.R")

out_dir <- "DataAudits/AgeCalibration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("AGE CALIBRATION: REPEATED-RECORD INDEPENDENCE AUDIT\n")
cat("============================================================\n")

## Production calibration record counts and unique individuals.
cal_support <- calibration_data %>%
  summarise(
    n_records = n(),
    n_unique_individuals = n_distinct(ID),
    records_per_individual_mean = n() / n_distinct(ID),
    records_per_individual_max = max(table(ID))
  )

cat("\nCurrent calibration support:\n")
print(cal_support, width = Inf)

support_site <- calibration_data %>%
  group_by(Species, Island) %>%
  summarise(
    n_records = n(),
    n_unique_individuals = n_distinct(ID),
    .groups = "drop"
  )

cat("\nCalibration support by species/site:\n")
print(support_site, n = Inf, width = Inf)

## One observation per individual. Use the earliest eligible record so
## this does not preferentially select a later/larger measurement.
calibration_one <- calibration_data %>%
  arrange(ID, CapDate) %>%
  group_by(ID) %>%
  slice(1L) %>%
  ungroup()

stopifnot(
  nrow(calibration_one) ==
    n_distinct(calibration_data$ID)
)

fit_min_one <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.05,
  data = calibration_one,
  method = "fn"
)

fit_max_one <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.95,
  data = calibration_one,
  method = "fn"
)

predict_bounds <- function(
    svl_value,
    sex_value,
    fit_min,
    fit_max
) {

  if(is.na(svl_value)) {
    return(c(min_age = 0, max_age = 40))
  }

  nd <- data.frame(
    svl = svl_value,
    sex_f = factor(
      sex_value,
      levels = c("F", "M")
    )
  )

  p_min <- as.numeric(
    exp(predict(fit_min, newdata = nd))
  )

  p_max <- as.numeric(
    exp(predict(fit_max, newdata = nd))
  )

  p_max <- min(p_max, 40)
  p_min <- min(p_min, p_max - 1)

  if(svl_value < 25) {
    p_max <- min(p_max, 4)
  }

  min_age <- max(0, p_min)
  max_age <- max(min_age + 1, p_max)

  c(min_age = min_age, max_age = max_age)
}

female_one <- t(
  vapply(
    growth_info$first_svl,
    predict_bounds,
    sex_value = "F",
    fit_min = fit_min_one,
    fit_max = fit_max_one,
    FUN.VALUE = c(min_age = 0, max_age = 0)
  )
)

male_one <- t(
  vapply(
    growth_info$first_svl,
    predict_bounds,
    sex_value = "M",
    fit_min = fit_min_one,
    fit_max = fit_max_one,
    FUN.VALUE = c(min_age = 0, max_age = 0)
  )
)

comparison <- growth_info %>%
  transmute(
    ID,
    ID_raw,
    island,
    species,
    sex_state,
    first_svl,
    current_min_F = pred_min_age_F,
    one_min_F = female_one[, "min_age"],
    current_max_F = pred_max_age_F,
    one_max_F = female_one[, "max_age"],
    current_min_M = pred_min_age_M,
    one_min_M = male_one[, "min_age"],
    current_max_M = pred_max_age_M,
    one_max_M = male_one[, "max_age"]
  ) %>%
  mutate(
    diff_min_F = one_min_F - current_min_F,
    diff_max_F = one_max_F - current_max_F,
    diff_min_M = one_min_M - current_min_M,
    diff_max_M = one_max_M - current_max_M,
    max_abs_diff = pmax(
      abs(diff_min_F),
      abs(diff_max_F),
      abs(diff_min_M),
      abs(diff_max_M),
      na.rm = TRUE
    )
  )

summary_all <- comparison %>%
  group_by(species) %>%
  summarise(
    n_individuals = n(),
    median_max_abs_diff_years =
      median(max_abs_diff, na.rm = TRUE),
    q95_max_abs_diff_years =
      quantile(max_abs_diff, 0.95, na.rm = TRUE, names = FALSE),
    max_abs_diff_years =
      max(max_abs_diff, na.rm = TRUE),
    n_gt_0_5yr =
      sum(max_abs_diff > 0.5, na.rm = TRUE),
    n_gt_1yr =
      sum(max_abs_diff > 1, na.rm = TRUE),
    n_gt_2yr =
      sum(max_abs_diff > 2, na.rm = TRUE),
    prop_gt_1yr =
      mean(max_abs_diff > 1, na.rm = TRUE),
    .groups = "drop"
  )

summary_total <- comparison %>%
  summarise(
    species = "ALL",
    n_individuals = n(),
    median_max_abs_diff_years =
      median(max_abs_diff, na.rm = TRUE),
    q95_max_abs_diff_years =
      quantile(max_abs_diff, 0.95, na.rm = TRUE, names = FALSE),
    max_abs_diff_years =
      max(max_abs_diff, na.rm = TRUE),
    n_gt_0_5yr =
      sum(max_abs_diff > 0.5, na.rm = TRUE),
    n_gt_1yr =
      sum(max_abs_diff > 1, na.rm = TRUE),
    n_gt_2yr =
      sum(max_abs_diff > 2, na.rm = TRUE),
    prop_gt_1yr =
      mean(max_abs_diff > 1, na.rm = TRUE)
  )

summary_out <- bind_rows(summary_all, summary_total)

cat("\nCurrent repeated-record vs one-record-per-individual age-bound differences:\n")
print(summary_out, n = Inf, width = Inf)

cat("\nProduction coefficients (all eligible records):\n")
cat("5th percentile:\n")
print(coef(fit_min_log))
cat("95th percentile:\n")
print(coef(fit_max_log))

cat("\nOne-record-per-individual coefficients:\n")
cat("5th percentile:\n")
print(coef(fit_min_one))
cat("95th percentile:\n")
print(coef(fit_max_one))

## SVL support relative to the production calibration range.
cal_svl_min <- min(calibration_data$svl, na.rm = TRUE)
cal_svl_max <- max(calibration_data$svl, na.rm = TRUE)

target_support <- growth_info %>%
  filter(!is.na(first_svl)) %>%
  mutate(
    outside_calibration_svl_range =
      first_svl < cal_svl_min |
      first_svl > cal_svl_max
  ) %>%
  group_by(species) %>%
  summarise(
    n_individuals_with_svl = n(),
    min_first_svl = min(first_svl),
    q05_first_svl =
      quantile(first_svl, 0.05, names = FALSE),
    median_first_svl = median(first_svl),
    q95_first_svl =
      quantile(first_svl, 0.95, names = FALSE),
    max_first_svl = max(first_svl),
    n_outside_calibration_svl_range =
      sum(outside_calibration_svl_range),
    prop_outside_calibration_svl_range =
      mean(outside_calibration_svl_range),
    .groups = "drop"
  )

cat("\nTarget SVL support relative to production calibration range [",
    cal_svl_min, ", ", cal_svl_max, "]:\n", sep = "")
print(target_support, n = Inf, width = Inf)

write_csv(
  cal_support,
  file.path(out_dir, "Calibration_RepeatedRecord_Support.csv")
)

write_csv(
  support_site,
  file.path(out_dir, "Calibration_UniqueIndividuals_BySpeciesSite.csv")
)

write_csv(
  comparison,
  file.path(out_dir, "Repeated_vs_OneRecord_AgeBounds.csv")
)

write_csv(
  summary_out,
  file.path(out_dir, "Repeated_vs_OneRecord_AgeBounds_Summary.csv")
)

write_csv(
  target_support,
  file.path(out_dir, "TargetSpecies_SVL_Extrapolation_Audit.csv")
)

cat("\nInterpretation guide:\n")
cat("- Small differences indicate that repeated known-age captures are not materially driving the age envelope.\n")
cat("- Large figginsi differences would make the cross-species age-transfer assumption more consequential.\n")
cat("- Individuals outside the calibration SVL range are true morphometric extrapolations and should be quantified.\n")

cat("\nSaved to ", out_dir, "\n", sep = "")
