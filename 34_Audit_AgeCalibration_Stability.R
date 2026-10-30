## ============================================================
## 34_Audit_AgeCalibration_Stability.R
##
## Audit the SVL -> age quantile-regression calibration used to
## construct latent birth intervals.
##
## The current DataPrep fit can emit "Solution may be nonunique"
## from the default rq() Barrodale-Roberts solver. This script
## compares its actual predicted age bounds with the Frisch-Newton
## solution without changing the fitted survival model.
## ============================================================

library(tidyverse)
library(quantreg)

source("DataPrep_Canonical_LatentSex_V5.R")

out_dir <- "DataAudits/AgeCalibration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("SVL -> AGE CALIBRATION STABILITY AUDIT\n")
cat("============================================================\n")

support_by_site <- calibration_data %>%
  count(
    Species,
    Island,
    name = "n_calibration_records"
  ) %>%
  arrange(
    Species,
    desc(n_calibration_records)
  )

print(support_by_site, n = Inf)

write_csv(
  support_by_site,
  file.path(
    out_dir,
    "CalibrationSupport_BySpeciesSite.csv"
  )
)

fit_min_fn <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.05,
  data = calibration_data,
  method = "fn"
)

fit_max_fn <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.95,
  data = calibration_data,
  method = "fn"
)

predict_bounds_from_fits <- function(
    svl_value,
    sex_value,
    fit_min,
    fit_max
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
      levels = c("F", "M")
    )
  )

  p_min <- as.numeric(
    exp(
      predict(
        fit_min,
        newdata = nd
      )
    )
  )

  p_max <- as.numeric(
    exp(
      predict(
        fit_max,
        newdata = nd
      )
    )
  )

  p_max <- min(p_max, 40)
  p_min <- min(p_min, p_max - 1)

  if(svl_value < 25) {
    p_max <- min(p_max, 4)
  }

  min_age <- max(0, p_min)
  max_age <- max(min_age + 1, p_max)

  c(
    min_age = min_age,
    max_age = max_age
  )
}

female_fn <- t(
  vapply(
    growth_info$first_svl,
    predict_bounds_from_fits,
    sex_value = "F",
    fit_min = fit_min_fn,
    fit_max = fit_max_fn,
    FUN.VALUE = c(
      min_age = 0,
      max_age = 0
    )
  )
)

male_fn <- t(
  vapply(
    growth_info$first_svl,
    predict_bounds_from_fits,
    sex_value = "M",
    fit_min = fit_min_fn,
    fit_max = fit_max_fn,
    FUN.VALUE = c(
      min_age = 0,
      max_age = 0
    )
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
    fn_min_F = female_fn[, "min_age"],
    current_max_F = pred_max_age_F,
    fn_max_F = female_fn[, "max_age"],
    current_min_M = pred_min_age_M,
    fn_min_M = male_fn[, "min_age"],
    current_max_M = pred_max_age_M,
    fn_max_M = male_fn[, "max_age"]
  ) %>%
  mutate(
    diff_min_F = fn_min_F - current_min_F,
    diff_max_F = fn_max_F - current_max_F,
    diff_min_M = fn_min_M - current_min_M,
    diff_max_M = fn_max_M - current_max_M,
    max_abs_diff = pmax(
      abs(diff_min_F),
      abs(diff_max_F),
      abs(diff_min_M),
      abs(diff_max_M),
      na.rm = TRUE
    )
  )

write_csv(
  comparison,
  file.path(
    out_dir,
    "CurrentBR_vs_FN_AgeBounds.csv"
  )
)

summary_diff <- comparison %>%
  summarise(
    n_individuals = n(),
    median_max_abs_diff_years =
      median(max_abs_diff, na.rm = TRUE),
    q95_max_abs_diff_years =
      quantile(
        max_abs_diff,
        0.95,
        na.rm = TRUE,
        names = FALSE
      ),
    max_abs_diff_years =
      max(max_abs_diff, na.rm = TRUE),
    n_gt_0_5yr =
      sum(max_abs_diff > 0.5, na.rm = TRUE),
    n_gt_1yr =
      sum(max_abs_diff > 1, na.rm = TRUE)
  )

print(summary_diff)

write_csv(
  summary_diff,
  file.path(
    out_dir,
    "CurrentBR_vs_FN_AgeBounds_Summary.csv"
  )
)

cat("\nCurrent/default rq coefficients (5th percentile):\n")
print(coef(fit_min_log))
cat("\nFrisch-Newton coefficients (5th percentile):\n")
print(coef(fit_min_fn))

cat("\nCurrent/default rq coefficients (95th percentile):\n")
print(coef(fit_max_log))
cat("\nFrisch-Newton coefficients (95th percentile):\n")
print(coef(fit_max_fn))

cat("\nInterpretation:\n")
cat(
  "- If prediction differences are effectively zero, the nonunique-solution ",
  "warning is computational rather than biologically material.\n"
)
cat(
  "- If many individuals move by >0.5-1 year, freeze a single stable ",
  "calibration method before rerunning long MCMC.\n"
)
cat(
  "- Calibration support by species/site should also be reported because ",
  "age information is not equally distributed among populations.\n"
)

cat("\nSaved to ", out_dir, "\n", sep = "")
