## ============================================================
## 36_Audit_AgeCalibration_Weighting_and_Site.R
##
## Final targeted audit of the SVL -> age calibration.
##
## Questions:
## 1. Does the pooled calibration change if each biological
##    individual receives equal total weight despite repeat captures?
## 2. Does it change if calibration sites also receive equal total
##    weight?
## 3. How much site-specific age-at-size heterogeneity is present
##    within the known-age inornata calibration data?
##
## AUDIT ONLY. Does not change production DataPrep.
## ============================================================

library(tidyverse)
library(quantreg)

source("DataPrep_Canonical_LatentSex_V5.R")

out_dir <- "DataAudits/AgeCalibration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("AGE CALIBRATION: WEIGHTING + SITE-HETEROGENEITY AUDIT\n")
cat("============================================================\n")

## ------------------------------------------------------------
## Weights
## ------------------------------------------------------------

cal_w <- calibration_data %>%
  add_count(ID, name = "n_records_id") %>%
  group_by(Island) %>%
  mutate(
    n_unique_ids_site = n_distinct(ID)
  ) %>%
  ungroup() %>%
  mutate(
    ## Each individual contributes total weight 1.
    w_individual =
      1 / n_records_id,

    ## Each site contributes total weight 1 and individuals within
    ## site contribute equally.
    w_site_individual =
      1 / (
        n_records_id *
        n_unique_ids_site
      ),

    ## Constant rescaling has no effect on rq coefficient estimates,
    ## but keeps numerical weight magnitudes convenient.
    w_individual =
      w_individual /
      mean(w_individual),

    w_site_individual =
      w_site_individual /
      mean(w_site_individual)
  )

weight_audit <- cal_w %>%
  group_by(Island) %>%
  summarise(
    n_records = n(),
    n_unique_individuals = n_distinct(ID),
    total_unweighted = n(),
    total_individual_weight = sum(w_individual),
    total_site_individual_weight = sum(w_site_individual),
    .groups = "drop"
  )

cat("\nWeight totals by calibration site:\n")
print(weight_audit, n = Inf, width = Inf)

## ------------------------------------------------------------
## Weighted pooled quantile regressions
## ------------------------------------------------------------

fit_min_iw <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.05,
  data = cal_w,
  weights = w_individual,
  method = "fn"
)

fit_max_iw <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.95,
  data = cal_w,
  weights = w_individual,
  method = "fn"
)

fit_min_sw <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.05,
  data = cal_w,
  weights = w_site_individual,
  method = "fn"
)

fit_max_sw <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.95,
  data = cal_w,
  weights = w_site_individual,
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

get_bounds <- function(fit_min, fit_max, sex_value) {
  t(
    vapply(
      growth_info$first_svl,
      predict_bounds,
      sex_value = sex_value,
      fit_min = fit_min,
      fit_max = fit_max,
      FUN.VALUE = c(min_age = 0, max_age = 0)
    )
  )
}

f_iw <- get_bounds(fit_min_iw, fit_max_iw, "F")
m_iw <- get_bounds(fit_min_iw, fit_max_iw, "M")
f_sw <- get_bounds(fit_min_sw, fit_max_sw, "F")
m_sw <- get_bounds(fit_min_sw, fit_max_sw, "M")

comparison <- growth_info %>%
  transmute(
    ID,
    ID_raw,
    island,
    species,
    first_svl,

    current_min_F = pred_min_age_F,
    current_max_F = pred_max_age_F,
    current_min_M = pred_min_age_M,
    current_max_M = pred_max_age_M,

    iw_min_F = f_iw[, "min_age"],
    iw_max_F = f_iw[, "max_age"],
    iw_min_M = m_iw[, "min_age"],
    iw_max_M = m_iw[, "max_age"],

    sw_min_F = f_sw[, "min_age"],
    sw_max_F = f_sw[, "max_age"],
    sw_min_M = m_sw[, "min_age"],
    sw_max_M = m_sw[, "max_age"]
  ) %>%
  mutate(
    individual_weight_max_abs_diff =
      pmax(
        abs(iw_min_F - current_min_F),
        abs(iw_max_F - current_max_F),
        abs(iw_min_M - current_min_M),
        abs(iw_max_M - current_max_M),
        na.rm = TRUE
      ),

    site_individual_weight_max_abs_diff =
      pmax(
        abs(sw_min_F - current_min_F),
        abs(sw_max_F - current_max_F),
        abs(sw_min_M - current_min_M),
        abs(sw_max_M - current_max_M),
        na.rm = TRUE
      )
  )

summarise_diff <- function(data, col, scheme) {
  x <- data[[col]]

  tibble(
    weighting_scheme = scheme,
    n_individuals = nrow(data),
    median_max_abs_diff_years = median(x, na.rm = TRUE),
    q95_max_abs_diff_years =
      quantile(x, 0.95, na.rm = TRUE, names = FALSE),
    max_abs_diff_years = max(x, na.rm = TRUE),
    n_gt_0_5yr = sum(x > 0.5, na.rm = TRUE),
    n_gt_1yr = sum(x > 1, na.rm = TRUE),
    n_gt_2yr = sum(x > 2, na.rm = TRUE),
    prop_gt_1yr = mean(x > 1, na.rm = TRUE)
  )
}

summary_weighting <- bind_rows(
  comparison %>%
    group_by(species) %>%
    group_modify(
      ~ bind_rows(
        summarise_diff(
          .x,
          "individual_weight_max_abs_diff",
          "equal total weight per individual"
        ),
        summarise_diff(
          .x,
          "site_individual_weight_max_abs_diff",
          "equal site + individual weight"
        )
      )
    ) %>%
    ungroup(),

  bind_rows(
    summarise_diff(
      comparison,
      "individual_weight_max_abs_diff",
      "equal total weight per individual"
    ),
    summarise_diff(
      comparison,
      "site_individual_weight_max_abs_diff",
      "equal site + individual weight"
    )
  ) %>%
    mutate(species = "ALL", .before = 1)
)

cat("\nProduction vs weighted-calibration age-bound differences:\n")
print(summary_weighting, n = Inf, width = Inf)

cat("\nProduction pooled coefficients:\n")
cat("5th percentile:\n")
print(coef(fit_min_log))
cat("95th percentile:\n")
print(coef(fit_max_log))

cat("\nEqual-individual-weight coefficients:\n")
cat("5th percentile:\n")
print(coef(fit_min_iw))
cat("95th percentile:\n")
print(coef(fit_max_iw))

cat("\nEqual-site + equal-individual-weight coefficients:\n")
cat("5th percentile:\n")
print(coef(fit_min_sw))
cat("95th percentile:\n")
print(coef(fit_max_sw))

## ------------------------------------------------------------
## Site heterogeneity within known-age inornata data
##
## Add site-specific intercept and SVL slope. This is not proposed
## as the production calibration; it quantifies whether age-at-size
## differs materially among the three calibration populations.
## ------------------------------------------------------------

cal_w <- cal_w %>%
  mutate(
    Island_f = factor(
      Island,
      levels = c("LeafNorth", "UCay", "FFRC")
    )
  )

fit_min_site <- rq(
  log_age ~ svl + I(svl^2) + sex_f + Island_f + svl:Island_f,
  tau = 0.05,
  data = cal_w,
  weights = w_individual,
  method = "fn"
)

fit_max_site <- rq(
  log_age ~ svl + I(svl^2) + sex_f + Island_f + svl:Island_f,
  tau = 0.95,
  data = cal_w,
  weights = w_individual,
  method = "fn"
)

grid <- expand_grid(
  svl = c(15, 20, 25, 30, 35, 40),
  sex_f = factor(c("F", "M"), levels = c("F", "M")),
  Island_f = factor(
    c("LeafNorth", "UCay", "FFRC"),
    levels = c("LeafNorth", "UCay", "FFRC")
  )
) %>%
  mutate(
    age_q05 = exp(
      predict(fit_min_site, newdata = .)
    ),
    age_q95 = exp(
      predict(fit_max_site, newdata = .)
    )
  )

site_spread <- grid %>%
  group_by(svl, sex_f) %>%
  summarise(
    min_site_q05 = min(age_q05),
    max_site_q05 = max(age_q05),
    spread_q05_years = max_site_q05 - min_site_q05,
    min_site_q95 = min(age_q95),
    max_site_q95 = max(age_q95),
    spread_q95_years = max_site_q95 - min_site_q95,
    .groups = "drop"
  )

cat("\nSite-specific known-age predictions at representative SVLs:\n")
print(grid, n = Inf, width = Inf)

cat("\nAcross-site spread in predicted age quantiles:\n")
print(site_spread, n = Inf, width = Inf)

write_csv(
  weight_audit,
  file.path(out_dir, "Calibration_WeightTotals_BySite.csv")
)

write_csv(
  comparison,
  file.path(out_dir, "Production_vs_Weighted_AgeBounds.csv")
)

write_csv(
  summary_weighting,
  file.path(out_dir, "Production_vs_Weighted_AgeBounds_Summary.csv")
)

write_csv(
  grid,
  file.path(out_dir, "SiteSpecific_AgeAtSVL_Predictions.csv")
)

write_csv(
  site_spread,
  file.path(out_dir, "SiteSpecific_AgeAtSVL_Spread.csv")
)

cat("\nInterpretation guide:\n")
cat("- Large individual-weight differences mean repeat captures materially influence the pooled age envelope.\n")
cat("- Large site-balanced differences mean calibration-site imbalance also matters.\n")
cat("- Large across-site age-at-SVL spreads show that a single pooled SVL clock does not represent all calibration populations equally.\n")
cat("- These results should be used to choose a conservative production age treatment before final MCMC.\n")

cat("\nSaved to ", out_dir, "\n", sep = "")
