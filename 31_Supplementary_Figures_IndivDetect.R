## ============================================================
## 31_Supplementary_Figures_IndivDetect.R
##
## Supplementary figures for the final individual-detection-RE model.
## Produces compact diagnostic / supporting plots from the final posterior
## and frozen covariate histories.
## ============================================================

library(tidyverse)

set.seed(20261027)

out_dir <- "Figures/Supplementary"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_figure <- function(plot, stem, width = 7.2, height = 5.2) {
  ggsave(
    file.path(out_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in"
  )

  ggsave(
    file.path(out_dir, paste0(stem, ".tiff")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw"
  )
}

extract_draw_matrix <- function(x) {
  if (inherits(x, "mcmc") || inherits(x, "mcmc.list")) return(as.matrix(x))
  if (is.matrix(x) || is.data.frame(x)) return(as.matrix(x))
  if (!is.list(x)) return(NULL)

  for (nm in c("samples", "combined_samples", "posterior_samples", "draws", "mcmc")) {
    if (!is.null(x[[nm]])) {
      ans <- extract_draw_matrix(x[[nm]])
      if (!is.null(ans) && ncol(ans) > 0) return(ans)
    }
  }

  NULL
}

read_chain_draws <- function(files) {
  mats <- map(files, ~ extract_draw_matrix(readRDS(.x)))
  if (any(map_lgl(mats, is.null))) stop("Could not extract posterior from one or more chain files.")
  common <- reduce(map(mats, colnames), intersect)
  bind_rows(map(mats, ~ as_tibble(.x[, common, drop = FALSE]))) %>%
    as.matrix()
}

read_model_draws <- function(combined_candidates, chain_pattern) {
  combined_file <- combined_candidates[file.exists(combined_candidates)][1]
  if (!is.na(combined_file)) {
    ans <- extract_draw_matrix(readRDS(combined_file))
    if (!is.null(ans)) return(ans)
  }

  files <- sort(Sys.glob(chain_pattern))
  files <- files[!grepl("SMOKE|smoke", files)]
  if (length(files) < 3) return(NULL)
  read_chain_draws(files[1:3])
}

posterior_summary <- function(draws, parameters) {
  keep <- intersect(parameters, colnames(draws))
  map_dfr(
    keep,
    function(parameter) {
      x <- draws[, parameter]
      tibble(
        parameter = parameter,
        median = median(x, na.rm = TRUE),
        lower = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
        upper = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
      )
    }
  )
}

primary <- read_model_draws(
  c(
    "Model_TimeVarying_Feeding_Matthew_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
    "ModelRuns/TimeVarying_Feeding_Matthew_IndivDetect_LONG/Model_TimeVarying_Feeding_Matthew_IndivDetect_LONG_3CHAINS_COMPLETE.rds"
  ),
  "ModelRuns/TimeVarying_Feeding_Matthew_IndivDetect_LONG/*CHAIN*.rds"
)

if (is.null(primary)) stop("Could not find final primary posterior.")

base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey30"),
    axis.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

## ============================================================
## FIGURE S1. WITHIN-SITE DENSITY HISTORIES
## ============================================================

density_file <- "DataAudits/ModelCovariates/01_Density_SiteYearModelTable.csv"

if (file.exists(density_file)) {
  density <- read_csv(density_file, show_col_types = FALSE) %>%
    mutate(
      observed = density_source == "observed_survey_year"
    )

  fig_s1 <- ggplot(
    density,
    aes(x = calendar_year, y = density_within_z_model)
  ) +
    geom_hline(yintercept = 0, colour = "grey75") +
    geom_line(linewidth = 0.7) +
    geom_point(
      data = density %>% filter(observed),
      size = 1.6
    ) +
    facet_wrap(~ island, ncol = 3, scales = "free_x") +
    labs(
      x = "Calendar year",
      y = "Within-site density (SD)",
      title = "Observed and interpolated within-site density histories",
      subtitle = "Points are observed survey years; lines between observations are interpolated. Outside observed support the model covariate is the site mean (0)."
    ) +
    base_theme

  save_figure(fig_s1, "Figure_S1_DensityHistories", width = 7.6, height = 7.2)
}

## ============================================================
## FIGURE S2. CORE SURVIVAL COEFFICIENTS
## ============================================================

survival_labels <- c(
  "betaSEX[1]" = "Male effect on age-10 mortality",
  "betaSEX[2]" = "Male effect on senescence rate",
  "betaDENS_a" = "Density effect on age-10 mortality",
  "betaDENS_b" = "Density effect on senescence rate",
  "betaCOHORT" = "Birth-cohort effect",
  "betaFEED" = "High/regular provisioning",
  "betaMATTHEW" = "Post-Matthew FFRC period"
)

survival_post <- posterior_summary(primary, names(survival_labels)) %>%
  mutate(
    label = recode(parameter, !!!survival_labels),
    label = factor(label, levels = rev(unname(survival_labels)))
  )

fig_s2 <- ggplot(survival_post, aes(x = median, y = label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18) +
  geom_point(size = 2.5) +
  labs(
    x = "Posterior coefficient (model scale)",
    y = NULL,
    title = "Core survival-model coefficients",
    subtitle = "Posterior median and 95% credible interval."
  ) +
  base_theme

save_figure(fig_s2, "Figure_S2_CoreSurvivalCoefficients", width = 7.4, height = 4.8)

## ============================================================
## FIGURE S3. DETECTION MODEL PARAMETERS
## ============================================================

detection_labels <- c(
  "alpha_p" = "Detection intercept",
  "betaP_SEX" = "Male detection effect",
  "betaP_AGE" = "Age detection effect",
  "sigma_island_p" = "Site detection SD",
  "sigma_survey_p" = "Survey-event detection SD",
  "sigma_indiv_p" = "Individual detection SD"
)

detection_post <- posterior_summary(primary, names(detection_labels)) %>%
  mutate(
    label = recode(parameter, !!!detection_labels),
    label = factor(label, levels = rev(unname(detection_labels)))
  )

fig_s3 <- ggplot(detection_post, aes(x = median, y = label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18) +
  geom_point(size = 2.5) +
  labs(
    x = "Posterior parameter value",
    y = NULL,
    title = "Detection-model parameters",
    subtitle = "The persistent individual detection SD is estimated alongside site and survey-event heterogeneity."
  ) +
  base_theme

save_figure(fig_s3, "Figure_S3_DetectionParameters", width = 7.4, height = 4.6)

## ============================================================
## FIGURE S4. MAGNITUDE OF INDIVIDUAL DETECTION HETEROGENEITY
## ============================================================

if ("sigma_indiv_p" %in% colnames(primary)) {
  indiv_or <- tibble(
    odds_ratio = exp(primary[, "sigma_indiv_p"])
  )

  fig_s4 <- ggplot(indiv_or, aes(x = odds_ratio)) +
    geom_density(fill = "grey80", colour = "black", linewidth = 0.7) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    labs(
      x = "Detection odds ratio for +1 SD individual effect",
      y = "Posterior density",
      title = "Persistent heterogeneity in individual detectability",
      subtitle = "A one-SD increase in the individual random effect multiplies detection odds by exp(sigma_indiv_p)."
    ) +
    base_theme

  save_figure(fig_s4, "Figure_S4_IndividualDetectionHeterogeneity", width = 6.8, height = 4.6)
}

## ============================================================
## FIGURE S5. LATENT SEX POSTERIOR PROBABILITIES
## ============================================================

sex_cols <- grep("^sex\\[[0-9]+\\]$", colnames(primary), value = TRUE)

if (length(sex_cols) > 0) {
  sex_prob <- tibble(
    node = sex_cols,
    p_male = colMeans(primary[, sex_cols, drop = FALSE])
  )

  fig_s5 <- ggplot(sex_prob, aes(x = p_male)) +
    geom_histogram(binwidth = 0.05, boundary = 0, closed = "left") +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
    labs(
      x = "Posterior probability male",
      y = "Latent-sex individuals",
      title = "Posterior resolution of latent sex",
      subtitle = "Only sex nodes monitored as latent in the final model are shown."
    ) +
    base_theme

  save_figure(fig_s5, "Figure_S5_LatentSexProbabilities", width = 6.8, height = 4.6)
}

## ============================================================
## FIGURE S6. CORE PARAMETER ROBUSTNESS ACROSS SENSITIVITIES
## Added automatically once all sensitivity outputs exist.
## ============================================================

model_specs <- tribble(
  ~model, ~combined_1, ~combined_2, ~chain_pattern,
  "Primary",
  "Model_TimeVarying_Feeding_Matthew_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/TimeVarying_Feeding_Matthew_IndivDetect_LONG/Model_TimeVarying_Feeding_Matthew_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/TimeVarying_Feeding_Matthew_IndivDetect_LONG/*CHAIN*.rds",
  "Any provisioning",
  "Model_Sensitivity_ANY_FEED_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_ANY_FEED_IndivDetect_LONG/Model_Sensitivity_ANY_FEED_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_ANY_FEED_IndivDetect_LONG/*CHAIN*.rds",
  "Acute Matthew",
  "Model_Sensitivity_MATTHEW_ACUTE_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_MATTHEW_ACUTE_IndivDetect_LONG/Model_Sensitivity_MATTHEW_ACUTE_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_MATTHEW_ACUTE_IndivDetect_LONG/*CHAIN*.rds",
  "Partial density",
  "Model_Sensitivity_DENSITY_PARTIAL_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_DENSITY_PARTIAL_IndivDetect_LONG/Model_Sensitivity_DENSITY_PARTIAL_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_DENSITY_PARTIAL_IndivDetect_LONG/*CHAIN*.rds"
)

robust_parameters <- c(
  "betaSEX[1]",
  "betaSEX[2]",
  "betaDENS_a",
  "betaDENS_b",
  "betaCOHORT",
  "betaFEED",
  "sigma_indiv_p"
)

robust <- pmap_dfr(
  model_specs,
  function(model, combined_1, combined_2, chain_pattern) {
    draws <- read_model_draws(
      c(combined_1, combined_2),
      chain_pattern
    )
    if (is.null(draws)) return(tibble())
    posterior_summary(draws, robust_parameters) %>%
      mutate(model = model)
  }
)

if (n_distinct(robust$model) > 1) {
  robust_labels <- c(
    "betaSEX[1]" = "Male: age-10 mortality",
    "betaSEX[2]" = "Male: senescence",
    "betaDENS_a" = "Density: age-10 mortality",
    "betaDENS_b" = "Density: senescence",
    "betaCOHORT" = "Birth cohort",
    "betaFEED" = "Provisioning",
    "sigma_indiv_p" = "Individual detection SD"
  )

  robust <- robust %>%
    mutate(
      parameter_label = recode(parameter, !!!robust_labels),
      model = factor(model, levels = model_specs$model)
    )

  fig_s6 <- ggplot(
    robust,
    aes(x = median, y = model)
  ) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey65") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18) +
    geom_point(size = 2.0) +
    facet_wrap(~ parameter_label, scales = "free_x", ncol = 2) +
    labs(
      x = "Posterior parameter value",
      y = NULL,
      title = "Robustness of core conclusions across pre-specified sensitivities",
      subtitle = "The individual detection random effect is retained in every model."
    ) +
    base_theme

  save_figure(fig_s6, "Figure_S6_SensitivityParameterRobustness", width = 7.6, height = 8.0)

  write_csv(
    robust,
    file.path(out_dir, "Figure_S6_SensitivityParameterRobustness_values.csv")
  )
}

cat("\n============================================================\n")
cat("SUPPLEMENTARY FIGURES COMPLETE\n")
cat("============================================================\n")
cat("Output directory: ", out_dir, "\n", sep = "")
