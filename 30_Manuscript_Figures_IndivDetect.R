## ============================================================
## 30_Manuscript_Figures_IndivDetect.R
##
## Main manuscript figures for the final individual-detection-RE model.
##
## Outputs:
##   Figure_1_ExposureHistory
##   Figure_2_SexMortalityTrajectories
##   Figure_3_ProvisioningSensitivity
##   Figure_4_DensityMortalityTrajectories
##
## Sensitivity estimates are added automatically when their completed
## chain/combined outputs are available.
## ============================================================

library(tidyverse)
library(coda)

set.seed(20261027)

out_dir <- "Figures/Manuscript"
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
  if (inherits(x, "mcmc.list")) {
    return(
      do.call(
        rbind,
        lapply(
          x,
          as.matrix
        )
      )
    )
  }

  if (inherits(x, "mcmc")) {
    return(as.matrix(x))
  }

  if (is.matrix(x) || is.data.frame(x)) {
    return(as.matrix(x))
  }

  if (!is.list(x)) {
    return(NULL)
  }

  preferred <- c(
    "samples",
    "combined_samples",
    "posterior_samples",
    "draws",
    "mcmc"
  )

  for (nm in preferred) {
    if (!is.null(x[[nm]])) {
      ans <- extract_draw_matrix(x[[nm]])
      if (!is.null(ans) && ncol(ans) > 0) {
        return(ans)
      }
    }
  }

  NULL
}

read_chain_draws <- function(files) {
  mats <- map(
    files,
    function(f) {
      obj <- readRDS(f)
      ans <- extract_draw_matrix(obj)
      if (is.null(ans)) {
        stop("Could not find posterior samples in: ", f)
      }
      ans
    }
  )

  common <- reduce(
    map(mats, colnames),
    intersect
  )

  if (length(common) == 0) {
    stop("No common posterior columns across chain files.")
  }

  bind_rows(
    map(
      mats,
      ~ as_tibble(.x[, common, drop = FALSE])
    )
  ) %>%
    as.matrix()
}

find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

read_model_draws <- function(
    combined_candidates,
    chain_pattern
) {
  combined_file <- find_first_existing(combined_candidates)

  if (!is.null(combined_file)) {
    obj <- readRDS(combined_file)
    ans <- extract_draw_matrix(obj)
    if (!is.null(ans) && ncol(ans) > 0) {
      message("Using combined posterior: ", combined_file)
      return(ans)
    }
  }

  chain_files <- Sys.glob(chain_pattern)
  chain_files <- chain_files[!grepl("SMOKE|smoke", chain_files)]
  chain_files <- sort(chain_files)

  if (length(chain_files) < 3) {
    return(NULL)
  }

  message("Using chain posterior files: ", paste(basename(chain_files), collapse = ", "))
  read_chain_draws(chain_files[1:3])
}

q_summary <- function(x) {
  tibble(
    median = median(x, na.rm = TRUE),
    lower = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    upper = quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )
}

require_cols <- function(draws, cols) {
  missing <- setdiff(cols, colnames(draws))
  if (length(missing) > 0) {
    stop(
      "Posterior is missing required parameter(s): ",
      paste(missing, collapse = ", ")
    )
  }
}

primary_draws <- read_model_draws(
  combined_candidates = c(
    "Model_TimeVarying_Feeding_Matthew_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
    "ModelRuns/TimeVarying_Feeding_Matthew_IndivDetect_LONG/Model_TimeVarying_Feeding_Matthew_IndivDetect_LONG_3CHAINS_COMPLETE.rds"
  ),
  chain_pattern = "ModelRuns/TimeVarying_Feeding_Matthew_IndivDetect_LONG/*CHAIN*.rds"
)

if (is.null(primary_draws)) {
  stop("Could not find the final primary individual-detection posterior.")
}

require_cols(
  primary_draws,
  c(
    "log_a10",
    "log_b",
    "betaSEX[1]",
    "betaSEX[2]",
    "betaDENS_a",
    "betaDENS_b",
    "betaFEED"
  )
)

plot_draws <- primary_draws[
  sample(
    seq_len(nrow(primary_draws)),
    min(2500, nrow(primary_draws))
  ),
  ,
  drop = FALSE
]

base_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey30"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

## ============================================================
## FIGURE 1. SITE-YEAR PROVISIONING HISTORY + MONITORING
## ============================================================

sites <- c(
  "Bitter",
  "WhiteB",
  "Gaulin",
  "LeafSouth",
  "NAdder",
  "Noddy",
  "LeafNorth",
  "UCay",
  "FFRC"
)

species <- c(
  Bitter = "C. figginsi",
  WhiteB = "C. figginsi",
  Gaulin = "C. figginsi",
  LeafSouth = "C. figginsi",
  NAdder = "C. figginsi",
  Noddy = "C. figginsi",
  LeafNorth = "C. inornata",
  UCay = "C. inornata",
  FFRC = "C. inornata"
)

exposure_intervals <- tribble(
  ~Site, ~start_year, ~end_year, ~state,
  "Bitter", 1980, 2005, "None",
  "Bitter", 2006, 2019, "High / regular",
  "WhiteB", 1980, 1997, "None",
  "WhiteB", 1998, 2019, "Low / intermittent",
  "Gaulin", 1980, 2019, "Low / intermittent",
  "LeafSouth", 1980, 2011, "None",
  "LeafSouth", 2012, 2019, "High / regular",
  "NAdder", 1980, 2019, "None",
  "Noddy", 1980, 2019, "None",
  "LeafNorth", 1980, 2019, "High / regular",
  "UCay", 1980, 2019, "High / regular",
  "FFRC", 1980, 2019, "None"
)

exposure_years <- exposure_intervals %>%
  rowwise() %>%
  mutate(calendar_year = list(seq.int(start_year, end_year))) %>%
  unnest(calendar_year) %>%
  ungroup() %>%
  mutate(
    Site = factor(Site, levels = rev(sites)),
    state = factor(
      state,
      levels = c("None", "Low / intermittent", "High / regular")
    )
  )

monitor_file <- "DataAudits/ModelCovariates/01_MNA_ObservedSurveyYears.csv"

monitor_years <- if (file.exists(monitor_file)) {
  read_csv(monitor_file, show_col_types = FALSE) %>%
    transmute(
      Site = factor(island, levels = rev(sites)),
      calendar_year
    ) %>%
    filter(!is.na(Site))
} else {
  tibble(
    Site = factor(character(), levels = rev(sites)),
    calendar_year = integer()
  )
}

site_labels <- setNames(
  paste0(sites, "  (", species[sites], ")"),
  sites
)

fig1 <- ggplot(
  exposure_years,
  aes(x = calendar_year, y = Site, fill = state)
) +
  geom_tile(height = 0.72, width = 1.02) +
  geom_point(
    data = monitor_years,
    aes(x = calendar_year, y = Site),
    inherit.aes = FALSE,
    shape = 124,
    size = 3.0,
    colour = "black"
  ) +
  scale_y_discrete(labels = site_labels) +
  scale_fill_brewer(
    palette = "YlOrRd",
    direction = 1,
    name = "Provisioning history"
  ) +
  scale_x_continuous(
    breaks = seq(1980, 2020, 5),
    limits = c(1979.5, 2019.5),
    expand = c(0, 0)
  ) +
  labs(
    x = "Calendar year",
    y = NULL,
    title = "Study populations and site-level provisioning history",
    subtitle = "Black ticks mark years with field monitoring; exposure is a site-by-calendar-year history, not an individual treatment.",
    caption = "Primary exposure contrasts high/regular provisioning against all other periods."
  ) +
  base_theme +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 9)
  )

save_figure(fig1, "Figure_1_ExposureHistory", width = 7.6, height = 5.2)

## ============================================================
## FIGURE 2. SEX-SPECIFIC GOMPERTZ MORTALITY TRAJECTORIES
## Reference: density = 0, cohort = 0, no high provisioning,
## no Matthew effect, average site/time random effects.
## ============================================================

ages <- seq(2, 50, by = 0.5)

sex_curve <- map_dfr(
  c(Female = 0, Male = 1),
  function(sex_value) {
    log_a10_i <-
      plot_draws[, "log_a10"] +
      plot_draws[, "betaSEX[1]"] * sex_value

    log_b_i <-
      plot_draws[, "log_b"] +
      plot_draws[, "betaSEX[2]"] * sex_value

    a10_i <- exp(log_a10_i)
    b_i <- exp(log_b_i)

    map_dfr(
      ages,
      function(age) {
        hazard <- a10_i * exp(b_i * (age - 10))
        tibble(
          age = age,
          median = median(hazard),
          lower = quantile(hazard, 0.025, names = FALSE),
          upper = quantile(hazard, 0.975, names = FALSE),
          sex = if_else(sex_value == 0, "Female", "Male")
        )
      }
    )
  }
)

fig2 <- ggplot(
  sex_curve,
  aes(x = age, y = median, colour = sex, fill = sex)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.16,
    colour = NA
  ) +
  geom_line(linewidth = 1.05) +
  scale_y_log10() +
  labs(
    x = "Age (years)",
    y = "Annual mortality hazard (log scale)",
    colour = "Sex",
    fill = "Sex",
    title = "Sex-specific actuarial ageing trajectories",
    subtitle = "Posterior median and 95% credible interval at reference density and calendar conditions."
  ) +
  base_theme +
  theme(legend.position = "bottom")

save_figure(fig2, "Figure_2_SexMortalityTrajectories")

## ============================================================
## FIGURE 3. PROVISIONING ASSOCIATION + PRE-SPECIFIED SENSITIVITIES
## ============================================================

model_specs <- tribble(
  ~label, ~combined_1, ~combined_2, ~chain_pattern,
  "Primary: high / regular",
  "Model_TimeVarying_Feeding_Matthew_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/TimeVarying_Feeding_Matthew_IndivDetect_LONG/Model_TimeVarying_Feeding_Matthew_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/TimeVarying_Feeding_Matthew_IndivDetect_LONG/*CHAIN*.rds",
  "Any provisioning vs none",
  "Model_Sensitivity_ANY_FEED_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_ANY_FEED_IndivDetect_LONG/Model_Sensitivity_ANY_FEED_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_ANY_FEED_IndivDetect_LONG/*CHAIN*.rds",
  "Acute Matthew window",
  "Model_Sensitivity_MATTHEW_ACUTE_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_MATTHEW_ACUTE_IndivDetect_LONG/Model_Sensitivity_MATTHEW_ACUTE_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_MATTHEW_ACUTE_IndivDetect_LONG/*CHAIN*.rds",
  "Partial density support",
  "Model_Sensitivity_DENSITY_PARTIAL_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_DENSITY_PARTIAL_IndivDetect_LONG/Model_Sensitivity_DENSITY_PARTIAL_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  "ModelRuns/Sensitivity_DENSITY_PARTIAL_IndivDetect_LONG/*CHAIN*.rds"
)

feeding_results <- pmap_dfr(
  model_specs,
  function(label, combined_1, combined_2, chain_pattern) {
    draws <- read_model_draws(
      combined_candidates = c(combined_1, combined_2),
      chain_pattern = chain_pattern
    )

    if (is.null(draws) || !("betaFEED" %in% colnames(draws))) {
      return(tibble())
    }

    hr <- exp(draws[, "betaFEED"])
    s <- q_summary(hr)

    s %>%
      mutate(
        label = label,
        prob_hr_gt_1 = mean(hr > 1),
        estimate_text = sprintf(
          "%.2f [%.2f, %.2f]",
          median,
          lower,
          upper
        )
      )
  }
)

feeding_results <- feeding_results %>%
  mutate(
    label = factor(label, levels = rev(model_specs$label))
  )

fig3 <- ggplot(
  feeding_results,
  aes(x = median, y = label)
) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey45") +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.18,
    linewidth = 0.7
  ) +
  geom_point(size = 2.7) +
  geom_text(
    aes(label = estimate_text, x = upper),
    hjust = -0.08,
    size = 3.2
  ) +
  scale_x_log10() +
  coord_cartesian(clip = "off") +
  labs(
    x = "Mortality hazard ratio associated with provisioning",
    y = NULL,
    title = "Provisioning association across the primary model and sensitivities",
    subtitle = "Posterior median and 95% credible interval; HR > 1 indicates higher mortality hazard during the specified provisioning exposure."
  ) +
  base_theme +
  theme(
    plot.margin = margin(5.5, 95, 5.5, 5.5)
  )

save_figure(fig3, "Figure_3_ProvisioningSensitivity", width = 7.6, height = 4.4)

write_csv(
  feeding_results %>%
    mutate(label = as.character(label)),
  file.path(out_dir, "Figure_3_ProvisioningSensitivity_values.csv")
)

## ============================================================
## FIGURE 4. WITHIN-SITE COHORT DENSITY AND PACE OF LIFE
## Female reference trajectory; density is the within-site z score
## used by the model.
## ============================================================

density_levels <- c(
  "Lower density (-1 SD)" = -1,
  "Site mean (0 SD)" = 0,
  "Higher density (+1 SD)" = 1
)

density_curve <- map_dfr(
  density_levels,
  function(density_value) {
    log_a10_i <-
      plot_draws[, "log_a10"] +
      plot_draws[, "betaDENS_a"] * density_value

    log_b_i <-
      plot_draws[, "log_b"] +
      plot_draws[, "betaDENS_b"] * density_value

    a10_i <- exp(log_a10_i)
    b_i <- exp(log_b_i)

    map_dfr(
      ages,
      function(age) {
        hazard <- a10_i * exp(b_i * (age - 10))
        tibble(
          age = age,
          median = median(hazard),
          lower = quantile(hazard, 0.025, names = FALSE),
          upper = quantile(hazard, 0.975, names = FALSE),
          density = names(density_levels)[density_levels == density_value]
        )
      }
    )
  }
)

density_curve <- density_curve %>%
  mutate(
    density = factor(
      density,
      levels = names(density_levels)
    )
  )

fig4 <- ggplot(
  density_curve,
  aes(x = age, y = median, colour = density, fill = density)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.10,
    colour = NA
  ) +
  geom_line(linewidth = 1.0) +
  scale_y_log10() +
  labs(
    x = "Age (years)",
    y = "Annual mortality hazard (log scale)",
    colour = "Birth-period density",
    fill = "Birth-period density",
    title = "Within-site cohort density is associated with mortality pace",
    subtitle = "Female reference trajectories: higher supported birth-period density shifts age-10 mortality upward but senescence acceleration downward."
  ) +
  base_theme +
  theme(legend.position = "bottom")

save_figure(fig4, "Figure_4_DensityMortalityTrajectories", width = 7.4, height = 5.2)

cat("\n============================================================\n")
cat("MAIN MANUSCRIPT FIGURES COMPLETE\n")
cat("============================================================\n")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Sensitivity rows available for Figure 3: ", nrow(feeding_results), " / 4\n", sep = "")
