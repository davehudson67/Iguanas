## ============================================================
## 51_Final_Robustness_Synthesis.R
##
## Final synthesis of the frozen audited primary model,
## pre-specified sensitivities, targeted two-site exclusion,
## and all leave-one-site-out fits.
##
## Outputs:
##   DataAudits/Final_Robustness/Final_Robustness_Summary.csv
##   DataAudits/Final_Robustness/Final_Robustness_KeyEffects_Long.csv
##   Figures/Final_Robustness/Final_Robustness_KeyEffects.pdf
##   Figures/Final_Robustness/Final_Robustness_KeyEffects.png
## ============================================================

library(tidyverse)

dir.create(
  "DataAudits/Final_Robustness",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Figures/Final_Robustness",
  recursive = TRUE,
  showWarnings = FALSE
)

primary_file <- "Model_Final_IndivDetect_Audited_3CHAINS_COMPLETE.rds"

sensitivity_files <- c(
  ANY_FEED =
    "Model_Sensitivity_ANY_FEED_AUDITED_LONG_3CHAINS_COMPLETE.rds",
  DENSITY_PARTIAL =
    "Model_Sensitivity_DENSITY_PARTIAL_AUDITED_LONG_3CHAINS_COMPLETE.rds",
  MATTHEW_ACUTE =
    "Model_Sensitivity_MATTHEW_ACUTE_AUDITED_LONG_3CHAINS_COMPLETE.rds"
)

two_site_file <-
  "Model_Sensitivity_Exclude_LeafNorth_U_Cay_3CHAINS_COMPLETE.rds"

if(!file.exists(two_site_file)) {
  two_site_file <-
    "Model_Sensitivity_Exclude_LeafNorth_UCay_3CHAINS_COMPLETE.rds"
}

loo_sites <- c(
  "bitter",
  "ffrc",
  "gaulin",
  "leafnorth",
  "leafsouth",
  "nadder",
  "noddy",
  "ucay",
  "whiteb"
)

loo_files <- setNames(
  paste0(
    "Model_LOO_",
    loo_sites,
    "_3CHAINS_COMPLETE.rds"
  ),
  loo_sites
)

all_required <- c(
  primary_file,
  unname(sensitivity_files),
  two_site_file,
  unname(loo_files)
)

missing_files <- all_required[!file.exists(all_required)]

if(length(missing_files) > 0) {
  stop(
    "Missing required combined model file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

read_fit <- function(file, scenario, class, omitted_site = NA_character_) {
  x <- readRDS(file)

  if(is.null(x$summary)) {
    stop("No $summary element in: ", file)
  }

  if(is.null(x$calendar_effect_summary)) {
    stop("No $calendar_effect_summary element in: ", file)
  }

  max_rhat <- if(
    is.null(x$summary_all) ||
      !"Rhat" %in% names(x$summary_all)
  ) {
    NA_real_
  } else {
    suppressWarnings(
      max(
        x$summary_all$Rhat,
        na.rm = TRUE
      )
    )
  }

  list(
    object = x,
    scenario = scenario,
    class = class,
    omitted_site = omitted_site,
    max_rhat = max_rhat
  )
}

fits <- list(
  read_fit(
    primary_file,
    "Primary",
    "Primary"
  ),
  read_fit(
    sensitivity_files[["ANY_FEED"]],
    "Sensitivity: any provisioning",
    "Sensitivity"
  ),
  read_fit(
    sensitivity_files[["DENSITY_PARTIAL"]],
    "Sensitivity: partial density support",
    "Sensitivity"
  ),
  read_fit(
    sensitivity_files[["MATTHEW_ACUTE"]],
    "Sensitivity: acute Matthew period",
    "Sensitivity"
  ),
  read_fit(
    two_site_file,
    "Exclude Leaf North + U Cay",
    "Targeted exclusion"
  )
)

loo_labels <- c(
  bitter = "Bitter",
  ffrc = "FFRC",
  gaulin = "Gaulin",
  leafnorth = "Leaf North",
  leafsouth = "Leaf South",
  nadder = "N Adder",
  noddy = "Noddy",
  ucay = "U Cay",
  whiteb = "White Bay"
)

for(site in loo_sites) {
  fits[[length(fits) + 1L]] <- read_fit(
    loo_files[[site]],
    paste0("LOO: ", loo_labels[[site]]),
    "Leave-one-site-out",
    loo_labels[[site]]
  )
}

extract_parameter <- function(fit, parameter) {
  z <- fit$object$summary %>%
    filter(.data$parameter == !!parameter)

  if(nrow(z) != 1) {
    stop(
      "Expected exactly one row for ",
      parameter,
      " in ",
      fit$scenario,
      "; found ",
      nrow(z)
    )
  }

  tibble(
    scenario = fit$scenario,
    class = fit$class,
    omitted_site = fit$omitted_site,
    parameter = parameter,
    median = z[["50%"]],
    low = z[["2.5%"]],
    high = z[["97.5%"]],
    Rhat = z$Rhat,
    max_Rhat = fit$max_rhat
  )
}

key_parameters <- c(
  "betaSEX[1]",
  "betaSEX[2]",
  "betaDENS_a",
  "betaDENS_b",
  "betaCOHORT"
)

parameter_long <- map_dfr(
  fits,
  function(fit) {
    map_dfr(
      key_parameters,
      function(p) extract_parameter(fit, p)
    )
  }
)

extract_calendar <- function(fit) {
  z <- fit$object$calendar_effect_summary

  if(nrow(z) < 2) {
    stop(
      "Expected at least two calendar-effect rows in ",
      fit$scenario
    )
  }

  tibble(
    scenario = fit$scenario,
    class = fit$class,
    omitted_site = fit$omitted_site,
    effect = c("Provisioning HR", "Matthew HR"),
    median = z$median_HR[1:2],
    low = z$low_HR[1:2],
    high = z$high_HR[1:2],
    Pr_gt_1 = z$Pr_HR_gt_1[1:2],
    max_Rhat = fit$max_rhat
  )
}

calendar_long <- map_dfr(
  fits,
  extract_calendar
)

summary_table <- parameter_long %>%
  select(
    scenario,
    class,
    omitted_site,
    parameter,
    median,
    low,
    high,
    Rhat,
    max_Rhat
  ) %>%
  pivot_wider(
    names_from = parameter,
    values_from = c(median, low, high, Rhat),
    names_glue = "{parameter}_{.value}"
  ) %>%
  left_join(
    calendar_long %>%
      select(
        scenario,
        effect,
        median,
        low,
        high,
        Pr_gt_1
      ) %>%
      pivot_wider(
        names_from = effect,
        values_from = c(median, low, high, Pr_gt_1),
        names_glue = "{effect}_{.value}"
      ),
    by = "scenario"
  )

write_csv(
  summary_table,
  "DataAudits/Final_Robustness/Final_Robustness_Summary.csv"
)

effect_labels <- c(
  "betaSEX[1]" =
    "Male effect on mortality at age 10 (log HR)",
  "betaSEX[2]" =
    "Male effect on Gompertz senescence rate (log ratio)",
  "betaDENS_a" =
    "Density effect on age-10 mortality (log HR)",
  "betaDENS_b" =
    "Density effect on senescence rate (log ratio)",
  "betaCOHORT" =
    "Birth-cohort effect on age-10 mortality"
)

plot_parameter_data <- parameter_long %>%
  mutate(
    effect = recode(
      parameter,
      !!!effect_labels
    ),
    scale = "Zero = no effect"
  ) %>%
  select(
    scenario,
    class,
    effect,
    scale,
    median,
    low,
    high
  )

plot_calendar_data <- calendar_long %>%
  mutate(
    scale = "One = no effect"
  ) %>%
  select(
    scenario,
    class,
    effect,
    scale,
    median,
    low,
    high
  )

key_effects_long <- bind_rows(
  plot_parameter_data,
  plot_calendar_data
)

write_csv(
  key_effects_long,
  "DataAudits/Final_Robustness/Final_Robustness_KeyEffects_Long.csv"
)

scenario_levels <- rev(
  c(
    "Primary",
    "Sensitivity: any provisioning",
    "Sensitivity: partial density support",
    "Sensitivity: acute Matthew period",
    "Exclude Leaf North + U Cay",
    paste0(
      "LOO: ",
      unname(loo_labels)
    )
  )
)

key_effects_long <- key_effects_long %>%
  mutate(
    scenario = factor(
      scenario,
      levels = scenario_levels
    ),
    class = factor(
      class,
      levels = c(
        "Primary",
        "Sensitivity",
        "Targeted exclusion",
        "Leave-one-site-out"
      )
    )
  )

plot_effects <- c(
  "Male effect on mortality at age 10 (log HR)",
  "Male effect on Gompertz senescence rate (log ratio)",
  "Density effect on age-10 mortality (log HR)",
  "Density effect on senescence rate (log ratio)",
  "Provisioning HR"
)

plot_data <- key_effects_long %>%
  filter(effect %in% plot_effects) %>%
  mutate(
    effect = factor(
      effect,
      levels = plot_effects
    )
  )

p <- ggplot(
  plot_data,
  aes(
    x = median,
    y = scenario,
    shape = class
  )
) +
  geom_vline(
    data = tibble(
      effect = factor(
        plot_effects[1:4],
        levels = plot_effects
      ),
      reference = 0
    ),
    aes(xintercept = reference),
    linewidth = 0.35,
    linetype = 2,
    inherit.aes = FALSE
  ) +
  geom_vline(
    data = tibble(
      effect = factor(
        "Provisioning HR",
        levels = plot_effects
      ),
      reference = 1
    ),
    aes(xintercept = reference),
    linewidth = 0.35,
    linetype = 2,
    inherit.aes = FALSE
  ) +
  geom_errorbarh(
    aes(
      xmin = low,
      xmax = high
    ),
    height = 0,
    linewidth = 0.5
  ) +
  geom_point(
    size = 2.0
  ) +
  facet_wrap(
    ~ effect,
    scales = "free_x",
    ncol = 2
  ) +
  labs(
    x = "Posterior median and 95% credible interval",
    y = NULL,
    shape = NULL,
    title = "Robustness of key survival-model effects",
    subtitle =
      "Frozen primary model, pre-specified sensitivities, targeted exclusion, and leave-one-site-out fits"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  "Figures/Final_Robustness/Final_Robustness_KeyEffects.pdf",
  p,
  width = 10.5,
  height = 12,
  units = "in"
)

ggsave(
  "Figures/Final_Robustness/Final_Robustness_KeyEffects.png",
  p,
  width = 10.5,
  height = 12,
  units = "in",
  dpi = 300
)

cat("\n============================================================\n")
cat("FINAL ROBUSTNESS SYNTHESIS COMPLETE\n")
cat("============================================================\n")
cat("Scenarios:", length(fits), "\n")
cat("Primary + sensitivities + targeted exclusion + 9 LOO fits.\n\n")

cat("Male senescence betaSEX[2] across scenarios:\n")
print(
  parameter_long %>%
    filter(parameter == "betaSEX[2]") %>%
    select(
      scenario,
      median,
      low,
      high,
      Rhat
    ),
  n = Inf
)

cat("\nProvisioning HR across scenarios:\n")
print(
  calendar_long %>%
    filter(effect == "Provisioning HR") %>%
    select(
      scenario,
      median,
      low,
      high,
      Pr_gt_1
    ),
  n = Inf
)

cat("\nMaximum monitored Rhat by scenario:\n")
print(
  parameter_long %>%
    distinct(
      scenario,
      max_Rhat
    ),
  n = Inf
)

cat("\nSaved:\n")
cat("  DataAudits/Final_Robustness/Final_Robustness_Summary.csv\n")
cat("  DataAudits/Final_Robustness/Final_Robustness_KeyEffects_Long.csv\n")
cat("  Figures/Final_Robustness/Final_Robustness_KeyEffects.pdf\n")
cat("  Figures/Final_Robustness/Final_Robustness_KeyEffects.png\n")
cat("============================================================\n")
