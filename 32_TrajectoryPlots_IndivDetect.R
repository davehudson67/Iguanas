## ============================================================
## 32_TrajectoryPlots_IndivDetect.R
##
## Standardised mortality-hazard and age-specific annual survival
## trajectories from the final individual-detection-RE model.
##
## IMPORTANT INTERPRETATION:
## Provisioning is represented exactly as fitted: a CURRENT site-year
## multiplicative effect on mortality hazard. The two plotted exposure
## scenarios hold that fitted indicator fixed across age only to show the
## model-implied contrast. They do NOT imply cumulative carry-over,
## permanent treatment classes, or individual causal effects.
##
## Reference conditions:
##   density = 0 (within-site mean)
##   cohort = 0 (centred reference)
##   Matthew = 0
##   island mortality RE = 0
##   site-period effect = 0
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(coda)
})

MODEL_FILE <- "Model_TimeVarying_Feeding_Matthew_IndivDetect_LONG_3CHAINS_COMPLETE.rds"

OUT_DIR <- "Figures/Manuscript"
dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

N_DRAWS <- 4000L
AGE_GRID <- seq(
  0,
  40,
  by = 0.1
)

## ============================================================
## 1. LOAD POSTERIOR DRAWS
## ============================================================

fit_obj <- readRDS(
  MODEL_FILE
)

if(
  is.null(
    fit_obj$samples
  )
) {
  stop(
    "Combined model object does not contain fit_obj$samples."
  )
}

if(
  inherits(
    fit_obj$samples,
    "mcmc.list"
  )
) {

  draws_all <- do.call(
    rbind,
    lapply(
      fit_obj$samples,
      as.matrix
    )
  )

} else if(
  inherits(
    fit_obj$samples,
    "mcmc"
  )
) {

  draws_all <- as.matrix(
    fit_obj$samples
  )

} else {

  draws_all <- as.matrix(
    fit_obj$samples
  )
}

required_parameters <- c(
  "log_a10",
  "log_b",
  "betaSEX[1]",
  "betaSEX[2]",
  "betaFEED"
)

missing_parameters <- setdiff(
  required_parameters,
  colnames(
    draws_all
  )
)

if(
  length(
    missing_parameters
  ) >
    0
) {
  stop(
    "Posterior is missing: ",
    paste(
      missing_parameters,
      collapse = ", "
    )
  )
}

cat(
  "Posterior matrix: ",
  nrow(draws_all),
  " draws x ",
  ncol(draws_all),
  " parameters\n",
  sep = ""
)

set.seed(
  20261028
)

if(
  nrow(
    draws_all
  ) >
    N_DRAWS
) {

  keep <- sample(
    seq_len(
      nrow(
        draws_all
      )
    ),
    N_DRAWS
  )

  draws <- draws_all[
    keep,
    ,
    drop = FALSE
  ]

} else {

  draws <- draws_all
}

cat(
  "Using ",
  nrow(draws),
  " posterior draws for trajectories.\n",
  sep = ""
)


## ============================================================
## 2. SCENARIOS
##
## sex:
##   0 = female
##   1 = male
##
## provisioning:
##   0 = other site-year periods
##   1 = high/regular provisioning site-year periods
## ============================================================

scenarios <- crossing(
  sex = c(
    "Female",
    "Male"
  ),
  provisioning = c(
    "Other site-year periods",
    "High/regular provisioning periods"
  )
) %>%
  mutate(
    sex_value =
      if_else(
        sex ==
          "Male",
        1,
        0
      ),

    feed_value =
      if_else(
        provisioning ==
          "High/regular provisioning periods",
        1,
        0
      )
  )


## ============================================================
## 3. EXACT GOMPERTZ HAZARD + ONE-YEAR SURVIVAL
##
## Instantaneous hazard at age a:
##
##   h(a) =
##     a10 * exp[b * (a - 10)] * exp(betaFEED * feed)
##
## Integrated hazard over the NEXT one-year interval:
##
##   H(a, a+1) =
##     exp(betaFEED * feed) * (a10 / b) *
##     {exp[b * (a + 1 - 10)] - exp[b * (a - 10)]}
##
## Annual survival at age a:
##
##   S_1yr(a) = exp[-H(a, a+1)]
##
## This is faithful to the fitted CURRENT site-year hazard multiplier
## and does not require an assumed lifelong provisioning history.
## ============================================================

trajectory_one_scenario <- function(
    sex_value,
    feed_value,
    sex_label,
    feed_label
) {

  log_a10_i <-
    draws[
      ,
      "log_a10"
    ] +
    draws[
      ,
      "betaSEX[1]"
    ] *
    sex_value

  log_b_i <-
    draws[
      ,
      "log_b"
    ] +
    draws[
      ,
      "betaSEX[2]"
    ] *
    sex_value

  a10_i <- exp(
    log_a10_i
  )

  b_i <- exp(
    log_b_i
  )

  feed_multiplier <- exp(
    draws[
      ,
      "betaFEED"
    ] *
    feed_value
  )

  hazard_matrix <- vapply(
    AGE_GRID,
    function(age) {

      a10_i *
        exp(
          b_i *
            (
              age -
                10
            )
        ) *
        feed_multiplier
    },
    numeric(
      length(
        a10_i
      )
    )
  )

  survival_matrix <- vapply(
    AGE_GRID,
    function(age) {

      one_year_hazard <-
        feed_multiplier *
        (
          a10_i /
            b_i
        ) *
        (
          exp(
            b_i *
              (
                age +
                  1 -
                  10
              )
          ) -
            exp(
              b_i *
                (
                  age -
                    10
                )
            )
        )

      exp(
        -one_year_hazard
      )
    },
    numeric(
      length(
        a10_i
      )
    )
  )

  hazard_summary <- tibble(
    age =
      AGE_GRID,

    median =
      apply(
        hazard_matrix,
        2,
        median
      ),

    lower =
      apply(
        hazard_matrix,
        2,
        quantile,
        probs =
          0.025
      ),

    upper =
      apply(
        hazard_matrix,
        2,
        quantile,
        probs =
          0.975
      ),

    metric =
      "Mortality hazard",

    sex =
      sex_label,

    provisioning =
      feed_label
  )

  survival_summary <- tibble(
    age =
      AGE_GRID,

    median =
      apply(
        survival_matrix,
        2,
        median
      ),

    lower =
      apply(
        survival_matrix,
        2,
        quantile,
        probs =
          0.025
      ),

    upper =
      apply(
        survival_matrix,
        2,
        quantile,
        probs =
          0.975
      ),

    metric =
      "Annual survival probability",

    sex =
      sex_label,

    provisioning =
      feed_label
  )

  bind_rows(
    hazard_summary,
    survival_summary
  )
}


trajectories <- pmap_dfr(
  scenarios,
  function(
      sex,
      provisioning,
      sex_value,
      feed_value
  ) {

    trajectory_one_scenario(
      sex_value =
        sex_value,

      feed_value =
        feed_value,

      sex_label =
        sex,

      feed_label =
        provisioning
    )
  }
) %>%
  mutate(
    metric =
      factor(
        metric,
        levels =
          c(
            "Mortality hazard",
            "Annual survival probability"
          )
      ),

    sex =
      factor(
        sex,
        levels =
          c(
            "Female",
            "Male"
          )
      ),

    provisioning =
      factor(
        provisioning,
        levels =
          c(
            "Other site-year periods",
            "High/regular provisioning periods"
          )
      )
  )


## ============================================================
## 4. SAVE PLOTTED VALUES
## ============================================================

write_csv(
  trajectories,
  file.path(
    OUT_DIR,
    "Figure_Trajectories_model_implied_values.csv"
  )
)


## ============================================================
## 5. COMBINED MANUSCRIPT FIGURE
##
## Rows = mortality hazard / survival
## Columns = female / male
## Colour = fitted provisioning contrast
##
## Mortality is deliberately on the natural scale so the Gompertz
## curvature remains visually apparent.
## ============================================================

plot_data <- trajectories %>%
  mutate(
    lower_plot =
      if_else(
        metric ==
          "Annual survival probability",
        pmax(
          0,
          lower
        ),
        lower
      ),

    upper_plot =
      if_else(
        metric ==
          "Annual survival probability",
        pmin(
          1,
          upper
        ),
        upper
      )
  )

trajectory_figure <- ggplot(
  plot_data,
  aes(
    x =
      age,
    y =
      median,
    colour =
      provisioning,
    fill =
      provisioning
  )
) +
  geom_ribbon(
    aes(
      ymin =
        lower_plot,
      ymax =
        upper_plot
    ),
    alpha =
      0.14,
    colour =
      NA
  ) +
  geom_line(
    linewidth =
      1
  ) +
  facet_grid(
    metric ~ sex,
    scales =
      "free_y"
  ) +
  scale_colour_brewer(
    palette =
      "Dark2",
    name =
      "Site-year provisioning state"
  ) +
  scale_fill_brewer(
    palette =
      "Dark2",
    name =
      "Site-year provisioning state"
  ) +
  scale_x_continuous(
    breaks =
      seq(
        0,
        40,
        5
      )
  ) +
  labs(
    x =
      "Age (years)",

    y =
      NULL,

    title =
      "Mortality hazard and annual survival across age",

    subtitle =
      NULL,

    caption =
      paste0(
        "Provisioning is applied as the current site-year hazard multiplier fitted by the model. ",
        "Annual survival is survival over the next year at each age; no lifelong exposure history is assumed. ",
        "Density and cohort are held at reference values, with Matthew and random effects set to zero."
      )
  ) +
  theme_classic(
    base_size =
      11
  ) +
  theme(
    plot.title =
      element_text(
        face =
          "bold"
      ),

    plot.subtitle =
      element_text(
        colour =
          "grey30"
      ),

    strip.background =
      element_blank(),

    strip.text =
      element_text(
        face =
          "bold",
        size =
          11
      ),

    legend.position =
      "bottom",

    legend.title =
      element_text(
        face =
          "bold"
      ),

    plot.caption =
      element_text(
        hjust =
          0,
        size =
          8.5,
        colour =
          "grey35"
      )
  )


ggsave(
  file.path(
    OUT_DIR,
    "Figure_Mortality_Survival_Trajectories.pdf"
  ),
  trajectory_figure,
  width =
    8.2,
  height =
    7.2,
  units =
    "in"
)

ggsave(
  file.path(
    OUT_DIR,
    "Figure_Mortality_Survival_Trajectories.tiff"
  ),
  trajectory_figure,
  width =
    8.2,
  height =
    7.2,
  units =
    "in",
  dpi =
    600,
  compression =
    "lzw"
)


## ============================================================
## 6. SIMPLE CHECKS
## ============================================================

age10_check <- trajectories %>%
  filter(
    metric ==
      "Mortality hazard",
    abs(
      age -
        10
    ) <
      1e-8
  ) %>%
  select(
    sex,
    provisioning,
    median,
    lower,
    upper
  )

cat(
  "\n============================================================\n"
)

cat(
  "TRAJECTORY FIGURE COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "Posterior matrix: ",
  nrow(draws_all),
  " x ",
  ncol(draws_all),
  "\n",
  sep = ""
)

cat(
  "Plot draws: ",
  nrow(draws),
  "\n",
  sep = ""
)

cat(
  "\nAge-10 mortality hazard check:\n"
)

print(
  age10_check,
  n =
    Inf
)

cat(
  "\nSaved:\n"
)

cat(
  "  ",
  file.path(
    OUT_DIR,
    "Figure_Mortality_Survival_Trajectories.pdf"
  ),
  "\n",
  sep = ""
)

cat(
  "  ",
  file.path(
    OUT_DIR,
    "Figure_Mortality_Survival_Trajectories.tiff"
  ),
  "\n",
  sep = ""
)

cat(
  "  ",
  file.path(
    OUT_DIR,
    "Figure_Trajectories_model_implied_values.csv"
  ),
  "\n",
  sep = ""
)
