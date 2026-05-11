## ============================================================
## PLOT ISLAND-SPECIFIC SURVIVAL + MORTALITY (hazard) CURVES
## with SEX differences + FED vs UNFED
library(dplyr)
library(tidyr)
library(ggplot2)
library(MCMCvis)

run_B_test <- readRDS("Model_simplePooled_samples.rds")

## ----------------------------
## 0) Grab posterior draws as a matrix
## ----------------------------
# NIMBLE runMCMC outputs a coda mcmc.list; convert to matrix of draws
samps <- as.matrix(run_B_test$samples)

## ----------------------------
## 1) Island metadata (your supplied scaffold)
## ----------------------------
island_names <- ordered_islands
fed_flag <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
island_info <- tibble(
  island = 1:n_islands,
  island_name = island_names,
  fed = fed_flag,
  fed_label = if_else(fed, "Fed islands", "Unfed islands")
)

## ----------------------------
## 2) Extract island effects + build amult/bmult by sex
## Model:
##   log(amult) = log(a) + betaSEX[1]*sex + betaISLAND_a[island]
##   log(bmult) = log(b) + betaSEX[2]*sex + betaISLAND_b[island]
## ----------------------------
a_base <- samps[, "a"]
b_base <- samps[, "b"]
betaS1 <- samps[, "betaSEX[1]"]
betaS2 <- samps[, "betaSEX[2]"]

betaIsland_a <- sapply(1:n_islands, \(j) samps[, paste0("betaISLAND_a[", j, "]")])
betaIsland_b <- sapply(1:n_islands, \(j) samps[, paste0("betaISLAND_b[", j, "]")])

# amult and bmult posterior draws per island (rows=draws, cols=islands)
amult_f <- sweep(exp(betaIsland_a), 1, a_base, "*")                   # Female sex=0
amult_m <- sweep(exp(betaIsland_a + betaS1), 1, a_base, "*")          # Male sex=1

bmult_f <- sweep(exp(betaIsland_b), 1, b_base, "*")                   # Female sex=0
bmult_m <- sweep(exp(betaIsland_b + betaS2), 1, b_base, "*")          # Male sex=1

## ----------------------------
## 3) Gompertz functions (age = years since birth)
## Survival: S(t) = exp(-(a/b)*(exp(b t)-1)) for b>0
## Hazard:   h(t) = a * exp(b t)
## ----------------------------
S_gomp <- function(a, b, t) {
  # safe for tiny b: fall back to exponential approx when b is extremely small
  eps <- 1e-8
  ifelse(b < eps, exp(-a * t), exp(-(a / b) * (exp(b * t) - 1)))
}

h_gomp <- function(a, b, t) a * exp(b * t)

# quantile summariser across posterior draws
qsum <- function(draws_by_t) {
  apply(draws_by_t, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
}

## ----------------------------
## 4) Create plotting grids + helper to build tidy data
## ----------------------------
tgrid <- seq(0, 45, by = 0.25)   # adjust to taste (0..40 if you prefer)

make_curve_df <- function(amult_draws, bmult_draws, sex_label, curve = c("survival","hazard")) {
  curve <- match.arg(curve)
  out <- vector("list", n_islands)
  
  for (j in 1:n_islands) {
    if (curve == "survival") {
      mat <- sapply(tgrid, \(tt) S_gomp(amult_draws[, j], bmult_draws[, j], tt))
    } else {
      mat <- sapply(tgrid, \(tt) h_gomp(amult_draws[, j], bmult_draws[, j], tt))
    }
    qs <- qsum(mat)  # 3 x length(tgrid)
    
    out[[j]] <- tibble(
      age = tgrid,
      lo  = qs[1, ],
      mid = qs[2, ],
      hi  = qs[3, ],
      island = j,
      sex = sex_label
    )
  }
  
  bind_rows(out) %>%
    left_join(island_info, by = "island") %>%
    mutate(curve = curve)
}

df_surv <- bind_rows(
  make_curve_df(amult_f, bmult_f, "Female", "survival"),
  make_curve_df(amult_m, bmult_m, "Male",   "survival")
)

df_haz <- bind_rows(
  make_curve_df(amult_f, bmult_f, "Female", "hazard"),
  make_curve_df(amult_m, bmult_m, "Male",   "hazard")
)

## ----------------------------
## 5) Plot: ISLAND-SPECIFIC SURVIVAL CURVES (fed vs unfed colours)
## ----------------------------
p_surv <- ggplot(df_surv, aes(age, mid, colour = fed_label, fill = fed_label)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  facet_grid(sex ~ island_name) +
  labs(x = "Age (years)", y = "Survival S(t)", colour = NULL, fill = NULL,
       title = "Island-specific survival curves (Gompertz), by sex",
       subtitle = "Lines = posterior median; ribbons = 95% credible intervals") +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 8),
    strip.text.x = element_text(size = 8)
  )

print(p_surv)

## ----------------------------
## 6) Plot: ISLAND-SPECIFIC HAZARD (mortality) CURVES (fed vs unfed colours)
## ----------------------------
p_haz <- ggplot(df_haz, aes(age, mid, colour = fed_label, fill = fed_label)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  facet_grid(sex ~ island_name) +
  #scale_y_log10() +
  labs(x = "Age (years)", y = "Hazard h(t) (log scale)", colour = NULL, fill = NULL,
       title = "Island-specific mortality (hazard) curves (Gompertz), by sex",
       subtitle = "Lines = posterior median; ribbons = 95% credible intervals; y-axis is log10") +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 8),
    strip.text.x = element_text(size = 8)
  ) +
  coord_cartesian(ylim = c(0, 5))

print(p_haz)

## ----------------------------
## 7) A “group summary” plot (Fed vs Unfed only)
##    (still separated by sex), pooling islands visually
## ----------------------------
# This summarises across islands by taking the median-of-medians across islands within each group,
# and similarly for lo/hi (a simple visual summary; not a formal pooled posterior).
df_surv_group <- df_surv %>%
  group_by(sex, fed_label, age) %>%
  summarise(
    mid = median(mid, na.rm = TRUE),
    lo  = median(lo,  na.rm = TRUE),
    hi  = median(hi,  na.rm = TRUE),
    .groups = "drop"
  )

p_surv_group <- ggplot(df_surv_group, aes(age, mid, colour = fed_label, fill = fed_label)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  facet_wrap(~ sex) +
  labs(x = "Age (years)", y = "Survival S(t)", colour = NULL, fill = NULL,
       title = "Fed vs unfed summary survival curves (visual grouping)",
       subtitle = "Median of island-level summaries within each group") +
  theme_bw() +
  theme(legend.position = "top")

print(p_surv_group)

## ============================================================
## Fed vs Unfed summary *mortality (hazard)* curves
## ============================================================

df_haz_group <- df_haz %>%
  group_by(sex, fed_label, age) %>%
  summarise(
    mid = median(mid, na.rm = TRUE),
    lo  = median(lo,  na.rm = TRUE),
    hi  = median(hi,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lo = pmax(0, lo),
    hi = pmax(lo, hi)
  )

p_haz_group <- ggplot(df_haz_group, aes(age, mid, colour = fed_label, fill = fed_label)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  facet_wrap(~ sex) +
  labs(x = "Age (years)", y = "Hazard h(t)",
       colour = NULL, fill = NULL,
       title = "Fed vs unfed summary mortality (hazard) curves (raw scale)",
       subtitle = "Median of island-level summaries within each group") +
  coord_cartesian(ylim = c(0, 5)) +
  theme_bw() +
  theme(legend.position = "top")

print(p_haz_group)

#==============================================================================
## ----------------------------
## 7 - IMPROVED) Statistically rigorous Group Summary (Fed vs Unfed)
## ----------------------------

# Identify which columns (islands) belong to which group
idx_fed <- which(fed_flag == TRUE)
idx_unfed <- which(fed_flag == FALSE)

# Helper function to average curves across a specific group of islands 
# within each MCMC draw, THEN calculate the quantiles.
make_pooled_curve_df <- function(amult_list, bmult_list, idx_group, group_label, sex_label, curve = c("survival", "hazard")) {
  curve <- match.arg(curve)
  
  # Evaluate curve for ALL islands in the group, resulting in an array: [MCMC draws, Timepoints, Islands]
  mat_list <- lapply(idx_group, function(j) {
    if (curve == "survival") {
      sapply(tgrid, \(tt) S_gomp(amult_list[, j], bmult_list[, j], tt))
    } else {
      sapply(tgrid, \(tt) h_gomp(amult_list[, j], bmult_list[, j], tt))
    }
  })
  
  # Convert list of matrices into a 3D array: dim = c(n_draws, n_timepoints, n_islands_in_group)
  arr <- array(unlist(mat_list), dim = c(nrow(amult_list), length(tgrid), length(idx_group)))
  
  # Calculate the mean curve across the chosen islands for EACH MCMC draw
  # Resulting matrix is [n_draws, n_timepoints]
  mean_curve_mat <- apply(arr, c(1, 2), mean)
  
  # Now take the quantiles over the MCMC draws
  qs <- qsum(mean_curve_mat)
  
  tibble(
    age = tgrid,
    lo  = qs[1, ],
    mid = qs[2, ],
    hi  = qs[3, ],
    fed_label = group_label,
    sex = sex_label
  )
}

# Build the true posterior group summaries
df_surv_group_true <- bind_rows(
  make_pooled_curve_df(amult_f, bmult_f, idx_fed,   "Fed islands",   "Female", "survival"),
  make_pooled_curve_df(amult_f, bmult_f, idx_unfed, "Unfed islands", "Female", "survival"),
  make_pooled_curve_df(amult_m, bmult_m, idx_fed,   "Fed islands",   "Male",   "survival"),
  make_pooled_curve_df(amult_m, bmult_m, idx_unfed, "Unfed islands", "Male",   "survival")
)

# Plot the mathematically rigorous summary
p_surv_group <- ggplot(df_surv_group_true, aes(age, mid, colour = fed_label, fill = fed_label)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  facet_wrap(~ sex) +
  labs(x = "Age (years)", y = "Survival S(t)", colour = NULL, fill = NULL,
       title = "True Bayesian Average Survival Curves (Fed vs Unfed)",
       subtitle = "Averaged across islands within each MCMC step") +
  theme_bw() +
  theme(legend.position = "top")

print(p_surv_group)

## ------------------------------------------------------------
## 7b - IMPROVED) Statistically rigorous Group Summary for HAZARD
## ------------------------------------------------------------

# Build the true posterior group summaries for the hazard curves
df_haz_group_true <- bind_rows(
  make_pooled_curve_df(amult_f, bmult_f, idx_fed,   "Fed islands",   "Female", "hazard"),
  make_pooled_curve_df(amult_f, bmult_f, idx_unfed, "Unfed islands", "Female", "hazard"),
  make_pooled_curve_df(amult_m, bmult_m, idx_fed,   "Fed islands",   "Male",   "hazard"),
  make_pooled_curve_df(amult_m, bmult_m, idx_unfed, "Unfed islands", "Male",   "hazard")
)

# Optional: Set a logical y-axis limit for the raw scale plot. 
# Gompertz hazards explode at old age, which can squish the interesting early-life 
# differences at the bottom of the plot. We cap it to look at the relevant range.
ymax_haz <- quantile(df_haz_group_true$hi, probs = 0.95, na.rm = TRUE)

# Plot the mathematically rigorous summary on raw scale
p_haz_group_true <- ggplot(df_haz_group_true, aes(age, mid, colour = fed_label, fill = fed_label)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  facet_wrap(~ sex) +
  labs(x = "Age (years)", y = "Hazard h(t)",
       colour = NULL, fill = NULL,
       title = "True Bayesian Average Mortality (Hazard) Curves (Fed vs Unfed)",
       subtitle = "Averaged across islands within each MCMC step") +
  coord_cartesian(ylim = c(0, 1)) + # adjust this 1 to ymax_haz if you want it dynamic
  theme_bw() +
  theme(legend.position = "top")

print(p_haz_group_true)
