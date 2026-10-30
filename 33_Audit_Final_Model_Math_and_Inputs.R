## ============================================================
## 33_Audit_Final_Model_Math_and_Inputs.R
##
## Deterministic audit of indexing and mathematics used by the
## final Iguana survival-recapture model. No MCMC is fitted.
## ============================================================

library(tidyverse)

source("DataPrep_Canonical_LatentSex_V5.R")

cov_file <- "Data/Iguana_ModelCovariates_Through2019.rds"
exp_file <- "Data/Iguana_SurveyOccasion_and_CalendarExposure_Audit.rds"

stopifnot(file.exists(cov_file), file.exists(exp_file))

cov_obj <- readRDS(cov_file)
exp_obj <- readRDS(exp_file)

cat("\n============================================================\n")
cat("FINAL MODEL MATH / INPUT AUDIT\n")
cat("============================================================\n")

## ------------------------------------------------------------
## 1. Core dimensions / site order
## ------------------------------------------------------------

stopifnot(
  identical(cov_obj$ordered_islands, ordered_islands),
  cov_obj$Origin == Origin,
  cov_obj$tMax == tMax,
  all(dim(cov_obj$density_matrix) == c(n_islands, tMax)),
  all(dim(exp_obj$feeding_high_matrix) == c(n_islands, tMax)),
  all(dim(exp_obj$feeding_any_matrix) == c(n_islands, tMax)),
  all(dim(exp_obj$matthew_post_matrix) == c(n_islands, tMax)),
  all(dim(exp_obj$matthew_acute_matrix) == c(n_islands, tMax))
)

cat("PASS: dimensions and site ordering are aligned.\n")

## ------------------------------------------------------------
## 2. Birth-time -> density calendar-year index
##
## Model time is:
##   t = calendar_time - Origin + 1
##
## Therefore:
##   birth calendar year = floor(tB + Origin - 1)
##   density matrix index = birth year - Origin + 1
##                        = floor(tB)
## ------------------------------------------------------------

tB_test <- rowMeans(cintB_adj)

birth_year_test <- floor(
  tB_test +
    Origin -
    1
)

expected_idx <-
  birth_year_test -
  Origin +
  1

model_idx <- floor(tB_test)

stopifnot(
  all(expected_idx == model_idx),
  all(model_idx >= 1),
  all(model_idx <= tMax)
)

cat("PASS: floor(tB) is exactly the birth-calendar-year density index.\n")

## ------------------------------------------------------------
## 3. Recompute density-support proportions from continuous
##    half-open birth support [lower, upper).
## ------------------------------------------------------------

density_table <- cov_obj$density_table

support_recalc <- map_dfr(
  seq_len(nind),
  function(i) {

    low_time <-
      cintB_adj[i, 1] +
      Origin -
      1

    high_time <-
      cintB_adj[i, 2] +
      Origin -
      1

    low_year <- floor(low_time)
    high_year <- ceiling(high_time) - 1L

    possible_years <- seq.int(
      low_year,
      high_year
    )

    support <- density_table %>%
      filter(
        island == growth_info$island[i],
        calendar_year %in% possible_years
      )

    tibble(
      ID = i,
      n_possible_birth_years = length(possible_years),
      n_density_supported_years = sum(support$density_supported),
      prop_density_supported =
        n_density_supported_years /
        n_possible_birth_years
    )
  }
)

saved_support <- cov_obj$birth_density_support %>%
  arrange(ID)

support_recalc <- support_recalc %>%
  arrange(ID)

stopifnot(
  identical(saved_support$ID, support_recalc$ID)
)

max_support_diff <- max(
  abs(
    saved_support$prop_density_supported -
      support_recalc$prop_density_supported
  )
)

cat(
  "Max saved-vs-recomputed density-support difference: ",
  signif(max_support_diff, 6),
  "\n",
  sep = ""
)

if(max_support_diff > 1e-12) {
  stop(
    "Density-support object is stale or inconsistent. ",
    "Rerun 05_Build_Density_and_Exposure_Histories.R."
  )
}

cat("PASS: density-support calculation matches continuous birth-time support.\n")

## ------------------------------------------------------------
## 4. Exposure matrices match their site-year tables
## ------------------------------------------------------------

for(r in seq_len(nrow(exp_obj$feeding_site_year))) {

  jj <- match(
    exp_obj$feeding_site_year$Site[r],
    ordered_islands
  )

  tt <- exp_obj$feeding_site_year$year_idx[r]

  expected_high <- replace_na(
    exp_obj$feeding_site_year$feeding_high[r],
    0L
  )

  expected_any <- replace_na(
    exp_obj$feeding_site_year$feeding_any[r],
    0L
  )

  stopifnot(
    exp_obj$feeding_high_matrix[jj, tt] == expected_high,
    exp_obj$feeding_any_matrix[jj, tt] == expected_any
  )
}

for(r in seq_len(nrow(exp_obj$matthew_site_year))) {

  jj <- match(
    exp_obj$matthew_site_year$Site[r],
    ordered_islands
  )

  tt <- exp_obj$matthew_site_year$year_idx[r]

  stopifnot(
    exp_obj$matthew_post_matrix[jj, tt] ==
      exp_obj$matthew_site_year$matthew_post2016[r],
    exp_obj$matthew_acute_matrix[jj, tt] ==
      exp_obj$matthew_site_year$matthew_acute_2017_2018[r]
  )
}

cat("PASS: feeding and Matthew matrices match their site-year definitions.\n")

## Explicit Matthew gate.
ffrc <- match("FFRC", ordered_islands)

year_to_idx <- function(year) {
  year - Origin + 1L
}

stopifnot(
  exp_obj$matthew_acute_matrix[ffrc, year_to_idx(2016)] == 0L,
  exp_obj$matthew_acute_matrix[ffrc, year_to_idx(2017)] == 1L,
  exp_obj$matthew_acute_matrix[ffrc, year_to_idx(2018)] == 1L,
  exp_obj$matthew_acute_matrix[ffrc, year_to_idx(2019)] == 0L,
  exp_obj$matthew_post_matrix[ffrc, year_to_idx(2016)] == 0L,
  exp_obj$matthew_post_matrix[ffrc, year_to_idx(2017)] == 1L,
  exp_obj$matthew_post_matrix[ffrc, year_to_idx(2018)] == 1L,
  exp_obj$matthew_post_matrix[ffrc, year_to_idx(2019)] == 1L
)

cat("PASS: FFRC Matthew coding is primary=2017+, acute=2017-2018.\n")

## ------------------------------------------------------------
## 5. Gompertz integrated-hazard algebra
##
## Base hazard:
##   h(a) = a10 * exp[b(a - 10)] * multiplier
##
## Exact segment integral:
##   H(a1,a2) = multiplier * (a10/b) *
##              {exp[b(a2-10)] - exp[b(a1-10)]}
## ------------------------------------------------------------

segment_exact <- function(a10, b, a1, a2, multiplier) {
  multiplier *
    (a10 / b) *
    (
      exp(b * (a2 - 10)) -
        exp(b * (a1 - 10))
    )
}

segment_numeric <- function(a10, b, a1, a2, multiplier) {
  integrate(
    function(age) {
      a10 *
        exp(
          b *
            (age - 10)
        ) *
        multiplier
    },
    lower = a1,
    upper = a2,
    rel.tol = 1e-12
  )$value
}

tests <- tribble(
  ~a10, ~b, ~a1, ~a2, ~multiplier,
  0.03, 0.08,  5.2,  6.0, 1.8,
  0.02, 0.05, 10.0, 11.0, 0.7,
  0.04, 0.12, 20.0, 21.0, 2.5
) %>%
  rowwise() %>%
  mutate(
    exact = segment_exact(
      a10,
      b,
      a1,
      a2,
      multiplier
    ),
    numeric = segment_numeric(
      a10,
      b,
      a1,
      a2,
      multiplier
    ),
    abs_diff = abs(exact - numeric)
  ) %>%
  ungroup()

print(tests)

stopifnot(
  max(tests$abs_diff) < 1e-9
)

cat("PASS: analytic Gompertz segment integral matches numerical integration.\n")

## ------------------------------------------------------------
## 6. Left-truncation density normalization for a constant segment
## ------------------------------------------------------------

a10 <- 0.03
b <- 0.08
L <- 7
mult <- 1.5

H_from_L <- function(x) {
  segment_exact(
    a10,
    b,
    L,
    x,
    mult
  )
}

f_cond <- function(x) {
  ifelse(
    x < L,
    0,
    a10 *
      exp(
        b *
          (x - 10)
      ) *
      mult *
      exp(
        -H_from_L(x)
      )
  )
}

mass <- integrate(
  f_cond,
  lower = L,
  upper = 250,
  rel.tol = 1e-10
)$value

cat(
  "Conditional Gompertz density mass from entry onward: ",
  signif(mass, 10),
  "\n",
  sep = ""
)

stopifnot(
  abs(mass - 1) < 1e-7
)

cat("PASS: left-truncated density is normalized.\n")

## ------------------------------------------------------------
## 7. Proportional-hazard exposure check
## ------------------------------------------------------------

beta <- 0.63
ages <- c(5, 10, 20, 30)

h0 <- 0.03 * exp(0.07 * (ages - 10))
h1 <- h0 * exp(beta)

stopifnot(
  max(
    abs(
      h1 / h0 -
        exp(beta)
    )
  ) <
    1e-12
)

cat("PASS: provisioning/Matthew coefficients are age-invariant hazard ratios.\n")

cat("\n============================================================\n")
cat("ALL DETERMINISTIC MODEL-MATH AUDITS PASSED\n")
cat("============================================================\n")
