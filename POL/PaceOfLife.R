## ============================================================
## MODEL_PaceOfLife.R
##
## Gompertz pace-of-life model
##
## log(a):
##   sex + cohort density + 5-year temporal cohort effect
##   + residual island effect
##
## log(b):
##   sex + feeding + species + residual island effect
##
## Feeding is an island/block-level treatment.
## Missing historical cohort densities use island mean density.
## ============================================================

source("DataPrep.R")

source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

library(nimble)
library(MCMCvis)
library(mcmcplots)
library(tidyverse)


## ============================================================
## 1. SPECIES + FEEDING
## ============================================================

## Species coding:
## 0 = figginsi
## 1 = inornata

species_island <- as.integer(
  ordered_islands %in% c("FFRC", "UCay", "U Cay")
)

names(species_island) <- ordered_islands

species <- species_island[island_idx]


## Feeding coding:
## 0 = none
## 1 = any feeding

feeding_by_island <- growth_info %>%
  mutate(feeding = tolower(trimws(feeding))) %>%
  group_by(island) %>%
  summarise(
    feed = as.integer(any(feeding != "none")),
    .groups = "drop"
  )

feed_island <- feeding_by_island$feed[
  match(
    ordered_islands,
    feeding_by_island$island
  )
]

feed <- feed_island[island_idx]


## Check island structure

island_table <- tibble(
  island_idx = seq_len(n_islands),
  island = ordered_islands,
  species = ifelse(
    species_island == 1,
    "inornata",
    "figginsi"
  ),
  feeding = ifelse(
    feed_island == 1,
    "fed",
    "unfed"
  )
)

print(island_table)

stopifnot(
  length(species) == nind,
  length(feed) == nind,
  !anyNA(species),
  !anyNA(feed)
)


## ============================================================
## 2. COHORT DENSITY
##
## Density = MNA / usable habitat area
##
## Observed island-year density where available.
## Missing historical cohort years use island mean density.
## ============================================================

habitat <- read_csv(
  "Data/Iguana_islands_habitat.csv",
  show_col_types = FALSE
) %>%
  transmute(
    island = islands,
    usable_m2 = Useable
  )


## ------------------------------------------------------------
## 2a. Observed life-history spans
## ------------------------------------------------------------

first_last <- growth_info %>%
  transmute(
    ID,
    island,
    first_year = tF_cal - Origin + 1,
    last_year  = tL_cal - Origin + 1
  )


## ------------------------------------------------------------
## 2b. Minimum number alive by island-year
## ------------------------------------------------------------

mna_long <- first_last %>%
  rowwise() %>%
  mutate(
    year_idx = list(
      seq(
        first_year,
        last_year
      )
    )
  ) %>%
  unnest(year_idx) %>%
  ungroup() %>%
  count(
    island,
    year_idx,
    name = "MNA"
  ) %>%
  left_join(
    habitat,
    by = "island"
  ) %>%
  mutate(
    density_raw =
      MNA /
      (usable_m2 / 10000)
  )


## ------------------------------------------------------------
## 2c. Restrict observed density to surveyed island-years
## ------------------------------------------------------------

survey_long <- actual_surveys %>%
  transmute(
    island = Island,
    year_idx = yr_adj
  )


density_observed <- survey_long %>%
  left_join(
    mna_long %>%
      select(
        island,
        year_idx,
        MNA,
        density_raw
      ),
    by = c(
      "island",
      "year_idx"
    )
  ) %>%
  mutate(
    MNA = replace_na(
      MNA,
      0L
    ),
    density_raw = replace_na(
      density_raw,
      0
    )
  )


## ------------------------------------------------------------
## 2d. Mean observed density for each island
## ------------------------------------------------------------

island_mean_density <- density_observed %>%
  group_by(island) %>%
  summarise(
    mean_density_raw =
      mean(
        density_raw,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

cat("\nMean observed density by island:\n")
print(island_mean_density)


## ------------------------------------------------------------
## 2e. Complete island x cohort density series
## ------------------------------------------------------------

density_complete <- expand_grid(
  island = ordered_islands,
  year_idx = 1:tMax
) %>%
  left_join(
    density_observed %>%
      select(
        island,
        year_idx,
        density_raw
      ),
    by = c(
      "island",
      "year_idx"
    )
  ) %>%
  left_join(
    island_mean_density,
    by = "island"
  ) %>%
  mutate(
    
    density_source =
      if_else(
        is.na(density_raw),
        "Island mean",
        "Observed cohort"
      ),
    
    density_filled_raw =
      coalesce(
        density_raw,
        mean_density_raw
      )
  )


cat("\nDensity source by island:\n")

print(
  density_complete %>%
    count(
      island,
      density_source
    ) %>%
    arrange(
      island,
      density_source
    )
)


## ------------------------------------------------------------
## 2f. Global standardisation after filling
##
## betaDENS_a = effect of 1 SD increase in log-density
## ------------------------------------------------------------

dens_mu <- mean(
  log1p(
    density_complete$density_filled_raw
  ),
  na.rm = TRUE
)

dens_sd <- sd(
  log1p(
    density_complete$density_filled_raw
  ),
  na.rm = TRUE
)

if(
  !is.finite(dens_sd) ||
  dens_sd <= 0
) {
  stop(
    "Invalid density standard deviation."
  )
}


density_complete <- density_complete %>%
  mutate(
    density_z =
      (
        log1p(
          density_filled_raw
        ) -
          dens_mu
      ) /
      dens_sd
  )


## ------------------------------------------------------------
## 2g. Island x year density matrix
## ------------------------------------------------------------

density_matrix <- matrix(
  NA_real_,
  nrow = n_islands,
  ncol = tMax,
  dimnames = list(
    ordered_islands,
    seq_len(tMax)
  )
)


for(
  r in seq_len(
    nrow(
      density_complete
    )
  )
) {
  
  j <- match(
    density_complete$island[r],
    ordered_islands
  )
  
  tt <-
    density_complete$year_idx[r]
  
  density_matrix[
    j,
    tt
  ] <-
    density_complete$density_z[r]
}


stopifnot(
  !anyNA(
    density_matrix
  ),
  nrow(
    density_matrix
  ) == n_islands,
  ncol(
    density_matrix
  ) == tMax
)


cat("\nDensity matrix summary:\n")

print(
  summary(
    as.vector(
      density_matrix
    )
  )
)


## ============================================================
## 3. TEMPORAL COHORT BLOCKS
##
## Five-year birth-cohort periods.
##
## We retain annual cohort index for density lookup,
## but temporal unexplained variation is modelled in broader
## 5-year blocks.
## ============================================================

block_width <- 5

n_time_blocks <- ceiling(
  tMax /
    block_width
)

cat(
  "\nTemporal structure:",
  n_time_blocks,
  "blocks of",
  block_width,
  "years\n"
)


## ============================================================
## 4. MODEL
## ============================================================

code_POL <- nimbleCode({
  
  for(i in 1:nind) {
    
    ## --------------------------------------------------------
    ## LATENT BIRTH
    ## --------------------------------------------------------
    
    tB[i] ~ dunif(
      cintB[i,1],
      cintB[i,2]
    )
    
    
    ## --------------------------------------------------------
    ## LEFT TRUNCATION
    ## --------------------------------------------------------
    
    L[i] <- max(
      0,
      island_start_time[i] -
        tB[i]
    )
    
    
    ## --------------------------------------------------------
    ## LATENT AGE AT DEATH
    ## --------------------------------------------------------
    
    tstar[i] ~ dGompertzLB(
      amult[i],
      bmult[i],
      lowerBound = L[i]
    )
    
    tD[i] <-
      tB[i] +
      tstar[i]
    
    censoredD[i] ~ dinterval(
      tD[i],
      cintD[i,1:2]
    )
    
    
    ## --------------------------------------------------------
    ## LATENT BIRTH COHORT
    ##
    ## Annual index retained for density lookup
    ## --------------------------------------------------------
    
    cohort_idx[i] <- max(
      1,
      min(
        tMax,
        floor(
          tB[i]
        ) + 1
      )
    )
    
    
    ## --------------------------------------------------------
    ## FIVE-YEAR TEMPORAL BLOCK
    ## --------------------------------------------------------
    
    time_block[i] <- max(
      1,
      min(
        n_time_blocks,
        floor(
          (
            cohort_idx[i] -
              1
          ) /
            block_width
        ) + 1
      )
    )
    
    
    ## --------------------------------------------------------
    ## COHORT DENSITY
    ## --------------------------------------------------------
    
    dens_cohort[i] <-
      density_matrix[
        island_idx[i],
        cohort_idx[i]
      ]
    
    
    ## ========================================================
    ## GOMPERTZ a
    ##
    ## Background mortality / pace
    ##
    ## sex
    ## cohort density
    ## temporal cohort block
    ## residual island deviation
    ## ========================================================
    
    log(
      amult[i]
    ) <-
      log(a) +
      betaSEX[1] *
      sex[i] +
      betaDENS_a *
      dens_cohort[i] +
      u_time[
        time_block[i]
      ] +
      island_a[
        island_idx[i]
      ]
    
    
    ## ========================================================
    ## GOMPERTZ b
    ##
    ## Actuarial ageing / senescence
    ##
    ## sex
    ## feeding treatment
    ## species
    ## residual island deviation
    ## ========================================================
    
    log(
      bmult[i]
    ) <-
      log(b) +
      betaSEX[2] *
      sex[i] +
      betaFEED_b *
      feed[i] +
      betaSPEC_b *
      species[i] +
      island_b[
        island_idx[i]
      ]
    
    
    ## ========================================================
    ## OBSERVATION / DETECTION COMPONENT
    ## ========================================================
    
    tB_year_idx[i] <- max(
      1,
      min(
        tMax,
        floor(
          tB[i]
        ) + 1
      )
    )
    
    
    nm_start[i] <- max(
      tB_year_idx[i],
      island_start_idx[i]
    )
    
    
    end_year[i] <- max(
      1,
      min(
        tMax,
        floor(
          tD[i]
        ) + 1
      )
    )
    
    
    nMpos[i] <-
      cum_effort[
        island_idx[i],
        end_year[i] + 1
      ] -
      cum_effort[
        island_idx[i],
        nm_start[i]
      ]
    
    
    nMpos_safe[i] <- max(
      y[i],
      nMpos[i]
    )
    
    
    pd_raw[i] <- exp(
      y[i] *
        log(
          mean.p +
            1.0E-10
        ) +
        (
          nMpos_safe[i] -
            y[i]
        ) *
        log(
          1 -
            mean.p +
            1.0E-10
        )
    )
    
    
    pd[i] <- min(
      0.999999,
      pd_raw[i]
    )
    
    
    dind[i] ~ dbern(
      pd[i]
    )
  }
  
  
  ## ==========================================================
  ## FIXED EFFECTS
  ## ==========================================================
  
  for(
    k in 1:2
  ) {
    
    betaSEX[k] ~ dnorm(
      0,
      sd = 1
    )
  }
  
  
  betaDENS_a ~ dnorm(
    0,
    sd = 1
  )
  
  
  betaFEED_b ~ dnorm(
    0,
    sd = 1
  )
  
  
  betaSPEC_b ~ dnorm(
    0,
    sd = 1
  )
  
  
  ## ==========================================================
  ## PARTIALLY POOLED RESIDUAL ISLAND EFFECTS
  ## ==========================================================
  
  sigma_island_a ~ dunif(
    0,
    2
  )
  
  sigma_island_b ~ dunif(
    0,
    2
  )
  
  
  for(
    r in 1:n_islands
  ) {
    
    island_a_raw[r] ~ dnorm(
      0,
      sd = 1
    )
    
    island_b_raw[r] ~ dnorm(
      0,
      sd = 1
    )
    
    
    island_a[r] <-
      sigma_island_a *
      island_a_raw[r]
    
    island_b[r] <-
      sigma_island_b *
      island_b_raw[r]
  }
  
  
  ## ==========================================================
  ## FIVE-YEAR TEMPORAL COHORT RANDOM EFFECT
  ##
  ## Independent zero-centred block deviations.
  ##
  ## Avoids cumulative RW1 drift and strong confounding with a.
  ## ==========================================================
  
  sigma_time ~ dunif(
    0,
    1
  )
  
  
  for(
    tt in 1:n_time_blocks
  ) {
    
    u_time[tt] ~ dnorm(
      0,
      sd = sigma_time
    )
  }
  
  
  ## ==========================================================
  ## BASELINE PARAMETERS
  ## ==========================================================
  
  a ~ dexp(1)
  
  b ~ dexp(1)
  
  mean.p ~ dunif(
    0,
    1
  )
})


## ============================================================
## 5. CONSTANTS + DATA
## ============================================================

consts_POL <- list(
  
  nind = nind,
  
  tMax = tMax,
  
  n_time_blocks =
    n_time_blocks,
  
  block_width =
    block_width,
  
  n_islands =
    n_islands,
  
  cintB =
    as.matrix(
      cintB_adj
    ),
  
  cintD =
    as.matrix(
      cintD_adj
    ),
  
  sex =
    as.numeric(
      sex
    ),
  
  feed =
    as.numeric(
      feed
    ),
  
  species =
    as.numeric(
      species
    ),
  
  island_idx =
    as.numeric(
      island_idx
    ),
  
  island_start_time =
    as.numeric(
      island_start_time
    ),
  
  island_start_idx =
    as.numeric(
      island_start_idx
    ),
  
  y =
    as.numeric(
      y
    )
)


data_POL <- list(
  
  censoredD =
    as.numeric(
      censoredD_vec
    ),
  
  dind =
    rep(
      1,
      nind
    ),
  
  cum_effort =
    as.matrix(
      cum_effort_for_model
    ),
  
  density_matrix =
    as.matrix(
      density_matrix
    )
)


## ============================================================
## 6. INITIAL VALUES
## ============================================================

make_inits_POL <- function(
    chain_id = 1
) {
  
  set.seed(
    1000 +
      chain_id
  )
  
  
  ## ----------------------------------------------------------
  ## Latent birth
  ## ----------------------------------------------------------
  
  tB_start <- runif(
    nind,
    min =
      cintB_adj[,1],
    max =
      cintB_adj[,2]
  )
  
  
  ## ----------------------------------------------------------
  ## Latent death
  ## ----------------------------------------------------------
  
  tD_start <- numeric(
    nind
  )
  
  
  for(
    i in seq_len(
      nind
    )
  ) {
    
    if(
      censoredD_vec[i] ==
      1
    ) {
      
      tD_start[i] <-
        runif(
          1,
          min =
            cintD_adj[i,1] +
            0.001,
          max =
            cintD_adj[i,2] -
            0.001
        )
      
    } else {
      
      tD_start[i] <-
        cintD_adj[i,2] +
        runif(
          1,
          1,
          5
        )
    }
    
    
    if(
      tD_start[i] <=
      tB_start[i]
    ) {
      
      tD_start[i] <-
        tB_start[i] +
        runif(
          1,
          1,
          5
        )
    }
  }
  
  
  tstar_start <- pmax(
    0.001,
    tD_start -
      tB_start
  )
  
  
  ## ----------------------------------------------------------
  ## Temporal random effects
  ## ----------------------------------------------------------
  
  u_time_start <- rnorm(
    n_time_blocks,
    0,
    0.05
  )
  
  
  ## ----------------------------------------------------------
  ## Return initials
  ## ----------------------------------------------------------
  
  list(
    
    tB =
      tB_start,
    
    tstar =
      tstar_start,
    
    a =
      exp(
        rnorm(
          1,
          log(0.002),
          0.25
        )
      ),
    
    b =
      exp(
        rnorm(
          1,
          log(0.15),
          0.20
        )
      ),
    
    mean.p =
      runif(
        1,
        0.20,
        0.35
      ),
    
    betaSEX =
      rnorm(
        2,
        0,
        0.10
      ),
    
    betaDENS_a =
      rnorm(
        1,
        0,
        0.10
      ),
    
    betaFEED_b =
      rnorm(
        1,
        0,
        0.10
      ),
    
    betaSPEC_b =
      rnorm(
        1,
        0,
        0.10
      ),
    
    sigma_island_a =
      runif(
        1,
        0.05,
        0.4
      ),
    
    sigma_island_b =
      runif(
        1,
        0.03,
        0.3
      ),
    
    island_a_raw =
      rnorm(
        n_islands,
        0,
        0.2
      ),
    
    island_b_raw =
      rnorm(
        n_islands,
        0,
        0.2
      ),
    
    sigma_time =
      runif(
        1,
        0.03,
        0.20
      ),
    
    u_time =
      u_time_start
  )
}


inits_POL <- list(
  make_inits_POL(1),
  make_inits_POL(2)
)


## ============================================================
## 7. BUILD MODEL
## ============================================================

model_POL <- nimbleModel(
  
  code =
    code_POL,
  
  constants =
    consts_POL,
  
  data =
    data_POL,
  
  inits =
    inits_POL[[1]],
  
  calculate =
    FALSE
)


lp <- model_POL$calculate()


cat(
  "\nInitial log probability:",
  lp,
  "\n"
)


if(
  !is.finite(
    lp
  )
) {
  
  stop(
    "Initial log probability is not finite."
  )
}


## ============================================================
## 8. CONFIGURE MCMC
## ============================================================

conf_POL <- configureMCMC(
  
  model_POL,
  
  monitors = c(
    
    "a",
    "b",
    "mean.p",
    
    "betaSEX",
    
    "betaDENS_a",
    "betaFEED_b",
    "betaSPEC_b",
    
    "sigma_island_a",
    "sigma_island_b",
    
    "island_a",
    "island_b",
    
    "sigma_time",
    "u_time"
  ),
  
  enableWAIC =
    TRUE
)


## ------------------------------------------------------------
## Baseline Gompertz parameters
## ------------------------------------------------------------

conf_POL$removeSamplers(
  c(
    "a",
    "b"
  )
)


conf_POL$addSampler(
  
  target =
    c(
      "a",
      "b"
    ),
  
  type =
    "AF_slice"
)


## ------------------------------------------------------------
## Hyperparameters
## ------------------------------------------------------------

for(
  v in c(
    "sigma_island_a",
    "sigma_island_b",
    "sigma_time"
  )
) {
  
  conf_POL$removeSamplers(
    v
  )
  
  conf_POL$addSampler(
    target = v,
    type = "slice"
  )
}


## ------------------------------------------------------------
## Fixed effects
## ------------------------------------------------------------

for(
  v in c(
    "betaSEX[1]",
    "betaSEX[2]",
    "betaDENS_a",
    "betaFEED_b",
    "betaSPEC_b"
  )
) {
  
  conf_POL$removeSamplers(
    v
  )
  
  conf_POL$addSampler(
    target = v,
    type = "slice"
  )
}


## ------------------------------------------------------------
## Latent birth + lifespan jointly
## ------------------------------------------------------------

for(
  i in seq_len(
    nind
  )
) {
  
  conf_POL$removeSamplers(
    c(
      paste0(
        "tB[",
        i,
        "]"
      ),
      paste0(
        "tstar[",
        i,
        "]"
      )
    )
  )
  
  
  conf_POL$addSampler(
    
    target =
      c(
        paste0(
          "tB[",
          i,
          "]"
        ),
        paste0(
          "tstar[",
          i,
          "]"
        )
      ),
    
    type =
      "AF_slice"
  )
}


## ============================================================
## 9. BUILD + COMPILE
## ============================================================

mcmc_POL <- buildMCMC(
  conf_POL
)


cModel_POL <- compileNimble(
  model_POL,
  resetFunctions = TRUE
)


cMCMC_POL <- compileNimble(
  mcmc_POL,
  project = model_POL
)


## ============================================================
## 10. TEST RUN
##
## Same test length as before.
## Do NOT go straight to 500k.
## ============================================================

run_POL_test2 <- runMCMC(
  
  cMCMC_POL,
  
  niter =
    50000,
  
  nburnin =
    14000,
  
  nchains =
    2,
  
  thin =
    10,
  
  inits =
    inits_POL,
  
  setSeed =
    TRUE,
  
  samplesAsCodaMCMC =
    TRUE,
  
  summary =
    TRUE,
  
  WAIC =
    TRUE,
  
  progressBar =
    TRUE
)


## ============================================================
## 11. SAVE
## ============================================================

saveRDS(
  run_POL_test2,
  "Model_PaceOfLife_5yr_TEST_samples.rds"
)


## ============================================================
## 12. MAIN DIAGNOSTICS
## ============================================================

main_pars <- c(
  
  "a",
  "b",
  
  "betaSEX",
  
  "betaDENS_a",
  "betaFEED_b",
  "betaSPEC_b",
  
  "sigma_island_a",
  "sigma_island_b",
  
  "sigma_time",
  
  "mean.p"
)


cat(
  "\nMAIN PARAMETER SUMMARY\n"
)


print(
  MCMCsummary(
    run_POL_test2$samples,
    params =
      main_pars
  )
)

## ============================================================
## 13. TEMPORAL EFFECTS
## ============================================================

cat(
  "\nTEMPORAL 5-YEAR COHORT EFFECTS\n"
)


print(
  MCMCsummary(
    run_POL_test2$samples,
    params =
      "u_time"
  )
)


## ============================================================
## 14. RESIDUAL ISLAND EFFECTS
## ============================================================

cat(
  "\nRESIDUAL ISLAND EFFECTS\n"
)


print(
  MCMCsummary(
    run_POL_test2$samples,
    params =
      c(
        "island_a",
        "island_b"
      )
  )
)


## ============================================================
## 15. WAIC
## ============================================================

cat(
  "\nWAIC\n"
)


print(
  run_POL_test2$WAIC
)


## ============================================================
## 16. TRACE PLOTS
## ============================================================

mcmcplot(
  run_POL_test2$samples,
  parms =
    main_pars
)


## ============================================================
## 17. FINAL LONG RUN
##
## Run only if:
##   a / b converge
##   feeding/species converge
##   sigma_time behaves
##   u_time effects no longer show chain-dependent drift
## ============================================================

# run_POL <- runMCMC(
#
#   cMCMC_POL,
#
#   niter = 500000,
#   nburnin = 100000,
#   nchains = 3,
#   thin = 10,
#
#   inits = inits_POL,
#
#   setSeed = TRUE,
#
#   samplesAsCodaMCMC = TRUE,
#
#   summary = TRUE,
#   WAIC = TRUE,
#
#   progressBar = TRUE
# )
#
#
# saveRDS(
#   run_POL,
#   "Model_PaceOfLife_5yr_samples.rds"
# )
#
#
# print(
#   MCMCsummary(
#     run_POL$samples
#   )
# )
#
#
# print(
#   run_POL$WAIC
# )