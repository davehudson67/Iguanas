## ============================================================
## IGUANA SURVIVAL MODEL: TWO-SPECIES SEPARATE ANALYSIS
## - Removes Island Random Effect to allow Feeding test
## - Binary covariates coded as 0 and 1
## - Fixes observation model NaNs
## ============================================================

library(data.table)
library(tidyverse)
library(nimble)

setwd("~/Iguanas")

## ---- 1. Load & Correct Global Data ----
load("igs_AllIslands_CleanDH_160226_obsStart.RData")
hab_area <- fread("Data/Iguana_islands_habitat.csv") %>% select(islands, Useable)

source("ModelComparison_FUNCTIONS.R")     
source("Distributions/Dist_GompertzLB.R") 
source("Distributions/Dist_Gompertz.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

# Correct FFRC figginsi -> inornata
problem_ids <- AllData %>% filter(Island == "FFRC", Species == "figginsi") %>% pull(ID)
if (length(problem_ids) > 0) {
  AllData <- AllData %>%
    mutate(Species = ifelse(ID %in% problem_ids & Island == "FFRC", "inornata", Species))
}

## ---- 2. Define the Master Modeling Function ----
run_species_model <- function(target_species, AllData, CH, cintB, cintD, hab_area) {
  
  message("\n============================================================")
  message("STARTING ANALYSIS FOR: ", target_species)
  message("============================================================\n")
  
  # -- A. Subset Data --
  full_id_df <- AllData %>% distinct(ID, .keep_all = TRUE)
  target_rows <- which(full_id_df$Species == target_species)
  
  sub_CH    <- CH[target_rows, ]
  sub_cintB <- cintB[target_rows, ]
  sub_cintD <- cintD[target_rows, ]
  
  sub_AllData <- AllData %>% filter(Species == target_species)
  id_df <- sub_AllData %>% distinct(ID, .keep_all = TRUE)
  
  nind <- nrow(id_df)
  tMax <- ncol(sub_CH)
  n_cohort <- tMax
  
  # -- B. Covariates (0 and 1 coding) --
  ordered_islands <- levels(as.factor(as.character(id_df$Island)))
  n_islands <- length(ordered_islands)
  
  # 0 = Female, 1 = Male
  sex <- as.numeric(as.factor(id_df$sex)) - 1
  # 0 = Unfed, 1 = Fed
  feeding <- ifelse(as.character(id_df$feeding) == "none", 0, 1)
  
  island_idx <- as.numeric(factor(id_df$Island, levels = ordered_islands))
  
  island_start_by_island <- id_df %>%
    distinct(Island, .keep_all = TRUE) %>%
    arrange(factor(Island, levels = ordered_islands)) %>%
    pull(IslandStart)
  IslandStart_ind <- island_start_by_island[island_idx]
  
  censoredD <- rep(2, nind)
  censoredD[sub_cintD[,1] > 0] <- 1
  y <- rowSums(sub_CH)
  dind <- rep(1, nind)
  last_occasion <- apply(sub_CH, 1, function(r) { w <- which(r > 0); if (length(w) == 0) 0L else max(w) })
  
  # -- C. Density Matrix --
  area_vec <- hab_area %>%
    mutate(islands = as.character(islands)) %>%
    right_join(tibble(islands = ordered_islands), by = "islands") %>%
    pull(Useable)
  
  id_df <- id_df %>% mutate(row_idx = 1:nind)
  
  first_last <- as.data.frame(sub_CH) %>%
    mutate(row_idx = 1:n()) %>%
    pivot_longer(-row_idx, names_to = "Year_Idx", values_to = "Detected") %>%
    mutate(Year_Idx = as.integer(str_remove(Year_Idx, "V"))) %>%
    left_join(id_df %>% select(row_idx, Island), by = "row_idx") %>%
    group_by(row_idx, Island) %>%
    summarise(
      First = suppressWarnings(min(Year_Idx[Detected == 1], na.rm = TRUE)),
      Last  = suppressWarnings(max(Year_Idx[Detected == 1], na.rm = TRUE)),
      anyDet = any(Detected == 1),
      .groups = "drop"
    ) %>%
    filter(anyDet)
  
  mna_df <- first_last %>%
    mutate(Year_Idx = map2(First, Last, seq)) %>%
    unnest(Year_Idx) %>%
    group_by(Island, Year_Idx) %>%
    summarise(MNA = n(), .groups = "drop") %>%
    mutate(island_idx = as.numeric(factor(Island, levels = ordered_islands))) %>%
    left_join(tibble(island_idx = 1:n_islands, Useable = area_vec), by = "island_idx") %>%
    mutate(Density = MNA / Useable) %>%
    filter(Year_Idx >= 1, Year_Idx <= tMax)
  
  dens_matrix <- matrix(0, nrow = n_islands, ncol = tMax)
  mna_scaled <- mna_df %>%
    group_by(island_idx) %>%
    mutate(dens_z = as.numeric(scale(log(Density)))) %>%  
    ungroup()
  mna_scaled$dens_z[!is.finite(mna_scaled$dens_z)] <- 0
  
  for (r in 1:nrow(mna_scaled)) {
    dens_matrix[mna_scaled$island_idx[r], mna_scaled$Year_Idx[r]] <- mna_scaled$dens_z[r]
  }
  
  # -- D. Effort Matrix --
  actual_surveys <- sub_AllData %>%
    filter(Status == 1) %>%
    distinct(Island, Year) %>%
    mutate(
      island_idx = as.numeric(factor(Island, levels = ordered_islands)),
      year_idx   = as.integer(Year)
    ) %>%
    filter(year_idx >= 1, year_idx <= tMax)
  
  effort_matrix <- matrix(0, nrow = n_islands, ncol = tMax)
  effort_matrix[as.matrix(actual_surveys[, c("island_idx", "year_idx")])] <- 1
  cum_effort <- t(apply(effort_matrix, 1, cumsum))
  cum_effort <- cbind(rep(0, nrow(cum_effort)), cum_effort)  
  cum_effort <- matrix(as.numeric(cum_effort), nrow = n_islands, ncol = tMax + 1)
  
  # -- E. Nimble Code (Species & Island RE removed) --
  code <- nimbleCode({
    for (i in 1:nind) {
      tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
      cohort_idx[i] <- max(1, min(n_cohort, floor(tB[i]) + 1))
      L[i] <- max(0, IslandStart_ind[i] - tB[i])
      tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
      tD[i] <- tB[i] + tstar[i]
      censoredD[i] ~ dinterval(tD[i], cintD[i, ])
      dens_entry[i] <- dens_matrix[ island[i], cohort_idx[i] ]
      
      log(amult[i]) <- log(a) +
        betaSEX[1]  * sex[i]        * zSEX[1] +
        betaFEED[1] * feeding[i]    * zFEED[1] +
        betaDENS[1] * dens_entry[i] * zDENS[1] +
        u_cohort_center[ cohort_idx[i] ]
      
      log(bmult[i]) <- log(b) +
        betaSEX[2]  * sex[i]        * zSEX[2] +
        betaFEED[2] * feeding[i]    * zFEED[2] +
        betaDENS[2] * dens_entry[i] * zDENS[2]
      
      tB_year_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
      island_start_year_idx[i] <- max(1, min(tMax, floor(IslandStart_ind[i]) + 1))
      nm_start[i] <- max(tB_year_idx[i], island_start_year_idx[i])
      end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
      
      nMpos[i] <- cum_effort[ island[i], end_year[i] + 1 ] - cum_effort[ island[i], nm_start[i] ]
      nMpos_safe[i] <- max(y[i], nMpos[i])
      
      pd_raw[i] <- exp(y[i] * log(mean.p + 1e-10) + (nMpos_safe[i] - y[i]) * log(1 - mean.p + 1e-10))
      pd[i] <- min(0.999999, pd_raw[i])
      dind[i] ~ dbern(pd[i])
    }
    
    u_cohort[1] ~ dnorm(0, sd = 10)
    for (t in 2:n_cohort) {
      u_cohort[t] ~ dnorm(u_cohort[t-1], sd = sigma_year_a)
    }
    u_bar <- mean(u_cohort[1:n_cohort])
    for (t in 1:n_cohort) {
      u_cohort_center[t] <- u_cohort[t] - u_bar
    }
    sigma_year_a ~ dunif(0, 5)
    
    for (k in 1:2) {
      betaSEX[k] ~ dnorm(0, sd = 1.5)
      betaFEED[k] ~ dnorm(0, sd = 1.5)
      betaDENS[k] ~ dnorm(0, sd = 1.5)
      
      zSEX[k] ~ dbern(0.5)
      zFEED[k] ~ dbern(0.5)
      zDENS[k] ~ dbern(0.5)
    }
    
    a ~ dexp(1)
    b ~ dexp(1)
    mean.p ~ dunif(0, 1)
  })
  
  # -- F. Constants, Data, Inits --
  consts <- list(
    nind = nind, tMax = tMax, n_cohort = n_cohort,
    sex = sex, feeding = feeding, island = island_idx,
    IslandStart_ind = IslandStart_ind
  )
  
  data_list <- list(
    y = as.numeric(y), cintB = as.matrix(sub_cintB), cintD = as.matrix(sub_cintD),
    censoredD = as.numeric(censoredD), dind = as.numeric(dind),
    dens_matrix = dens_matrix, cum_effort = cum_effort
  )
  
  initFn_fast <- function() {
    tiny <- 1e-6
    tBinit <- runif(nind, sub_cintB[,1], sub_cintB[,2])
    min_death <- pmax(tBinit, IslandStart_ind) + tiny
    min_death <- pmax(min_death, last_occasion + tiny)
    
    tDinit <- numeric(nind)
    id_right <- which(censoredD == 2)
    id_int   <- which(censoredD != 2)
    
    if (length(id_int)) {
      lo <- pmax(sub_cintD[id_int,1], min_death[id_int])
      hi <- sub_cintD[id_int,2]
      tDinit[id_int] <- runif(length(id_int), lo, hi)
    }
    if (length(id_right)) {
      tDprop <- sub_cintD[id_right,2] + rexp(length(id_right), rate = 5)
      tDinit[id_right] <- pmax(tDprop, min_death[id_right])
    }
    
    list(
      tB = tBinit, tD = tDinit, tstar = tDinit - tBinit,
      cohort_idx = pmax(1, pmin(n_cohort, floor(tBinit) + 1)), 
      a = rexp(1, rate = 10), b = rexp(1, rate = 10), mean.p = runif(1, 0.2, 0.8),
      betaSEX = rnorm(2, 0, 0.1), betaFEED = rnorm(2, 0, 0.1), betaDENS = rnorm(2, 0, 0.1),
      zSEX = c(0,0), zFEED = c(0,0), zDENS = c(0,0),
      u_cohort = rnorm(n_cohort, 0, 0.1), sigma_year_a = runif(1, 0.1, 1)
    )
  }
  
  inits_list <- list(initFn_fast(), initFn_fast())
  
  # -- G. Build & Run --
  model <- nimbleModel(code, constants = consts, data = data_list, inits = inits_list[[1]])
  cModel <- compileNimble(model)
  
  conf <- configureMCMC(model,
                        monitors = c("a","b","mean.p",
                                     "betaSEX","betaFEED","betaDENS",
                                     "zSEX","zFEED","zDENS",
                                     "sigma_year_a","u_cohort_center"),
                        enableWAIC = TRUE)
  
  conf$removeSamplers(c("a","b"))
  conf$addSampler(target = c("a","b"), type = "AF_slice")
  
  configureRJ(conf = conf,
              targetNodes = c("betaSEX","betaFEED","betaDENS"),
              indicatorNodes = c("zSEX","zFEED","zDENS"),
              control = list(mean = 0, scale = 0.5))
  
  mcmc  <- buildMCMC(conf)
  cMCMC <- compileNimble(mcmc, project = model)
  
  run <- runMCMC(cMCMC,
                 niter = 250000,
                 nburnin = 75000,
                 nchains = 2,
                 inits = inits_list,
                 setSeed = TRUE,
                 samplesAsCodaMCMC = TRUE,
                 summary = TRUE,
                 WAIC = TRUE,
                 progressBar = TRUE)
  
  # Save the output dynamically based on species name
  save_name <- paste0("iguana_survival_", target_species, ".rds")
  saveRDS(run, save_name)
  message("\nFinished ", target_species, ". Saved to ", save_name)
  
  return(run)
}

## ---- 3. Execute the Models ----

# Run for figginsi
run_figginsi <- run_species_model(target_species = "figginsi", 
                                  AllData = AllData, CH = CH, 
                                  cintB = cintB, cintD = cintD, hab_area = hab_area)
saveRDS(run_figginsi, "Run_figginsi.rds")

# Run for inornata
run_inornata <- run_species_model(target_species = "inornata", 
                                  AllData = AllData, CH = CH, 
                                  cintB = cintB, cintD = cintD, hab_area = hab_area)
saveRDS(run_inornata, "Run_inornata.rds")
run_inornata$summary

#------------------------------------------------------------------------------

