## ============================================================
## FINAL IGUANA SURVIVAL MODEL
## - Left-truncated Gompertz (Island-specific start)
## - RJMC Variable Selection: Sex, Feeding, Species, Density
## - INTERACTION: Species x Feeding (RJMC selected)
## - Time-Varying Density (Cohort-specific at entry)
## ============================================================
library(data.table)
library(tidyverse)
library(nimble)
library(mcmcplots)
library(bayesplot)
library(MCMCvis)

# Set working directory and load data
setwd("~/Iguanas")
load("igs_AllIslands_CleanDH_160226_obsStart.RData")
hab_area <- fread("Data/Iguana_islands_habitat.csv") %>% select(islands, Useable)

## Load custom distributions and functions
source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

## ------------------------------------------------------------
## 1) Build individual-level covariates
## ------------------------------------------------------------
id_df <- AllData %>% distinct(ID, .keep_all = TRUE)

# Numeric conversion for NIMBLE (0/1 coding)
sex <- as.numeric(as.factor(id_df$sex)) - 1
feeding <- ifelse(id_df$feeding == "none", 0, 1)
species <- ifelse(id_df$Species == "inornata", 1, 0)
island_idx <- as.numeric(as.factor(id_df$Island))
IslandStart <- id_df$IslandStart

## ------------------------------------------------------------
## 2) Time-varying density matrix (Fixed Size)
## ------------------------------------------------------------
n_islands <- length(levels(as.factor(id_df$Island)))
tMax <- ncol(CH)

# Calculate MNA (Minimum Number Alive)
mna_df <- as.data.frame(CH) %>%
  mutate(ID = 1:n()) %>%
  pivot_longer(-ID, names_to = "Year_Idx", values_to = "Detected") %>%
  mutate(Year_Idx = as.numeric(str_remove(Year_Idx, "V"))) %>%
  left_join(id_df %>% select(ID, Island), by = "ID") %>%
  group_by(ID, Island) %>%
  filter(any(Detected == 1)) %>%
  summarise(First = min(Year_Idx[Detected == 1]),
            Last = max(Year_Idx[Detected == 1]), .groups = "drop") %>%
  rowwise() %>%
  do(data.frame(ID = .$ID, Island = .$Island, Year_Idx = .$First:.$Last)) %>%
  group_by(Island, Year_Idx) %>%
  summarise(MNA = n(), .groups = "drop") %>%
  left_join(hab_area, by = c("Island" = "islands")) %>%
  mutate(Density = MNA / Useable)

# Pre-allocate matrix [Islands x Years] to prevent out-of-bounds errors
dens_matrix_fixed <- matrix(0, nrow = n_islands, ncol = tMax)

# Fill with Z-scored density data
mna_clean <- mna_df %>%
  group_by(Island) %>%
  mutate(dens_z = as.numeric(scale(log(Density)))) %>%
  ungroup() %>%
  filter(Year_Idx >= 1, Year_Idx <= tMax) %>%
  mutate(island_idx = as.numeric(as.factor(Island)))

for(i in 1:nrow(mna_clean)) {
  dens_matrix_fixed[mna_clean$island_idx[i], mna_clean$Year_Idx[i]] <- mna_clean$dens_z[i]
}

## Check sampling effort
# Identify years where surveys actually occurred (at least one capture)
survey_effort <- as.data.frame(CH) %>%
  mutate(ID = 1:n()) %>%
  pivot_longer(-ID, names_to = "Year_Idx", values_to = "Detected") %>%
  mutate(Year_Idx = as.numeric(str_remove(Year_Idx, "V"))) %>%
  left_join(id_df %>% select(ID, Island), by = "ID") %>%
  group_by(Island, Year_Idx) %>%
  summarise(Any_Captures = max(Detected), .groups = "drop") %>%
  filter(Any_Captures > 0)

# Calculate MNA
mna_plot_df <- as.data.frame(CH) %>%
  mutate(ID = 1:n()) %>%
  pivot_longer(-ID, names_to = "Year_Idx", values_to = "Detected") %>%
  mutate(Year_Idx = as.numeric(str_remove(Year_Idx, "V"))) %>%
  left_join(id_df %>% select(ID, Island), by = "ID") %>%
  group_by(ID, Island) %>%
  filter(any(Detected == 1)) %>%
  summarise(First = min(Year_Idx[Detected == 1]),
            Last = max(Year_Idx[Detected == 1]), .groups = "drop") %>%
  rowwise() %>%
  do(data.frame(ID = .$ID, Island = .$Island, Year_Idx = .$First:.$Last)) %>%
  group_by(Island, Year_Idx) %>%
  summarise(MNA = n(), .groups = "drop") %>%
  left_join(hab_area, by = c("Island" = "islands")) %>%
  mutate(Density = MNA / Useable) %>%
  # Flag years that were actually surveyed
  left_join(survey_effort, by = c("Island", "Year_Idx")) %>%
  mutate(Surveyed = ifelse(is.na(Any_Captures), "Interpolated", "Surveyed"))

# Plot Density over Time
ggplot(mna_plot_df, aes(x = Year_Idx, y = Density, group = Island)) +
  # Line shows the MNA trajectory
  geom_line(color = "grey70", linetype = "dashed") +
  # Points highlight actual survey years
  geom_point(aes(color = Surveyed, size = Surveyed)) +
  facet_wrap(~Island, scales = "free_y") +
  scale_color_manual(values = c("Interpolated" = "grey80", "Surveyed" = "firebrick")) +
  scale_size_manual(values = c("Interpolated" = 1, "Surveyed" = 2)) +
  theme_minimal() +
  labs(title = "Iguana Density (MNA / Useable Area) by Island",
       subtitle = "Red points indicate years with actual capture data; Dashed lines are interpolated MNA",
       x = "Year Index (Since 1955)",
       y = "Density (Iguanas per Useable Unit)") +
  theme(legend.position = "bottom")

## Create a binary effort matrix [n_islands x tMax]
n_islands <- length(levels(as.factor(id_df$Island)))
tMax <- ncol(CH)

# Identify actual survey events (Years where ANYONE was caught on that island)
actual_surveys <- AllData %>%
  filter(Status == 1) %>%
  distinct(Island, Year) %>%
  mutate(
    island_idx = as.numeric(as.factor(Island)),
    year_idx = Year
  ) %>%
  filter(year_idx >= 1, year_idx <= tMax)

# Pre-allocate a matrix of Zeros [Islands x Years]
effort_matrix <- matrix(0, nrow = n_islands, ncol = tMax)

# Fill with 1s using matrix-coordinate indexing (Fast and No NAs)
effort_matrix[as.matrix(actual_surveys[, c("island_idx", "year_idx")])] <- 1

# 5. Final Sanity Checks
stopifnot(ncol(effort_matrix) == ncol(CH))
stopifnot(nrow(effort_matrix) == n_islands)

# Calculate cumulative sum across years for each island
# effort_matrix is [n_islands x tMax]
cum_effort <- t(apply(effort_matrix, 1, cumsum))

# Add a column of Zeros at the beginning
cum_effort <- cbind(rep(0, nrow(cum_effort)), cum_effort)

# Ensure it is a plain numeric matrix
cum_effort <- matrix(as.numeric(cum_effort), 
                     nrow = nrow(cum_effort), 
                     ncol = ncol(cum_effort))

rm(actual_surveys)

## ------------------------------------------------------------
## 2) NIMBLE Model Code
## ------------------------------------------------------------
code <- nimbleCode({
  for (i in 1:nind) {
    # Latent birth time
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    
    # Left truncation
    L[i] <- max(0, IslandStart[i] - tB[i])
    
    # Survival
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    
    # Censoring
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    # Cohort Density Lookup
    entry_yr[i] <- max(1, min(tMax, trunc(tB[i])))
    idx_dens[i] <- dens_matrix[ island[i], entry_yr[i] ]
    
    # Linear Predictors (Strong Heredity)
    log(amult[i]) <- log(a) +
      betaSEX[1]     * sex[i]     * zSEX[1] +
      betaFEED[1]    * feeding[i] * zFEED[1] +
      betaSPEC[1]    * species[i] * zSPEC[1] +
      betaDENS[1]    * idx_dens[i] * zDENS[1] +
      betaSPEC_FEED[1] * (species[i] * feeding[i]) * (zSPEC[1] * zFEED[1] * zSPEC_FEED[1])
    
    log(bmult[i]) <- log(b) +
      betaSEX[2]     * sex[i]     * zSEX[2] +
      betaFEED[2]    * feeding[i] * zFEED[2] +
      betaSPEC[2]    * species[i] * zSPEC[2] +
      betaDENS[2]    * idx_dens[i] * zDENS[2] +
      betaSPEC_FEED[2] * (species[i] * feeding[i]) * (zSPEC[2] * zFEED[2] * zSPEC_FEED[2])
    
    # --- FIXED OBSERVATION WINDOW (Cumulative Logic) ---
    # nm: First year of risk
    # nM: Last year of risk
    nm[i] <- max(max(ceiling(tB[i]), IslandStart[island[i]]), 1)
    nM[i] <- min(min(floor(tD[i]), IslandEnd[island[i]]), tMax)
    
    # nMpos: Total survey years between nm and nM
    # Using the cumulative matrix: Effort(End) - Effort(Start-1)
    # Because cum_effort has a zero-column at [1], index j is at j+1
    nMpos[i] <- cum_effort[island[i], nM[i] + 1] - cum_effort[island[i], nm[i]]
    
    pd[i] <- exp(y[i] * log(mean.p + 1e-10) + (nMpos[i] - y[i]) * log(1 - mean.p + 1e-10))
    dind[i] ~ dbern(pd[i])
  }
  
  # Priors + RJMC
  for (k in 1:2) {
    betaSEX[k] ~ dnorm(0, sd = 1.5); betaFEED[k] ~ dnorm(0, sd = 1.5)
    betaSPEC[k] ~ dnorm(0, sd = 1.5); betaDENS[k] ~ dnorm(0, sd = 1.5)
    betaSPEC_FEED[k] ~ dnorm(0, sd = 1.5)
    zSEX[k] ~ dbern(0.5); zFEED[k] ~ dbern(0.5)
    zSPEC[k] ~ dbern(0.5); zDENS[k] ~ dbern(0.5); zSPEC_FEED[k] ~ dbern(0.5)
    z_eff_interaction[k] <- zSPEC[k] * zFEED[k] * zSPEC_FEED[k]
  }
  a ~ dexp(1); b ~ dexp(1); mean.p ~ dunif(0, 1)
})

## ------------------------------------------------------------
## 3) Setup Data, Constants, and Inits
## ------------------------------------------------------------
y <- rowSums(CH)

data_list <- list(
  y = as.numeric(unname(y)),
  cintB = as.matrix(unname(cintB)),
  cintD = as.matrix(unname(cintD)),
  censoredD = as.numeric(unname(censoredD)),
  dind = rep(1, nind),
  dens_matrix = dens_matrix_fixed,
  effort_matrix = effort_matrix # Strictly numeric, no NAs
)

consts <- list(
  nind = nind, 
  tMax = tMax, 
  sex = sex, 
  feeding = feeding, 
  species = species,
  island = island_idx, 
  IslandStart = IslandStart,
  IslandEnd = rep(69, tMax)
)

last_occasion <- apply(CH, 1, function(r) { w <- which(r > 0); if (length(w) == 0) 0L else max(w) })

initFn <- function(cintB, cintD, censoredD, IslandStart, last_occasion) {
  n <- nrow(cintB); tBinit <- runif(n, cintB[,1], cintB[,2])
  min_death <- pmax(tBinit + pmax(0, IslandStart - tBinit) + 1e-6, last_occasion + 1e-6)
  tDinit <- numeric(n); id_right <- which(censoredD == 2); id_int <- which(censoredD != 2)
  if (length(id_int)) tDinit[id_int] <- runif(length(id_int), pmax(cintD[id_int,1], min_death[id_int]), cintD[id_int,2])
  if (length(id_right)) tDinit[id_right] <- pmax(cintD[id_right,2] + rexp(length(id_right), 5), min_death[id_right])
  list(
    tB = tBinit, tD = tDinit, tstar = tDinit - tBinit, a = rexp(1, 10), b = rexp(1, 10), mean.p = runif(1, 0.2, 0.8),
    betaSEX = rnorm(2, 0, 0.1), betaFEED = rnorm(2, 0, 0.1), betaSPEC = rnorm(2, 0, 0.1), 
    betaDENS = rnorm(2, 0, 0.1), betaSPEC_FEED = rnorm(2, 0, 0.1),
    zSEX = c(0,0), zFEED = c(0,0), zSPEC = c(0,0), zDENS = c(0,0), zSPEC_FEED = c(0,0)
  )
}

inits_list <- list(initFn(cintB, cintD, censoredD, IslandStart, last_occasion),
                   initFn(cintB, cintD, censoredD, IslandStart, last_occasion))

## ------------------------------------------------------------
## 4) Compilation and Execution
## ------------------------------------------------------------
model <- nimbleModel(code, constants = consts, data = data_list, inits = inits_list)
cModel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c("a", "b", "mean.p", "betaSEX", "betaFEED", 
                                          "betaSPEC", "betaDENS", "betaSPEC_FEED",
                                          "zSEX", "zFEED", "zSPEC", "zDENS", "zSPEC_FEED",
                                          "z_eff_interaction"))

conf$removeSamplers(c("a", "b"))
conf$addSampler(target = c("a", "b"), type = "AF_slice")

configureRJ(conf = conf,
            targetNodes = c("betaSEX", "betaFEED", "betaSPEC", "betaDENS", "betaSPEC_FEED"),
            indicatorNodes = c("zSEX", "zFEED", "zSPEC", "zDENS", "zSPEC_FEED"),
            control = list(mean = 0, scale = 0.5))

mcmc <- buildMCMC(conf)
cMCMC <- compileNimble(mcmc, project = model)

run <- runMCMC(cMCMC, niter = 250000, nburnin = 75000, nchains = 2,
               samplesAsCodaMCMC = TRUE, summary = TRUE, progressBar = TRUE)

saveRDS(run, "iguana_finalInteraction_samples.rds")
run <- readRDS("iguana_finalInteraction_samples.rds")

p <- mcmc_trace(run$samples, regex_pars = "beta") +
  theme_minimal() +
  labs(title = "Iguana Survival Model: Beta Traceplots")
ggsave("Iguana_Betas_Trace.png", plot = p, width = 12, height = 10, dpi = 300)

# Create a clean summary table for the main effects
stats_table <- MCMCsummary(run$samples, 
                           params = c('a', 'b', 'betaSEX', 'betaFEED', 'betaSPEC', 'betaDENS',
                                      'betaSPEC_FEED', 'mean.p'), 
                           ISB = TRUE, 
                           round = 3)
stats_table

# Plot the Betas for Initial Mortality [1] and Senescence [2]
p_effects <- mcmc_intervals(run$samples, regex_pars = "beta") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Posterior Estimates of Covariate Effects",
       subtitle = "Values > 0 increase mortality risk; Values < 0 decrease it")

p_effects

# Extract means of the 'z' parameters
pip_results <- as.data.frame(run$summary$all.chains) %>%
  rownames_to_column("Parameter") %>%
  filter(str_detect(Parameter, "^z")) %>%
  select(Parameter, PIP = Mean) %>%
  arrange(desc(PIP))

pip_results

# Function to calculate Gompertz Survival
gompertz_surv <- function(t, a, b) { exp(-(a/b) * (exp(b*t) - 1)) }

# Extract means for the 4 groups (assuming Female, Average Density)
# Note: Adjust these based on your specific 'run$summary' row names
m <- run$summary$all.chains[, "Mean"]

plot_data <- expand.grid(
  Age = 0:50,
  Species = c("figginsi", "inornata"),
  Feeding = c("None", "Fed")
) %>%
  mutate(
    # Logic for a (Initial Mortality)
    a_val = case_when(
      Species == "figginsi" ~ m["a"],
      Species == "inornata" ~ m["a"] * exp(m["betaSPEC[1]"])
    ),
    # Logic for b (Senescence + Interaction)
    b_val = case_when(
      Species == "figginsi" & Feeding == "None" ~ m["b"],
      Species == "figginsi" & Feeding == "Fed"  ~ m["b"] * exp(m["betaFEED[2]"]),
      Species == "inornata" & Feeding == "None" ~ m["b"] * exp(m["betaSPEC[2]"]),
      Species == "inornata" & Feeding == "Fed"  ~ m["b"] * exp(m["betaSPEC[2]"] + m["betaFEED[2]"] + m["betaSPEC_FEED[2]"])
    ),
    Survival = gompertz_surv(Age, a_val, b_val)
  )

ggplot(plot_data, aes(x = Age, y = Survival, color = Feeding, linetype = Species)) +
  geom_line(size = 1) +
  theme_bw() +
  labs(y = "Probability of Survival", x = "Age (Years)") +
  scale_color_manual(values = c("None" = "black", "Fed" = "blue"))

# 1. Extract means from your results summary
# (Ensure these names match your run$summary$all.chains exactly)
m <- run$summary$all.chains[, "Mean"]

# 2. Define the Gompertz Hazard Function: h(t) = a * exp(b * t)
gompertz_hazard <- function(t, a, b) { a * exp(b * t) }

# 3. Create the plotting data frame
hazard_data <- expand.grid(
  Age = 0:50,
  Species = c("figginsi", "inornata"),
  Feeding = c("None", "Supplemental")
) %>%
  mutate(
    # Species index for logic
    spec_idx = ifelse(Species == "inornata", 1, 0),
    feed_idx = ifelse(Feeding == "Supplemental", 1, 0),
    
    # Calculate group-specific a (Initial Mortality)
    # Note: betaSPEC[1] was significant in your results
    a_eff = m["a"] * exp(m["betaSPEC[1]"] * spec_idx),
    
    # Calculate group-specific b (Senescence + Interaction)
    # Logic: b_base + spec_effect + feed_effect + interaction
    b_eff = m["b"] * exp(
      m["betaSPEC[2]"] * spec_idx + 
        m["betaFEED[2]"] * feed_idx + 
        m["betaSPEC_FEED[2]"] * (spec_idx * feed_idx)
    ),
    
    # Calculate Hazard
    Hazard = gompertz_hazard(Age, a_eff, b_eff)
  )

# 4. Plot 1: Linear Scale (Shows the 'explosion' of risk at old age)
ggplot(hazard_data, aes(x = Age, y = Hazard, color = Feeding, linetype = Species)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("None" = "black", "Supplemental" = "red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Iguana Mortality Curves (Hazard Function)",
       subtitle = "Instantaneous risk of death by age",
       x = "Age (Years)",
       y = "Hazard Rate (Risk of death per year)")

# 5. Plot 2: Log Scale (The 'Gompertz Law' Plot)
# This is the best way to see differences in the 'Aging Rate' (Slope)
ggplot(hazard_data, aes(x = Age, y = Hazard, color = Feeding, linetype = Species)) +
  geom_line(size = 1.2) +
  scale_y_log10() + # Log scale turns Gompertz into straight lines
  scale_color_manual(values = c("None" = "black", "Supplemental" = "red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Log-Mortality Curves",
       subtitle = "Parallel lines = same aging rate; Different slopes = different senescence",
       x = "Age (Years)",
       y = "Log Hazard Rate")
