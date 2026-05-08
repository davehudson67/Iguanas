library(nimble)
library(MCMCvis)

# Custom Distributions (Ensure these files are in your working directory)
source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

# --- 1. DECADAL MAPPING ---
# Map years (1 to tMax_adj) to 10-year chunks
decade_lookup <- floor((1:tMax_adj - 1) / 10) + 1
n_decades <- max(decade_lookup)

# --- 2. CONSTANTS & DATA LISTS ---
consts <- list(
  nind = nind, 
  tMax = tMax_adj, 
  n_islands = n_islands,
  n_decades = n_decades,
  sex = sex_vec, 
  island = island_idx,
  y_val = as.numeric(y_vec), 
  is_start_vec = as.numeric(is_start_vec_adj), 
  IslandStart_ind = is_start_vec_adj
)

data <- list(
  y = as.numeric(y_vec), 
  cintB = cintB_adj, 
  cintD = cintD_adj, 
  censoredD = censoredD_vec,
  cum_effort = cum_effort_adj,
  decade_lookup = decade_lookup # Must be in data for dynamic indexing
)

# --- 3. ROBUST INITIAL VALUES ---
init_Func <- function() {
  tBinit <- numeric(nind)
  tDinit <- numeric(nind)
  for(i in 1:nind) {
    tBinit[i] <- runif(1, cintB_adj[i, 1], cintB_adj[i, 2])
    d_low <- max(cintD_adj[i, 1], tBinit[i] + 0.1)
    d_hi  <- if(censoredD_vec[i] == 1) cintD_adj[i, 2] else d_low + 10
    tDinit[i] <- runif(1, d_low, max(d_low + 0.5, d_hi))
  }
  list(
    tB = tBinit, tstar = tDinit - tBinit,
    log_a = rnorm(n_islands, -8, 0.5),
    log_b = rnorm(n_islands, -2.3, 0.2),
    u_decade = c(NA, rnorm(n_decades - 1, 0, 0.1)),
    sigma_decade = 0.1,
    mean.p = 0.2,
    betaSEX = c(0, 0)
  )
}

# --- 4. NIMBLE MODEL CODE ---
code_IslandDecadal <- nimbleCode({
  for (i in 1:nind) {
    ## Latent Birth
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    cohort_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    L[i] <- max(0, IslandStart_ind[i] - tB[i])
    
    ## Latent Survival (Gompertz)
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    censoredD[i] ~ dinterval(tD[i], cintD[i, ])
    
    ## Linear Predictors
    # log(a): Island Baseline + Shared Decadal Noise + Sex
    log(amult[i]) <- log_a[island[i]] + 
      u_decade[decade_lookup[cohort_idx[i]]] + 
      betaSEX[1] * sex[i]
    
    # log(b): Island Senescence + Sex
    log(bmult[i]) <- log_b[island[i]] + 
      betaSEX[2] * sex[i]
    
    ## Observation Model
    nm_start[i] <- max(cohort_idx[i], is_start_vec[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    nMpos_raw[i] <- max(0, cum_effort[island[i], end_year[i] + 1] - 
                          cum_effort[island[i], nm_start[i]])
    nMpos_safe[i] <- max(y_val[i], nMpos_raw[i])
    y[i] ~ dbin(mean.p, nMpos_safe[i])
  }
  
  ## Shared Decadal Random Walk (The "Regional Noise Filter")
  u_decade[1] <- 0 
  for (d in 2:n_decades) {
    u_decade[d] ~ dnorm(u_decade[d-1], sd = sigma_decade)
  }
  sigma_decade ~ dexp(1)
  
  ## Island-Specific Parameters (8 Islands)
  for (j in 1:n_islands) {
    log_a[j] ~ dnorm(-8, sd = 4)
    log_b[j] ~ dnorm(-2, sd = 2)
  }
  
  ## Global Sex Effect
  for (k in 1:2) {
    betaSEX[k] ~ dnorm(0, sd = 1.5)
  }
  
  mean.p ~ dunif(0, 1)
})

# --- 5. BUILD AND COMPILE ---
model <- nimbleModel(code_IslandDecadal, constants = consts, data = data, inits = init_Func())
cModel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c("log_a", "log_b", "mean.p", "betaSEX", "sigma_decade", "u_decade"))

# --- 6. SAMPLER OPTIMIZATION ---

# Block Individuals (Diagonal Ridge Fix)
cat("Blocking individual birth/lifespan pairs...\n")
for (i in 1:nind) {
  conf$removeSamplers(paste0('tB[', i, ']'))
  conf$removeSamplers(paste0('tstar[', i, ']'))
  conf$addSampler(target = c(paste0('tB[', i, ']'), paste0('tstar[', i, ']')), type = "RW_block")
}

# Block Island Parameters (Gompertz Seesaw Fix)
for (j in 1:n_islands) {
  conf$removeSamplers(c(paste0("log_a[", j, "]"), paste0("log_b[", j, "]")))
  conf$addSampler(target = c(paste0("log_a[", j, "]"), paste0("log_b[", j, "]")), type = "AF_slice")
}

# --- 7. EXECUTION ---
cMCMC <- compileNimble(buildMCMC(conf), project = model)
results <- runMCMC(cMCMC, niter = 100000, nburnin = 40000, nchains = 2, summary = TRUE)

# --- 8. OUTPUT ---
# This gives you a separate a and b for every island
MCMCsummary(results$samples, params = c("log_a", "log_b", "betaSEX", "sigma_decade"))