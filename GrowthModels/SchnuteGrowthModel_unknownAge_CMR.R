# Load necessary library
library(nimble)
library(tidyverse)
library(lubridate)

# Load data
igs <- read_csv("Data/LeafCayThru2019_CleanOct2020.csv")
igs$age <- as.numeric(igs$age)

# Select known age individuals
igsKA <- igs %>%
  filter(!is.na(age), !is.na(svl)) %>%
  mutate(animal_id = as.numeric(as.factor(animal_id))) %>%
  mutate(date = mdy(date)) %>%
  mutate(year = year(date)) %>%
  group_by(animal_id) %>%
  distinct(year, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(Year = year - 1983)

min(igsKA$year)
max(igsKA$year)

N <- length(unique(igsKA$animal_id))
nT_ints <- max(igsKA$year) - min(igsKA$year) + 2
deltaT <- rep(1, nT_ints)

# Data matrices
svl_matrix <- matrix(NA, nrow = N, ncol = nT_ints)
y <- matrix(0, nrow = N, ncol = nT_ints)
age_matrix <- matrix(NA, nrow = N, ncol = nT_ints)

# Determine the first capture times
f <- apply(!is.na(svl_matrix), 1, which.max)
max(f)

# Fill the matrix with SVL measurements
for (i in 1:nrow(igsKA)) {
  id <- igsKA$animal_id[i]
  year <- igsKA$Year[i]
  svl <- igsKA$svl[i]
  age <- igsKA$age[i]
  svl_matrix[id, year] <- svl
  age_matrix[id, year] <- age
  y[id, year] <- 1
}

Lr <- svl_matrix
mean_svl <- mean(igs$svl[igs$age > 20], na.rm = TRUE)
sd_svl <- sd(igs$svl, na.rm = TRUE)

# Define the nimble model code
schnute_code <- nimbleCode({
  ## Growth model
  for(i in 1:nind){
    #SVL at hatching size
    Hs[i] <- mu_Hs + eps_Hs[i]
    eps_Hs[i] ~ dnorm(0, tau_Hs)
    #maximum size
    y2[i] <- mu_y2 + eps_y2[i]
    eps_y2[i] ~ dnorm(0, tau_y2)
    #characteristic grwth rate
    log(k[i]) <- mu_k + eps_k[i]
    eps_k[i] ~ dnorm(0, tau_k)
    
    # first capture
    Lr[i, f[i]] ~ dnorm(Lr_hat[i, f[i]], tau_Lr)         # SVL - observed
    Lr_hat[i, f[i]] ~ dunif(5, 40)   # expected size of animal at unknown age at first capture
    
  }
  
  #Priors Hsize
  mu_Hs ~ dunif(0, 20)
  tau_Hs <- 1 / sd_Hs^2
  sd_Hs ~ dunif(2, 3) 
  
  #Priors for y2
  mu_y2 ~ dunif(25, 80)
  tau_y2 <- 1 / sd_y2^2
  sd_y2 ~ dunif(0, 5)
  
  #Priors for K
  mu_k ~ dunif(0, -20)
  tau_k <- 1 / sd_k^2
  sd_k ~ dunif(0,2)
  
  tau_Lr <- 1 / sd_Lr^2
  sd_Lr ~ dunif(0, 5)

  for(i in 1:nind){
    for(t in (f[i] + 1):nT_ints){
      
      Lr[i,t] ~ dnorm(Lr_hat[i,t], tau_Lr)
      #Lr_hat[i,t] <- Lr_hat[i, t-1] * exp(-k[i] * deltaT[t-1]) + (y2[i]-Hs[i]*exp(-k[i]*(40-0))) * (1-exp(-k[i]*deltaT[t-1]))/(1-exp(-k[i]*(40-0)))
      Lr_hat[i,t] <- (Lr_hat[i, t-1] * exp(-k[i] * deltaT[t-1])) + ((y2[i] - (Hs[i] * exp(-k[i] * (40 - 0)))) * (1 - exp(-k[i] * deltaT[t-1])) / (1 - exp(-k[i] * (40 - 0))))
      
    }
    
    Linf[i] <- (y2[i] - Hs[i] * exp(-k[i] * (40 - 0))) / 1 - exp(-k[i] * (40 - 0))
    t0[i] <- 1 / k[i] * log((y2[i] - Hs[i]) / (y2[i] - Hs[i] * exp(-k[i] * (40 - 0))))
    
    # standardize SVL
    for(t in (f[i]):nT_ints){
      SVL_st[i,t] <- (Lr_hat[i,t] - mean.svl) / sd.svl
    }
  }
})

# Constants
constants <- list(
  nind = N,
  nT_ints = nT_ints,
  f = f,
  deltaT = deltaT,
  mean.svl = mean_svl,
  sd.svl = sd_svl
)

# Data for nimble
dataList <- list(
  Lr = svl_matrix
)

# Initial values
inits <- list(
  mu_Hs = 10,
  sd_Hs = 2.5,
  mu_y2 = 55,
  sd_y2 = 2.5,
  mu_k = -6.5,
  sd_k = 1,
  sd_Lr = 2.5
)

# Create and compile the nimble model
schnute_model <- nimbleModel(code = schnute_code, constants = constants, data = dataList, inits = inits)

C_schnute_model <- compileNimble(schnute_model)

## check config
config <- configureMCMC(C_schnute_model, monitors = c("Linf", "k", "Hs"), thin = 1)

#Build the model
built <- buildMCMC(config)
cBuilt <- compileNimble(built)
config$printMonitors()

#Run the model
system.time(run <- runMCMC(cBuilt, 
                           niter = 5000,
                           nburnin = 1000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run$samples, "growth_model_samples.rds")
samples <- readRDS("growth_model_samples.rds")

plot(run$samples)
run$summary

# Extract samples from the MCMC output
samples_matrix <- as.matrix(samples)

# Define a function for the Schnute growth model prediction
schnute_growth <- function(age, Hs, Linf, k) {
  (Hs * exp(-k * age)) + ((Linf - (Hs * exp(-k * 40))) * (1 - exp(-k * age)) / (1 - exp(-k * 40)))
}

# Generate age range for predictions
age_range <- seq(min(age_matrix, na.rm = TRUE), max(age_matrix, na.rm = TRUE), length.out = 100)

# Initialize an array to store the predicted SVL
predicted_svl <- matrix(NA, nrow = nrow(samples_matrix), ncol = length(age_range))

# Dynamically reference columns using grep
Hs_col <- grep("Hs\\[", colnames(samples_matrix))
Linf_col <- grep("Linf\\[", colnames(samples_matrix))
k_col <- grep("k\\[", colnames(samples_matrix))

for (i in 1:nrow(samples_matrix)) {
  Hs <- samples_matrix[i, Hs_col[1]]
  Linf <- samples_matrix[i, Linf_col[1]]
  k <- samples_matrix[i, k_col[1]]
  
  predicted_svl[i, ] <- schnute_growth(age_range, Hs, Linf, k)
}

# Calculate mean, lower and upper quantiles for the predictions
mean_pred <- apply(predicted_svl, 2, mean)
lower_pred <- apply(predicted_svl, 2, quantile, probs = 0.025)
upper_pred <- apply(predicted_svl, 2, quantile, probs = 0.975)

# Prepare data for plotting
plot_data <- data.frame(
  age = age_range,
  mean_svl = mean_pred,
  lower_svl = lower_pred,
  upper_svl = upper_pred
)

# Prepare the recorded data for plotting
recorded_data <- data.frame(
  age = as.vector(t(age_matrix)),
  svl = as.vector(t(svl_matrix))
)

# Remove rows with NA values from the recorded data
recorded_data <- recorded_data[!is.na(recorded_data$svl), ]

# Plot the predicted growth curve with CIs and recorded data
ggplot(plot_data, aes(x = age, y = mean_svl)) +
  geom_line(linewidth = 1, color = "blue") +
  geom_ribbon(aes(ymin = lower_svl, ymax = upper_svl), alpha = 0.2, fill = "blue") +
  geom_point(data = recorded_data, aes(x = age, y = svl), color = "red") +
  labs(title = "Predicted Growth Curve with Recorded Data", x = "Age", y = "SVL") +
  theme_minimal()
