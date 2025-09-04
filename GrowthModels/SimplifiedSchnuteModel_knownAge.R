# Install and load required packages
library(nimble)
library(tidyverse)
library(lubridate)
library(coda)

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
  mutate(Year = year - 1983) %>%
  mutate(birthYr = Year - round(min(age))) %>%
  ungroup() %>%
  select(animal_id, sex, age, svl) %>%
  arrange(animal_id)

# Prepare the data
N_individuals <- length(unique(igsKA$animal_id))
N_observations <- igsKA %>%
  group_by(animal_id) %>%
  summarise(count = n()) %>%
  pull(count)
max_obs <- max(N_observations)

# Initialize matrices with NA
data_list <- list(
  svl = matrix(NA, nrow = N_individuals, ncol = max_obs),
  age = matrix(NA, nrow = N_individuals, ncol = max_obs)
)

# Fill the matrices with data
for (i in 1:N_individuals) {
  id_data <- igsKA[igsKA$animal_id == unique(igsKA$animal_id)[i], ]
  n_obs <- nrow(id_data)  # Number of observations for this individual
  
  # Assign the observed data to the corresponding positions in the matrix
  data_list$svl[i, 1:n_obs] <- id_data$svl
  data_list$age[i, 1:n_obs] <- id_data$age
}

# Constants for the NIMBLE model
constants_list <- list(
  N_individuals = N_individuals,
  N_observations = N_observations,
  min_age = min(igsKA$age),
  max_age = max(igsKA$age)
)

# Define the NIMBLE model code
schnute_code <- nimbleCode({
  
  # Priors for the Schnute growth model parameters
  y1 ~ dnorm(4, sd = 10)
  y2 ~ dnorm(40, sd = 10)
  theta ~ dnorm(0.5, sd = 10)
  
  # Random effects for each individual
  for (i in 1:N_individuals) {
    RE_y[i] ~ dnorm(0, sd = sigma_y1)
  
    for (j in 1:N_observations[i]) {
      mu[i, j] <- (y1 + RE_y[i]) ^ (1 - theta) + ((y2 + RE_y[i]) ^ (1 - theta) - (y1 + RE_y[i]) ^ (1 - theta)) * (age[i, j] - min_age) / (max_age - min_age) 
      
      svl[i, j] ~ dnorm(mu[i, j] ^ (1 / (1 - theta)), sd = sigma_obs)
    }
  }
  
  # Priors for the random effects standard deviations
  sigma_y1 ~ dunif(0, 10)
  sigma_y2 ~ dunif(0, 10)
  sigma_obs ~ dunif(0, 10)
})

# Initial values for the MCMC
inits_list <- list(
  y1 = 8.8,
  y2 = 43.5,
  theta = 0.5,
  sigma_y1 = 1,
  sigma_y2 = 1,
  sigma_obs = 1,
  RE_y = rnorm(N_individuals, 0, 1)
)

# Parameters to monitor
parameters <- c("y1", "y2", "theta", "sigma_y1", "sigma_y2", "sigma_obs")

# Build the NIMBLE model
schnute_model <- nimbleModel(code = schnute_code, 
                             data = data_list, 
                             constants = constants_list, 
                             inits = inits_list)

#schnute_model$initializeInfo()

# Compile the model
compiled_schnute_model <- compileNimble(schnute_model)

# Configure the MCMC
schnute_mcmc <- buildMCMC(schnute_model, monitors = parameters)

# Compile the MCMC
compiled_schnute_mcmc <- compileNimble(schnute_mcmc, project = schnute_model)

# Run the MCMC
mcmc_samples <- runMCMC(compiled_schnute_mcmc, 
                        niter = 50000, 
                        nburnin = 12000, 
                        thin = 10,
                        nchains = 2,
                        samplesAsCodaMCMC = TRUE,
                        summary = TRUE)

saveRDS(mcmc_samples, "ThetaSchnuteModelRun.rds")

# Summary of MCMC samples
mcmc_samples$summary
plot(mcmc_samples$samples)

# Convert to data frame for easier manipulation
samples_df <- as.data.frame(as.matrix(mcmc_samples$samples))

# Define age range
age_sequence <- seq(from = min(igsKA$age), to = 75, length.out = 100)

# Sample a subset of posterior samples for computational efficiency
# Adjust the number as needed
n_samples <- nrow(samples_df)  # Number of posterior samples to use
sample_indices <- sample(1:nrow(samples_df), n_samples)

# Initialize a matrix to store predictions
predictions <- matrix(NA, nrow = nrow(samples_df), ncol = length(age_sequence))

# Loop over sampled posterior draws
for (i in 1:n_samples) {
  # Extract parameter values for this sample
  y1 <- samples_df$y1[sample_indices[i]]
  y2 <- samples_df$y2[sample_indices[i]]
  theta <- samples_df$theta[sample_indices[i]]
  
  # Compute predicted SVL for each age
  mu <- (y1)^(1 - theta) + ((y2)^(1 - theta) - (y1)^(1 - theta)) * 
    (age_sequence - min(igsKA$age)) / (max(igsKA$age) - min(igsKA$age))
  
  svl_pred <- mu^(1 / (1 - theta))
  
  # Store predictions
  predictions[i, ] <- svl_pred
}

# Convert predictions to data frame
predictions_df <- as.data.frame(predictions)
colnames(predictions_df) <- paste0("Age_", round(age_sequence, 2))

# Convert to long format for easier summarization
predictions_long <- predictions_df %>%
  pivot_longer(cols = everything(), names_to = "Age", values_to = "SVL") %>%
  mutate(Age = as.numeric(sub("Age_", "", Age)))

# Summarize predictions
prediction_summary <- predictions_long %>%
  group_by(Age) %>%
  summarize(
    mean_SVL = mean(SVL),
    median_SVL = median(SVL),
    lower_CI = quantile(SVL, 0.025),
    upper_CI = quantile(SVL, 0.975)
  )

# Plot using ggplot2
ggplot(prediction_summary, aes(x = Age)) +
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = mean_SVL), color = "blue", linewidth = 1) +
  labs(
    title = "Predicted Population-Level Growth Curve",
    x = "Age",
    y = "Snout-Vent Length (SVL)"
  ) +
  #theme_minimal() +
  geom_point(data = igsKA, aes(x = age, y = svl), color = "black", alpha = 0.5, size = 2)
