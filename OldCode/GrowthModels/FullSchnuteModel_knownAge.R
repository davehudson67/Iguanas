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

# Define the NIMBLE model code with the full Schnute growth model
schnute_code_full <- nimbleCode({
  
  # Priors for the parameters
  a ~ dnorm(0, sd = 1)   # Prior for the growth rate parameter a
  b ~ dnorm(0, sd = 1)   # Prior for the shape parameter b
  Y1 ~ dnorm(4, sd = 10)
  Y2 ~ dnorm(40, sd = 10)

  # Random effects for each individual (if applicable)
  for (i in 1:N_individuals) {
    RE_Y[i] ~ dnorm(0, sd = sigma_Y)
    
    for (j in 1:N_observations[i]) {
      
      # Full Schnute model equation
      mu[i, j] <- (Y1 + RE_Y[i] + ((Y2 + RE_Y[i] - (Y1 + RE_Y[i])) / (1 - exp(-a * (max_age - min_age)))) * (1 - exp(-a * (age[i, j] - min_age))))^(1 / b)
      
      # Observed SVL follows a normal distribution around the expected value
      svl[i,j] ~ dnorm(mu[i,j], sd = sigma_obs)
    }
  }
  
  # Priors for the random effects and observation error
  sigma_Y ~ dunif(0, 10)
  sigma_obs ~ dunif(0, 10)
})

# Define initial values for MCMC sampling
inits_list_full <- list(
  a = 0.1,
  b = 0.5,
  Y1 = 8,
  Y2 = 43,
  sigma_Y = 1,
  sigma_obs = 1
)

constants_list <- list(N_individuals = N_individuals,
                       N_observations = N_observations,
                       min_age = min(igsKA$age),
                       max_age = max(igsKA$age))

# Parameters to monitor during MCMC
parameters_full <- c("a", "b", "Y1", "Y2", "sigma_Y", "sigma_obs")

# Build and compile the model
schnute_model_full <- nimbleModel(code = schnute_code_full, 
                                  data = data_list, 
                                  constants = constants_list, 
                                  inits = inits_list_full)

compiled_schnute_model_full <- compileNimble(schnute_model_full)

# Configure and compile the MCMC
schnute_mcmc_full <- buildMCMC(schnute_model_full, monitors = parameters_full)

compiled_schnute_mcmc_full <- compileNimble(schnute_mcmc_full, project = schnute_model_full)

# Run the MCMC
mcmc_samples_full <- runMCMC(compiled_schnute_mcmc_full, 
                             niter = 50000, 
                             nburnin = 12000, 
                             thin = 10,
                             nchains = 2,
                             summary = TRUE,
                             samplesAsCodaMCMC = TRUE)

# Summary of MCMC samples
mcmc_samples_full$summary
plot(mcmc_samples_full$samples)

# Convert MCMC samples to a matrix for further analysis
samples_matrix_full <- as.matrix(mcmc_samples_full$samples)


# Create a sequence of ages to plot the predicted growth curve
age_seq <- seq(min(constants_list$min_age), max(constants_list$max_age), length.out = 100)

# Initialize a matrix to store the predictions from each MCMC sample
n_samples <- nrow(samples_matrix_full)
predicted_svl_full <- matrix(NA, nrow = n_samples, ncol = length(age_seq))

# Loop over each sample to calculate the predicted SVL at each age
for (s in 1:n_samples) {
  a <- samples_matrix_full[s, "a"]
  b <- samples_matrix_full[s, "b"]
  Y1 <- samples_matrix_full[s, "Y1"]
  Y2 <- samples_matrix_full[s, "Y2"]
  
  exp_a_T1 <- exp(-a * (max_age - min_age))
  
  for (j in 1:length(age_seq)) {
    exp_a_t <- exp(-a * (age_seq[j] - min_age))
    mu <- (Y1 + (Y2 - Y1) / (1 - exp_a_T1) * (1 - exp_a_t))^(1 / b)
    predicted_svl_full[s, j] <- mu
  }
}

# Calculate the mean and 95% credible intervals across all samples
mean_pred_svl_full <- apply(predicted_svl_full, 2, mean)
lower_cred_int_full <- apply(predicted_svl_full, 2, quantile, probs = 0.025)
upper_cred_int_full <- apply(predicted_svl_full, 2, quantile, probs = 0.975)

# Plot the predicted growth curve with credible intervals
ggplot(prediction_summary, aes(x = Age)) +
  geom_ribbon(aes(ymin = lower_cred_int_full, ymax = upper_cred_int_full), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = mean_pred_svl_full), color = "blue", linewidth = 1) +
  geom_point(data = igsKA, aes(x = age, y = svl), color = "black", alpha = 0.5, size = 2) +
  labs(
    title = "Predicted Population-Level Growth Curve with Actual Data (Full Schnute Model)",
    x = "Age",
    y = "Snout-Vent Length (SVL)"
  ) +
  theme_minimal()
