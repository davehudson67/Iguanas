# Load required libraries
library(reshape2)
library(tidyverse)
library(coda)

out <- readRDS("ig_modelrun_out.rds")

mcmc_samples <- as.matrix(out$samples)

# Extract relevant posterior samples
mu_Hs_post <- mcmc_samples[ , 'mu_Hs']
mu_y2_post <- mcmc_samples[ , 'mu_y2']
mu_k_post <- mcmc_samples[ , 'mu_k']
maleEff_y2_post <- mcmc_samples[ , 'maleEff_y2']
maleEff_k_post <- mcmc_samples[ , 'maleEff_k']

# Simulate growth using Schnute function for a matrix of inputs
schnute_growth_vectorized <- function(t, Hs, y2, k) {
  Hs * exp(-k * t) + (y2 - Hs * exp(-k * 55)) * (1 - exp(-k * t)) / (1 - exp(-k * 55))
}

# Age sequence for prediction
age <- seq(0, 55, length.out = 100)

# Create matrices for male and female parameters from the MCMC samples
Hs_post <- mu_Hs_post  # No sex effect for Hs, so this is the same for males and females
y2_male <- mu_y2_post + maleEff_y2_post  # Male max size
k_male <- exp(mu_k_post + maleEff_k_post)  # Male growth rate

y2_female <- mu_y2_post  # Female max size
k_female <- exp(mu_k_post)  # Female growth rate

# Calculate SVL for all ages and all iterations simultaneously using outer
svl_male <- outer(Hs_post, age, FUN = function(Hs, t) schnute_growth_vectorized(t, Hs, y2_male, k_male))
svl_female <- outer(Hs_post, age, FUN = function(Hs, t) schnute_growth_vectorized(t, Hs, y2_female, k_female))

# Convert results to data frames for plotting
predicted_male_df <- melt(svl_male)
predicted_female_df <- melt(svl_female)

# Add age, sex, and iteration information to the data frames
colnames(predicted_male_df) <- c("iteration", "age_index", "svl")
colnames(predicted_female_df) <- c("iteration", "age_index", "svl")

# Add age (based on age index)
predicted_male_df$age <- age[predicted_male_df$age_index]
predicted_female_df$age <- age[predicted_female_df$age_index]

# Add sex information
predicted_male_df$sex <- "male"
predicted_female_df$sex <- "female"

# Combine male and female data into one data frame
predicted_df <- rbind(predicted_male_df, predicted_female_df)

# Clean up the data frame
predicted_df <- predicted_df[, c("iteration", "age", "svl", "sex")]

# Summarize predicted data for males and females
predicted_summary <- predicted_df %>%
  group_by(age, sex) %>%
  summarise(
    median_svl = median(svl),
    lower_svl = quantile(svl, 0.025),
    upper_svl = quantile(svl, 0.975)
  )

# known age observed data points
observed_svl <- igsKA %>%
  select(animal_id, age, sex, svl)

summary(as.factor(igsKA$sex))

# Plot using ggplot2 (code from before remains the same)
ggplot() +
  geom_ribbon(data = predicted_summary, aes(x = age, ymin = lower_svl, ymax = upper_svl, fill = sex), alpha = 0.3) +
  geom_line(data = predicted_summary, aes(x = age, y = median_svl, color = sex), size = 1) +
  geom_point(data = observed_svl, aes(x = age, y = svl, color = sex), size = 2, alpha = 0.8) +
  labs(x = "Age", y = "SVL", title = "Predicted SVL Growth Curves with Observed Data (Males and Females)") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) +  
  scale_fill_manual(values = c("lightblue", "pink")) +
  theme(legend.position = "right")
