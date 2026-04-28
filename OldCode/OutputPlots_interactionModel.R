library(tidyverse)
library(data.table)
library(coda)
library(patchwork)
library(bayestestR)

## Load custom distributions
source("Distributions/Dist_Gompertz.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

# Load your new run
# run <- readRDS("iguana_survival_heredity_results.rds")

################################################################################
# 1) Inclusion Probabilities (Updated for all 5 variables)
################################################################################

samples <- as.matrix(run$samples)

# Use regex to find all z indicators
z_cols <- grep("^z", colnames(samples), value = TRUE)
# Exclude the 'z_eff' parameters if you just want the raw indicators
z_cols <- z_cols[!grepl("z_eff", z_cols)]

z_probs <- colMeans(samples[, z_cols] != 0) %>%
  enframe(name = "parameter", value = "Inclusion_Prob") %>%
  mutate(
    # Parse names like zSPEC_FEED[1]
    variable = str_extract(parameter, "(?<=z)[A-Z_]+"),
    type = ifelse(str_detect(parameter, "\\[1\\]"), "a (Initial)", "b (Senescence)")
  )

incl_plot <- ggplot(z_probs, aes(x = variable, y = Inclusion_Prob, color = variable)) +
  geom_point(size = 3) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  facet_wrap(~type) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  labs(title = "Variable Inclusion Probabilities", y = "Probability", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

incl_plot

################################################################################
# 2) Posterior Density Plots (Conditional on Inclusion)
################################################################################

# We want to plot the betas only for iterations where they were "active" (z=1)
beta_cols <- grep("^beta", colnames(samples), value = TRUE)

# Create a long format of betas filtered by their corresponding z-indicator
post_betas <- map_df(beta_cols, function(b_name) {
  z_name <- str_replace(b_name, "beta", "z")
  
  # Extract samples where this specific variable was included
  valid_samples <- samples[samples[, z_name] == 1, b_name]
  
  tibble(
    Estimate = valid_samples,
    Parameter = b_name,
    Variable = str_extract(b_name, "(?<=beta)[A-Z_]+"),
    Type = ifelse(str_detect(b_name, "\\[1\\]"), "a", "b")
  )
})

dens_plot <- ggplot(post_betas, aes(x = Estimate, fill = Variable)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, color = "black", linetype = "dotted") +
  facet_grid(Variable ~ Type, scales = "free") +
  theme_bw() +
  labs(title = "Conditional Posterior Densities (z=1)") +
  theme(legend.position = "none")

dens_plot

################################################################################
# 3) Trajectories (Updated for Interaction & Density)
################################################################################

# Sub-sample for speed
post_sub <- as.data.frame(samples) %>% sample_n(5000)

# Define prediction grid (Density set to 0 = average)
newdata <- expand.grid(
  t = 0:50, 
  sex = 0:1, 
  feeding = 0:1, 
  species = 0:1, 
  dens = 0 
)

pars <- list()

for(k in 1:2) { # k=1 for 'a', k=2 for 'b'
  
  # Extract relevant columns for this parameter index [k]
  b_sex  <- post_sub[[paste0("betaSEX[", k, "]")]]
  b_feed <- post_sub[[paste0("betaFEED[", k, "]")]]
  b_spec <- post_sub[[paste0("betaSPEC[", k, "]")]]
  b_dens <- post_sub[[paste0("betaDENS[", k, "]")]]
  b_int  <- post_sub[[paste0("betaSPEC_FEED[", k, "]")]]
  
  # Extract indicators
  z_sex  <- post_sub[[paste0("zSEX[", k, "]")]]
  z_feed <- post_sub[[paste0("zFEED[", k, "]")]]
  z_spec <- post_sub[[paste0("zSPEC[", k, "]")]]
  z_dens <- post_sub[[paste0("zDENS[", k, "]")]]
  z_int  <- post_sub[[paste0("zSPEC_FEED[", k, "]")]]
  
  # Intercept (a or b)
  intercept <- log(post_sub[[ifelse(k==1, "a", "b")]])
  
  # Matrix multiplication for the linear predictor
  # Note: Interaction logic follows Strong Heredity: (zSPEC * zFEED * zINT)
  lin_pred <- matrix(NA, nrow = nrow(post_sub), ncol = nrow(newdata))
  
  for(j in 1:nrow(newdata)) {
    row <- newdata[j,]
    lin_pred[,j] <- intercept +
      (b_sex  * z_sex  * row$sex) +
      (b_feed * z_feed * row$feeding) +
      (b_spec * z_spec * row$species) +
      (b_dens * z_dens * row$dens) +
      (b_int  * (z_spec * z_feed * z_int) * (row$species * row$feeding))
  }
  
  pars[[k]] <- exp(lin_pred)
}

# Calculate Survival and Mortality
# pars[[1]] is amult, pars[[2]] is bmult
postSurv <- matrix(NA, nrow = nrow(newdata), ncol = nrow(post_sub))
postMort <- matrix(NA, nrow = nrow(newdata), ncol = nrow(post_sub))

for(j in 1:nrow(newdata)) {
  t_val <- newdata$t[j]
  a_vec <- pars[[1]][,j]
  b_vec <- pars[[2]][,j]
  
  # Gompertz formulas
  postSurv[j,] <- exp(-(a_vec/b_vec) * (exp(b_vec * t_val) - 1))
  postMort[j,] <- a_vec * exp(b_vec * t_val)
}

# Summarize for plotting
summarize_post <- function(mat, newdata) {
  as_tibble(t(apply(mat, 1, function(x) quantile(x, c(0.025, 0.5, 0.975))))) %>%
    set_names(c("LCI", "Median", "UCI")) %>%
    bind_cols(newdata) %>%
    mutate(
      Feeding = ifelse(feeding == 1, "Fed", "None"),
      Sex = ifelse(sex == 1, "Male", "Female"),
      Species = ifelse(species == 1, "inornata", "figginsi")
    )
}

surv_df <- summarize_post(postSurv, newdata)
mort_df <- summarize_post(postMort, newdata)

# Plotting
p_surv <- ggplot(surv_df, aes(x = t, y = Median, color = Feeding, linetype = Species)) +
  geom_line(size = 1) +
  facet_wrap(~Sex) +
  theme_bw() +
  labs(y = "Survival Probability", x = "Age")

p_mort <- ggplot(mort_df, aes(x = t, y = Median, color = Feeding, linetype = Species)) +
  geom_line(size = 1) +
  scale_y_log10() +
  facet_wrap(~Sex) +
  theme_bw() +
  labs(y = "Log Mortality Rate", x = "Age")

(p_surv / p_mort) + plot_layout(guides = "collect")

# 1. Mortality Plot (Linear Scale)
p_mort_linear <- ggplot(mort_df, aes(x = t, y = Median, color = Feeding, fill = Feeding)) +
  # Add ribbons for 95% Credible Intervals
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.1, color = NA) +
  geom_line(aes(linetype = Species), size = 1) +
  facet_wrap(Species~Sex, scales = "free") +
  theme_bw() +
  # Use coord_cartesian to zoom in if the 'explosion' at age 50 is too high
  # coord_cartesian(ylim = c(0, 0.5)) + 
  labs(
    title = "Instantaneous Mortality Rate (Hazard)",
    subtitle = "Linear scale shows the exponential increase in risk with age",
    y = "Hazard (Risk of death per year)", 
    x = "Age (Years)"
  ) +
  scale_color_manual(values = c("None" = "black", "Fed" = "#2596be")) +
  scale_fill_manual(values = c("None" = "black", "Fed" = "#2596be")) +
  theme(strip.background = element_blank(), text = element_text(size = 14))

p_mort_linear

# 2. Combine with Survival Plot
p_final <- p_surv / p_mort_linear + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

#------------------------------------------------------------------------------

# 1. Flatten MCMC chains into a single data frame
# as.matrix(run$samples) automatically merges all chains
post_df <- as.data.frame(as.matrix(run$samples))

# 2. Clean the column names 
# NIMBLE names like "betaDENS[1]" can be tricky for some functions.
# We'll change them to "betaDENS_1" for easier handling.
colnames(post_df) <- gsub("\\[|\\]", "_", colnames(post_df))
colnames(post_df) <- gsub("__", "_", colnames(post_df)) # Clean up double underscores

# 3. Select the parameters you are suspicious of
# We'll look at the initial mortality effects ([1])
plot_cols <- c("betaDENS_1_", "betaFEED_1_", "betaSPEC_1_")

# 4. Thin the data significantly
# ggpairs is very slow with 100k+ points. 1,000 points is enough to see correlations.
post_thinned <- post_df %>%
  select(all_of(plot_cols)) %>%
  slice_sample(n = 1000)

# 5. Create the plot
p_pairs <- ggpairs(post_thinned,
                   lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
                   diag = list(continuous = wrap("densityDiag", fill = "blue", alpha = 0.3)),
                   upper = list(continuous = wrap("cor", size = 5))) +
  theme_minimal() +
  labs(title = "Correlation between Density, Feeding, and Species Effects",
       subtitle = "Diagonal shows multimodality; Scatter shows parameter trade-offs")

# 6. Save to file (Safe for servers/headless environments)
ggsave("MCMC_Pairs_Plot.png", p_pairs, width = 10, height = 10, dpi = 300)

# Display the plot if you are in a local RStudio session
# print(p_pairs)

#-------------------------------------------------------------------------------#

# Find the individuals with the longest capture histories
influential_check <- AllData %>%
  group_by(ID, Species, Island, feeding) %>%
  summarise(
    Years_Observed = max(Year) - min(Year),
    First_Year = min(Year),
    Last_Year = max(Year),
    Known_Death = any(Status == 2),
    n_captures = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Years_Observed))

# Look specifically at the 'Fed' inornata (likely UCay)
fed_inornata_senescence <- influential_check %>%
  filter(Species == "inornata", feeding != "none") %>%
  head(20)

print(fed_inornata_senescence)

# Look for 'Found Dead' individuals and their estimated 'study age'
deaths_check <- AllData %>%
  group_by(ID) %>%
  filter(any(Status == 2)) %>%
  summarise(
    Species = first(Species),
    Island = first(Island),
    Feeding = first(feeding),
    Years_Since_First_Cap = max(Year) - min(Year),
    Death_Year = max(Year)
  ) %>%
  arrange(Species, Feeding, desc(Years_Since_First_Cap))

# Compare death timing between species
ggplot(deaths_check, aes(x = Years_Since_First_Cap, fill = Feeding)) +
  geom_histogram(binwidth = 2) +
  facet_wrap(~Species) +
  theme_minimal() +
  labs(title = "Years from First Capture to Death",
       x = "Years in Study", y = "Number of Deaths")

AllData %>%
  group_by(Island) %>%
  summarise(
    Study_Duration = max(Year) - min(Year),
    Avg_Indiv_Span = mean(reframe(group_by(cur_data(), ID), s = max(Year)-min(Year))$s)
  )

#------------------------------------------------------------------------------#
library(patchwork)

# 1) Identify influential deaths and "Ghost Deaths"
inornata_death_diagnostics <- AllData %>%
  filter(Species == "inornata") %>%
  group_by(ID) %>%
  summarise(
    Island = first(Island),
    Feeding = first(feeding),
    Has_Death_Record = any(Status == 2),
    # Count how many records have a valid Year
    Valid_Year_Count = sum(!is.na(Year)),
    # Calculate span only for those with valid years
    Span = if(Valid_Year_Count > 1) max(Year, na.rm = TRUE) - min(Year, na.rm = TRUE) else 0,
    # Identify if the death record itself is missing a year
    Ghost_Death = any(Status == 2 & is.na(Year)),
    .groups = "drop"
  ) %>%
  filter(Has_Death_Record == TRUE)

# 2) Compare the "Death Profile" between Fed and None
death_profile_summary <- inornata_death_diagnostics %>%
  group_by(Feeding) %>%
  summarise(
    Total_Deaths = n(),
    Avg_Span_Before_Death = mean(Span, na.rm = TRUE),
    Ghost_Death_Count = sum(Ghost_Death),
    Short_Life_Deaths = sum(Span < 5 & !Ghost_Death), # Died within 5 years of first cap
    .groups = "drop"
  )

print(death_profile_summary)

# 2) Identify "Dead on Arrival" (DOA) individuals
# These are iguanas caught for the first time and found dead in the same year.
# These are HIGHLY influential for the 'a' (initial mortality) parameter.
doa_check <- AllData %>%
  group_by(ID) %>%
  summarise(
    Species = first(Species),
    Feeding = first(feeding),
    Is_DOA = all(Status == 2) || (n() == 1 & any(Status == 2)),
    .groups = "drop"
  ) %>%
  group_by(Species, Feeding) %>%
  summarise(
    Total_Indiv = n(),
    DOA_Count = sum(Is_DOA),
    Percent_DOA = (DOA_Count / Total_Indiv) * 100
  )

print(doa_check)


table(AllData$Status)

AllData %>%
  filter(if_any(everything(), ~ grepl("dead", .x, ignore.case = TRUE))) %>%
  select(ID, Island, Year, contains("dead"), everything()) %>%
  head()

# Check the 'Vanishing' profile for inornata
vanishing_inornata <- AllData %>%
  filter(Species == "inornata") %>%
  group_by(ID, feeding) %>%
  summarise(
    Last_Year = max(Year),
    Span = max(Year) - min(Year),
    .groups = "drop"
  ) %>%
  # Study ended around Year 64 (2019)
  mutate(Years_Missing = 64 - Last_Year)

# Visualize: Are Fed iguanas 'missing' for longer?
ggplot(vanishing_inornata, aes(x = Years_Missing, fill = feeding)) +
  geom_histogram(binwidth = 2, position = "dodge") +
  theme_minimal() +
  labs(title = "Inornata: Years since last capture",
       subtitle = "The model treats 'Years Missing' as 'Likely Dead'",
       x = "Years since last seen", y = "Count of Individuals")

AllData %>% filter(Status == 2) %>% select(ID, Species, Island, Year)
