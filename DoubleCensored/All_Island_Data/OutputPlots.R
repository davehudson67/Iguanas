library(tidyverse)
library(data.table)
library(coda)
library(patchwork)
library(bayestestR)

## load custom distributions
source("Distributions/Dist_Gompertz.R")
source("../NIMBLE_Distributions/Dist_GompertzNim.R")

run <- readRDS("DoubleCensored/All_Island_Data/Samples/G_AllIslands_SexFeeding2.rds")
run$summary

################################################################################
#
#  Inclusion Probability
#
################################################################################

## extract samples
samples <- as.matrix(run$samples)

## Marginal probabilities of inclusion for each variable
zNames <- model$expandNodeNames(c('zSEX', 'zFEED'))
zCols <- which(colnames(samples) %in% zNames)
binary <- as.data.table((samples[, zCols] != 0) + 0)
res <- binary[ , .N, by=names(binary)]
res <- res[order(N, decreasing = T)]
res <- res[, prob := N/dim(samples)[1]]
res
res

samples <- as.data.frame(samples)

z_indicators <- samples %>%
  select(all_of(zCols)) %>%
  colSums()

z_indicators <- data.frame(z_indicators/sum(res$N))

z_indicators$z <- rownames(z_indicators)
colnames(z_indicators) <- c("Inclusion_Prob", "z")
z_indicators$variable <- rep(c("Feeding", "Sex"), each = 2)
z_indicators$parameter <- rep(c("a", "b"), times = 2)

incl <- ggplot(z_indicators, aes(x = variable, y = Inclusion_Prob)) +
  geom_point(aes(colour = variable)) + 
  geom_hline(yintercept = 0.5, colour = "red") +
  scale_y_continuous(limits = c(0,1)) + 
  labs(x = "Variable", y = "Inclusion Probability") +
  facet_grid(~parameter) +
  theme(text = element_text(size = 15), legend.position = "none")
incl

################################################################################
#
#  Posterior Desnity PLots
#
################################################################################

post <- as.matrix(run$samples)
post <- sample_n(as.data.frame(post), 10000)
postdf <- as.data.frame(post)

colnames(postdf) 
colnames(postdf) <- c("a", "b", "betaFEEDa", "betaFEEDb", "betaSEXa", "betaSEXb",
                      "mean.p", "zFEEDa", "zFEEDb", "zSEXa", "zSEXb") 

# Find global x-axis limits across all distributions
global_x_range <- range(postdf$betaFEEDa, postdf$betaFEEDb, postdf$betaSEXa, postdf$betaSEXb)

# Prepare the data for both betaFEEDa and betaFEEDb in a long format
betaFEEDa <- postdf %>%
  filter(zFEEDa == 1) %>%
  select(betaFEEDa) %>%
  mutate(parameter = "a") %>%
  mutate(variable = "Feeding") %>%
  rename(Estimate = betaFEEDa)

betaFEEDb <- postdf %>%
  filter(zFEEDb == 1) %>%
  select(betaFEEDb) %>%
  mutate(parameter = "b") %>%
  mutate(variable = "Feeding") %>%
  rename(Estimate = betaFEEDb)

# Combine the two datasets into one
combined_data <- bind_rows(betaFEEDa, betaFEEDb)

# Compute credible intervals separately for each parameter
ci_hdi_a <- ci(betaFEEDa$Estimate, method = "HDI", ci = 0.95)
ci_hdi_b <- ci(betaFEEDb$Estimate, method = "HDI", ci = 0.95)

# Store the CI values in a dataframe with matching parameter names
ci_df <- data.frame(
  parameter = c("a", "b"),
  CI_low = c(ci_hdi_a$CI_low, ci_hdi_b$CI_low),
  CI_high = c(ci_hdi_a$CI_high, ci_hdi_b$CI_high)
)

# Plot both density distributions in the same plot, with color by parameter
#feed <- ggplot(combined_data, aes(x = betaFEED, fill = parameter, color = parameter)) +
#  geom_density(alpha = 0.3) +  # Transparent fill for better visibility
#  geom_vline(data = ci_df, aes(xintercept = CI_low, color = parameter), linetype = "dashed") +
#  geom_vline(data = ci_df, aes(xintercept = CI_high, color = parameter), linetype = "dashed") +
#  geom_vline(xintercept = 0, colour = "red") +
#  labs(x = "Beta Coefficient", y = "Density", title = "Density Plots of betaFEEDa and betaFEEDb") +
#  theme_minimal() +
#  scale_fill_manual(values = c("betaFEEDa" = "#2596be", "betaFEEDb" = "#e28743")) +  # Custom colors for fill
#  scale_color_manual(values = c("betaFEEDa" = "#2596be", "betaFEEDb" = "#e28743"))   # Custom colors for lines

# Feed plot
feed <- ggplot(combined_data, aes(x = Estimate, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.3) +  
  geom_vline(data = ci_df, aes(xintercept = CI_low, color = parameter), linetype = "dashed") +
  geom_vline(data = ci_df, aes(xintercept = CI_high, color = parameter), linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "red") +
  labs(x = "Beta Coefficient", y = "Density", title = "Feeding (none)") +
  theme_bw() +
  scale_fill_manual(values = c("a" = "#2596be", "b" = "#2596be")) +
  scale_color_manual(values = c("a" = "#2596be", "b" = "#2596be")) +
  coord_cartesian(xlim = global_x_range) +  # Align the x-axis across both plots
  facet_grid(~parameter) +
  theme(legend.position = "none")


# Prepare the data for both betaSEXa and betaSEXb in a long format
betaSEXa <- postdf %>%
  filter(zSEXa == 1) %>%
  select(betaSEXa) %>%
  mutate(parameter = "a") %>%
  mutate(variable = "Sex (male)") %>%
  rename(Estimate = betaSEXa)  # Add identifier

betaSEXb <- postdf %>%
  filter(zSEXb == 1) %>%
  select(betaSEXb) %>%
  mutate(parameter = "b")  %>%
  mutate(variable = "Sex (male)") %>%
  rename(Estimate = betaSEXb)

# Combine the two datasets into one
combined_data <- bind_rows(betaSEXa, betaSEXb)

# Compute credible intervals separately for each parameter
ci_hdi_a <- ci(betaSEXa$Estimate, method = "HDI", ci = 0.95)
ci_hdi_b <- ci(betaSEXb$Estimate, method = "HDI", ci = 0.95)

# Store the CI values in a dataframe with matching parameter names
ci_df <- data.frame(
  parameter = c("a", "b"),
  CI_low = c(ci_hdi_a$CI_low, ci_hdi_b$CI_low),
  CI_high = c(ci_hdi_a$CI_high, ci_hdi_b$CI_high)
)

# Sex plot
sex <- ggplot(combined_data, aes(x = Estimate, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.3) +  
  geom_vline(data = ci_df, aes(xintercept = CI_low, color = parameter), linetype = "dashed") +
  geom_vline(data = ci_df, aes(xintercept = CI_high, color = parameter), linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "red") +
  labs(x = "Beta Coefficient", y = "Density", title = "Sex (male)") +
  theme_bw() +
  scale_fill_manual(values = c("a" = "#2596be", "b" = "#2596be")) +
  scale_color_manual(values = c("a" = "#2596be", "b" = "#2596be")) +
  coord_cartesian(xlim = global_x_range) +  # Align the x-axis across both plots
  facet_grid(~parameter) +
  theme(legend.position = "none")

# Combine plots
feed / sex

################################################################################
#
# Trajectories
#
################################################################################

## load posterior samples
post <- as.matrix(run$samples)
post <- sample_n(as.data.frame(post), 10000)
post <- as.matrix(post)

####################################################################################################################################################################################
## create points to predict to for different groups
newdata <- expand.grid(t = 0:50, sex = 0:1, feeding = 0:1)

## loop over Gompertz parameters
parnms <- c("a", "b")
pars <- list()

for(i in 1:2) {
  ## extract samples for linear predictor for e.g. a
  ## making sure they line up with newdata

  parPost <- post[, match(paste0("beta", c("SEX", "FEED"), "[", i, "]"), colnames(post))]
  
  ## extract inclusion samples
  zpost <- post[, match(paste0("z", c("SEX", "FEED"), "[", i, "]"), colnames(post))]
  
  ## now multiply by relevant z indexes
  ## (expanding zpost matrix to make elementwise
  ## multiplication straightforward)
  zpost <- cbind(
    zpost[, paste0("zSEX[", i, "]")], ## Sex
    zpost[, paste0("zFEED[", i, "]")] ## Feeding
  )
  
  ## posterior samples for a linear predictor
  ## using matrix multiplication because it's faster
  ## (first brackets uses ELEMENTWISE multiplication (*)
  ## and then we use MATRIX multiplication e.g. %*%)
  parLinpred <- (parPost * zpost) %*% t(newdata[, -1])
  
  ## add intercept
  parLinpred <- parLinpred + log(post[, parnms[i]])
  
  ## transformed parameters
  pars[[i]] <- exp(parLinpred)
}

## for efficient sampling of the Gompertz collapse down to large matrix
pars <- purrr::map(pars, t) %>%
  purrr::map(as.vector) %>%
  reduce(cbind) %>%
  {cbind(rep(newdata$t, nrow(zpost)), .)}

## Siler sampling (calculate on the log-scale to try and avoid rounding errors)
postSurv <- pGompertz(pars[, 1], pars[, 2], pars[, 3], lower.tail = FALSE, log = TRUE)
postDens <- dGompertz(pars[, 1], pars[, 2], pars[, 3], log = TRUE)
postMort <- exp(postDens - postSurv)
postSurv <- exp(postSurv)

## posterior predictive summaries
postSurv <- postSurv %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025),#, na.rm = TRUE),
      Median = median(x),
      UCI = quantile(x, probs = 0.975))#, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  cbind(newdata[, 1:3]) %>%
  mutate(Feeding = ifelse(feeding == 1, "some", "none")) %>%
  mutate(Sex = ifelse(sex == 1, "Male", "Female")) %>%
  select(-sex, -feeding)

postMort <- postMort %>%
  matrix(nrow = nrow(newdata)) %>%
  apply(1, function(x) {
    c(LCI = quantile(x, probs = 0.025),#, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      UCI = quantile(x, probs = 0.975))#, na.rm = TRUE))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  cbind(newdata[, 1:3]) %>%
  mutate(Feeding = ifelse(feeding == 1, "some", "none")) %>%
  mutate(Sex = ifelse(sex == 1, "Male", "Female")) %>%
  select(-sex, -feeding)

## mortality plot
p1 <- mutate(postMort, Feeding = factor(Feeding)) %>%
  mutate(Time = t) %>%
  ggplot(aes(x = Time)) +
  geom_line(aes(y = Median, colour = Feeding)) +
  #geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = Feeding), alpha = 0.2) +
  facet_grid(~ Sex) +
  ylab("Mortality rate") +
  xlab("Time (yrs)") +
  labs(colour = "Feeding", fill = "Feeding", xlab = "Time(yrs)") +
  #coord_cartesian(ylim = c(0, 0.25), xlim = c(0,12)) +
  theme_bw() +
  theme(text = element_text(size = 15), 
        #axis.ticks = element_line(colour = "grey"),
        #panel.border = element_rect(fill = NA),
        #panel.background = element_blank(),
        #panel.grid.major = element_line(colour = "grey", size = 0.1), 
        #panel.grid.minor = element_line(colour = "grey", size = 0.05),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank()) 
#  scale_x_continuous(limits = c(0, 12), breaks = c(0, 4, 8, 12))
p1

## survival plot
p2 <- mutate(postSurv, Feeding = factor(Feeding)) %>%
  mutate(Time = t) %>%
  ggplot(aes(x = Time)) +
  geom_line(aes(y = Median, colour = Feeding)) +
  #geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = Inbreeding), alpha = 0.2) +
  facet_grid(~ Sex) +
  #labs(title = "Survival", subtitle = "tD ~ Sex * Infection + Inbreeding:Sex + Inbreeding:Infection") +
  scale_y_continuous(limits = c(0, 1)) +
  ylab("Survival probability") +
  xlab("Time (yrs)") +
  #ggtitle("b.") +
  labs(colour = "Feeding", fill = "Feeding") +
  theme_bw() +
  theme(text = element_text(size = 15),
        #axis.ticks = element_line(colour = "grey"),
        #panel.border = element_rect(fill = NA),
        #panel.background = element_blank(),
        #panel.grid.major = element_line(colour = "grey", size = 0.1), 
        #panel.grid.minor = element_line(colour = "grey", size = 0.05),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank())

p <- p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p

ggsave("predictivePosteriors_SexInfectionInbrCAT.pdf", p, width = 10, height = 5)


