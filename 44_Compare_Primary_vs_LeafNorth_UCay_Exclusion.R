## ============================================================
## 44_Compare_Primary_vs_LeafNorth_UCay_Exclusion.R
## ============================================================
library(tidyverse)
library(coda)

primary_file <- "Model_Final_IndivDetect_Audited_3CHAINS_COMPLETE.rds"
sens_file <- "Model_Sensitivity_Exclude_LeafNorth_UCay_3CHAINS_COMPLETE.rds"
if(!file.exists(primary_file)) stop("Run 38_Combine_Final_IndivDetect_Audited_Chains.R first.")
if(!file.exists(sens_file)) stop("Run 42_Combine_Sensitivity_Exclude_LeafNorth_UCay_Chains.R first.")

primary <- readRDS(primary_file)
sens <- readRDS(sens_file)
out_dir <- "DataAudits/Sensitivity_Exclude_LeafNorth_UCay"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summarise_param <- function(samples, parameter, analysis) {
  x <- as.matrix(samples)[, parameter]
  tibble(analysis = analysis, parameter = parameter, median = median(x), low95 = unname(quantile(x, 0.025)), high95 = unname(quantile(x, 0.975)), Pr_positive = mean(x > 0))
}

comparison <- bind_rows(
  summarise_param(primary$samples, "betaSEX[1]", "Primary 9-site"),
  summarise_param(primary$samples, "betaSEX[2]", "Primary 9-site"),
  summarise_param(sens$samples, "betaSEX[1]", "Exclude LeafNorth + UCay"),
  summarise_param(sens$samples, "betaSEX[2]", "Exclude LeafNorth + UCay")
) %>%
  mutate(effect = recode(parameter, `betaSEX[1]` = "Male effect on age-10 mortality", `betaSEX[2]` = "Male effect on senescence"))

write_csv(comparison, file.path(out_dir, "Primary_vs_exclusion_sex_effects.csv"))
print(comparison, n = Inf)

b_primary <- as.matrix(primary$samples)[, "betaSEX[2]"]
b_sens <- as.matrix(sens$samples)[, "betaSEX[2]"]
cat("\nKEY QUESTION\n")
cat("Primary P(male senescence effect > 0): ", round(mean(b_primary > 0), 4), "\n", sep = "")
cat("Excluding Leaf North + U Cay P(male senescence effect > 0): ", round(mean(b_sens > 0), 4), "\n", sep = "")
cat("Primary median betaSEX[2]: ", round(median(b_primary), 3), "\n", sep = "")
cat("Exclusion median betaSEX[2]: ", round(median(b_sens), 3), "\n", sep = "")
