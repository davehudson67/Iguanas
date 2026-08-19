## ============================================================
## FINAL CORRECTED PIPELINE: SVL-Tightened with 30-Year Cap
## ============================================================
library(tidyverse)
library(lubridate)
library(quantreg)
library(data.table)

rm(list = ls())

# --- 1. LOAD RAW DATA ---
island_files <- list(
  Bitter = "DoubleCensored/All_Island_Data/Data/BitterGuanaCay_clean.csv",
  FFRC   = "DoubleCensored/All_Island_Data/Data/FlatRockThru2019_CleanOct2020.csv",
  Gaulin = "DoubleCensored/All_Island_Data/Data/GaulinCay_clean.csv",
  Leaf   = "DoubleCensored/All_Island_Data/Data/LeafCayThru2019_CleanOct2020.csv",
  Noddy  = "DoubleCensored/All_Island_Data/Data/NoddyCCay_Clean.csv",
  NAdder = "DoubleCensored/All_Island_Data/Data/NorthAdderly_clean.csv",
  UCay   = "DoubleCensored/All_Island_Data/Data/UcayThru2019_CleanOct2020.csv",
  WhiteB = "DoubleCensored/All_Island_Data/Data/WhiteBay_clean.csv"
)

AllData_raw <- map2_df(island_files, names(island_files), ~{
  read_csv(.x, col_types = cols(.default = "c")) %>% mutate(Island = .y)
})

# --- 2. CLEANING & STATUS ---
AllData <- AllData_raw %>%
  mutate(
    CapDate = parse_date_time(date, orders = c("dmy", "mdy", "ymd")),
    CalYear = year(CapDate),
    svl = as.numeric(svl),
    age = as.numeric(age),
    sex = case_when(grepl("^f", sex, ignore.case = TRUE) ~ "F",
                    grepl("^m", sex, ignore.case = TRUE) ~ "M",
                    TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(CapDate), !is.na(animal_id))

# Detect Status
dead_keywords <- "dead|deceased|died|carcass"
row_text <- do.call(paste, c(AllData, sep = " "))
AllData$Status <- if_else(grepl(dead_keywords, row_text, ignore.case = TRUE), 2, 1)

AllData <- AllData %>%
  mutate(ID_raw = paste0(Island, animal_id)) %>%
  group_by(ID_raw) %>%
  mutate(sex = (na.omit(sex))[1]) %>% 
  filter(!is.na(sex)) %>%
  filter(n_distinct(feeding, na.rm = TRUE) <= 1) %>%
  filter(!(Island == "Bitter" & feeding == "none")) %>%
  ungroup()

AllData$ID <- as.numeric(as.factor(AllData$ID_raw))
nind <- length(unique(AllData$ID))

# --- 3. GROWTH CALIBRATION (SVL -> Age) ---
calibration_data <- AllData %>%
  filter(!is.na(age), !is.na(svl), age > 0) %>%
  # IMPORTANT: If your data has "known ages" > 30, cap them to 30 for the model
  mutate(age = pmin(age, 30)) %>% 
  mutate(log_age = log(age), sex_f = as.factor(sex))

fit_min_log <- rq(log_age ~ svl + I(svl^2) + sex_f, tau = 0.05, data = calibration_data)
fit_max_log <- rq(log_age ~ svl + I(svl^2) + sex_f, tau = 0.95, data = calibration_data)

# --- 4. INDIVIDUAL LIFE HISTORIES ---
growth_info <- AllData %>%
  group_by(ID) %>%
  summarise(
    island = Island[1], sex = sex[1], feeding = feeding[1],
    first_svl = svl[!is.na(svl)][1],
    first_svl_year = CalYear[!is.na(svl)][1],
    # Cap known age at 30
    known_age = pmin(age[!is.na(age)][1], 30),
    known_birth_cal = if_else(!is.na(known_age), CalYear[!is.na(age)][1] - known_age, NA_real_),
    total_growth = if(any(!is.na(svl))) max(svl, na.rm=T) - min(svl, na.rm=T) else 0,
    tF_cal = min(CalYear),
    tL_cal = max(CalYear[Status == 1]), 
    tD_cal = if(any(Status == 2)) max(CalYear[Status == 2]) else NA,
    .groups = "drop"
  ) %>% arrange(ID)

# --- 5. PREDICT BIRTH WINDOWS WITH HARD BIOLOGICAL CAP ---
growth_info <- growth_info %>%
  mutate(sex_f = as.factor(sex)) %>%
  mutate(
    # Raw predictions from regression
    p_min_raw = if_else(!is.na(first_svl), exp(predict(fit_min_log, newdata = list(svl = first_svl, sex_f = sex_f))), 0),
    p_max_raw = if_else(!is.na(first_svl), exp(predict(fit_max_log, newdata = list(svl = first_svl, sex_f = sex_f))), 30)
  ) %>%
  mutate(
    # APPLY HARD BIOLOGICAL CAP (Max age = 30)
    p_max = pmin(p_max_raw, 30),
    p_min = pmin(p_min_raw, p_max - 1)
  ) %>%
  mutate(
    # Logic Gates
    p_max = if_else(!is.na(first_svl) & first_svl < 25, pmin(p_max, 4), p_max),
    p_max = if_else(total_growth > 5, pmin(p_max, 15), p_max)
  ) %>%
  mutate(pred_min_age = pmax(0, p_min), pred_max_age = pmax(pred_min_age + 1, p_max))

# --- 6. DYNAMIC TIMELINE ALIGNMENT ---
cintB_cal <- matrix(NA, nrow = nind, ncol = 2)
for(i in 1:nind) {
  if(!is.na(growth_info$known_birth_cal[i])) {
    cintB_cal[i, ] <- growth_info$known_birth_cal[i] + c(-1, 0)
  } else {
    ref_yr <- if(!is.na(growth_info$first_svl_year[i])) growth_info$first_svl_year[i] else growth_info$tF_cal[i]
    cintB_cal[i, 1] <- floor(ref_yr - growth_info$pred_max_age[i])
    cintB_cal[i, 2] <- floor(ref_yr - growth_info$pred_min_age[i])
  }
}

# Now Origin will be ~1950 (1980 - 30)
Origin <- min(cintB_cal[, 1], na.rm = TRUE)

# Shift relative to Origin
cintB_adj <- cintB_cal - Origin + 1
cintD_adj <- matrix(NA, nrow = nind, ncol = 2)
for(i in 1:nind) {
  cintD_adj[i, 1] <- growth_info$tL_cal[i] - Origin + 1
  cintD_adj[i, 2] <- if_else(!is.na(growth_info$tD_cal[i]), growth_info$tD_cal[i] - Origin + 1, cintD_adj[i, 1])
}

# Final logical check for cintB
for(i in 1:nind) {
  tF_adj <- growth_info$tF_cal[i] - Origin + 1
  cintB_adj[i, 2] <- min(cintB_adj[i, 2], tF_adj)
  if(cintB_adj[i, 1] >= cintB_adj[i, 2]) cintB_adj[i, 1] <- cintB_adj[i, 2] - 1
  cintB_adj[i, 1] <- pmax(1, cintB_adj[i, 1])
}

# --- 7. CAPTURE HISTORY & EFFORT ---
tMax_adj <- max(AllData$CalYear) - Origin + 1
CH <- matrix(0, nrow = nind, ncol = tMax_adj)
for(r in 1:nrow(AllData)) {
  CH[AllData$ID[r], AllData$CalYear[r] - Origin + 1] <- 1
}

id_df_final <- growth_info %>%
  mutate(sex_num = if_else(sex == "M", 1, 0),
         feed_num = if_else(feeding == "none", 0, 1))

ordered_islands <- levels(as.factor(id_df_final$island))
n_islands <- length(ordered_islands)
island_idx <- as.numeric(factor(id_df_final$island, levels = ordered_islands))

effort_matrix <- matrix(0, nrow = n_islands, ncol = tMax_adj)
actual_surveys <- AllData %>% 
  distinct(Island, CalYear) %>% 
  mutate(idx = as.numeric(factor(Island, levels = ordered_islands)),
         yr_adj = CalYear - Origin + 1)
for(r in 1:nrow(actual_surveys)) { effort_matrix[actual_surveys$idx[r], actual_surveys$yr_adj[r]] <- 1 }
cum_effort_adj <- t(apply(effort_matrix, 1, cumsum))
cum_effort_adj <- cbind(rep(0, n_islands), cum_effort_adj)

island_starts_cal <- AllData %>% group_by(Island) %>% summarise(S = min(CalYear), .groups="drop")
is_start_vec_adj <- island_starts_cal$S[match(id_df_final$island, island_starts_cal$Island)] - Origin + 1

# --- 8. DIAGNOSTICS ---
message("Origin Year (Model Year 1): ", Origin)
message("Average Birth Window: ", round(mean(cintB_adj[,2]-cintB_adj[,1]), 2), " years")

# --- 1. PREPARE VECTORS ---
y_vec <- as.numeric(rowSums(CH))
sex_vec <- id_df_final$sex_num
feed_vec <- id_df_final$feed_num

# --- 2. DEFINE CENSORING STATUS DIRECTLY ---
# Logic: If cintD[i,1] == cintD[i,2], it's right-censored (Value 2)
# If cintD[i,1] < cintD[i,2], it's interval-censored (Value 1)
censoredD_vec <- ifelse(cintD_adj[,1] == cintD_adj[,2], 2, 1)

# --- 3. CLEANUP ENVIRONMENT ---
# Keep only what is strictly necessary for the NIMBLE run
rm(list = setdiff(ls(), c("nind", "tMax_adj", "CH", "cintB_adj", "cintD_adj", 
                          "id_df_final", "island_idx", "n_islands", 
                          "cum_effort_adj", "is_start_vec_adj", "Origin",
                          "y_vec", "sex_vec", "feed_vec", "censoredD_vec")))

message("Setup verified. Ready to run NIMBLE.")
message("Cleanup complete. Model ready for ", nind, " individuals over ", tMax_adj, " years.")
