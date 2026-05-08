## ============================================================
## FINAL CORRECTED DATA-PREP PIPELINE
## SVL-tightened birth windows + 40-year cap
## 2-column death censoring matrix for NIMBLE dinterval()
## ============================================================
library(tidyverse)
library(lubridate)
library(quantreg)
library(data.table)
library(nimble)

rm(list = ls())


## ============================================================
## 1. LOAD RAW DATA
## ============================================================

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
  read_csv(.x, col_types = cols(.default = "c")) %>%
    mutate(Island = .y)
})


## ============================================================
## 2. CLEANING, DATE PARSING, STATUS
## ============================================================

AllData <- AllData_raw %>%
  mutate(
    CapDate = parse_date_time(date, orders = c("dmy", "mdy", "ymd")),
    CalYear = year(CapDate),
    
    svl = as.numeric(svl),
    age = as.numeric(age),
    
    sex = case_when(
      grepl("^f", sex, ignore.case = TRUE) ~ "F",
      grepl("^m", sex, ignore.case = TRUE) ~ "M",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(CapDate), !is.na(CalYear), !is.na(animal_id))


## Detect rows that indicate dead animals
dead_keywords <- "dead|deceased|died|carcass"

row_text <- do.call(
  paste,
  c(AllData, sep = " ")
)

AllData$Status <- if_else(
  grepl(dead_keywords, row_text, ignore.case = TRUE),
  2,   # dead / found dead
  1    # alive capture / live observation
)

## Individual ID and filtering
AllData <- AllData %>%
  mutate(ID_raw = paste0(Island, "_", animal_id)) %>%
  group_by(ID_raw) %>%
  mutate(
    sex = {
      sx <- na.omit(sex)
      if(length(sx) == 0) NA_character_ else sx[1]
    }
  ) %>%
  ungroup() %>%
  filter(!is.na(sex)) %>%
  group_by(ID_raw) %>%
  filter(n_distinct(feeding, na.rm = TRUE) <= 1) %>%
  ungroup() %>%
  filter(!(Island == "Bitter" & feeding == "none"))


## Numeric ID
AllData <- AllData %>%
  mutate(ID = as.numeric(as.factor(ID_raw)))

nind <- length(unique(AllData$ID))

## ------------------------------------------------------------
## REMOVE POST-DEATH RECORDS
## ------------------------------------------------------------
death_years <- AllData %>%
  filter(Status == 2) %>%
  group_by(ID_raw) %>%
  summarise(
    first_dead_year = min(CalYear, na.rm = TRUE),
    .groups = "drop"
  )

post_death_records <- AllData %>%
  left_join(death_years, by = "ID_raw") %>%
  filter(
    !is.na(first_dead_year),
    CalYear > first_dead_year
  )

if(nrow(post_death_records) > 0) {
  message("Removing ", nrow(post_death_records), " post-death records.")
  
  print(
    post_death_records %>%
      select(
        ID_raw,
        Island,
        animal_id,
        date,
        CapDate,
        CalYear,
        Status,
        first_dead_year,
        sex,
        age,
        svl,
        feeding
      ) %>%
      arrange(ID_raw, CalYear, CapDate),
    n = Inf
  )
}

AllData <- AllData %>%
  left_join(death_years, by = "ID_raw") %>%
  filter(
    is.na(first_dead_year) |
      CalYear <= first_dead_year
  ) %>%
  select(-first_dead_year)


## ============================================================
## 3. GROWTH CALIBRATION: SVL -> AGE
## ============================================================

calibration_data <- AllData %>%
  filter(!is.na(age), !is.na(svl), age > 0) %>%
  mutate(
    age = pmin(age, 30),
    log_age = log(age),
    sex_f = as.factor(sex)
  )

if(nrow(calibration_data) == 0) {
  stop("No usable calibration data for SVL-age quantile regression.")
}

fit_min_log <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.05,
  data = calibration_data
)

fit_max_log <- rq(
  log_age ~ svl + I(svl^2) + sex_f,
  tau = 0.95,
  data = calibration_data
)


## ============================================================
## 4. INDIVIDUAL LIFE HISTORIES
## ============================================================

growth_info <- AllData %>%
  arrange(ID, CalYear, CapDate) %>%
  group_by(ID) %>%
  summarise(
    ID_raw = ID_raw[1],
    island = Island[1],
    sex = sex[1],
    feeding = feeding[1],
    
    first_svl = {
      x <- svl[!is.na(svl)]
      if(length(x) == 0) NA_real_ else x[1]
    },
    
    first_svl_year = {
      x <- CalYear[!is.na(svl)]
      if(length(x) == 0) NA_real_ else x[1]
    },
    
    known_age = {
      x <- age[!is.na(age)]
      if(length(x) == 0) NA_real_ else pmin(x[1], 40)
    },
    
    known_age_year = {
      x <- CalYear[!is.na(age)]
      if(length(x) == 0) NA_real_ else x[1]
    },
    
    total_growth = if(any(!is.na(svl))) {
      max(svl, na.rm = TRUE) - min(svl, na.rm = TRUE)
    } else {
      0
    },
    
    tF_cal = min(CalYear, na.rm = TRUE),
    
    ## Last confirmed alive capture/observation
    tL_cal = {
      x <- CalYear[Status == 1]
      if(length(x) == 0) NA_real_ else max(x, na.rm = TRUE)
    },
    
    ## Found-dead year, if any
    tD_cal = {
      x <- CalYear[Status == 2]
      if(length(x) == 0) NA_real_ else max(x, na.rm = TRUE)
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    known_birth_cal = if_else(
      !is.na(known_age) & !is.na(known_age_year),
      known_age_year - known_age,
      NA_real_
    )
  ) %>%
  arrange(ID)


## Safety checks
if(any(is.na(growth_info$tF_cal))) {
  stop("Some individuals have missing first capture year.")
}

if(any(is.na(growth_info$tL_cal))) {
  bad <- growth_info %>% filter(is.na(tL_cal))
  print(bad)
  stop("Some individuals have no live observation, so tL_cal is missing.")
}

if(any(!is.na(growth_info$tD_cal) & growth_info$tD_cal <= growth_info$tL_cal)) {
  bad <- growth_info %>%
    filter(!is.na(tD_cal), tD_cal <= tL_cal)
  
  print(bad)
  
  stop(
    "Some individuals have tD_cal <= tL_cal. ",
    "For interval censoring you need last_alive < found_dead. ",
    "Check whether dead rows are also being counted as live observations."
  )
}


## ============================================================
## 5. PREDICT BIRTH WINDOWS WITH HARD BIOLOGICAL CAP
## ============================================================

growth_info <- growth_info %>%
  mutate(sex_f = as.factor(sex)) %>%
  mutate(
    p_min_raw = if_else(
      !is.na(first_svl),
      as.numeric(exp(predict(
        fit_min_log,
        newdata = data.frame(
          svl = first_svl,
          sex_f = sex_f
        )
      ))),
      0
    ),
    
    p_max_raw = if_else(
      !is.na(first_svl),
      as.numeric(exp(predict(
        fit_max_log,
        newdata = data.frame(
          svl = first_svl,
          sex_f = sex_f
        )
      ))),
      40
    )
  ) %>%
  mutate(
    ## Hard biological cap
    p_max = pmin(p_max_raw, 40),
    p_min = pmin(p_min_raw, p_max - 1),
    
    ## Logic gates
    p_max = if_else(
      !is.na(first_svl) & first_svl < 25,
      pmin(p_max, 4),
      p_max
    ),
    
    p_max = if_else(
      total_growth > 5,
      pmin(p_max, 15),
      p_max
    ),
    
    pred_min_age = pmax(0, p_min),
    pred_max_age = pmax(pred_min_age + 1, p_max)
  )


## ============================================================
## 6. BUILD CALENDAR-SCALE BIRTH WINDOWS
## ============================================================

cintB_cal <- matrix(NA_real_, nrow = nind, ncol = 2)

for(i in seq_len(nind)) {
  
  if(!is.na(growth_info$known_birth_cal[i])) {
    
    ## Known age gives a tight birth interval
    ## Birth is allowed in the year before / up to estimated birth year.
    cintB_cal[i, ] <- growth_info$known_birth_cal[i] + c(-1, 0)
    
  } else {
    
    ## Otherwise infer from first SVL, or first capture if first SVL missing
    ref_yr <- if(!is.na(growth_info$first_svl_year[i])) {
      growth_info$first_svl_year[i]
    } else {
      growth_info$tF_cal[i]
    }
    
    cintB_cal[i, 1] <- floor(ref_yr - growth_info$pred_max_age[i])
    cintB_cal[i, 2] <- floor(ref_yr - growth_info$pred_min_age[i])
  }
  
  if(any(is.na(cintB_cal[i, ]))) {
    stop("Missing birth interval for individual ID ", growth_info$ID[i])
  }
  
  if(cintB_cal[i, 1] >= cintB_cal[i, 2]) {
    cintB_cal[i, 1] <- cintB_cal[i, 2] - 1
  }
}


## ============================================================
## 7. DEFINE DYNAMIC TIME ORIGIN
## ============================================================
## Model year 1 corresponds to Origin.
## This keeps all model times positive.

Origin <- min(cintB_cal[, 1], na.rm = TRUE)


## ============================================================
## 8. SHIFT BIRTH WINDOWS TO MODEL TIME
## ============================================================

cintB_adj <- cintB_cal - Origin + 1


## Tighten birth windows so birth must occur before first capture.
for(i in seq_len(nind)) {
  
  tF_adj <- growth_info$tF_cal[i] - Origin + 1
  
  ## Birth upper bound cannot be after first capture
  cintB_adj[i, 2] <- min(cintB_adj[i, 2], tF_adj)
  
  ## Birth lower bound should be at least model year 1
  cintB_adj[i, 1] <- max(1, cintB_adj[i, 1])
  
  ## Ensure positive-width birth interval
  if(cintB_adj[i, 1] >= cintB_adj[i, 2]) {
    cintB_adj[i, 1] <- max(1, cintB_adj[i, 2] - 1)
  }
  
  ## Final rescue if first capture is at model year 1
  if(cintB_adj[i, 1] >= cintB_adj[i, 2]) {
    cintB_adj[i, 1] <- 1
    cintB_adj[i, 2] <- 2
  }
  
  if(cintB_adj[i, 1] >= cintB_adj[i, 2]) {
    stop("Invalid adjusted birth interval for individual ID ", growth_info$ID[i])
  }
}


## ============================================================
## 9. BUILD 2-COLUMN DEATH CENSORING MATRIX
## ============================================================
## NIMBLE interpretation:
##
##   censoredD_vec[i] ~ dinterval(tD[i], cintD_adj[i, 1:2])
##
## With 2 cutpoints:
##
##   censoredD_vec = 0: tD <= cintD[,1]
##   censoredD_vec = 1: cintD[,1] < tD <= cintD[,2]
##   censoredD_vec = 2: tD > cintD[,2]

cintD_cal <- matrix(NA_real_, nrow = nind, ncol = 2)
censoredD_vec <- integer(nind)

for(i in seq_len(nind)) {
  
  if(is.na(growth_info$tD_cal[i])) {
    
    ## Right-censored: death is after last known alive
    cintD_cal[i, 1] <- 0
    cintD_cal[i, 2] <- growth_info$tL_cal[i]
    censoredD_vec[i] <- 2
    
  } else {
    
    ## Interval-censored: last alive < death <= found dead
    cintD_cal[i, 1] <- growth_info$tL_cal[i]
    cintD_cal[i, 2] <- growth_info$tD_cal[i]
    censoredD_vec[i] <- 1
  }
}


## Shift death censoring matrix to model time
cintD_adj <- cintD_cal

cintD_adj[, 1] <- ifelse(
  cintD_cal[, 1] == 0,
  0,
  cintD_cal[, 1] - Origin + 1
)

cintD_adj[, 2] <- cintD_cal[, 2] - Origin + 1


## Safety checks for death censoring
if(any(is.na(cintD_adj))) {
  stop("NA values found in cintD_adj.")
}

if(any(cintD_adj[, 2] <= 0)) {
  bad <- which(cintD_adj[, 2] <= 0)
  print(bad)
  stop("Some death/censoring upper bounds are <= 0.")
}

if(any(censoredD_vec == 1 & cintD_adj[, 2] <= cintD_adj[, 1])) {
  bad <- which(censoredD_vec == 1 & cintD_adj[, 2] <= cintD_adj[, 1])
  
  print(
    cbind(
      i = bad,
      ID = growth_info$ID[bad],
      ID_raw = growth_info$ID_raw[bad],
      island = growth_info$island[bad],
      tL_cal = cintD_cal[bad, 1],
      tD_cal = cintD_cal[bad, 2],
      tL_adj = cintD_adj[bad, 1],
      tD_adj = cintD_adj[bad, 2]
    )
  )
  
  stop("Invalid interval-censored death windows found: upper <= lower.")
}

if(any(censoredD_vec == 2 & cintD_adj[, 1] != 0)) {
  stop("Right-censored individuals should have cintD_adj[,1] == 0.")
}

if(any(!censoredD_vec %in% c(1, 2))) {
  stop("censoredD_vec should only contain 1 or 2.")
}

# cleanup

rm(AllData_raw, calibration_data, cintB_cal, cintD_cal, fit_max_log, fit_min_log, post_death_records)
rm(row_text, dead_keywords, tF_adj)
rm(death_years)
rm(island_files)

## ============================================================
## 10. INDIVIDUAL COVARIATE DATAFRAME
## ============================================================

sex <- growth_info %>%
  mutate(sex = if_else(sex == "M", 1, 0)) %>%
  pull(sex)

feed <- growth_info %>%
  mutate(feeding = if_else(feeding == "none", 0, 1)) %>%
  pull(feeding)

ordered_islands <- levels(as.factor(growth_info$island))
n_islands <- length(ordered_islands)

island_idx <- as.numeric(
  factor(growth_info$island, levels = ordered_islands)
)

## ============================================================
## 11. CAPTURE HISTORY MATRIX
## ============================================================

tMax <- max(AllData$CalYear, na.rm = TRUE) - Origin + 1

CH <- matrix(0, nrow = nind, ncol = tMax)

for(r in seq_len(nrow(AllData))) {
  if(AllData$Status[r] == 1){
    i <- AllData$ID[r]
    t <- AllData$CalYear[r] - Origin + 1
    CH[i, t] <- 1
  } else {
    i <- AllData$ID[r]
    t <- AllData$CalYear[r] - Origin + 1
    CH[i, t] <- 0
  }
}

## Number of observed years/captures per animal
y <- as.numeric(rowSums(CH))


## ============================================================
## 12. ISLAND-LEVEL EFFORT MATRIX
## ============================================================

effort_matrix <- matrix(0, nrow = n_islands, ncol = tMax)

actual_surveys <- AllData %>%
  distinct(Island, CalYear) %>%
  mutate(
    idx = as.numeric(factor(Island, levels = ordered_islands)),
    yr_adj = CalYear - Origin + 1) %>%
  filter(!is.na(idx), !is.na(yr_adj), yr_adj >= 1, yr_adj <= tMax)

for(r in seq_len(nrow(actual_surveys))) {
  effort_matrix[
    actual_surveys$idx[r],
    actual_surveys$yr_adj[r]
  ] <- 1
}


## Cumulative island effort
cum_effort_adj <- t(apply(effort_matrix, 1, cumsum))

## Island start year on adjusted time scale
island_starts_cal <- AllData %>%
  group_by(Island) %>%
  summarise(S = min(CalYear, na.rm = TRUE), .groups = "drop")

island_start <- island_starts_cal$S[match(growth_info$island, island_starts_cal$Island)] - Origin + 1

## ============================================================
## 13. FINAL MODEL VECTORS
## ============================================================

tF <- growth_info$tF_cal - Origin + 1
tL <- growth_info$tL_cal - Origin + 1
tD <- ifelse(is.na(growth_info$tD_cal), NA, growth_info$tD_cal - Origin + 1)

## ============================================================
## 14. INITIAL VALUES FOR LATENT BIRTH AND DEATH
## ============================================================
## Model uses:
##   tB[i]     = latent birth time
##   tstar[i] = lifespan / age at death
##   tD[i]     = tB[i] + tstar[i]

tB_init <- rowMeans(cintB_adj)
tD_init <- numeric(nind)

for(i in seq_len(nind)) {
  if(censoredD_vec[i] == 1) {
    ## interval-censored death: initialise inside interval
    tD_init[i] <- mean(cintD_adj[i, ])
  } else if(censoredD_vec[i] == 2) {
    ## right-censored: initialise after censoring point
    tD_init[i] <- cintD_adj[i, 2] + 1
  } else {
    stop("Unexpected censoring category for individual ", i)
  }
  
  ## Ensure death is after birth
  if(tD_init[i] <= tB_init[i]) {
    tD_init[i] <- tB_init[i] + 1
  }
}

tstar_init <- pmax(1, tD_init - tB_init)

## ============================================================
## 15. FINAL SANITY CHECKS
## ============================================================

stopifnot(
  nrow(cintB_adj) == nind,
  ncol(cintB_adj) == 2,
  nrow(cintD_adj) == nind,
  ncol(cintD_adj) == 2,
  
  length(censoredD_vec) == nind,
  length(sex) == nind,
  length(feed) == nind,
  length(island_idx) == nind,
  length(y) == nind,
  length(tB_init) == nind,
  length(tstar_init) == nind,
  
  all(cintB_adj[, 1] < cintB_adj[, 2]),
  all(cintD_adj[, 2] > 0),
  all(censoredD_vec %in% c(1, 2)),
  all(tB_init > 0),
  all(tstar_init > 0)
)

## Interval-censored only: lower must be below upper
stopifnot(
  all(cintD_adj[censoredD_vec == 1, 1] < cintD_adj[censoredD_vec == 1, 2])
)

## Right-censored only: first cutpoint should be zero
stopifnot(
  all(cintD_adj[censoredD_vec == 2, 1] == 0)
)

## ------------------------------------------------------------
## Inspect birth/death conflicts
## ------------------------------------------------------------

bad_birth_death <- which(cintB_adj[, 2] >= cintD_adj[, 2])

bad_birth_death_df <- tibble(
  i = bad_birth_death,
  ID = growth_info$ID_raw[bad_birth_death],
  island = growth_info$island[bad_birth_death],
  sex = growth_info$sex[bad_birth_death],
  feeding = growth_info$feeding[bad_birth_death],
  last_seen = growth_info$tL_cal[bad_birth_death] - Origin + 1,
  censor_cat = censoredD_vec[bad_birth_death],
  birth_low = cintB_adj[bad_birth_death, 1],
  birth_high = cintB_adj[bad_birth_death, 2],
  death_low = cintD_adj[bad_birth_death, 1],
  death_high = cintD_adj[bad_birth_death, 2],
  birth_low_cal = cintB_adj[bad_birth_death, 1] + Origin - 1,
  birth_high_cal = cintB_adj[bad_birth_death, 2] + Origin - 1,
  death_low_cal = ifelse(
    cintD_adj[bad_birth_death, 1] == 0,
    NA,
    cintD_adj[bad_birth_death, 1] + Origin - 1
  ),
  death_high_cal = cintD_adj[bad_birth_death, 2] + Origin - 1
)

print(bad_birth_death_df, n = Inf)

cintB_adj[bad_birth_death, 2] <- cintB_adj[bad_birth_death, 2] - 1 

## Optional warning if birth upper bound is after death/censoring upper point
bad_birth_death <- which(cintB_adj[, 2] >= cintD_adj[, 2])
if(length(bad_birth_death) > 0) {
  warning(
    length(bad_birth_death),
    " individuals have birth upper bound >= death/censoring upper bound. ",
    "This may be okay for right-censored animals with very short histories, ",
    "but inspect bad_birth_death if model initialisation fails."
  )
}

## ============================================================
## 16. DIAGNOSTIC OUTPUT
## ============================================================

message("------------------------------------------------------------")
message("Setup verified.")
message("Origin Year / Model Year 1: ", Origin)
message("Individuals: ", nind)
message("Model duration: ", tMax, " years")
message("Number of islands: ", n_islands)
message("Average birth window width: ",
        round(mean(cintB_adj[, 2] - cintB_adj[, 1]), 2),
        " years")
message("Death censoring categories:")
print(table(censoredD_vec))
message("Birth interval summary:")
print(summary(cintB_adj))
message("Death censoring interval summary:")
print(summary(cintD_adj))
message("Initial lifespan summary:")
print(summary(tstar_init))
message("------------------------------------------------------------")


## ============================================================
## 17. FINAL PRE-NIMBLE FIXES
## ============================================================
## You adjusted some birth upper bounds above, so rebuild initial values.
tB_init <- rowMeans(cintB_adj)
tD_init <- numeric(nind)

for(i in seq_len(nind)) {
  
  if(censoredD_vec[i] == 1) {
    ## interval-censored death: initialise inside interval
    tD_init[i] <- mean(cintD_adj[i, ])
    
  } else if(censoredD_vec[i] == 2) {
    ## right-censored: initialise after censoring point
    tD_init[i] <- cintD_adj[i, 2] + 1
    
  } else {
    stop("Unexpected censoring category for individual ", i)
  }
  
  ## Ensure death is after birth
  if(tD_init[i] <= tB_init[i]) {
    tD_init[i] <- tB_init[i] + 1
  }
}

tstar_init <- pmax(1, tD_init - tB_init)

## ------------------------------------------------------------
## cum_effort needs a leading zero column
## ------------------------------------------------------------
## If effort_matrix has tMax columns, cum_effort_for_model has tMax + 1.
## Column 1 = cumulative effort before year 1.
## Column t + 1 = cumulative effort up to and including year t.

cum_effort_for_model <- cbind(0, cum_effort_adj)

stopifnot(
  nrow(cum_effort_for_model) == n_islands,
  ncol(cum_effort_for_model) == tMax + 1
)

## ------------------------------------------------------------
## Island start time and index
## ------------------------------------------------------------
## island_start is already in model-time occasion units.
## We keep one version for left truncation in continuous time,
## and one integer version for indexing the effort matrix.

island_start_time <- island_start

island_start_idx <- pmax(
  1,
  pmin(tMax, floor(island_start_time) + 1)
)

## ------------------------------------------------------------
## Final checks before NIMBLE
## ------------------------------------------------------------

stopifnot(
  length(sex) == nind,
  length(feed) == nind,
  length(island_idx) == nind,
  length(island_start_time) == nind,
  length(island_start_idx) == nind,
  length(y) == nind,
  all(island_idx >= 1),
  all(island_idx <= n_islands),
  all(island_start_idx >= 1),
  all(island_start_idx <= tMax),
  all(cintB_adj[, 1] < cintB_adj[, 2]),
  all(cintD_adj[, 2] > 0),
  all(censoredD_vec %in% c(1, 2)),
  all(tB_init >= cintB_adj[, 1]),
  all(tB_init <= cintB_adj[, 2]),
  all(tstar_init > 0)
)


## ============================================================
## 18. LOAD CUSTOM DISTRIBUTIONS
## ============================================================

source("ModelComparison_FUNCTIONS.R")
source("Distributions/Dist_GompertzLB.R")
source("Distributions/Dist_Gompertz.R")
source("Distributions/Dist_GompertzNim.R")

## ============================================================
## MODEL B: GOMPERTZ WITH PARTIAL POOLING ON ISLAND EFFECTS FOR b
## ============================================================

code_B <- nimbleCode({
  
  for (i in 1:nind) {
    
    ## Latent birth/entry time
    tB[i] ~ dunif(cintB[i, 1], cintB[i, 2])
    L[i] <- max(0, island_start_time[i] - tB[i])
    ## Left-truncated Gompertz age-at-death
    tstar[i] ~ dGompertzLB(amult[i], bmult[i], lowerBound = L[i])
    tD[i] <- tB[i] + tstar[i]
    censoredD[i] ~ dinterval(tD[i], cintD[i, 1:2])
    
    log(amult[i]) <- log(a) +
      betaSEX[1] * sex[i] +
      betaISLAND_a[island_idx[i]]
    
    log(bmult[i]) <- log(b) +
      betaSEX[2] * sex[i] +
      betaISLAND_b[island_idx[i]]
    

    tB_year_idx[i] <- max(1, min(tMax, floor(tB[i]) + 1))
    nm_start[i] <- max(tB_year_idx[i], island_start_idx[i])
    end_year[i] <- max(1, min(tMax, floor(tD[i]) + 1))
    nMpos[i] <- cum_effort[island_idx[i], end_year[i] + 1] - cum_effort[island_idx[i], nm_start[i]]
    nMpos_safe[i] <- max(y[i], nMpos[i])
    
    pd_raw[i] <- exp(y[i] * log(mean.p + 1.0E-10) + (nMpos_safe[i] - y[i]) * log(1 - mean.p + 1.0E-10))
    pd[i] <- min(0.999999, pd_raw[i])
    dind[i] ~ dbern(pd[i])
  }
  
  for (k in 1:2) {
    betaSEX[k]  ~ dnorm(0, 1)
  }
  
  for (r in 1:n_islands) {
    betaISLAND_a[r] ~ dnorm(0, 1)
  }
  
  sigma_b_island ~ dunif(0, 2)
  for (r in 1:n_islands) {
    betaISLAND_b_raw[r] ~ dnorm(0, 1)
    betaISLAND_b[r] <- sigma_b_island * betaISLAND_b_raw[r]
  }
  
  a ~ dexp(1)
  b ~ dexp(1)
  
  mean.p ~ dunif(0, 1)
})

consts <- list(
  nind = nind,
  tMax = tMax,
  n_islands = n_islands,
  cintB = as.matrix(cintB_adj),
  cintD = as.matrix(cintD_adj),
  sex = as.numeric(sex),
  island_idx = as.numeric(island_idx),
  island_start_time = as.numeric(island_start_time),
  island_start_idx = as.numeric(island_start_idx),
  y = as.numeric(y)
)

data_B <- list(
  censoredD = as.numeric(censoredD_vec),
  dind = rep(1, nind),
  cum_effort = as.matrix(cum_effort_for_model)
)

make_inits_B <- function(chain_id = 1) {
  
  ## Jitter tB inside each birth interval
  tB_start <- runif(
    nind,
    min = cintB_adj[, 1],
    max = cintB_adj[, 2]
  )
  
  ## Initialise tD consistently with censoring
  tD_start <- numeric(nind)
  
  for(i in seq_len(nind)) {
    
    if(censoredD_vec[i] == 1) {
      ## interval-censored: inside death interval
      tD_start[i] <- runif(
        1,
        min = cintD_adj[i, 1] + 0.001,
        max = cintD_adj[i, 2] - 0.001
      )
    }
    
    if(censoredD_vec[i] == 2) {
      ## right-censored: after censoring point
      tD_start[i] <- cintD_adj[i, 2] + runif(1, 1, 5)
    }
    
    ## Ensure death is after birth
    if(tD_start[i] <= tB_start[i]) {
      tD_start[i] <- tB_start[i] + runif(1, 1, 5)
    }
  }
  
  tstar_start <- pmax(0.001, tD_start - tB_start)
  
  list(
    tB = tB_start,
    tstar = tstar_start,
    
    ## Baseline Gompertz parameters
    a = exp(rnorm(1, log(0.01), 0.3)),
    b = exp(rnorm(1, log(0.05), 0.3)),
    
    ## Fixed effects
    betaSEX = rnorm(2, 0, 0.2),
    
    ## Island effects on a
    betaISLAND_a = rnorm(n_islands, 0, 0.2),
    
    ## Partial pooling on b
    sigma_b_island = runif(1, 0.05, 0.7),
    betaISLAND_b_raw = rnorm(n_islands, 0, 0.2),
    
    ## Detection
    mean.p = runif(1, 0.15, 0.4)
  )
}

inits_list_B <- list(
  make_inits_B(1),
  make_inits_B(2),
  make_inits_B(3)
)

## ============================================================
## BUILD MODEL B
## ============================================================

model_B <- nimbleModel(
  code = code_B,
  constants = consts,
  data = data_B,
  inits = inits_list_B,
  calculate = FALSE
)

initial_lp_B <- model_B$calculate()

message("Initial model B log probability:")
print(initial_lp_B)

## ============================================================
## CONFIGURE MCMC FOR MODEL B
## ============================================================

conf_B <- configureMCMC(
  model_B,
  monitors = c("a", "b", "mean.p", "betaSEX", "betaISLAND_a", "betaISLAND_b", "sigma_b_island"),
  monitors2 = c("tB", "tstar", "tD", "amult", "bmult"),
  enableWAIC = TRUE
)

## More robust samplers for baseline Gompertz parameters
conf_B$removeSamplers(c("a", "b"))
conf_B$addSampler(target = c("a", "b"), type = "AF_slice")

## Slice sampler for the partial-pooling scale
conf_B$removeSamplers("sigma_b_island")
conf_B$addSampler(target = "sigma_b_island", type = "slice")

## Slice samplers for latent continuous states
for (i in 1:nind) {
  conf_B$removeSamplers(paste0('tB[', i, ']'))
  conf_B$removeSamplers(paste0('tstar[', i, ']'))
  conf_B$addSampler(target = c(paste0('tB[', i, ']'), paste0('tstar[', i, ']')), type = "AF_slice")
}

mcmc_B <- buildMCMC(conf_B)

## ============================================================
## COMPILE MODEL B
## ============================================================

cModel_B <- compileNimble(model_B)

cMCMC_B <- compileNimble(mcmc_B, project = model_B)

run_B_test <- runMCMC(
  cMCMC_B,
  niter = 500000,
  nburnin = 90000,
  nchains = 3,
  inits = inits_list_B,
  setSeed = TRUE,
  samplesAsCodaMCMC = TRUE,
  summary = TRUE,
  WAIC = TRUE,
  progressBar = TRUE
)

saveRDS(run_B_test, "Model_simplePooled_samples.rds")

run_B_test$summary
run_B_test$WAIC

MCMCsummary(run_B_test$samples)
MCMCtrace(run_B_test$samples)
