## ============================================================
## DATA HANDLING PIPELINE (CLEANED + CONSISTENT)
## - Reads island CSVs
## - Harmonises columns + adds Island + builds unique ID
## - Cleans sex, defines Status (alive/dead), builds IslandStart
## - Removes individuals with >1 feeding regime (e.g., 6 Bitter dual-feeding)
## - Assigns Species (inornata vs figginsi)
## - Converts ID to numeric index for CH construction
## - Builds CH, censoring intervals (cintB, cintD), y, etc.
## - Saves image
## ============================================================

library(haven)
library(survival)
library(tidyverse)
library(lubridate)
library(nimble)
library(purrr)
library(mclust)
library(GGally)
library(boot)
library(parallel)
library(mcmcplots)
library(coda)
library(stringr)

rm(list = ls())

## -----------------------------
## 1) Load island datasets
## -----------------------------
Bitter <- read_csv("DoubleCensored/All_Island_Data/Data/BitterGuanaCay_clean.csv") %>%
  mutate(date = dmy(date),
         feeding = as.factor(feeding))

FFRC <- read_csv("DoubleCensored/All_Island_Data/Data/FlatRockThru2019_CleanOct2020.csv") %>%
  mutate(date = mdy(date),
         feeding = as.factor(feeding))

Gaulin <- read_csv("DoubleCensored/All_Island_Data/Data/GaulinCay_clean.csv") %>%
  mutate(date = dmy(date)) %>%
  filter(!(Pit == "4264176E36" & `body-mass` == 2070)) %>%
  filter(!(Pit == "442A0F6D30" & `body-mass` == 2460)) %>%
  mutate(feeding = as.factor(feeding))

Leaf <- read_csv("DoubleCensored/All_Island_Data/Data/LeafCayThru2019_CleanOct2020.csv") %>%
  mutate(date = mdy(date),
         feeding = as.factor(feeding))

Noddy <- read_csv("DoubleCensored/All_Island_Data/Data/NoddyCCay_Clean.csv") %>%
  mutate(date = dmy(date),
         feeding = as.factor(feeding))

NAdder <- read_csv("DoubleCensored/All_Island_Data/Data/NorthAdderly_clean.csv") %>%
  mutate(date = dmy(date),
         feeding = as.factor(feeding))

UCay <- read_csv("DoubleCensored/All_Island_Data/Data/UcayThru2019_CleanOct2020.csv") %>%
  mutate(date = mdy(date),
         feeding = as.factor(feeding))

WhiteB <- read_csv("DoubleCensored/All_Island_Data/Data/WhiteBay_clean.csv") %>%
  mutate(date = dmy(date),
         feeding = as.factor(feeding))

## Quick sanity checks (optional)
# lapply(list(Bitter, FFRC, Gaulin, Leaf, Noddy, NAdder, UCay, WhiteB),
#        function(df) list(range = range(df$date, na.rm = TRUE),
#                          feeding = summary(df$feeding)))

## -----------------------------
## 2) Harmonise columns & bind
## -----------------------------
island_names <- c("Bitter", "FFRC", "Gaulin", "Leaf", "Noddy", "NAdder", "UCay", "WhiteB")
island_list  <- mget(island_names)

## Add Island column
island_list <- lapply(names(island_list), function(nm) {
  df <- island_list[[nm]]
  df$Island <- nm
  df
})
names(island_list) <- island_names

## Common columns across all islands
common_cols <- Reduce(intersect, lapply(island_list, colnames))
print(common_cols)

## Columns required to keep a row
NACols <- c("CapDate", "animal_id")

AllData <- do.call(rbind, lapply(island_list, function(df) df[, common_cols, drop = FALSE])) %>%
  rename(CapDate = date) %>%
  filter(rowSums(is.na(.)) != ncol(.)) %>%  # remove fully-NA rows
  filter(complete.cases(select(., all_of(NACols)))) %>%
  mutate(
    sex = as.factor(sex),
    ID  = str_replace_all(paste(Island, animal_id), " ", "")
  )

rm(list = c("Bitter","FFRC","Gaulin","Leaf","Noddy","NAdder","UCay","WhiteB","island_list","common_cols"))

## -----------------------------
## 3) Clean/standardise sex
## -----------------------------
AllData <- AllData %>%
  mutate(
    sex = case_when(
      sex == "female" ~ "F",
      sex == "male"   ~ "M",
      TRUE            ~ as.character(sex)
    )
  ) %>%
  group_by(ID) %>%
  filter(any(sex %in% c("M","F"))) %>%
  ungroup()

get_mode <- function(x) {
  x <- na.omit(x)
  if (length(x) == 0) return(NA)
  tab <- table(x)
  names(tab)[which.max(tab)]
}

AllData <- AllData %>%
  arrange(ID) %>%
  group_by(ID) %>%
  mutate(sex = get_mode(sex)) %>%
  ungroup()

## Drop individuals with only 1 record and sex still NA (belt & braces)
AllData <- AllData %>%
  group_by(ID) %>%
  mutate(n_entries = n()) %>%
  ungroup() %>%
  filter(!(n_entries == 1 & is.na(sex))) %>%
  select(-n_entries)

## -----------------------------
## 4) Status and Year variables
## -----------------------------
## NOTE: This flags dead if any column contains "dead" (broad but matches your prior approach)
AllData <- AllData %>%
  mutate(
    Status = ifelse(
      apply(., 1, function(row) any(grepl("(?i)dead", row, perl = TRUE))),
      2, 1
    )
  ) %>%
  arrange(ID)

AllData$Year <- year(AllData$CapDate)

AllData <- AllData %>%
  group_by(Island) %>%
  mutate(IslandStart = min(Year, na.rm = TRUE)) %>%
  ungroup()

## -----------------------------
## 5) Remove individuals with >1 feeding regime (ignore NA)
## -----------------------------
## (This should catch the 6 Bitter dual-feeding individuals)
multi_feed <- AllData %>%
  group_by(ID) %>%
  summarise(n_feed_types = n_distinct(feeding, na.rm = TRUE), .groups = "drop") %>%
  filter(n_feed_types > 1) %>%
  pull(ID)

## Audit what will be removed
multi_feed_check <- AllData %>%
  filter(ID %in% multi_feed) %>%
  distinct(Island, ID, feeding) %>%
  arrange(Island, ID, feeding)

cat("\n--- Individuals with >1 feeding regime (to be removed) ---\n")
print(multi_feed_check)

cat("\n--- Counts by island for multi-feed individuals ---\n")
AllData %>%
  filter(ID %in% multi_feed) %>%
  distinct(ID, Island) %>%
  count(Island) %>%
  arrange(Island) %>%
  print()

## Remove them
AllData <- AllData %>%
  filter(!(ID %in% multi_feed))

## Verify now “one feeding per island” at the individual level
cat("\n--- Individual-level feeding counts after removal ---\n")
AllData %>%
  distinct(ID, .keep_all = TRUE) %>%
  count(Island, feeding) %>%
  arrange(Island, feeding) %>%
  print()

## ============================================================
## REMOVE Bitter individuals with feeding == "none"
## (drops whole individuals, not just the rows)
## ============================================================
bitter_none_ids <- AllData %>%
  filter(Island == "Bitter", feeding == "none") %>%
  distinct(ID) %>%
  pull(ID)

cat("\nBitter 'none' IDs to remove:", length(bitter_none_ids), "\n")

## optional audit table
AllData %>%
  filter(ID %in% bitter_none_ids) %>%
  distinct(Island, ID, feeding) %>%
  arrange(ID) %>%
  print(n = 200)

## remove those individuals completely
AllData <- AllData %>%
  filter(!(ID %in% bitter_none_ids))

## quick check
AllData %>%
  distinct(ID, .keep_all = TRUE) %>%
  filter(Island == "Bitter") %>%
  count(feeding) %>%
  print()

## -----------------------------
## 6) Species assignment
## -----------------------------
## Reference list for inornata IDs
inorn <- read_csv("Data/Iguana_ID_InornataSp.csv") %>%
  mutate(
    Island = fct_recode(as_factor(Island), FFRC = "FlatRock"),
    ID     = str_c(Island, animal_id, sep = "")
  )

sep_ids <- inorn %>% distinct(ID)

AllData <- AllData %>%
  mutate(
    Species = case_when(
      Island %in% c("FFRC", "UCay") ~ "inornata",
      ID %in% sep_ids$ID            ~ "inornata",
      TRUE                          ~ "figginsi"
    )
  )

## Check for any individuals labelled as both species (should be none)
species_conflicts <- AllData %>%
  group_by(ID) %>%
  summarise(n_species = n_distinct(Species), .groups = "drop") %>%
  filter(n_species > 1)

if (nrow(species_conflicts) > 0) {
  warning("Some IDs have multiple Species labels. Inspect species_conflicts.")
  print(species_conflicts)
}

## -----------------------------
## 7) Reindex time + ID for modelling
## -----------------------------
YearZ <- 1950

## Convert ID to numeric (model index)
AllData$ID <- as.numeric(as.factor(AllData$ID))
AllData <- arrange(AllData, ID)

## Put Year / IslandStart onto “since 1950” scale
AllData$Year       <- AllData$Year - YearZ
AllData$IslandStart <- AllData$IslandStart - YearZ

## Death time per individual (if any Status==2)
AllData <- AllData %>%
  mutate(Death = ifelse(Status == 2, Year, NA_real_)) %>%
  group_by(ID) %>%
  mutate(Death = ifelse(all(is.na(Death)), NA_real_, max(Death, na.rm = TRUE))) %>%
  ungroup()

saveRDS(AllData, "AllData.rds")

## -----------------------------
## 8) Birth times for known-aged individuals
## -----------------------------
AllData <- AllData %>%
  mutate(age = suppressWarnings(as.numeric(age))) %>%
  group_by(ID) %>%
  mutate(ageR = round(age, digits = 0)) %>%
  mutate(
    Birth = ifelse(
      all(is.na(Year)) | all(is.na(ageR)),
      NA_real_,
      min(Year, na.rm = TRUE) - min(ageR, na.rm = TRUE)
    )
  ) %>%
  mutate(
    LastSeen = case_when(
      any(Status == 1, na.rm = TRUE) ~ max(Year[Status == 1], na.rm = TRUE),
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

tKB <- AllData %>% distinct(ID, .keep_all = TRUE) %>% pull(Birth)
tKD <- AllData %>% distinct(ID, .keep_all = TRUE) %>% pull(Death)

## -----------------------------
## 9) Capture History (CH)
## -----------------------------
## CH columns: 1..(last_year_since_1950)
last_year_since_1950 <- max(AllData$Year, na.rm = TRUE)
nind <- length(unique(AllData$ID))

CH <- matrix(0, nrow = nind, ncol = last_year_since_1950)

alive_data <- AllData %>% filter(Status == 1)

## Fill CH with 1 where alive detections occurred
CH[as.matrix(alive_data[, c("ID", "Year")])] <- 1

## First / last detection occasion
tL <- apply(CH, 1, function(x) { w <- which(x == 1); if (length(w)==0) NA_integer_ else max(w) })
tF <- apply(CH, 1, function(x) { w <- which(x == 1); if (length(w)==0) NA_integer_ else min(w) })
names(tF) <- NULL

## Drop any individuals with no alive detections (shouldn’t happen, but safe)
no_det <- which(is.na(tF) | is.na(tL))
if (length(no_det) > 0) {
  warning("Dropping individuals with no alive detections: ", length(no_det))
  keep <- setdiff(1:nind, no_det)
  CH  <- CH[keep, , drop = FALSE]
  tL  <- tL[keep]
  tF  <- tF[keep]
  tKB <- tKB[keep]
  tKD <- tKD[keep]
  nind <- nrow(CH)
}

## -----------------------------
## 10) Censoring intervals
## -----------------------------
## Birth censoring: known births vs unknown (adult) births
cintB <- matrix(NA, nrow = nind, ncol = 2)

known_idx <- !is.na(tKB)
cintB[known_idx, 1] <- tKB[known_idx] - 1
cintB[known_idx, 2] <- tKB[known_idx]

unknown_idx <- is.na(tKB)
cintB[unknown_idx, 1] <- tF[unknown_idx] - 29
cintB[unknown_idx, 2] <- tF[unknown_idx]

## Death censoring: interval for known dead, right-censor for unknown
cintD <- cbind(tL, tKD)
cintD[is.na(tKD), 2] <- cintD[is.na(tKD), 1]
cintD[is.na(tKD), 1] <- 0
censoredD <- ifelse(!is.na(tKD), 1, 2)

tD    <- rep(NA, length(tKD))
tstar <- rep(NA, length(tKD))
dind  <- rep(1, length(tKD))

## Detections count
y <- rowSums(CH)
names(y) <- NULL

## Basic checks
stopifnot(all(tL >= tF))
if (any(!is.na(tKD))) {
  stopifnot(all(tKD[!is.na(tKD)] > tL[!is.na(tKD)]))
}

## -----------------------------
## 11) Save
## -----------------------------
save.image("igs_AllIslands_CleanDH_280426_obsStart.RData")
