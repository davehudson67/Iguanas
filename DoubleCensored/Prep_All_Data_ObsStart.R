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
rm(list=ls())

#source("ModelComparison_FUNCTIONS.R")

## Load data
#Allig <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_Alligatory_Cay.csv")
#Bitter <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_Bitter_Guana_Cay.csv")
#FFRC <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_FFRC_Allens.csv")
#Gaulin <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_Gaulin_Cay.csv")
#Leaf_A <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_Leaf_Cay_Allens.csv")
#Leaf <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_Leaf_Cay.csv")
#Noddy <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_Noddy_Cay.csv")
#NAdder <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_North_Adderly.csv")
#UCay <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_U_Cay_Allens.csv")
#WhiteB <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_White_Bay.csv")

#Allig <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_Alligatory_Cay.csv")
Bitter <- read_csv("DoubleCensored/All_Island_Data/Data/BitterGuanaCay_clean.csv") %>%
  mutate(date = dmy(date)) %>%
  mutate(feeding = as.factor(feeding))
summary(Bitter$feeding)
range(Bitter$date)
FFRC <- read_csv("DoubleCensored/All_Island_Data/Data/FlatRockThru2019_CleanOct2020.csv")%>%
  mutate(date = mdy(date))  %>%
  mutate(feeding = as.factor(feeding))
summary(FFRC$feeding)
range(FFRC$date)
Gaulin <- read_csv("DoubleCensored/All_Island_Data/Data/GaulinCay_clean.csv") %>%
  mutate(date = dmy(date)) %>%
  filter(!(Pit == "4264176E36" & `body-mass` == 2070)) %>%
  filter(!(Pit == "442A0F6D30" & `body-mass` == 2460))  %>%
  mutate(feeding = as.factor(feeding))
summary(Gaulin$feeding)
range(Gaulin$date)
#Leaf_A <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_Leaf_Cay_Allens.csv")
Leaf <- read_csv("DoubleCensored/All_Island_Data/Data/LeafCayThru2019_CleanOct2020.csv") %>%
  mutate(date = mdy(date))  %>%
  mutate(feeding = as.factor(feeding),
         species = "")
summary(Leaf$feeding)
range(Leaf$date)
Noddy <- read_csv("DoubleCensored/All_Island_Data/Data/NoddyCCay_Clean.csv")%>%
  mutate(date = dmy(date)) %>%
  mutate(feeding = as.factor(feeding))
summary(Noddy$feeding)
range(Noddy$date)
NAdder <- read_csv("DoubleCensored/All_Island_Data/Data/NorthAdderly_clean.csv")%>%
  mutate(date = dmy(date))  %>%
  mutate(feeding = as.factor(feeding))
summary(NAdder$feeding)
UCay <- read_csv("DoubleCensored/All_Island_Data/Data/UcayThru2019_CleanOct2020.csv") %>%
  mutate(date = mdy(date))  %>%
  mutate(feeding = as.factor(feeding))
summary(UCay$feeding)
WhiteB <- read_csv("DoubleCensored/All_Island_Data/Data/WhiteBay_clean.csv") %>%
  mutate(date = dmy(date))  %>%
  mutate(feeding = as.factor(feeding))
summary(WhiteB$feeding)

summary(Bitter$feeding)
summary(FFRC$feeding)
summary(Gaulin$feeding)
summary(Leaf$feeding)
summary(NAdder$feeding)
summary(Noddy$feeding)
summary(UCay$feeding)
summary(WhiteB$feeding)

## Consistent cols
#island_list <- list(Allig, Bitter, FFRC, Gaulin, Leaf_A, Leaf, Noddy, NAdder, UCay, WhiteB)
island_list <- list(Bitter, FFRC, Gaulin, Leaf, Noddy, NAdder, UCay, WhiteB)

## Find common column names across all dataframes
common_cols <- Reduce(intersect, lapply(island_list, colnames))
print(common_cols)

## Add Island identifier
#island_names <- c("Allig", "Bitter", "FFRC", "Gaulin", "Leaf_A", "Leaf", "Noddy", "NAdder", "UCay", "WhiteB")
island_names <- c("Bitter", "FFRC", "Gaulin", "Leaf", "Noddy", "NAdder", "UCay", "WhiteB")

## Get the dataframes as a list
island_list <- mget(island_names)

## Add the "Island" column to each dataframe
island_list <- lapply(names(island_list), function(name) {
  df <- island_list[[name]]
  df$Island <- name  # Assign dataframe name as a new column
  df
})

## Assign modified dataframes back to the global environment
names(island_list) <- island_names
list2env(island_list, envir = .GlobalEnv)

## Combine the dfs
common_cols <- Reduce(intersect, lapply(island_list, colnames))

## Remove rows where Capture Date or Pit # are NA
NACols <- c("CapDate", "animal_id")

## Keep only common columns and bind rows
AllData <- do.call(rbind, lapply(island_list, function(df) df[, common_cols, drop = FALSE])) %>%
  rename("CapDate" = "date") %>%
  #rename("Pit" = "animal_id") %>%
  filter(rowSums(is.na(.)) != ncol(.)) %>%  # Remove fully NA rows
  filter(complete.cases(select(., all_of(NACols)))) %>%
  mutate(sex = as.factor(sex)) %>%
  mutate(ID = paste(Island, animal_id)) %>%
  mutate(ID = str_replace_all(ID, " ", ""))
  #filter(Pit != "N/A") %>%
  #filter(Pit != "0666F-DFF0") %>% # remove this individual - only record is a death
  #filter(Pit != "No Pit" | Pit != "NO PIT" | Pit == "see notes")

## Sort out Sex
levels(AllData$sex)

rm(list = c("Bitter", "FFRC", "Gaulin", "Leaf", "NAdder", "Noddy", "UCay", "WhiteB", "island_list", "common_cols"))
## Identify individuals with multiple "Sex" values
#inconsistent_ids <- AllData %>%
#  group_by(Pit) %>%
#  summarise(unique_sex = n_distinct(Sex), .groups = "drop") %>%
#  filter(unique_sex > 1) %>%
#  pull(Pit)

# View affected individuals
#UKSex <- AllData %>%
#  filter(Pit %in% inconsistent_ids) %>%
#  arrange(Pit)

#AllData <- AllData %>%
#  mutate(
#    Sex_edit = substr(Sex, 1, 1),  # Extract first letter
#    Sex_edit = case_when(
#      Sex_edit == "m" ~ "M",  # Normalize lowercase 'm' to 'M'
#      Sex_edit == "j" ~ "J",  # Normalize lowercase 'j' to 'J'
#      TRUE ~ Sex_edit  # Keep other values unchanged
#    )) %>%
#  mutate(Sex_edit = case_when(Sex_edit == "J" ~ "Unknown", TRUE ~ Sex_edit))
AllData <- AllData %>%
  mutate(
    sex = case_when(
      sex == "female" ~ "F",  # Normalize lowercase 'm' to 'M'
      sex == "male" ~ "M",  # Normalize lowercase 'j' to 'J'
      TRUE ~ sex  # Keep other values unchanged
    )) %>%
  group_by(ID) %>%
  # Keep only those individuals that have at least one "M" or "F"
  filter(any(sex %in% c("M", "F"))) %>%
  ungroup()

levels(as.factor(AllData$sex))

# Helper function to compute the mode (most common value) of a vector.
get_mode <- function(x) {
  # Remove NA values.
  x <- na.omit(x)
  if(length(x) == 0) return(NA)
  tab <- table(x)
  # In case of a tie, which.max returns the first value encountered.
  mode_val <- names(tab)[which.max(tab)]
  return(mode_val)
}

AllData %>%
  group_by(ID) %>%
  summarise(n_sexes = n_distinct(sex)) %>%
  filter(n_sexes > 1)
AllData <- arrange(AllData, ID)

# Process the data frame: for each individual, compute the mode of 'sex'
# and then assign that mode to all rows for that individual.
AllData <- AllData %>%
  group_by(ID) %>%
  mutate(common_sex = get_mode(sex),
         # If any entry is missing or inconsistent, set it to the most common value.
         sex = common_sex) %>%
  ungroup() %>%
  select(-common_sex)

AllData <- AllData %>%
  group_by(ID) %>%
  mutate(n_entries = n()) %>%
  ungroup() %>%
  filter(!(n_entries == 1 & is.na(sex))) %>%
  select(-n_entries)

which(is.na(AllData$sex))
## Adjust other random sex entries
#affected_ids <- AllData %>%
#  filter(Sex_edit %in% c(".", "d")) %>%
#  pull(Pit)

#AllData %>%
#  filter(Pit %in% affected_ids) %>%
#  arrange(Pit)

#AllData$Sex_edit[AllData$Pit == "4B21115615"] <- "M"
#AllData$Sex_edit[AllData$Pit == "4265337A2D"] <- "F"
#AllData$Sex_edit[AllData$Pit == "45335D6A7F"] <- "F"

## Add found dead entries
AllData <- AllData %>%
  mutate(Status = ifelse(
    apply(., 1, function(row) any(grepl("(?i)dead", row, perl = TRUE))), 2, 1)
  ) 

summary(as.factor(AllData$Status))

AllData <- arrange(AllData, ID)

## Check dead entries
alive_ids <- AllData %>%
  filter(Status != "2") %>%
  pull(ID)

# Keep only rows where the individual has at least one "alive" entry
#AllData <- AllData %>%
#  filter(ID %in% alive_ids) %>%
#  arrange(ID)

## Set date
#AllData$CapDate <- mdy(AllData$CapDate)
AllData$Year <- year(AllData$CapDate)
AllData <- AllData %>%
  group_by(Island) %>%
  mutate(IslandStart = min(Year))

# Check feeding information - Summary of feeding
feeding <- AllData %>%
  distinct(ID, .keep_all = TRUE) %>%
  count(Island, feeding) %>%
  arrange(Island, feeding)

# Check for individuals who have multiple feeding types
multi_feed <- AllData %>%
  group_by(ID) %>%
  summarise(n_feed_types = n_distinct(feeding)) %>%
  filter(n_feed_types > 1) %>%
  pull(ID)

# Remove them from the dataset
AllData <- AllData %>%
  filter(!(ID %in% multi_feed))

## Prepare data
YearZ <- 1950
Obs_Start <- min(AllData$Year)
first_year <- min(AllData$Year) - 1
last_year <- max(AllData$Year)
tMax <- last_year - first_year
nind <- length(levels(as.factor(AllData$ID)))

# Separate species
inorn <- read_csv("Data/Iguana_ID_InornataSp.csv") %>%
  mutate(
    Island = fct_recode(as_factor(Island), FFRC = "FlatRock"),  # rename level
    ID     = str_c(Island, animal_id, sep = "")                 # join columns
  )

summary(inorn)

sep_ids <- inorn %>% distinct(ID)  # just in case

AllData <- AllData %>%
  mutate(
    Species = case_when(
      # Force known inornata islands
      Island %in% c("FFRC", "UCay") ~ "inornata",
      # Use the reference list for others
      ID %in% sep_ids$ID ~ "inornata",
      # Default to figginsi
      TRUE ~ "figginsi"
    )
  )

#Check for individuals with more than one species label
species_conflicts <- AllData %>%
  group_by(ID) %>%
  summarise(n_species = n_distinct(Species), .groups = "drop") %>%
  filter(n_species > 1)

summary(as.factor(AllData$Species))

# Create a summary table of individuals per Island and Species
species_dist <- AllData %>%
  group_by(Island, Species) %>%
  summarise(n_individuals = n_distinct(ID), .groups = "drop") %>%
  # Pivot wider to see the comparison clearly
  pivot_wider(names_from = Species, values_from = n_individuals, values_fill = 0)

print(species_dist)

AllData$ID <- as.numeric(as.factor(AllData$ID))
AllData <- arrange(AllData, ID)
AllData$Year <- AllData$Year - YearZ
AllData$IslandStart <- AllData$IslandStart - YearZ

AllData <- AllData %>%
  mutate(Death = ifelse(Status == "2", Year, NA)) %>%  # Assign Year for deaths, NA otherwise
  group_by(ID) %>%
  mutate(Death = ifelse(all(is.na(Death)), NA, max(Death, na.rm = TRUE))) %>%
  ungroup()

AllData <- as.data.frame(AllData)
saveRDS(AllData, "AllData.rds")

## Add Birth date for known aged individuals
AllData <- AllData %>%
  mutate(age = suppressWarnings(as.numeric(age))) %>%
  group_by(ID) %>%
  mutate(ageR = round(age, digits = 0)) %>%  # Ensure age is rounded correctly
  mutate(
    Birth = ifelse(
      all(is.na(Year)) | all(is.na(ageR)), NA,  # If everything is NA, return NA
      min(Year, na.rm = TRUE) - min(ageR, na.rm = TRUE)
    )
  ) %>%
  mutate(LastSeen = case_when(any(Status == 1, na.rm = TRUE) ~ max(Year[Status == 1], na.rm = TRUE), TRUE ~ NA_real_)) %>%
  ungroup()

tKB <- AllData %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(Birth)

tKD <- AllData %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(Death)

## Create Capture History
CH <- matrix(NA, nind, last_year - YearZ)
AllData <- arrange(AllData, ID)

# Initialize with 0s
CH <- matrix(0, nrow = nind, ncol = last_year - YearZ)

# Filter for alive records
alive_data <- AllData %>% filter(Status == 1)

# Use a two-column matrix to index CH directly
# This is essentially: CH[row_index, col_index] <- 1
CH[as.matrix(alive_data[, c("ID", "Year")])] <- 1

## extract last alive time
tL <- apply(CH, 1, function(x) max(which(x == 1)))

## extract first alive time
tF <- apply(CH, 1, function(x) min(which(x == 1)))
names(tF) <- NULL

## define censoring times for birth time
# Initialize a blank matrix
cintB <- matrix(NA, nrow = nind, ncol = 2)

# Handle Known Births (Hatchlings and Aged Juveniles)
known_idx <- !is.na(tKB)
cintB[known_idx, 1] <- tKB[known_idx] - 1
cintB[known_idx, 2] <- tKB[known_idx]

# Handle Unknown Births (Adults)
# Use the first capture year minus 30 years
unknown_idx <- is.na(tKB)
cintB[unknown_idx, 1] <- tF[unknown_idx] - 29
cintB[unknown_idx, 2] <- tF[unknown_idx]

min(cintB) # check >0)
colnames(cintB) <- NULL

## define censoring matrices for death time
cintD <- cbind(tL, tKD)
cintD[is.na(tKD), 2] <- cintD[is.na(tKD), 1]
cintD[is.na(tKD), 1] <- 0
colnames(cintD) <- NULL
censoredD <- ifelse(!is.na(tKD), 1, 2)
tD <- rep(NA, length(tKD))
tstar <- rep(NA, length(tKD))

## some checks
stopifnot(all(tL >= tF))
stopifnot(all(tKD[!is.na(tKD)] > tL[!is.na(tKD)]))

## set up dummy variables
dind <- rep(1, length(tKD))

## extract number of captures
y <- apply(CH, 1, sum)
names(y) <- NULL

## set up nind
nind <- length(y)

# Clean up
rm(feeding)
rm(alive_ids)
rm(island_names)
rm(multi_feed)
rm(NACols)
rm(tag_index)
rm(tMax)
rm(inorn)
rm(sep_ids)
rm(alive_data)
rm(species_conflicts)
rm(species_dist)
rm(unknown_idx)
rm(known_idx)

save.image("igs_AllIslands_CleanDH_230226_obsStart.RData")
