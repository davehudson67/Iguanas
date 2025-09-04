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

source("ModelComparison_FUNCTIONS.R")

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
FFRC <- read_csv("DoubleCensored/All_Island_Data/Data/FlatRockThru2019_CleanOct2020.csv")%>%
  mutate(date = mdy(date))  %>%
  mutate(feeding = as.factor(feeding))
summary(FFRC$feeding)
Gaulin <- read_csv("DoubleCensored/All_Island_Data/Data/GaulinCay_clean.csv") %>%
  mutate(date = dmy(date)) %>%
  filter(!(Pit == "4264176E36" & `body-mass` == 2070)) %>%
  filter(!(Pit == "442A0F6D30" & `body-mass` == 2460))  %>%
  mutate(feeding = as.factor(feeding))
summary(Gaulin$feeding)
#Leaf_A <- read_csv("DoubleCensored/All_Island_Data/Cfigginsi_Leaf_Cay_Allens.csv")
Leaf <- read_csv("DoubleCensored/All_Island_Data/Data/LeafCayThru2019_CleanOct2020.csv") %>%
  mutate(date = mdy(date))  %>%
  mutate(feeding = as.factor(feeding))
summary(Leaf$feeding)
Noddy <- read_csv("DoubleCensored/All_Island_Data/Data/NoddyCCay_Clean.csv")%>%
  mutate(date = dmy(date)) %>%
  mutate(feeding = as.factor(feeding))
summary(Noddy$feeding)
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
AllData <- AllData %>%
  filter(ID %in% alive_ids) %>%
  arrange(ID)

## Set date
#AllData$CapDate <- mdy(AllData$CapDate)
AllData$Year <- year(AllData$CapDate)

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
Obs_Start <- min(AllData$Year) - 1
first_year <- min(AllData$Year) - 1
last_year <- max(AllData$Year)
tMax <- last_year - first_year
nind <- length(levels(as.factor(AllData$ID)))
AllData$ID <- as.numeric(as.factor(AllData$ID))
AllData <- arrange(AllData, ID)
AllData$Year <- AllData$Year - first_year

AllData <- AllData %>%
  mutate(Death = ifelse(Status == "2", Year, NA)) %>%  # Assign Year for deaths, NA otherwise
  group_by(ID) %>%
  mutate(Death = ifelse(all(is.na(Death)), NA, max(Death, na.rm = TRUE))) %>%
  ungroup()

AllData <- as.data.frame(AllData)
saveRDS(AllData, "AllData.rds")

## Add Birth date for known aged individuals
#AllData <- AllData %>%
#  mutate(age = as.numeric(age)) %>%
#  mutate(Year = Year - first_year) %>%
#  group_by(ID) %>%
#  mutate(ageR = round(age)) %>%
#  mutate(Birth = min(Year) - min(ageR)) %>%
#  mutate(LastSeen = max(Year)) %>%
#  #distinct(tag, Year, .keep_all = TRUE) %>%
#  #mutate(captures = n()) %>%
#  #distinct(animal_id, .keep_all = TRUE) %>%
#  #mutate(death = NA) %>%
#  ungroup() #%>%
##mutate(Max_captures = max(Year)) %>%
##mutate(n = 1)

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

###################################################
#AllData <- arrange(AllData, ID, age)

#igsAll <- igsAll %>%
#  mutate(Year = Year - first_year) %>%
#  group_by(ID) %>%
#  mutate(LastSeen = max(Year)) %>%
#  mutate(Birth = Birth - first_year) %>%
#  #distinct(animal_id, Year, .keep_all = TRUE) %>%
  #mutate(captures = n()) %>%
  #distinct(animal_id, .keep_all = TRUE) %>%
  #mutate(death = NA) %>%
  #mutate(Max_captures = max(Year)) %>%
#  droplevels() %>%
#  ungroup()
###################################################

## Create Capture History
CH <- matrix(NA, nind, last_year - first_year)
#AllData$Year <- AllData$Year - first_year
AllData <- arrange(AllData, ID)

## Fill CH matrix using a loop (only when Status == 1)
for (i in 1:nrow(AllData)) {
  if (AllData$Status[i] == 1) {  # Only process if the individual is alive
    tag_index <- AllData$ID[i]
    year_index <- AllData$Year[i]
    
    # Fill the CH matrix
    CH[tag_index, year_index] <- 1
  }
}

## Add births
#filtered <- filter(AllData, !is.na(AllData$Birth), .preserve = TRUE) %>%
#  distinct(ID, .keep_all = TRUE)
#for(i in 1:nrow(filtered)){
#  tag <- AllData$ID[i]
#  birth <- AllData$Birth[i]
#  CH[tag, birth] <- 1
#}
#rm(filtered)

tKB <- AllData %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(Birth)

tKD <- AllData %>%
  distinct(ID, .keep_all = TRUE) %>%
  pull(Death)

######################################################
#tKB <- igsAll %>%
#  distinct(animal_id, .keep_all = TRUE) %>%
#  pull(Birth)
######################################################

CH[is.na(CH)] <- 0

## extract last alive time
tL <- apply(CH, 1, function(x) max(which(x == 1)))

## extract first alive time
tF <- apply(CH, 1, function(x) min(which(x == 1)))
names(tF) <- NULL

## define censoring times for birth time
#cintB <- cbind(tKB - 1, tKB)
#cintB[is.na(tKB), 2] <- tF[is.na(tKB)] - 25
#cintB[cintB < 0] <- 1
#cintB[is.na(tKB), 1] <- 0
#colnames(cintB) <- NULL
#censoredB <- rep(1, nrow(cintB))
#tB <- rep(NA, length(tKB))

cintB <- cbind(tF - 1, tF)
cintB[is.na(tKB), 1] <- tF[is.na(tKB)] - 25
#cintB[cintB < 0] <- 1
#cintB[is.na(tKB), 2] <- 0
colnames(cintB) <- NULL
censoredB <- rep(1, nrow(cintB))
tB <- rep(NA, length(tKB))

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
save.image("igs_AllIslands_CleanDH_290825_obsStart.RData")
