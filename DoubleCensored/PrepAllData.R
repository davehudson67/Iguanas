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
rm(list=ls())

source("ModelComparison_FUNCTIONS.R")

## Load data
leafClay <- read_csv("Data/LeafCayThru2019_CleanOct2020.csv")
flatRock <- read_csv("Data/FlatRockThru2019_CleanOct2020.csv")
Ucay <- read_csv("Data/UcayThru2019_CleanOct2020.csv")

leafClay <- select(leafClay, animal_id, sex, age, date, svl, feeding) %>%
  mutate(animal_id = paste0("LC_", animal_id),
         island = "leafClay")

flatRock <- select(flatRock, animal_id, sex, age, date, svl, feeding) %>%
  mutate(animal_id = paste0("FR_", animal_id),
         island = "flatRock")

Ucay <- select(Ucay, animal_id, sex, age, date, svl, feeding) %>%
  mutate(animal_id = paste0("UC_", animal_id),
         island = "Ucay")

igs <- rbind(leafClay, flatRock, Ucay)

## Set date
igs$date <- mdy(igs$date)
igs$Year <- year(igs$date)
igs$age[igs$age == "2.7?"] <- 2.7

## sort data
igsAll <- igs %>%
  filter(!is.na(age), .preserve = TRUE) %>%
  mutate(age = as.numeric(age)) %>%
  #mutate(Year = Year - first_year) %>%
  mutate(ageR = round(age)) %>%
  group_by(animal_id) %>%
  mutate(Birth = min(Year) - min(ageR)) %>%
  #mutate(LastSeen = max(Year)) %>%
  #distinct(animal_id, Year, .keep_all = TRUE) %>%
  #mutate(captures = n()) %>%
  #distinct(animal_id, .keep_all = TRUE) %>%
  #mutate(death = NA) %>%
  #ungroup() %>%
  #mutate(Max_captures = max(Year)) %>%
  mutate(n = 1) %>%
  droplevels() %>%
  ungroup()

first_year <- min(igsAll$Birth)
last_year <- max(igs$Year)
tMax <- last_year - first_year

igsAll <- igsAll %>%
  mutate(Year = Year - first_year) %>%
  group_by(animal_id) %>%
  mutate(LastSeen = max(Year)) %>%
  mutate(Birth = Birth - first_year) %>%
  #distinct(animal_id, Year, .keep_all = TRUE) %>%
  #mutate(captures = n()) %>%
  #distinct(animal_id, .keep_all = TRUE) %>%
  #mutate(death = NA) %>%
  #mutate(Max_captures = max(Year)) %>%
  droplevels() %>%
  ungroup()

nind <- length(levels(as.factor(igsAll$animal_id)))
igsAll$tag <- as.numeric(as.factor(igsAll$animal_id))
igsAll <- arrange(igsAll, tag)

CH <- matrix(NA, nind, last_year - first_year)

## Fill CH matrix using a loop
for (i in 1:nrow(igsAll)) {
  tag_index <- igsAll$tag[i]
  year_index <- igsAll$Year[i]
  
  # Fill the CH matrix
  CH[tag_index, year_index] <- igsAll$n[i]
}

## Add births
#filtered <- filter(igsAll, !is.na(igsAll$age), .preserve = TRUE)
#for(i in 1:nrow(filtered)){
#  tag <- igsAll$tag[i]
#  birth <- igsAll$Birth[i]
#  CH[tag, birth] <- 1
#}
#rm(filtered)

tKD <- rep(NA, nrow(CH))
tKB <- igsAll %>%
  distinct(animal_id, .keep_all = TRUE) %>%
  pull(Birth)

CH[is.na(CH)] <- 0

## extract last alive time
tL <- apply(CH, 1, function(x) max(which(x == 1)))
names(tL) <- NULL

## extract first alive time
tF <- apply(CH, 1, function(x) min(which(x == 1)))
names(tF) <- NULL

## define censoring times for birth time
cintB <- cbind(tKB - 1, tKB)
cintB[is.na(tKB), 1] <- tF[is.na(tKB)] - 25
cintB[cintB < 0] <- 1
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

save.image("igs_ready_AllIslands_KA.RData")
