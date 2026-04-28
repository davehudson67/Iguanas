library(tidyverse)
library(lubridate)

igs <- read_csv("Data/LeafCayThru2019_CleanOct2020.csv")

igs$n <- 1
igs$date <- mdy(igs$date)
igs$Year <- year(igs$date)

create_capture_history <- function(data, first_year, last_year) {
  # Check for multiple records in the same year for the same animal
  #browser()
  duplicates <- data %>%
    group_by(animal_id, Year) %>%
    filter(n() > 1) %>%
    distinct(animal_id)
  
  # Warning for animals with multiple entries
  if (nrow(duplicates) > 0) {
    warning("The following individuals have multiple entries for a given year: ",
            paste(duplicates[[2]], duplicates[[1]], collapse = ", "))
  }
  
  # Remove duplicates, keeping only the first record for each animal-year pair
  data <- data %>%
    distinct(animal_id, Year, .keep_all = TRUE)
  
  # Create the capture history
  result <- data %>%
    group_by(animal_id) %>%
    complete(Year = seq(first_year, last_year, by = 1)) %>%
    arrange(animal_id, Year) %>%
    mutate(capture_history = ifelse(!is.na(n), 1, 0)) %>%
    mutate(capture_history = ifelse(cumsum(capture_history) == 0, NA, capture_history)) %>%
    select(animal_id, Year, capture_history) %>%
    pivot_wider(names_from = Year, values_from = capture_history, names_prefix = "Y_")
  
  return(list(CH = result, duplicates = duplicates))
}

out <- create_capture_history(igs, 1970, 2019)

# Access CH and duplicates
CH <- out$CH
duplicates <- out$duplicates

igs_firstEncounters <- distinct(igs, animal_id, .keep_all = TRUE)

length(which(is.na(igs_firstEncounters$age)))




igs <- igs %>%
  group_by(animal_id) %>%
  summarise(n = n(), .groups = 'drop')
max(igs$n)
