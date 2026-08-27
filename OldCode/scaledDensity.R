# Make sure you have the tidyverse packages installed
# install.packages("tidyverse")
library(dplyr)
library(tidyr)

# --- 1. Create a Sample Dataset ---
# This mimics your raw data. Your real data frame should be named 'AllData'.
# It MUST have columns for the animal ID, the year (on your 0-64 scale), and the island.
AllData <- data.frame(
  animal_id = c('A', 'B', 'C', 'A', 'D', 'E', 'F', 'G', 'C', 'H'),
  Year = c(24, 24, 25, 25, 25, 26, 26, 26, 25, 27), # On your 0-64 scale
  Island = c('IslandA', 'IslandA', 'IslandA', 'IslandA', 'IslandB', 'IslandB', 'IslandB', 'IslandB', 'IslandB', 'IslandA')
)

# Your real data also has a complete range of years, let's say 0-64
all_years <- 0:64
all_islands <- unique(as.factor(AllData$Island))
levels(all_islands)

# --- 2. Calculate Unique Individuals Per Year and Island ---
# This creates a "long" table of the raw counts.
density_long <- AllData %>%
  group_by(Year, Island) %>%
  summarise(
    n_unique = n_distinct(ID),
    .groups = 'drop' # Good practice to ungroup after summarising
  )

# --- 3. Complete the Grid to Fill in Missing Zeros ---
# This is the CRITICAL step. It ensures that if an island was not surveyed in
# a given year, its observed density is correctly recorded as 0, not NA.
density_complete <- density_long %>%
  # 'expand' creates all combinations, then we join the counts back in
  tidyr::expand(Year = all_years, Island = all_islands) %>%
  left_join(density_long, by = c("Year", "Island")) %>%
  # Replace the NA's that were just created with 0
  mutate(n_unique = ifelse(is.na(n_unique), 0, n_unique))


# --- 4. Pivot from Long to Wide Matrix Format ---
# The result is a data frame with years as rows and islands as columns.
density_wide_df <- density_complete %>%
  pivot_wider(
    names_from = Island,
    values_from = n_unique
  ) %>%
  arrange(Year) # VERY IMPORTANT: Ensure rows are sorted by year.

# Let's inspect the result. Note the explicit zeros for years with no surveys.
 print(head(density_wide_df))
#   Year IslandA IslandB
#   <dbl>   <dbl>   <dbl>
# 1     0       0       0
# 2     1       0       0
# ...
# 24    24      2       0
# 25    25      2       2
# 26    26      0       3

# --- 5. Convert to a Matrix and Scale ---

# First, get the final column names in the correct order. You will need this
# later to create the 'island_index' for your model.
island_names_in_order <- colnames(density_wide_df)[-1] # Exclude the 'Year' column

# Convert the data frame to a numeric matrix (stripping the Year column)
density_raw_matrix <- as.matrix(density_wide_df[, -1])

# Finally, scale the entire matrix.
scaled_density_matrix <- scale(density_raw_matrix)

# This is the final object that goes into your NIMBLE 'consts' list.
# It should have 65 rows (for years 0-64) and n_islands columns.
cat("Dimensions of the final matrix:", dim(scaled_density_matrix), "\n")
cat("First few rows of the scaled matrix:\n")
print(head(scaled_density_matrix))
