## ============================================================
## MNA / COHORT-DENSITY AUDIT
##
## PURPOSE
##
## Check whether early MNA-based density is likely to contain
## a monitoring start-up / ascertainment artefact.
##
## Specifically:
##
## 1. When did monitoring begin on each island?
## 2. How many animals were known alive in each island-year?
## 3. How many entered the marked population for the first time?
## 4. What proportion of MNA consisted of newly encountered animals?
## 5. Does MNA systematically ramp up after monitoring begins?
## 6. Which density years are directly represented by survey years?
## 7. Which birth cohorts currently receive island-mean imputed density?
##
## IMPORTANT:
##
## This is an audit only.
## It does NOT modify density_matrix or any model object.
## ============================================================


library(tidyverse)


cat("\n")
cat("============================================================\n")
cat("MNA / COHORT-DENSITY AUDIT\n")
cat("============================================================\n")


## ============================================================
## 1. HABITAT AREA
## ============================================================

habitat_audit <- read_csv(
  "Data/Iguana_islands_habitat.csv",
  show_col_types = FALSE
) %>%
  transmute(
    island = islands,
    usable_m2 = Useable,
    usable_ha = usable_m2 / 10000
  )


cat("\nIsland habitat areas:\n")

print(
  habitat_audit,
  n = Inf
)


## ============================================================
## 2. INDIVIDUAL FIRST / LAST KNOWN YEARS
##
## These are calendar years.
## ============================================================

life_years <- growth_info %>%
  transmute(
    ID,
    ID_raw,
    island,
    
    first_year =
      tF_cal,
    
    last_year =
      tL_cal
  )


stopifnot(
  !anyNA(life_years$first_year),
  !anyNA(life_years$last_year),
  
  all(
    life_years$last_year >=
      life_years$first_year
  )
)


## ============================================================
## 3. ACTUAL SURVEY YEARS
##
## Keep number of survey records as well as survey presence.
## ============================================================

survey_years_audit <- actual_surveys %>%
  transmute(
    island = Island,
    
    year_idx =
      yr_adj,
    
    calendar_year =
      yr_adj +
      Origin -
      1
  ) %>%
  count(
    island,
    year_idx,
    calendar_year,
    name = "n_survey_records"
  )


cat("\nActual survey-year coverage:\n")

survey_years_audit %>%
  group_by(island) %>%
  summarise(
    first_survey =
      min(calendar_year),
    
    last_survey =
      max(calendar_year),
    
    n_survey_years =
      n_distinct(calendar_year),
    
    n_survey_records =
      sum(n_survey_records),
    
    .groups = "drop"
  ) %>%
  arrange(first_survey) %>%
  print(n = Inf)


## ============================================================
## 4. RECONSTRUCT MNA
##
## An individual contributes to MNA from first known year
## through last known year inclusive.
##
## This is the same conceptual construction currently used
## for the cohort-density covariate.
## ============================================================

mna_calendar <- life_years %>%
  rowwise() %>%
  mutate(
    calendar_year = list(
      seq(
        first_year,
        last_year
      )
    )
  ) %>%
  unnest(calendar_year) %>%
  ungroup() %>%
  count(
    island,
    calendar_year,
    name = "MNA"
  )


## ============================================================
## 5. NEW ENTRANTS INTO THE MARKED / KNOWN POPULATION
##
## These are individuals whose FIRST recorded year is this year.
##
## A high value early in monitoring can indicate accumulation
## of previously unmarked animals rather than true population
## growth.
## ============================================================

new_entries <- life_years %>%
  count(
    island,
    first_year,
    name = "n_new"
  ) %>%
  rename(
    calendar_year =
      first_year
  )


## ============================================================
## 6. LAST OBSERVATION COUNTS
##
## Useful descriptively when looking at turnover.
## ============================================================

last_entries <- life_years %>%
  count(
    island,
    last_year,
    name = "n_last"
  ) %>%
  rename(
    calendar_year =
      last_year
  )


## ============================================================
## 7. COMPLETE ISLAND-YEAR AUDIT TABLE
##
## Start from observed survey years.
## ============================================================

mna_survey_audit <- survey_years_audit %>%
  left_join(
    mna_calendar,
    by = c(
      "island",
      "calendar_year"
    )
  ) %>%
  left_join(
    new_entries,
    by = c(
      "island",
      "calendar_year"
    )
  ) %>%
  left_join(
    last_entries,
    by = c(
      "island",
      "calendar_year"
    )
  ) %>%
  left_join(
    habitat_audit,
    by = "island"
  ) %>%
  mutate(
    MNA =
      replace_na(
        MNA,
        0L
      ),
    
    n_new =
      replace_na(
        n_new,
        0L
      ),
    
    n_last =
      replace_na(
        n_last,
        0L
      ),
    
    density_raw =
      MNA /
      usable_ha
  ) %>%
  group_by(island) %>%
  arrange(
    calendar_year,
    .by_group = TRUE
  ) %>%
  mutate(
    survey_number =
      row_number(),
    
    years_since_first_survey =
      calendar_year -
      min(calendar_year),
    
    cumulative_new =
      cumsum(n_new),
    
    n_previously_known =
      pmax(
        0,
        MNA -
          n_new
      ),
    
    entrant_fraction =
      if_else(
        MNA > 0,
        n_new / MNA,
        NA_real_
      ),
    
    annual_MNA_change =
      MNA -
      lag(MNA),
    
    annual_MNA_ratio =
      if_else(
        lag(MNA) > 0,
        MNA / lag(MNA),
        NA_real_
      )
  ) %>%
  ungroup()


## ============================================================
## 8. BASIC MONITORING SUMMARY BY ISLAND
## ============================================================

monitoring_summary <- life_years %>%
  group_by(island) %>%
  summarise(
    first_capture_year =
      min(first_year),
    
    last_capture_year =
      max(last_year),
    
    n_individuals =
      n(),
    
    .groups = "drop"
  ) %>%
  left_join(
    survey_years_audit %>%
      group_by(island) %>%
      summarise(
        first_survey_year =
          min(calendar_year),
        
        last_survey_year =
          max(calendar_year),
        
        n_survey_years =
          n_distinct(calendar_year),
        
        .groups = "drop"
      ),
    by = "island"
  ) %>%
  mutate(
    capture_minus_survey_start =
      first_capture_year -
      first_survey_year
  )


cat("\n")
cat("============================================================\n")
cat("MONITORING START SUMMARY\n")
cat("============================================================\n")

print(
  monitoring_summary,
  n = Inf
)


## ============================================================
## 9. FIRST TEN SURVEY YEARS FOR EACH ISLAND
##
## This is probably the most useful table to inspect first.
## ============================================================

early_year_audit <- mna_survey_audit %>%
  filter(
    survey_number <= 10
  ) %>%
  select(
    island,
    survey_number,
    calendar_year,
    n_survey_records,
    MNA,
    n_new,
    n_previously_known,
    entrant_fraction,
    cumulative_new,
    density_raw,
    annual_MNA_change,
    annual_MNA_ratio
  )


cat("\n")
cat("============================================================\n")
cat("FIRST 10 SURVEY YEARS BY ISLAND\n")
cat("============================================================\n")

print(
  early_year_audit,
  n = Inf
)


## ============================================================
## 10. START-UP SUMMARY
##
## Compare first 3, first 5 and later survey years.
##
## These cut-offs are descriptive only.
## We are NOT yet declaring the first 3/5 years invalid.
## ============================================================

startup_summary <- mna_survey_audit %>%
  mutate(
    monitoring_phase = case_when(
      survey_number <= 3 ~
        "Survey years 1-3",
      
      survey_number <= 5 ~
        "Survey years 4-5",
      
      TRUE ~
        "Survey year 6+"
    )
  ) %>%
  group_by(
    island,
    monitoring_phase
  ) %>%
  summarise(
    n_years =
      n(),
    
    mean_MNA =
      mean(
        MNA,
        na.rm = TRUE
      ),
    
    median_MNA =
      median(
        MNA,
        na.rm = TRUE
      ),
    
    mean_density =
      mean(
        density_raw,
        na.rm = TRUE
      ),
    
    median_density =
      median(
        density_raw,
        na.rm = TRUE
      ),
    
    mean_new =
      mean(
        n_new,
        na.rm = TRUE
      ),
    
    mean_entrant_fraction =
      mean(
        entrant_fraction,
        na.rm = TRUE
      ),
    
    median_entrant_fraction =
      median(
        entrant_fraction,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("START-UP PHASE SUMMARY\n")
cat("============================================================\n")

print(
  startup_summary,
  n = Inf
)


## ============================================================
## 11. EARLY VERSUS ESTABLISHED MONITORING
##
## A simple island-level summary.
##
## "Established" is survey year 6 onward ONLY for diagnosis.
## This is not yet a modelling decision.
## ============================================================

early_late_summary <- mna_survey_audit %>%
  group_by(island) %>%
  summarise(
    MNA_first_year =
      first(MNA),
    
    MNA_year5 =
      if(
        n() >= 5
      ) {
        MNA[5]
      } else {
        NA_real_
      },
    
    median_MNA_first3 =
      median(
        MNA[
          survey_number <= 3
        ],
        na.rm = TRUE
      ),
    
    median_MNA_year6plus =
      if(
        any(
          survey_number >= 6
        )
      ) {
        median(
          MNA[
            survey_number >= 6
          ],
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    
    mean_entry_fraction_first3 =
      mean(
        entrant_fraction[
          survey_number <= 3
        ],
        na.rm = TRUE
      ),
    
    mean_entry_fraction_year6plus =
      if(
        any(
          survey_number >= 6
        )
      ) {
        mean(
          entrant_fraction[
            survey_number >= 6
          ],
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    
    first3_to_later_MNA_ratio =
      median_MNA_first3 /
      median_MNA_year6plus,
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("EARLY VERSUS ESTABLISHED MONITORING\n")
cat("============================================================\n")

print(
  early_late_summary,
  n = Inf
)


## ============================================================
## 12. ZERO / VERY LOW MNA SURVEY YEARS
## ============================================================

low_mna_years <- mna_survey_audit %>%
  filter(
    MNA == 0 |
      MNA <= 5
  ) %>%
  select(
    island,
    survey_number,
    calendar_year,
    MNA,
    n_new,
    entrant_fraction,
    density_raw
  )


cat("\n")
cat("============================================================\n")
cat("SURVEY YEARS WITH MNA <= 5\n")
cat("============================================================\n")

print(
  low_mna_years,
  n = Inf
)


## ============================================================
## 13. CHECK FOR AN ASCERTAINMENT RAMP
##
## Within each island, calculate correlation between MNA and
## years since monitoring began during the first 10 survey years.
##
## IMPORTANT:
## A positive correlation does not prove ascertainment bias;
## real populations can increase too.
##
## It is simply a useful flag.
## ============================================================

startup_correlations <- mna_survey_audit %>%
  filter(
    survey_number <= 10
  ) %>%
  group_by(island) %>%
  summarise(
    n_years =
      n(),
    
    cor_MNA_time =
      if(
        n() >= 3 &&
        sd(MNA) > 0
      ) {
        cor(
          years_since_first_survey,
          MNA
        )
      } else {
        NA_real_
      },
    
    cor_density_time =
      if(
        n() >= 3 &&
        sd(density_raw) > 0
      ) {
        cor(
          years_since_first_survey,
          density_raw
        )
      } else {
        NA_real_
      },
    
    cor_entry_fraction_time =
      if(
        sum(
          is.finite(
            entrant_fraction
          )
        ) >= 3
      ) {
        cor(
          years_since_first_survey[
            is.finite(
              entrant_fraction
            )
          ],
          
          entrant_fraction[
            is.finite(
              entrant_fraction
            )
          ]
        )
      } else {
        NA_real_
      },
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("FIRST-10-YEAR START-UP CORRELATIONS\n")
cat("============================================================\n")

print(
  startup_correlations,
  n = Inf
)


## ============================================================
## 14. REBUILD CURRENT DENSITY SOURCE TABLE
##
## This reproduces the current logic:
##
## observed survey year -> MNA density
## otherwise            -> island mean density
## ============================================================

density_observed_audit <- survey_years_audit %>%
  select(
    island,
    year_idx,
    calendar_year
  ) %>%
  distinct() %>%
  left_join(
    mna_calendar %>%
      mutate(
        year_idx =
          calendar_year -
          Origin +
          1
      ) %>%
      select(
        island,
        year_idx,
        MNA
      ),
    by = c(
      "island",
      "year_idx"
    )
  ) %>%
  left_join(
    habitat_audit,
    by = "island"
  ) %>%
  mutate(
    MNA =
      replace_na(
        MNA,
        0L
      ),
    
    density_raw =
      MNA /
      usable_ha
  )


island_mean_density_audit <- density_observed_audit %>%
  group_by(island) %>%
  summarise(
    mean_density_raw =
      mean(
        density_raw,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


density_source_audit <- expand_grid(
  island = ordered_islands,
  year_idx = seq_len(tMax)
) %>%
  mutate(
    calendar_year =
      year_idx +
      Origin -
      1
  ) %>%
  left_join(
    density_observed_audit %>%
      select(
        island,
        year_idx,
        density_raw
      ),
    by = c(
      "island",
      "year_idx"
    )
  ) %>%
  left_join(
    island_mean_density_audit,
    by = "island"
  ) %>%
  mutate(
    density_source =
      if_else(
        is.na(density_raw),
        "Island mean imputation",
        "Observed survey-year MNA"
      ),
    
    density_used =
      coalesce(
        density_raw,
        mean_density_raw
      )
  )


## ============================================================
## 15. OBSERVED VERSUS IMPUTED DENSITY YEARS
## ============================================================

density_source_summary <- density_source_audit %>%
  group_by(
    island,
    density_source
  ) %>%
  summarise(
    n_years =
      n(),
    
    first_year =
      min(calendar_year),
    
    last_year =
      max(calendar_year),
    
    .groups = "drop"
  )


cat("\n")
cat("============================================================\n")
cat("DENSITY SOURCE BY ISLAND\n")
cat("============================================================\n")

print(
  density_source_summary,
  n = Inf
)


## ============================================================
## 16. BIRTH-MIDPOINT DENSITY SOURCE
##
## Descriptive only.
##
## The actual model keeps birth time latent, so an individual's
## cohort year can move during MCMC.
## ============================================================

birth_mid_year_idx <- floor(
  rowMeans(
    cintB_adj
  )
)

birth_mid_year_idx <- pmax(
  1,
  pmin(
    tMax,
    birth_mid_year_idx
  )
)


birth_mid_density_source <- tibble(
  ID = growth_info$ID_raw,
  island = growth_info$island,
  
  birth_mid_idx =
    birth_mid_year_idx,
  
  birth_mid_calendar =
    birth_mid_year_idx +
    Origin -
    1
) %>%
  left_join(
    density_source_audit %>%
      select(
        island,
        year_idx,
        density_source,
        density_used
      ),
    by = c(
      "island",
      "birth_mid_idx" = "year_idx"
    )
  )


cat("\n")
cat("============================================================\n")
cat("DENSITY SOURCE AT BIRTH-INTERVAL MIDPOINT\n")
cat("============================================================\n")

birth_mid_density_source %>%
  count(
    island,
    density_source
  ) %>%
  group_by(island) %>%
  mutate(
    proportion =
      n /
      sum(n)
  ) %>%
  ungroup() %>%
  print(n = Inf)


## ============================================================
## 17. FULL BIRTH-WINDOW EXPOSURE TO OBSERVED / IMPUTED DENSITY
##
## More informative than midpoint alone.
##
## For every animal, calculate what proportion of candidate
## integer cohort years in its birth interval use observed MNA
## versus island-mean imputation.
## ============================================================

birth_density_window_audit <- vector(
  "list",
  nind
)


for(i in seq_len(nind)) {
  
  possible_years <- seq(
    floor(
      cintB_adj[i, 1]
    ),
    floor(
      cintB_adj[i, 2]
    )
  )
  
  possible_years <- unique(
    pmax(
      1,
      pmin(
        tMax,
        possible_years
      )
    )
  )
  
  
  src <- density_source_audit %>%
    filter(
      island ==
        growth_info$island[i],
      
      year_idx %in%
        possible_years
    )
  
  
  birth_density_window_audit[[i]] <- tibble(
    ID =
      growth_info$ID_raw[i],
    
    island =
      growth_info$island[i],
    
    n_possible_years =
      length(
        possible_years
      ),
    
    n_observed_density =
      sum(
        src$density_source ==
          "Observed survey-year MNA"
      ),
    
    n_imputed_density =
      sum(
        src$density_source ==
          "Island mean imputation"
      ),
    
    prop_observed_density =
      n_observed_density /
      n_possible_years,
    
    prop_imputed_density =
      n_imputed_density /
      n_possible_years
  )
}


birth_density_window_audit <- bind_rows(
  birth_density_window_audit
)


cat("\n")
cat("============================================================\n")
cat("BIRTH-WINDOW DENSITY-SOURCE SUMMARY\n")
cat("============================================================\n")


print(
  summary(
    birth_density_window_audit$
      prop_observed_density
  )
)


cat("\nIndividuals with NO directly observed density year ")
cat("within their possible birth interval:\n")


birth_density_window_audit %>%
  filter(
    prop_observed_density == 0
  ) %>%
  count(island) %>%
  print(n = Inf)


cat("\nIndividuals whose entire possible birth interval ")
cat("has observed annual density:\n")


birth_density_window_audit %>%
  filter(
    prop_observed_density == 1
  ) %>%
  count(island) %>%
  print(n = Inf)


## ============================================================
## 18. VISUAL AUDIT 1:
## MNA THROUGH TIME
## ============================================================

p_mna <- ggplot(
  mna_survey_audit,
  aes(
    x = calendar_year,
    y = MNA
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(
    ~ island,
    scales = "free_y"
  ) +
  labs(
    x = "Calendar year",
    y = "Minimum number alive",
    title = "Minimum number alive through monitoring",
    subtitle = "Survey years only"
  ) +
  theme_bw()


print(
  p_mna
)


## ============================================================
## 19. VISUAL AUDIT 2:
## NEW ENTRANTS AS PROPORTION OF MNA
## ============================================================

p_entry <- ggplot(
  mna_survey_audit,
  aes(
    x = calendar_year,
    y = entrant_fraction
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(
    ~ island
  ) +
  labs(
    x = "Calendar year",
    y = "Newly encountered / MNA",
    title = "Accumulation of animals into the known population",
    subtitle = "High early values may indicate monitoring start-up effects"
  ) +
  theme_bw()


print(
  p_entry
)


## ============================================================
## 20. VISUAL AUDIT 3:
## MNA VERSUS SURVEY NUMBER
##
## Aligns islands by monitoring start rather than calendar year.
## ============================================================

p_startup <- ggplot(
  mna_survey_audit,
  aes(
    x = survey_number,
    y = MNA
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(
    ~ island,
    scales = "free_y"
  ) +
  labs(
    x = "Survey year number since monitoring began",
    y = "Minimum number alive",
    title = "MNA relative to monitoring start"
  ) +
  theme_bw()


print(
  p_startup
)


## ============================================================
## 21. VISUAL AUDIT 4:
## RAW MNA-BASED DENSITY
## ============================================================

p_density <- ggplot(
  mna_survey_audit,
  aes(
    x = calendar_year,
    y = density_raw
  )
) +
  geom_line() +
  geom_point() +
  facet_wrap(
    ~ island,
    scales = "free_y"
  ) +
  labs(
    x = "Calendar year",
    y = "MNA per usable hectare",
    title = "Raw MNA-based density through time"
  ) +
  theme_bw()


print(
  p_density
)


## ============================================================
## 22. SAVE AUDIT TABLES
## ============================================================

write_csv(
  mna_survey_audit,
  "Audit_MNA_island_year.csv"
)


write_csv(
  monitoring_summary,
  "Audit_MNA_monitoring_summary.csv"
)


write_csv(
  early_year_audit,
  "Audit_MNA_first10_years.csv"
)


write_csv(
  startup_summary,
  "Audit_MNA_startup_phases.csv"
)


write_csv(
  early_late_summary,
  "Audit_MNA_early_vs_established.csv"
)


write_csv(
  density_source_summary,
  "Audit_Density_source_by_island.csv"
)


write_csv(
  birth_density_window_audit,
  "Audit_Birth_window_density_source.csv"
)


## ============================================================
## 23. SAVE FIGURES
## ============================================================

ggsave(
  "Audit_MNA_through_time.png",
  p_mna,
  width = 12,
  height = 8,
  dpi = 300
)


ggsave(
  "Audit_MNA_new_entrant_fraction.png",
  p_entry,
  width = 12,
  height = 8,
  dpi = 300
)


ggsave(
  "Audit_MNA_since_monitoring_start.png",
  p_startup,
  width = 12,
  height = 8,
  dpi = 300
)


ggsave(
  "Audit_MNA_density_through_time.png",
  p_density,
  width = 12,
  height = 8,
  dpi = 300
)


cat("\n")
cat("============================================================\n")
cat("MNA / DENSITY AUDIT COMPLETE\n")
cat("============================================================\n")