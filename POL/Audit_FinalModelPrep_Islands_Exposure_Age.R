library(tidyverse)
library(mgcv)

## ============================================================
## FINAL-MODEL PREPARATION AUDIT
## ============================================================

source("DataPrep.R")

out_dir <- "DataAudits/FinalModelPrep"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 1. Exact island labels ------------------------------------------------

island_audit <- AllData %>%
  mutate(Island = as.character(Island)) %>%
  group_by(Island) %>%
  summarise(
    n_records = n(),
    n_individuals = n_distinct(ID),
    first_year = min(CalYear, na.rm = TRUE),
    last_year = max(CalYear, na.rm = TRUE),
    feeding_labels = paste(sort(unique(as.character(feeding))), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Island)

cat("\n============================================================\n")
cat("EXACT ISLAND LABELS\n")
cat("============================================================\n")
print(island_audit, n = Inf)

write_csv(
  island_audit,
  file.path(out_dir, "Audit_IslandLabels.csv")
)

## 2. Leaf North / South check ------------------------------------------

leaf_records <- AllData %>%
  mutate(Island = as.character(Island)) %>%
  filter(str_detect(str_to_lower(Island), "leaf"))

leaf_summary <- leaf_records %>%
  group_by(Island) %>%
  summarise(
    n_records = n(),
    n_individuals = n_distinct(ID),
    first_year = min(CalYear, na.rm = TRUE),
    last_year = max(CalYear, na.rm = TRUE),
    .groups = "drop"
  )

leaf_ids <- leaf_records %>%
  distinct(Island, ID, ID_raw, animal_id) %>%
  arrange(Island, ID)

cat("\n============================================================\n")
cat("LEAF CAY LABEL CHECK\n")
cat("============================================================\n")
print(leaf_summary, n = Inf)
cat("Distinct Leaf labels:", n_distinct(leaf_records$Island), "\n")

if(n_distinct(leaf_records$Island) < 2) {
  cat("FLAG: Leaf North and Leaf South are not obviously separated in Island.\n")
}

write_csv(
  leaf_summary,
  file.path(out_dir, "Audit_LeafLabels.csv")
)

write_csv(
  leaf_ids,
  file.path(out_dir, "Audit_LeafIDs.csv")
)

## 3. Known-age support -------------------------------------------------

known_people <- growth_info %>%
  as_tibble() %>%
  mutate(
    island = as.character(island),
    has_known_age = !is.na(known_age) & !is.na(known_age_year)
  )

known_by_island <- known_people %>%
  group_by(island) %>%
  summarise(
    n_individuals = n(),
    n_known_age = sum(has_known_age),
    proportion_known = mean(has_known_age),
    .groups = "drop"
  ) %>%
  arrange(desc(n_known_age))

cat("\n============================================================\n")
cat("KNOWN-AGE SUPPORT BY ISLAND\n")
cat("============================================================\n")
print(known_by_island, n = Inf)

write_csv(
  known_by_island,
  file.path(out_dir, "Audit_KnownAge_ByIsland.csv")
)

## 4. Known-age capture-level SVL audit --------------------------------

known_ref <- known_people %>%
  filter(has_known_age) %>%
  select(ID, island, known_age, known_age_year)

svl_known <- AllData %>%
  inner_join(known_ref, by = "ID") %>%
  mutate(
    svl_num = suppressWarnings(as.numeric(svl)),
    age_at_capture = known_age + CalYear - known_age_year,
    island_age = as.character(island),
    ID_f = factor(ID),
    sex_f = factor(sex),
    calendar_decade = (CalYear - 2000) / 10
  ) %>%
  filter(
    is.finite(svl_num),
    is.finite(age_at_capture),
    age_at_capture >= 0
  )

svl_support <- svl_known %>%
  group_by(island_age) %>%
  summarise(
    n_individuals = n_distinct(ID),
    n_records = n(),
    first_year = min(CalYear),
    last_year = max(CalYear),
    min_age = min(age_at_capture),
    max_age = max(age_at_capture),
    .groups = "drop"
  ) %>%
  arrange(desc(n_individuals))

cat("\n============================================================\n")
cat("KNOWN-AGE CAPTURE SUPPORT\n")
cat("============================================================\n")
print(svl_support, n = Inf)

write_csv(
  svl_support,
  file.path(out_dir, "Audit_KnownAge_CaptureSupport.csv")
)

supported <- svl_support %>%
  filter(n_individuals >= 20) %>%
  pull(island_age)

svl_fit_data <- svl_known %>%
  filter(island_age %in% supported) %>%
  mutate(island_f = factor(island_age))

if(n_distinct(svl_fit_data$island_f) >= 2) {

  m_common <- gam(
    svl_num ~ s(age_at_capture, k = 7) + sex_f + island_f + s(ID_f, bs = "re"),
    data = svl_fit_data,
    method = "REML"
  )

  m_island_age <- gam(
    svl_num ~ island_f + sex_f +
      s(age_at_capture, by = island_f, k = 6) +
      s(ID_f, bs = "re"),
    data = svl_fit_data,
    method = "REML"
  )

  m_time <- gam(
    svl_num ~ s(age_at_capture, k = 7) + sex_f + island_f +
      calendar_decade + s(ID_f, bs = "re"),
    data = svl_fit_data,
    method = "REML"
  )

  model_compare <- tibble(
    model = c(
      "common_age_smooth",
      "island_specific_age_smooths",
      "common_age_plus_calendar_time"
    ),
    AIC = c(AIC(m_common), AIC(m_island_age), AIC(m_time)),
    deviance_explained = c(
      summary(m_common)$dev.expl,
      summary(m_island_age)$dev.expl,
      summary(m_time)$dev.expl
    )
  ) %>%
    arrange(AIC)

  cat("\n============================================================\n")
  cat("KNOWN-AGE SVL MODEL COMPARISON\n")
  cat("============================================================\n")
  print(model_compare, n = Inf)

  cat("\n============================================================\n")
  cat("SVL + CALENDAR-TIME MODEL\n")
  cat("============================================================\n")
  print(summary(m_time))

  write_csv(
    model_compare,
    file.path(out_dir, "Audit_KnownAge_SVL_ModelComparison.csv")
  )

  saveRDS(
    list(
      common = m_common,
      island_age = m_island_age,
      time = m_time
    ),
    file.path(out_dir, "Audit_KnownAge_SVL_Models.rds")
  )
}

## 5. Provisional tourism/provisioning history --------------------------

provisioning_history <- tribble(
  ~site_name, ~period_start, ~period_end, ~state, ~candidate_exposure, ~notes,
  "Leaf Cay North", 1980, 2018, "present and increasing", "high/regular", "Already present in 1980 (~20 visitors/day); approximate increase thereafter.",
  "Leaf Cay North", 2019, 2023, "marked reduction", "low/intermittent", "Visitation reported to drop sharply 2019-2023.",
  "U Cay", 1980, 2018, "present and increasing", "high/regular", "Already present in 1980 (~5 visitors/day); approximate increase thereafter.",
  "U Cay", 2019, 2023, "marked reduction", "low/intermittent", "Visitation reported to drop sharply 2019-2023.",
  "Bitter Guana", 1980, 2005, "no organised tourism described", "none", "Organised tour boats reported to begin in 2006.",
  "Bitter Guana", 2006, 2019, "intense daily feeding", "high/regular", "Tour boats began in 2006; intense daily feeding thereafter.",
  "Gaulin", 1980, 2009, "visited/fed at lower intensity", "low/intermittent", "Longstanding limited feeding, mainly west beach/cave.",
  "Gaulin", 2010, 2019, "declining activity", "low/intermittent", "Activity declined as Bitter Guana became more popular.",
  "White Bay", 1980, 1997, "feeding not documented", "none", "Feeding thought to begin in late 1990s.",
  "White Bay", 1998, 2009, "limited/local feeding", "low/intermittent", "Local visitors; not mass tourism.",
  "White Bay", 2010, 2019, "waning feeding", "low/intermittent", "Feeding began to wane around 2010.",
  "Leaf Cay South", 1980, 2011, "feeding not documented", "none", "Feeding reported to begin around 2012.",
  "Leaf Cay South", 2012, 2019, "intense daily feeding", "high/regular", "Intense daily feeding thereafter.",
  "Noddy", 1980, 2019, "not visited", "none", "Reported not visited.",
  "North Adderly", 1980, 2019, "not visited", "none", "Reported not visited.",
  "FFRC", 1980, 2019, "no provisioning reported", "none", "No tourism/provisioning exposure described."
)

write_csv(
  provisioning_history,
  file.path(out_dir, "PROVISIONAL_TourismProvisioning_History.csv")
)

cat("\n============================================================\n")
cat("PROVISIONAL TOURISM / PROVISIONING HISTORY\n")
cat("============================================================\n")
print(provisioning_history, n = Inf)

## 6. FFRC disturbance history -----------------------------------------

disturbance_history <- tibble(
  site_name = "FFRC",
  event_year = 2016,
  event = "Hurricane Matthew",
  primary_candidate = "post-2016 regime indicator",
  sensitivity_candidate = "acute 2016-2018 window"
)

write_csv(
  disturbance_history,
  file.path(out_dir, "PROVISIONAL_Disturbance_History.csv")
)

cat("\n============================================================\n")
cat("AUDIT COMPLETE\n")
cat("============================================================\n")
cat("Outputs:", out_dir, "\n")
cat("Next: resolve Leaf North/South mapping before final exposure coding.\n")
cat("============================================================\n")

