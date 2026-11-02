## ============================================================
## 40_Audit_Male_Disappearance_LeafNorth_UCay.R
##
## Targeted descriptive sensitivity requested by collaborators:
## did large/older males disappear unusually often from Leaf North
## and U Cay in the mid-2000s, and was the same pattern present at
## the other populations?
##
## This is deliberately descriptive. It does not label unexplained
## disappearance as mortality or human removal.
## ============================================================
library(tidyverse)
library(lubridate)
Sys.unsetenv("IGUANA_EXCLUDE_SITES")
source("DataPrep_Canonical_LatentSex_V5.R")

out_dir <- "DataAudits/MaleRemoval_Sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

live <- AllData %>%
  filter(Status == 1L | capture_then_death) %>%
  arrange(ID_raw, CapDate)

id_terminal <- growth_info %>%
  transmute(ID_raw, Island = island, sex = sex_clean, last_live_time = tL_cal, known_death = !is.na(tD_cal), birth_mid = rowMeans(cintB_cal)) %>%
  left_join(live %>% group_by(ID_raw) %>% slice_max(CapDate, n = 1, with_ties = FALSE) %>% ungroup() %>% transmute(ID_raw, last_date = CapDate, last_year = CalYear, last_svl = svl, last_mass = body_mass, last_comments = comments, permanent_removal = permanent_removal), by = "ID_raw") %>%
  mutate(est_age_last = as.numeric(decimal_date(last_date)) - birth_mid)

site_years <- live %>% distinct(Island, CalYear)
site_last_year <- site_years %>% group_by(Island) %>% summarise(site_end = max(CalYear), .groups = "drop")

id_terminal <- id_terminal %>%
  left_join(site_last_year, by = "Island") %>%
  mutate(years_of_monitoring_after = site_end - last_year, unexplained_terminal = !known_death & !permanent_removal & years_of_monitoring_after >= 2)

live_year <- live %>%
  distinct(ID_raw, Island, CalYear) %>%
  left_join(id_terminal %>% select(ID_raw, sex, last_year, unexplained_terminal, est_age_last, last_svl), by = "ID_raw") %>%
  mutate(terminal_this_year = unexplained_terminal & CalYear == last_year)

terminal_rate <- live_year %>%
  filter(sex %in% c("F", "M")) %>%
  group_by(Island, CalYear, sex) %>%
  summarise(n_seen = n_distinct(ID_raw), n_terminal = n_distinct(ID_raw[terminal_this_year]), terminal_fraction = n_terminal / n_seen, .groups = "drop")

write_csv(id_terminal, file.path(out_dir, "Individual_terminal_histories.csv"))
write_csv(terminal_rate, file.path(out_dir, "Terminal_fraction_by_site_year_sex.csv"))

focus <- terminal_rate %>%
  filter(CalYear >= 2000, CalYear <= 2010)
write_csv(focus, file.path(out_dir, "Terminal_fraction_2000_2010.csv"))

male_terminal <- id_terminal %>%
  filter(sex == "M", unexplained_terminal) %>%
  arrange(Island, last_year, desc(est_age_last))
write_csv(male_terminal, file.path(out_dir, "Male_unexplained_terminal_records.csv"))

period_summary <- id_terminal %>%
  filter(sex %in% c("F", "M")) %>%
  mutate(period = case_when(last_year >= 2003 & last_year <= 2007 ~ "2003-2007", last_year >= 2000 & last_year <= 2010 ~ "Other 2000-2010", TRUE ~ "Outside 2000-2010")) %>%
  group_by(Island, sex, period) %>%
  summarise(n_individuals = n(), n_unexplained_terminal = sum(unexplained_terminal, na.rm = TRUE), median_age_terminal = median(est_age_last[unexplained_terminal], na.rm = TRUE), median_svl_terminal = median(last_svl[unexplained_terminal], na.rm = TRUE), .groups = "drop") %>%
  mutate(median_age_terminal = if_else(is.nan(median_age_terminal), NA_real_, median_age_terminal), median_svl_terminal = if_else(is.nan(median_svl_terminal), NA_real_, median_svl_terminal))
write_csv(period_summary, file.path(out_dir, "Period_summary_site_sex.csv"))

p1 <- terminal_rate %>%
  filter(CalYear <= 2017, sex %in% c("F", "M")) %>%
  ggplot(aes(CalYear, terminal_fraction, linetype = sex, group = sex)) +
  geom_line() + geom_point() + facet_wrap(~Island, scales = "free_x") +
  labs(x = "Calendar year", y = "Fraction of animals seen that year never recaptured", linetype = "Sex", title = "Terminal-disappearance fraction by site, year and sex") +
  theme_bw()
ggsave(file.path(out_dir, "Terminal_fraction_by_site_year_sex.png"), p1, width = 12, height = 8, dpi = 300)

p2 <- male_terminal %>%
  filter(last_year >= 1998, last_year <= 2012) %>%
  ggplot(aes(last_year, est_age_last)) +
  geom_point(aes(size = last_svl), alpha = 0.65) + facet_wrap(~Island) +
  labs(x = "Year of last capture", y = "Estimated age at last capture (birth-window midpoint)", size = "Last SVL", title = "Unexplained terminal male captures around the mid-2000s") +
  theme_bw()
ggsave(file.path(out_dir, "Male_terminal_age_size_1998_2012.png"), p2, width = 12, height = 8, dpi = 300)

cat("\nMale-removal descriptive audit complete.\n")
cat("Primary comparison files:\n")
cat("  ", file.path(out_dir, "Terminal_fraction_2000_2010.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "Male_unexplained_terminal_records.csv"), "\n", sep = "")
cat("  ", file.path(out_dir, "Terminal_fraction_by_site_year_sex.png"), "\n", sep = "")
