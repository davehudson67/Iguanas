# Empirical age support, not estimated survival or an at-risk life table.
# Canonical preparation clears the workspace; source before defining objects.
Sys.unsetenv("IGUANA_EXCLUDE_SITES")
source("DataPrep_Canonical_LatentSex_V5.R")
library(tidyverse)
out <- "DataAudits/Final_AgeSupport"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
stopifnot(nrow(growth_info)==4152L, nrow(cintB_cal)==nrow(growth_info),
          identical(as.integer(growth_info$ID), seq_len(nrow(growth_info))))
ind <- growth_info %>%
  transmute(ID, island, species, sex=coalesce(as.character(sex_clean),"Unresolved"),
            known_age=known_age_consistent, first=tF_cal, last=tL_cal,
            death=tD_cal, birth_low=cintB_cal[,1], birth_high=cintB_cal[,2])
live <- AllData %>%
  filter(alive_encounter) %>%
  transmute(ID, time=lubridate::decimal_date(CapDate)) %>%
  distinct() %>%
  left_join(ind, by="ID") %>%
  mutate(age_low=pmax(0,time-birth_high),
         age_high=pmax(0,time-birth_low),
         age_mid=(age_low+age_high)/2)
stopifnot(all(is.finite(live$age_mid)))
summ <- function(d) {
  d %>% summarise(individuals=n_distinct(ID), live_dates=n(),
    age_min=min(age_low), age_max=max(age_high),
    midpoint_q05=quantile(age_mid,.05), midpoint_median=median(age_mid),
    midpoint_q95=quantile(age_mid,.95), midpoint_max=max(age_mid),
    .groups="drop")
}
write_csv(summ(live %>% group_by(species,island,sex,known_age)),
          file.path(out,"AgeSupport_ByPopulationSexKnownAge.csv"))
write_csv(summ(live %>% group_by(sex,known_age)),
          file.path(out,"AgeSupport_BySexKnownAge.csv"))
# Distinct observed individuals in each five-year band; possible vs midpoint.
# Bands refer to observed live encounters only, not follow-up time at risk.
bands <- map_dfr(seq(0,80,5), function(a) {
  live %>% group_by(species,island,sex,known_age) %>%
    summarise(possible=n_distinct(ID[age_low < a+5 & age_high >= a]),
              midpoint=n_distinct(ID[age_mid >= a & age_mid < a+5]),
              .groups="drop") %>% mutate(age_start=a,age_end=a+5)
})
write_csv(bands,file.path(out,"AgeSupport_ObservedAgeBands.csv"))
# Common display range: 5th–95th percentiles of one first/last midpoint
# per individual, separately for resolved females and males.
# This transparent display convention is not a statistical identifiability threshold.
limits <- ind %>% filter(sex %in% c("F","M")) %>%
  mutate(first_mid=pmax(0,first-(birth_low+birth_high)/2),
         last_mid=pmax(0,last-(birth_low+birth_high)/2)) %>%
  group_by(sex) %>%
  summarise(n=n(), lower=quantile(first_mid,.05),
            upper=quantile(last_mid,.95), .groups="drop")
stopifnot(nrow(limits)==2L)
range <- tibble(age_min=max(0,ceiling(max(limits$lower))),
                age_max=floor(min(limits$upper)))
stopifnot(range$age_max>range$age_min)
write_csv(limits,file.path(out,"AgeSupport_DisplayRangeBySex.csv"))
write_csv(range,file.path(out,"AgeSupport_DisplayRange.csv"))
print(limits)
print(range)
cat("Age support audit complete. Display range is descriptive, not a proof of support at every age.\n")
