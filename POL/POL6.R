source("DataPrep.R")


## ============================================================
## BIRTH-TIME / SVL AUDIT
## ============================================================

library(tidyverse)

cat("\n")
cat("============================================================\n")
cat("1. RELEVANT DataPrep.R SOURCE CODE\n")
cat("============================================================\n")

prep_lines <- readLines(
  "DataPrep.R",
  warn = FALSE
)

hits <- grep(
  "svl|birth|cintB|ageR|quantreg|rq\\(",
  prep_lines,
  ignore.case = TRUE
)

show_lines <- sort(
  unique(
    unlist(
      lapply(
        hits,
        function(i) {
          seq(
            max(1, i - 4),
            min(length(prep_lines), i + 7)
          )
        }
      )
    )
  )
)

cat(
  paste0(
    sprintf("%4d", show_lines),
    ": ",
    prep_lines[show_lines],
    collapse = "\n"
  ),
  "\n"
)


## ============================================================
## 2. RAW BIRTH INFORMATION BY INDIVIDUAL
## ============================================================

birth_raw <- AllData %>%
  group_by(ID) %>%
  summarise(
    Birth_raw = {
      z <- suppressWarnings(
        as.numeric(Birth)
      )
      
      z <- z[
        is.finite(z) &
          z > 0
      ]
      
      if(length(z)) z[1] else NA_real_
    },
    .groups = "drop"
  )


## ============================================================
## 3. SVL INFORMATION BY INDIVIDUAL
## ============================================================

svl_raw <- AllData %>%
  mutate(
    svl_num = suppressWarnings(
      as.numeric(svl)
    )
  ) %>%
  arrange(
    ID,
    Year,
    CapDate
  ) %>%
  group_by(ID) %>%
  summarise(
    n_svl = sum(
      is.finite(svl_num)
    ),
    
    first_svl = {
      z <- svl_num[
        is.finite(svl_num)
      ]
      
      if(length(z)) z[1] else NA_real_
    },
    
    min_svl = {
      z <- svl_num[
        is.finite(svl_num)
      ]
      
      if(length(z)) min(z) else NA_real_
    },
    
    .groups = "drop"
  )


## ============================================================
## 4. MODEL BIRTH WINDOWS
## ============================================================

stopifnot(
  nrow(growth_info) == nrow(cintB_adj)
)

birth_audit <- growth_info %>%
  transmute(
    ID,
    island,
    first_capture_year = tF_cal,
    last_capture_year = tL_cal
  ) %>%
  mutate(
    birth_lower_rel = cintB_adj[,1],
    birth_upper_rel = cintB_adj[,2],
    
    birth_lower_cal = Origin + birth_lower_rel,
    birth_upper_cal = Origin + birth_upper_rel,
    
    birth_mid_cal = (
      birth_lower_cal +
        birth_upper_cal
    ) / 2,
    
    birth_window_width = (
      birth_upper_cal -
        birth_lower_cal
    )
  ) %>%
  left_join(
    birth_raw,
    by = "ID"
  ) %>%
  left_join(
    svl_raw,
    by = "ID"
  ) %>%
  mutate(
    birth_known_raw = !is.na(Birth_raw),
    has_svl = n_svl > 0
  )


## ============================================================
## 5. BASIC COUNTS
## ============================================================

cat("\n")
cat("============================================================\n")
cat("BIRTH INFORMATION AUDIT\n")
cat("============================================================\n")

cat(
  "Individuals:",
  nrow(birth_audit),
  "\n"
)

cat(
  "Raw known birth:",
  sum(
    birth_audit$birth_known_raw
  ),
  "\n"
)

cat(
  "Raw unknown birth:",
  sum(
    !birth_audit$birth_known_raw
  ),
  "\n"
)

cat(
  "Unknown birth WITH SVL:",
  sum(
    !birth_audit$birth_known_raw &
      birth_audit$has_svl
  ),
  "\n"
)

cat(
  "Unknown birth WITHOUT SVL:",
  sum(
    !birth_audit$birth_known_raw &
      !birth_audit$has_svl
  ),
  "\n"
)


## ============================================================
## 6. BIRTH-WINDOW WIDTHS
## ============================================================

cat("\nBirth-window widths:\n")

print(
  birth_audit %>%
    mutate(
      group = case_when(
        birth_known_raw ~
          "Known birth",
        
        !birth_known_raw & has_svl ~
          "Unknown birth + SVL",
        
        TRUE ~
          "Unknown birth + no SVL"
      )
    ) %>%
    group_by(group) %>%
    summarise(
      n = n(),
      min_width = min(
        birth_window_width,
        na.rm = TRUE
      ),
      median_width = median(
        birth_window_width,
        na.rm = TRUE
      ),
      mean_width = mean(
        birth_window_width,
        na.rm = TRUE
      ),
      max_width = max(
        birth_window_width,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
)


## ============================================================
## 7. CHECK UNKNOWN-BIRTH ANIMALS WITHOUT SVL
## ============================================================

no_svl_births <- birth_audit %>%
  filter(
    !birth_known_raw,
    !has_svl
  )

cat(
  "\nUnknown-birth individuals without SVL:",
  nrow(no_svl_births),
  "\n"
)

if(nrow(no_svl_births) > 0) {
  print(
    no_svl_births,
    n = Inf
  )
}


## ============================================================
## 8. CHECK FOR INVALID BIRTH WINDOWS
## ============================================================

cat("\nInvalid windows:\n")

cat(
  "NA endpoints:",
  sum(
    !is.finite(
      birth_audit$birth_lower_cal
    ) |
      !is.finite(
        birth_audit$birth_upper_cal
      )
  ),
  "\n"
)

cat(
  "Upper <= lower:",
  sum(
    birth_audit$birth_upper_cal <=
      birth_audit$birth_lower_cal,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Birth midpoint after first capture:",
  sum(
    birth_audit$birth_mid_cal >
      birth_audit$first_capture_year,
    na.rm = TRUE
  ),
  "\n"
)


## ============================================================
## 9. REFERENCE-AGE CANDIDATE
##
## Use midpoint of each individual's observed age range.
## This is ONLY used to choose a fixed parameterisation pivot.
## Birth remains latent in the survival model.
## ============================================================

birth_audit <- birth_audit %>%
  mutate(
    observed_age_mid = (
      (
        first_capture_year +
          last_capture_year
      ) / 2
    ) -
      birth_mid_cal
  )

age_ref_all <- median(
  birth_audit$observed_age_mid[
    is.finite(
      birth_audit$observed_age_mid
    ) &
      birth_audit$observed_age_mid >= 0
  ],
  na.rm = TRUE
)

age_ref_all_round <- round(
  age_ref_all
)

cat(
  "\nMedian midpoint age represented in model:",
  round(age_ref_all, 2),
  "years\n"
)

cat(
  "Candidate reference age:",
  age_ref_all_round,
  "years\n"
)


## ============================================================
## 10. KNOWN-BIRTH-ONLY CHECK
##
## Useful independent comparison: does the candidate reference
## age look similar when based only on actual known birth years?
## ============================================================

known_ref <- birth_audit %>%
  filter(
    birth_known_raw
  ) %>%
  mutate(
    known_observed_age_mid = (
      (
        first_capture_year +
          last_capture_year
      ) / 2
    ) -
      Birth_raw
  )

age_ref_known <- median(
  known_ref$known_observed_age_mid[
    is.finite(
      known_ref$known_observed_age_mid
    ) &
      known_ref$known_observed_age_mid >= 0
  ],
  na.rm = TRUE
)

cat(
  "Median midpoint age from KNOWN births only:",
  round(age_ref_known, 2),
  "years\n"
)

cat(
  "\nDifference between all-individual and known-birth medians:",
  round(
    age_ref_all -
      age_ref_known,
    2
  ),
  "years\n"
)