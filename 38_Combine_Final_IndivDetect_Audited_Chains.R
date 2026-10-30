## ============================================================
## 11_Combine_TimeVarying_Feeding_Matthew_LONG_Chains_V3.R
##
## Combine independently run LONG chains, calculate convergence,
## posterior summaries, feeding/Matthew hazard ratios, and latent
## sex posterior probabilities.
## ============================================================


library(tidyverse)
library(coda)
library(MCMCvis)


run_dir <-
  "ModelRuns/Final_IndivDetect_Audited"


chain_ids <- c(
  1L,
  2L,
  3L
)


chain_files <- file.path(
  run_dir,
  paste0(
    "Model_Final_IndivDetect_Audited_CHAIN",
    chain_ids,
    ".rds"
  )
)


missing_files <- chain_files[
  !file.exists(
    chain_files
  )
]


if(
  length(
    missing_files
  ) >
    0
) {
  stop(
    "Missing chain file(s):\n",
    paste(
      missing_files,
      collapse =
        "\n"
    )
  )
}


chain_objects <- lapply(
  chain_files,
  readRDS
)


## Refuse to combine chains generated from different code or inputs.
if(
  any(
    vapply(
      chain_objects,
      function(x) is.null(x$analysis_hashes),
      logical(1)
    )
  )
) {
  stop(
    "One or more chain files lack analysis_hashes. ",
    "Do not combine legacy and audited chains."
  )
}

reference_hashes <- chain_objects[[1]]$analysis_hashes

for(k in seq_along(chain_objects)) {
  if(
    !identical(
      chain_objects[[k]]$analysis_hashes,
      reference_hashes
    )
  ) {
    stop(
      "Chain ",
      k,
      " was generated from different model/data inputs. ",
      "Refusing to combine incompatible chains."
    )
  }
}

cat("\nVerified identical model/data hashes across all chains.\n")
print(reference_hashes)


loaded_chain_ids <- vapply(
  chain_objects,
  function(x) {
    as.integer(
      x$chain_id
    )
  },
  integer(
    1
  )
)


if(
  !identical(
    loaded_chain_ids,
    chain_ids
  )
) {
  stop(
    "Chain-ID mismatch. Expected: ",
    paste(
      chain_ids,
      collapse = ", "
    ),
    "; loaded: ",
    paste(
      loaded_chain_ids,
      collapse = ", "
    )
  )
}


## Ensure all chains used identical run settings.
settings <- map_dfr(
  chain_objects,
  function(x) {

    tibble(
      chain_id =
        x$chain_id,

      niter =
        x$niter,

      nburnin =
        x$nburnin,

      thin =
        x$thin,

      seed =
        x$chain_seed
    )
  }
)


if(
  n_distinct(
    settings$niter
  ) !=
    1 |
  n_distinct(
    settings$nburnin
  ) !=
    1 |
  n_distinct(
    settings$thin
  ) !=
    1
) {
  stop(
    "Chains were not run with identical NITER/NBURN/THIN settings."
  )
}


extract_one_chain <- function(
    obj
) {

  x <- obj$fit$samples


  if(
    inherits(
      x,
      "mcmc.list"
    )
  ) {

    if(
      length(
        x
      ) !=
        1
    ) {
      stop(
        "Expected exactly one chain in each chain file."
      )
    }

    return(
      x[[
          1
  ]]
    )
  }


  if(
    inherits(
      x,
      "mcmc"
    )
  ) {
    return(
      x
    )
  }


  as.mcmc(
    x
  )
}


chain_samples <- lapply(
  chain_objects,
  extract_one_chain
)


reference_names <- colnames(
  chain_samples[[
      1
  ]]
)


for(k in seq_along(chain_samples)) {

  if(
    !identical(
      colnames(
        chain_samples[[
            k
  ]]
      ),
      reference_names
    )
  ) {
    stop(
      "Parameter columns differ between chains."
    )
  }
}


samples_combined <- mcmc.list(
  chain_samples
)


primary_parameters <- c(
  "a10_ref",
  "b_ref",
  "log_a10",
  "log_b",
  "betaSEX",
  "betaDENS_a",
  "betaDENS_b",
  "betaCOHORT",
  "betaFEED",
  "betaMATTHEW",
  "sigma_site_period",
  "sigma_island_a",
  "alpha_p",
  "betaP_SEX",
  "betaP_AGE",
  "sigma_island_p",
  "sigma_survey_p",
  "sigma_indiv_p"
)


summary_primary <- MCMCsummary(
  samples_combined,
  params =
    primary_parameters
) %>%
  as.data.frame() %>%
  rownames_to_column(
    "parameter"
  ) %>%
  as_tibble()


sample_matrix <- as.matrix(
  samples_combined
)


feeding_hr <- exp(
  sample_matrix[
    ,
    "betaFEED"
  ]
)


matthew_hr <- exp(
  sample_matrix[
    ,
    "betaMATTHEW"
  ]
)


calendar_effect_summary <- tibble(
  effect =
    c(
      "High/regular provisioning vs other",
      "FFRC post-2016 Matthew period"
    ),

  median_HR =
    c(
      median(
        feeding_hr
      ),
      median(
        matthew_hr
      )
    ),

  low_HR =
    c(
      unname(
        quantile(
          feeding_hr,
          0.025
        )
      ),
      unname(
        quantile(
          matthew_hr,
          0.025
        )
      )
    ),

  high_HR =
    c(
      unname(
        quantile(
          feeding_hr,
          0.975
        )
      ),
      unname(
        quantile(
          matthew_hr,
          0.975
        )
      )
    ),

  Pr_HR_gt_1 =
    c(
      mean(
        feeding_hr >
          1
      ),
      mean(
        matthew_hr >
          1
      )
    )
)


latent_lookup <-
  chain_objects[[
      1
  ]]$latent_lookup


latent_sex_summary <- latent_lookup


latent_sex_summary$Pr_male <- vapply(
  latent_sex_summary$ID,
  function(i) {

    node_name <- paste0(
      "sex[",
      i,
      "]"
    )

    mean(
      sample_matrix[
        ,
        node_name
      ]
    )
  },
  numeric(
    1
  )
)


latent_sex_summary <- latent_sex_summary %>%
  mutate(
    posterior_class =
      case_when(
        Pr_male >=
          0.90 ~
          "strong_male",

        Pr_male <=
          0.10 ~
          "strong_female",

        TRUE ~
          "uncertain"
      )
  ) %>%
  arrange(
    desc(
      abs(
        Pr_male -
          0.5
      )
    )
  )


## Full convergence table for all monitored nodes.
summary_all <- MCMCsummary(
  samples_combined
) %>%
  as.data.frame() %>%
  rownames_to_column(
    "parameter"
  ) %>%
  as_tibble()


convergence_flags <- summary_all %>%
  filter(
    !is.na(
      Rhat
    ),
    Rhat >
      1.05
  ) %>%
  arrange(
    desc(
      Rhat
    )
  )


cat("\n")
cat("============================================================\n")
cat("FINAL AUDITED INDIVID-DETECTION MODEL: COMBINED SUMMARY\n")
cat("============================================================\n")

print(
  settings,
  n = Inf
)


cat("\nPRIMARY PARAMETERS\n")
cat("------------------------------------------------------------\n")

print(
  summary_primary,
  n = Inf
)


cat("\nCALENDAR-EXPOSURE EFFECTS\n")
cat("------------------------------------------------------------\n")

print(
  calendar_effect_summary,
  n = Inf
)


cat("\nLATENT SEX\n")
cat("------------------------------------------------------------\n")
cat(
  "Strong male P>=0.90: ",
  sum(
    latent_sex_summary$Pr_male >=
      0.90
  ),
  "\n",
  sep =
    ""
)
cat(
  "Strong female P<=0.10: ",
  sum(
    latent_sex_summary$Pr_male <=
      0.10
  ),
  "\n",
  sep =
    ""
)
cat(
  "Remaining uncertain: ",
  sum(
    latent_sex_summary$Pr_male >
      0.10 &
      latent_sex_summary$Pr_male <
      0.90
  ),
  "\n",
  sep =
    ""
)


cat("\nCONVERGENCE: RHAT > 1.05\n")
cat("------------------------------------------------------------\n")

if(
  nrow(
    convergence_flags
  ) ==
    0
) {

  cat(
    "None.\n"
  )

} else {

  print(
    convergence_flags,
    n = Inf
  )
}


out_dir <-
  "DataAudits/Final_IndivDetect_Audited"


dir.create(
  out_dir,
  recursive =
    TRUE,
  showWarnings =
    FALSE
)


write_csv(
  settings,
  file.path(
    out_dir,
    "LongRun_ChainSettings.csv"
  )
)


write_csv(
  summary_primary,
  file.path(
    out_dir,
    "LongRun_PrimarySummary.csv"
  )
)


write_csv(
  summary_all,
  file.path(
    out_dir,
    "LongRun_AllMonitoredSummary.csv"
  )
)


write_csv(
  convergence_flags,
  file.path(
    out_dir,
    "LongRun_RhatAbove1.05.csv"
  )
)


write_csv(
  calendar_effect_summary,
  file.path(
    out_dir,
    "LongRun_CalendarEffects.csv"
  )
)


write_csv(
  latent_sex_summary,
  file.path(
    out_dir,
    "LongRun_LatentSexPosterior.csv"
  )
)


saveRDS(
  list(
    samples =
      samples_combined,

    summary =
      summary_primary,

    summary_all =
      summary_all,

    convergence_flags =
      convergence_flags,

    calendar_effect_summary =
      calendar_effect_summary,

    latent_sex =
      latent_sex_summary,

    chain_settings =
      settings,

    density_eligible =
      chain_objects[[
          1
  ]]$density_eligible,

    density_support =
      chain_objects[[
          1
  ]]$birth_density_support,

    feeding_high_matrix =
      chain_objects[[
          1
  ]]$feeding_high_matrix,

    matthew_post_matrix =
      chain_objects[[
          1
  ]]$matthew_post_matrix,

    period_ref_year =
      chain_objects[[
          1
  ]]$period_ref_year,

    visit_gap_days =
      chain_objects[[
          1
  ]]$visit_gap_days,

    survey_events_used =
      chain_objects[[
          1
  ]]$survey_events_used,

    survey_events_all =
      chain_objects[[
          1
  ]]$survey_events_all,

    Origin =
      chain_objects[[
          1
  ]]$Origin,

    tMax =
      chain_objects[[
          1
  ]]$tMax,

    ordered_islands =
      chain_objects[[
          1
  ]]$ordered_islands
  ),
  "Model_Final_IndivDetect_Audited_3CHAINS_COMPLETE.rds"
)


cat("\n")
cat("============================================================\n")
cat("FINAL AUDITED INDIVID-DETECTION COMBINATION COMPLETE\n")
cat("============================================================\n")
cat("Saved:\n")
cat("  Model_TimeVarying_Feeding_Matthew_LONG_3CHAINS_COMPLETE.rds\n")
cat("  DataAudits/Final_IndivDetect_Audited/*.csv\n")
cat("============================================================\n")
