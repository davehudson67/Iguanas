## ============================================================
## 26_Build_IndivDetect_Sensitivity_Files.R
##
## Build sensitivity-analysis scripts matching the revised primary
## model with persistent individual detection heterogeneity.
##
## Uses the user's LOCAL 12/13 sensitivity scripts as the source.
## Existing scripts/results are not modified.
## ============================================================

model_src <- "12_Model_Sensitivity_LONG_1CHAIN.R"
combine_src <- "13_Combine_Sensitivity_LONG_Chains.R"

model_out <- "26_Model_Sensitivity_IndivDetect_LONG_1CHAIN.R"
combine_out <- "27_Combine_Sensitivity_IndivDetect_LONG_Chains.R"
launcher_out <- "28_Launch_All_IndivDetect_Sensitivities.sh"
combine_all_out <- "29_Combine_All_IndivDetect_Sensitivities.sh"

stopifnot(
  file.exists(model_src),
  file.exists(combine_src)
)

replace_exact <- function(x, old, new, expected = 1L, label = "replacement") {

  pos <- gregexpr(
    old,
    x,
    fixed = TRUE
  )[[1]]

  n <- if(
    length(pos) == 1L &&
    pos[1] == -1L
  ) {
    0L
  } else {
    length(pos)
  }

  if(
    !is.null(expected) &&
    n != expected
  ) {
    stop(
      label,
      ": expected ",
      expected,
      " match(es); found ",
      n
    )
  }

  if(
    n == 0L
  ) {
    stop(
      label,
      ": no matches found"
    )
  }

  gsub(
    old,
    new,
    x,
    fixed = TRUE
  )
}


## ============================================================
## 1. BUILD REVISED SENSITIVITY MODEL
## ============================================================

txt <- paste(
  readLines(
    model_src,
    warn = FALSE
  ),
  collapse = "\n"
)

life_marker <- paste0(
  "  ## ----------------------------------------------------------\n",
  "  ## INDIVIDUAL LIFE HISTORIES\n",
  "  ## ----------------------------------------------------------"
)

indiv_re_block <- paste0(
  "  ## ----------------------------------------------------------\n",
  "  ## PERSISTENT INDIVIDUAL DETECTION HETEROGENEITY\n",
  "  ## ----------------------------------------------------------\n\n",
  "  sigma_indiv_p ~ dunif(\n",
  "    0,\n",
  "    3\n",
  "  )\n\n",
  "  for(i in 1:nind) {\n\n",
  "    indiv_p_raw[i] ~ dnorm(\n",
  "      0,\n",
  "      sd =\n",
  "        1\n",
  "    )\n",
  "  }\n\n\n",
  life_marker
)

txt <- replace_exact(
  txt,
  life_marker,
  indiv_re_block,
  expected = 1L,
  label = "individual RE block"
)

txt <- replace_exact(
  txt,
  paste0(
    "      survey_p[\n",
    "        obs_event[k]\n",
    "      ]"
  ),
  paste0(
    "      survey_p[\n",
    "        obs_event[k]\n",
    "      ] +\n",
    "      sigma_indiv_p *\n",
    "      indiv_p_raw[\n",
    "        obs_ind[k]\n",
    "      ]"
  ),
  expected = 1L,
  label = "detection predictor"
)

txt <- replace_exact(
  txt,
  paste0(
    "    survey_p_raw =\n",
    "      rnorm(\n",
    "        n_events,\n",
    "        0,\n",
    "        0.3\n",
    "      )\n",
    "  )"
  ),
  paste0(
    "    survey_p_raw =\n",
    "      rnorm(\n",
    "        n_events,\n",
    "        0,\n",
    "        0.3\n",
    "      ),\n\n",
    "    sigma_indiv_p =\n",
    "      runif(\n",
    "        1,\n",
    "        0.2,\n",
    "        0.8\n",
    "      ),\n\n",
    "    indiv_p_raw =\n",
    "      rnorm(\n",
    "        nind,\n",
    "        0,\n",
    "        0.3\n",
    "      )\n",
    "  )"
  ),
  expected = 1L,
  label = "individual RE initial values"
)

txt <- replace_exact(
  txt,
  paste0(
    '  "sigma_survey_p",\n',
    "  latent_sex_nodes"
  ),
  paste0(
    '  "sigma_survey_p",\n',
    '  "sigma_indiv_p",\n',
    "  latent_sex_nodes"
  ),
  expected = 1L,
  label = "sigma_indiv_p monitor"
)

txt <- replace_exact(
  txt,
  paste0(
    '  "sigma_island_p",\n',
    '  "sigma_survey_p"\n',
    ")) {"
  ),
  paste0(
    '  "sigma_island_p",\n',
    '  "sigma_survey_p",\n',
    '  "sigma_indiv_p"\n',
    ")) {"
  ),
  expected = 1L,
  label = "sigma_indiv_p slice sampler"
)

txt <- replace_exact(
  txt,
  paste0(
    "  ANY_FEED =\n",
    "    100000L,\n",
    "  MATTHEW_ACUTE =\n",
    "    200000L,\n",
    "  DENSITY_PARTIAL =\n",
    "    300000L"
  ),
  paste0(
    "  ANY_FEED =\n",
    "    500000L,\n",
    "  MATTHEW_ACUTE =\n",
    "    600000L,\n",
    "  DENSITY_PARTIAL =\n",
    "    700000L"
  ),
  expected = 1L,
  label = "seed offsets"
)

txt <- replace_exact(
  txt,
  paste0(
    '    "Sensitivity_",\n',
    "    MODEL_VARIANT,\n",
    '    "_LONG"'
  ),
  paste0(
    '    "Sensitivity_",\n',
    "    MODEL_VARIANT,\n",
    '    "_IndivDetect_LONG"'
  ),
  expected = 1L,
  label = "model output directory"
)

txt <- replace_exact(
  txt,
  paste0(
    '    "Model_Sensitivity_",\n',
    "    MODEL_VARIANT,\n",
    '    "_LONG_CHAIN",'
  ),
  paste0(
    '    "Model_Sensitivity_",\n',
    "    MODEL_VARIANT,\n",
    '    "_IndivDetect_LONG_CHAIN",'
  ),
  expected = 1L,
  label = "model chain filename"
)

writeLines(
  strsplit(
    txt,
    "\n",
    fixed = TRUE
  )[[1]],
  model_out
)


## ============================================================
## 2. BUILD REVISED COMBINER
## ============================================================

ctxt <- paste(
  readLines(
    combine_src,
    warn = FALSE
  ),
  collapse = "\n"
)

ctxt <- replace_exact(
  ctxt,
  paste0(
    '    "Sensitivity_",\n',
    "    MODEL_VARIANT,\n",
    '    "_LONG"'
  ),
  paste0(
    '    "Sensitivity_",\n',
    "    MODEL_VARIANT,\n",
    '    "_IndivDetect_LONG"'
  ),
  expected = 2L,
  label = "combiner run/audit directories"
)

ctxt <- replace_exact(
  ctxt,
  paste0(
    '    "Model_Sensitivity_",\n',
    "    MODEL_VARIANT,\n",
    '    "_LONG_CHAIN",'
  ),
  paste0(
    '    "Model_Sensitivity_",\n',
    "    MODEL_VARIANT,\n",
    '    "_IndivDetect_LONG_CHAIN",'
  ),
  expected = 1L,
  label = "combiner chain filename"
)

ctxt <- replace_exact(
  ctxt,
  paste0(
    '  "sigma_island_p",\n',
    '  "sigma_survey_p"\n',
    ")"
  ),
  paste0(
    '  "sigma_island_p",\n',
    '  "sigma_survey_p",\n',
    '  "sigma_indiv_p"\n',
    ")"
  ),
  expected = 1L,
  label = "combiner primary parameter list"
)

ctxt <- replace_exact(
  ctxt,
  "_LONG_3CHAINS_COMPLETE.rds",
  "_IndivDetect_LONG_3CHAINS_COMPLETE.rds",
  expected = 2L,
  label = "combined RDS names"
)

writeLines(
  strsplit(
    ctxt,
    "\n",
    fixed = TRUE
  )[[1]],
  combine_out
)


## ============================================================
## 3. LAUNCHER: NINE INDEPENDENT LONG CHAINS
## ============================================================

launcher <- c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  "",
  "cd ~/Iguanas",
  "mkdir -p logs",
  "",
  'variants=("ANY_FEED" "MATTHEW_ACUTE" "DENSITY_PARTIAL")',
  "",
  'for variant in "${variants[@]}"',
  "do",
  "  for chain in 1 2 3",
  "  do",
  '    session="iguana_sens_${variant,,}_c${chain}"',
  '    log="logs/sens_indivdetect_${variant,,}_chain${chain}.log"',
  "",
  '    if tmux has-session -t "$session" 2>/dev/null',
  "    then",
  '      echo "$session already exists; not starting duplicate."',
  "    else",
  '      tmux new-session -d -s "$session" \\',
  '        "cd ~/Iguanas && \\',
  '         MODEL_VARIANT=${variant} \\',
  '         CHAIN_ID=${chain} \\',
  '         NITER=150000 \\',
  '         NBURN=50000 \\',
  '         THIN=10 \\',
  '         Rscript 26_Model_Sensitivity_IndivDetect_LONG_1CHAIN.R \\',
  '         2>&1 | tee ${log}"',
  '      echo "Started $session -> $log"',
  "    fi",
  "  done",
  "done",
  "",
  'tmux ls 2>/dev/null | grep "iguana_sens_" || true'
)

writeLines(
  launcher,
  launcher_out
)

Sys.chmod(
  launcher_out,
  mode = "0755"
)


## ============================================================
## 4. COMBINE ALL THREE VARIANTS AFTER CHAINS FINISH
## ============================================================

combine_all <- c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  "",
  "cd ~/Iguanas",
  "mkdir -p logs",
  "",
  'variants=("ANY_FEED" "MATTHEW_ACUTE" "DENSITY_PARTIAL")',
  "",
  'for variant in "${variants[@]}"',
  "do",
  '  echo "============================================================"',
  '  echo "Combining ${variant}"',
  '  echo "============================================================"',
  "",
  '  MODEL_VARIANT=${variant} \\',
  '  Rscript 27_Combine_Sensitivity_IndivDetect_LONG_Chains.R \\',
  '  2>&1 | tee "logs/combine_sens_indivdetect_${variant,,}.log"',
  "done"
)

writeLines(
  combine_all,
  combine_all_out
)

Sys.chmod(
  combine_all_out,
  mode = "0755"
)


cat("\nCreated:\n")
cat("  ", model_out, "\n", sep = "")
cat("  ", combine_out, "\n", sep = "")
cat("  ", launcher_out, "\n", sep = "")
cat("  ", combine_all_out, "\n", sep = "")
cat("\nPrevious sensitivity scripts/results were not modified.\n")
