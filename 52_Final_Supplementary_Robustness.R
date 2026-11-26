# 52_Final_Supplementary_Robustness.R
# Post-processing only: no model fitting.
library(tidyverse)

root <- "DataAudits/Final_Robustness"
figdir <- "Figures/Final_Robustness"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
x <- read_csv(file.path(root, "Final_Robustness_Summary.csv"), show_col_types = FALSE)
stopifnot(nrow(x) == 14L, !anyDuplicated(x$scenario))
fmt <- function(p, digits = 3) {
  cols <- paste0(p, c("_median", "_low", "_high"))
  stopifnot(all(cols %in% names(x)))
  f <- paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)")
  sprintf(f, x[[cols[1]]], x[[cols[2]]], x[[cols[3]]])
}
tab <- tibble(
  Model = x$scenario,
  `Male age-10 log HR` = fmt("betaSEX[1]"),
  `Male log senescence-rate ratio` = fmt("betaSEX[2]"),
  `Density age-10 log HR` = fmt("betaDENS_a"),
  `Density log senescence-rate ratio` = fmt("betaDENS_b"),
  `Cohort log HR` = fmt("betaCOHORT"),
  `Provisioning HR` = fmt("Provisioning HR", 2),
  `P(provisioning HR > 1)` = sprintf("%.3f", x[["Provisioning HR_Pr_gt_1"]]),
  `Matthew HR` = fmt("Matthew HR", 2),
  `P(Matthew HR > 1)` = sprintf("%.3f", x[["Matthew HR_Pr_gt_1"]]),
  `Max Rhat` = sprintf("%.3f", x$max_Rhat),
  `Provisioning contrast` = if_else(x$scenario == "Sensitivity: any provisioning", "Any vs none", "High/regular vs other"),
  `Matthew contrast` = case_when(
    x$scenario == "LOO: FFRC" ~ "Not estimable: FFRC omitted",
    x$scenario == "Sensitivity: acute Matthew period" ~ "FFRC 2017-2018",
    TRUE ~ "FFRC post-2016"
  )
)
tab[tab$Model == "LOO: FFRC", c("Matthew HR", "P(Matthew HR > 1)")] <- NA_character_
write_csv(tab, file.path(root, "Supplementary_Robustness_Table.csv"), na = "NA")
writeLines(c(
  "Posterior medians (95% credible intervals). Coefficients are on the log scale; calendar exposure effects are hazard ratios.",
  "The any-provisioning and acute-Matthew analyses change the exposure contrasts, as identified in the table.",
  "The Matthew coefficient is not informed by exposure data when FFRC is omitted and is therefore not reported for that fit.",
  "Density coefficients apply to the density index used in each model; partial-support sensitivity changes eligibility.",
  "Single-site omission fits are influence analyses, not independent replications."
), file.path(root, "Supplementary_Robustness_Table_Notes.txt"))

d <- read_csv(file.path(root, "Final_Robustness_KeyEffects_Long.csv"), show_col_types = FALSE)
stopifnot(nrow(d) == 14L * 7L)
d <- d %>%
  mutate(
    scenario = factor(scenario, levels = rev(x$scenario)),
    effect = factor(effect, levels = unique(effect)),
    class = factor(class, levels = c("Primary", "Sensitivity", "Targeted exclusion", "Leave-one-site-out"))
  ) %>%
  filter(!(scenario == "LOO: FFRC" & effect == "Matthew HR"))
refs <- d %>%
  distinct(effect, scale) %>%
  mutate(reference = if_else(scale == "One = no effect", 1, 0))
p <- ggplot(d, aes(median, scenario, shape = class, colour = class)) +
  geom_vline(data = refs, aes(xintercept = reference), linetype = 2, linewidth = 0.4, colour = "grey55") +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = 0.12, linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~effect, scales = "free_x", ncol = 2, labeller = label_wrap_gen(35)) +
  scale_colour_manual(values = c("black", "#38678D", "#A33F32", "grey45")) +
  labs(x = "Posterior median and 95% credible interval", y = NULL, colour = NULL, shape = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(figdir, "Supplementary_Robustness.pdf"), p, width = 12, height = 15)
ggsave(file.path(figdir, "Supplementary_Robustness.png"), p, width = 12, height = 15, dpi = 300)
print(tab, width = Inf, n = Inf)
cat("\nSupplementary table and figure saved. No models refitted.\n")
