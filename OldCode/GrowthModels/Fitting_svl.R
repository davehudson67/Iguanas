library(dplyr)
library(nlme)
library(stringr)

# Load and clean data
AllData <- readRDS("AllData.rds") %>%
  filter(!is.na(svl))

# Identify all rows with “bad” svl or age
bad_rows <- AllData %>%
  filter(
    # detect anything other than digits or dot in the text
    str_detect(svl,  "[^0-9\\.]") |
      str_detect(age, "[^0-9\\.]")
  )

# See which IDs they are
bad_ids <- bad_rows %>% distinct(ID) %>% pull(ID)

# Extract all captures for those individuals
problem_igs <- AllData %>%
  filter(ID %in% bad_ids)

# Inspect
print(problem_igs)

# Correct svl
AllData <- AllData %>%
  mutate(
    svl = case_when(
      str_to_lower(svl) == "dead" ~ NA_real_,
      str_detect(svl, "/") ~ as.numeric(str_extract(svl, "^[^/]+")),
      TRUE ~ as.numeric(svl)
    )
  )

# Correct age
AllData <- AllData %>%
  mutate(
    age = case_when(
      str_detect(age, regex("^\\s*juv\\s*$", ignore_case = TRUE)) ~ NA_character_,
      TRUE ~ age),
    age = str_remove_all(age, "\\?"),
    age = as.numeric(age)
  )

#------------------------------------------------------------------------------#

Ig_svl <- AllData %>%
  group_by(ID) %>%
  arrange(CapDate) %>%
  filter(!is.na(svl)) %>%
  mutate(
    t0 = first(CapDate),
    time = as.numeric(difftime(CapDate, t0, units="days"))/365
  ) %>%
  select(-t0) %>%
  ungroup()

# Self‐start function for Von Bertalanffy
VBSS <- function(time, Asym, k, t0) {
  Asym * (1 - exp(-k * (time - t0)))
}

# You need good starting values: use your known‐age subset to get roughs
known <- filter(Ig_svl, !is.na(age)) %>%
  filter(!is.na(svl)) 

start_vals <- coef(nls(
  svl ~ VBSS(age, Asym, k, t0), data = known,
  start = list(Asym = max(known$svl), k = 0.3, t0 = -0.1)
))

# Fit nlme with a random effect on t0 by individual
fit <- nlme(
  svl    ~ VBSS(time, Asym, k, t0),
  data   = known,
  fixed  = Asym + k + t0 ~ 1,
  random = t0 ~ 1 | ID,
  start  = start_vals,
  control = nlmeControl(
    maxIter    = 100,    # allow more outer iterations (default 50)
    msMaxIter  = 200,    # allow more inner (mixed-effects) iterations
    pnlsTol    = 0.1,    # relax the tolerance on the profiled nls
    pnlsMaxIter = 50     # up from default 5
  )
)

summary(fit)

# Population (fixed) effects:
pop_Asym <- fixef(fit)["Asym"]   # ≈ 64.6 mm
pop_k    <- fixef(fit)["k"]      # ≈ 0.024 per year
pop_t0   <- fixef(fit)["t0"]     # ≈ −14.82 years

# Individual random effects on t0:
ind_t0   <- ranef(fit)$t0        # vector of length = # individuals

re_df <- ranef(fit) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  rename(t0_dev = t0)

# Build a df of true t0_i = pop_t0 + deviation
re_df <- re_df %>%
  mutate(t0_i = pop_t0 + t0_dev) %>%
  select(ID, t0_i)

# Left-join onto full capture data
Ig_svl_k <- known %>%
  left_join(re_df, by = "ID") %>%
  mutate(age_est = time - t0_i)

#------------------------------------------------------------------------------

# Extract known‐age records
known2 <- Ig_svl_k %>%
  filter(!is.na(age))

ggplot(known2, aes(x = age, y = age_est)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x     = "Known age (years)",
    y     = "Predicted age (years)",
    title = "Predicted vs. Known Age for Iguanas"
  ) +
  theme_minimal()
