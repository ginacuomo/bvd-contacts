# =============================================================================
# Ebola Contact Risk Analysis
# Statistical Analysis Plan Implementation
# =============================================================================
# Data: multi-study extraction in standard format
# Mapping: exposure_mapping.csv — maps raw study definitions to canonical levels
# Outcome: Secondary attack rate (SAR) by contact type
#
# Canonical exposure hierarchy (reference = "No/minimal contact"):
#   Class 0: No/minimal contact                          [reference]
#   Class 1: Indirect contact only
#   Class 2: Direct physical contact — no fluids and no nursing
#   Class 3: Nursing care — no body fluids
#   Class 4: Body fluid contact
#   Class 5: Handled corpse
#
# Two analyses:
#   1. Disease prevalence (SAR) by exposure level
#   2. Distribution of contacts across exposure levels
# =============================================================================

library(tidyverse)
library(readxl)
library(lme4)
library(broom.mixed)


# 1. Define the canonical levels -- based on Bower et al. -----------------

canonical_levels <- c(
  "No/minimal contact",                                 # Class 0 — reference
  "Indirect contact only",                              # Class 1
  "Direct physical contact — no fluids and no nursing", # Class 2
  "Nursing care — no body fluids",                      # Class 3
  "Body fluid contact",                                 # Class 4
  "Handled corpse"                                      # Class 5
)

reference_level    <- canonical_levels[1]
non_reference_levels <- canonical_levels[-1]


# 2. Read data ------------------------------------------------------------
raw <- read_excel("data/sar_extraction.xlsx", sheet = "data_me")

# TODO: fix the mapping in the files - this is currently my preliminary mapping
mapping_raw <- read_excel("data/sar_extraction.xlsx", sheet = "exposure_mapping")

# =============================================================================
# 2. DATA CLEANING
# =============================================================================

dat <- raw %>%
  mutate(
    numerator    = as.integer(numerator),
    denominator  = as.integer(denominator),
    year         = as.character(year),
    study_id     = paste(first_author, year_publication, location, sep = "_"),
    sar_observed = numerator / denominator,
    uninfected_contacts     = denominator - numerator
  ) %>%
  filter(
    !is.na(numerator),
    !is.na(denominator),
    denominator > 0,
    numerator <= denominator
  ) %>% filter(include == TRUE)

cat("Studies and contact definitions in dataset:\n")
dat %>% count(first_author, definition_contact_me) %>% print()

# =============================================================================
# 3. PARSE AND VALIDATE MAPPING TABLE
# =============================================================================
 
mapping_table <- mapping_raw %>%
  filter(is.na(include) | include == 1) %>%
  mutate(
    # Split semicolon-separated canonical levels into one row per level
    canonical_level = str_split(canonical_level, ";\\s*")
  ) %>%
  unnest(canonical_level) %>%
  mutate(
    canonical_level = str_trim(canonical_level),
    household       = as.integer(household)
  ) %>%
  # Classify rows: disaggregated (1 canonical level) vs combined (>1)
  group_by(first_author, definition_contact_me) %>%
  mutate(
    n_canonical = n(),
    is_combined = n_canonical > 1
  ) %>%
  ungroup()

# Validate: every canonical level must be in the defined hierarchy
unrecognised <- mapping_table %>%
  filter(!canonical_level %in% canonical_levels) %>%
  distinct(canonical_level)

if (nrow(unrecognised) > 0) {
  stop(
    "Unrecognised canonical levels in mapping table — check spelling:\n",
    paste(unrecognised$canonical_level, collapse = "\n")
  )
} else {
  cat("\nMapping table validated: all canonical levels recognised.\n")
}

# =============================================================================
# 4. JOIN MAPPING ONTO DATA
# =============================================================================
# Join on first_author + definition_contact_me so that identical strings
# used by different authors can map to different canonical levels if needed.

dat_mapped <- dat %>%
  left_join(
    mapping_table %>%
      select(first_author, definition_contact_me, canonical_level,
             household, is_combined, n_canonical),
    by = c("first_author", "definition_contact_me"),
    relationship = "many-to-many"
  )

# Flag and remove unmapped rows
unmapped <- dat_mapped %>%
  filter(is.na(canonical_level)) %>%
  distinct(first_author, definition_contact_me)

if (nrow(unmapped) > 0) {
  warning(
    "No mapping found for the following rows — excluded from analysis:\n",
    paste(unmapped$first_author, unmapped$definition_contact_me,
          sep = " | ", collapse = "\n")
  )
  dat_mapped <- dat_mapped %>% filter(!is.na(canonical_level))
}

# =============================================================================
# 5. COMPUTE REFERENCE FRACTIONS FROM DISAGGREGATED STUDIES
# =============================================================================

# (a) Contact level fractions — from rows mapping to exactly one canonical level
# TODO: fix this -- this is currently biased becasue the encoding is completely wrong
ref_fractions_contact <- dat_mapped %>%
  filter(!is_combined) %>%
  group_by(canonical_level) %>%
  summarise(total_n = sum(denominator), .groups = "drop") %>%
  mutate(fraction_contact = total_n / sum(total_n))

cat("\nReference contact level fractions (from disaggregated studies):\n")
print(ref_fractions_contact)

# =============================================================================
# 6. APPLY CONTACT LEVEL FRACTIONS
# =============================================================================

# TODO: figure out the best way to enforce integer values for numerators and denominators
#   when we have fractional encoding
dat_fractions <- dat_mapped %>%
  left_join(
    ref_fractions_contact %>% select(canonical_level, fraction_contact),
    by = "canonical_level"
  ) %>%
  mutate(
    fraction_contact = if_else(!is_combined, 1, fraction_contact)
  ) %>%
  # Renormalise within study x raw definition in case some levels have no
  # reference data (guards against fractions not summing to 1)
  group_by(study_id, definition_contact_me) %>%
  mutate(fraction_contact = fraction_contact / sum(fraction_contact)) %>%
  ungroup() %>%
  mutate(
    numerator_adj   = numerator   * fraction_contact,
    denominator_adj = denominator * fraction_contact,
    uninfected_contacts_adj    = denominator_adj - numerator_adj
  )

# =============================================================================
# 8. BUILD DESIGN MATRIX
# =============================================================================

design_matrix <- dat_fractions %>%
  group_by(study_id, canonical_level) %>%
  summarise(
    numerator_adj     = sum(numerator_adj),
    denominator_adj   = sum(denominator_adj),
    uninfected_contacts_adj      = sum(uninfected_contacts_adj),
    fraction_contact  = sum(fraction_contact),
    .groups           = "drop"
  ) %>%
  pivot_wider(
    names_from  = canonical_level,
    values_from = fraction_contact,
    values_fill = 0,
    id_cols     = c(study_id, numerator_adj, denominator_adj,
                    uninfected_contacts_adj)
  ) %>%
  select(-any_of(reference_level)) %>%
  # Enforce canonical column order
  select(
    study_id, numerator_adj, denominator_adj, uninfected_contacts_adj,
    any_of(non_reference_levels)
  )

# TODO: figure out why the values of the design matrix aren't between 0 and 1?
cat("\nDesign matrix:\n")
print(design_matrix)

# Validate: check all expected columns are present
missing_cols <- setdiff(non_reference_levels, names(design_matrix))
if (length(missing_cols) > 0) {
  warning(
    "The following canonical levels have no observations and are absent ",
    "from the design matrix:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# =============================================================================
# 9. ANALYSIS 1: DISEASE PREVALENCE BY EXPOSURE LEVEL
# =============================================================================
# Model: binomial GLM with logit link (fixed effects)
#   logit(p_ij) = intercept + beta_household * household + sum(beta_k * X_k)
#
# Note: with few studies the random effect is not identifiable.
# The mixed-effects model is provided commented out for use as studies accumulate.

X <- design_matrix %>%
  select(any_of(non_reference_levels)) %>%
  as.matrix()

X <- cbind(intercept = 1, X)

y_num  <- design_matrix$numerator_adj
y_uninfected <- design_matrix$uninfected_contacts_adj

# TODO: need to fix the non-integers before we get to this stage of the analysis
fit_fixed <- glm(
  cbind(y_num, y_uninfected) ~ X - 1,
  family = binomial(link = "logit")
)

cat("\n=== Fixed-effects logistic regression ===\n")
print(summary(fit_fixed))

# Odds ratios and 95% CIs
or_fixed <- broom::tidy(fit_fixed, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(term = str_remove(term, "^X")) %>%
  rename(parameter = term, odds_ratio = estimate,
         ci_lower = conf.low, ci_upper = conf.high)

cat("\nOdds ratios (reference: No/minimal contact):\n")
print(or_fixed)

# Predicted SAR at each canonical level, holding household at reference (0)
pred_grid <- tibble(canonical_level = factor(canonical_levels, levels = canonical_levels)) %>%
  mutate(household_encoded = 0)

# Build prediction matrix matching X structure
X_pred <- matrix(0, nrow = length(canonical_levels),
                 ncol = ncol(X), dimnames = list(NULL, colnames(X)))
X_pred[, "intercept"] <- 1

for (i in seq_along(canonical_levels)) {
  lvl <- canonical_levels[i]
  col <- paste0("X", lvl)
  if (col %in% colnames(X_pred)) X_pred[i, col] <- 1
}

pred_link <- X_pred %*% coef(fit_fixed)
pred_se   <- sqrt(diag(X_pred %*% vcov(fit_fixed) %*% t(X_pred)))

pred_grid <- pred_grid %>%
  mutate(
    predicted_sar = plogis(pred_link),
    ci_lower      = plogis(pred_link - 1.96 * pred_se),
    ci_upper      = plogis(pred_link + 1.96 * pred_se)
  )

cat("\nPredicted SAR by canonical exposure level (household = non-household):\n")
print(pred_grid %>% select(canonical_level, predicted_sar, ci_lower, ci_upper))

# --- Mixed-effects model (uncomment when >1 study available) ----------------
# fit_mixed <- glmer(
#   cbind(numerator_adj, uninfected_contacts_adj) ~ household_encoded +
#     `Indirect contact only` +
#     `Direct physical contact — no fluids and no nursing` +
#     `Nursing care — no body fluids` +
#     `Body fluid contact` +
#     `Handled corpse` +
#     (1 | study_id),
#   data   = design_matrix,
#   family = binomial(link = "logit")
# )
# summary(fit_mixed)

# =============================================================================
# 10. ANALYSIS 2: DISTRIBUTION OF CONTACTS ACROSS EXPOSURE LEVELS
# =============================================================================

exposure_dist <- dat_fractions %>%
  group_by(study_id, canonical_level) %>%
  summarise(n = sum(denominator_adj), .groups = "drop") %>%
  group_by(study_id) %>%
  mutate(
    total_n = sum(n),
    prop    = n / total_n,
    canonical_level = factor(canonical_level, levels = canonical_levels)
  ) %>%
  ungroup()

cat("\n=== Exposure distribution within each study ===\n")
print(exposure_dist)

exposure_summary <- exposure_dist %>%
  group_by(canonical_level) %>%
  summarise(
    n_studies      = n(),
    mean_prop      = weighted.mean(prop, total_n),
    total_contacts = sum(n),
    .groups        = "drop"
  ) %>%
  mutate(overall_prop = total_contacts / sum(total_contacts))

cat("\n=== Summary exposure distribution across studies ===\n")
print(exposure_summary)

# =============================================================================
# 11. JOINT INTERPRETATION: POPULATION-ATTRIBUTABLE SAR
# =============================================================================

pop_prev <- pred_grid %>%
  left_join(
    exposure_summary %>% select(canonical_level, overall_prop),
    by = "canonical_level"
  ) %>%
  mutate(contribution = predicted_sar * overall_prop)

# TODO: fix all of these pooled SARs -- this isn't what we wanted to do
cat("\n=== Population-attributable prevalence ===\n")
print(pop_prev %>% select(canonical_level, predicted_sar, overall_prop, contribution))
cat(sprintf("\nEstimated overall SAR: %.3f\n", sum(pop_prev$contribution, na.rm = TRUE)))

# =============================================================================
# 12. SENSITIVITY ANALYSIS: EXCLUDE FRACTIONAL ROWS
# =============================================================================

dat_disaggregated <- dat_fractions %>% filter(!is_combined)

# TODO: fix the mapping -- there should be far less disaggregated
cat(sprintf(
  "\nSensitivity: %d fully disaggregated rows out of %d total mapped rows\n",
  nrow(dat_disaggregated), nrow(dat_fractions)
))

# Refit on disaggregated rows only — same pipeline as above
# (uncomment and adapt once multiple studies are present)

# =============================================================================
# 13. VISUALISATION
# =============================================================================

# Enforce factor order for plotting
# TODO: fix pred_grid
pred_grid <- pred_grid %>%
  mutate(canonical_level = factor(canonical_level, levels = rev(canonical_levels)))

p1 <- ggplot(pred_grid, aes(x = canonical_level, y = predicted_sar)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  coord_flip() +
  labs(
    title    = "Secondary attack rate by exposure level",
    subtitle = "Reference: No/minimal contact",
    x        = NULL,
    y        = "Secondary attack rate (predicted)"
  ) +
  theme_minimal()

ggsave("sar_by_exposure_level.png", p1, width = 8, height = 5, dpi = 150)

p2 <- exposure_dist %>%
  mutate(canonical_level = factor(canonical_level, levels = rev(canonical_levels))) %>%
  ggplot(aes(x = canonical_level, y = prop, fill = study_id)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() +
  labs(
    title = "Distribution of contacts across exposure levels",
    x     = NULL,
    y     = "Proportion of total contacts",
    fill  = "Study"
  ) +
  theme_minimal()

ggsave("exposure_distribution.png", p2, width = 8, height = 5, dpi = 150)

cat("\nDone. Plots saved to working directory.\n")
