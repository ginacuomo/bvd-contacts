# ebola_contact_analysis.R
# Author: Gina Cuomo-Dannenburg
# Date: 2026-06-26
# Purpose: Implementation of the statistic analysis plan to understand the secondary
# attack rate by contact definition from data from previous Ebola outbreaks
# Input: Data extracted from review of literature 
# Input: Mapping of exposures from data into a shared exposure hierarchy defined below
# Outcome 1: Secondary attack rate (SAR) by contact type
# Outcome 2: Contacts per index case by contact type
#
# Canonical exposure hierarchy (reference = "No direct physical contact"):
#   Reference: No direct physical contact    [Classes 0 + 1 merged]
#   Class 2: Direct physical contact - no fluids and no nursing
#   Class 3: Nursing care - no body fluids
#   Class 4: Body fluid contact
#   Class 5: Handled corpse

# =============================================================================

library(tidyverse)
library(readxl)
library(lme4)
library(broom.mixed)


# 1. Define levels --------------------------------------------------------
# Classes 0 and 1 merged into a single reference "No direct physical contact"
# because both are sparsely reported and was propagating bias through the analysis

canonical_levels <- c(
  "No direct physical contact",                         # Classes 0+1 — reference
  "Direct physical contact - no fluids and no nursing", # Class 2
  "Nursing care - no body fluids",                      # Class 3
  "Body fluid contact",                                 # Class 4
  "Handled corpse"                                      # Class 5
)

canonical_levels_sensitivity <- c(
  "No/minimal contact",
  "Indirect contact only",
  "Direct physical contact - no fluids and no nursing",
  "Nursing care - no body fluids",
  "Body fluid contact",
  "Handled corpse"
)

reference_level      <- canonical_levels[1]
non_reference_levels <- canonical_levels[-1]


# 2. Import data and clean --------------------------

raw <- read_excel("data/sar_extraction.xlsx", sheet = "data_me")
mapping_raw <- read_excel("data/sar_extraction.xlsx", sheet = "exposure_mapping")

dat <- raw %>%
  mutate(numerator = as.integer(numerator),
         denominator = as.integer(denominator),
         num_index = as.integer(num_index),
         year = as.character(year),
         study_id = paste(first_author, year_publication, location, sep = "_"),
         sar_observed = numerator / denominator,
         uninfected_contacts = denominator - numerator) %>%
  filter(!is.na(numerator),
         !is.na(denominator),
         denominator > 0,
         numerator <= denominator,
         include == TRUE)

cat("Studies and contact definitions in dataset:\n")
dat %>%
  group_by(first_author) %>%
  summarise(n = length(unique(definition_contact_me))) %>%
  print()

# 3. Mapping table ---------------------------
mapping_long <- mapping_raw %>%
  pivot_longer(cols      = starts_with("exposure"),
               names_to  = "exposurenum",
               values_to = "exposure") %>%
  dplyr::select(-exposurenum) %>%       # explicit namespace avoids MASS::select conflict
  filter(!is.na(exposure)) %>%
  # Collapse Classes 0 and 1 into the merged reference category
  mutate(exposure = case_when(
    exposure %in% c("No/minimal contact", "Indirect contact only") ~
      "No direct physical contact",
    TRUE ~ exposure))

mapping_table <- mapping_long %>%
  # After collapsing 0+1, deduplicate in case both mapped to the same definition
  distinct(first_author, definition_contact_me, exposure, household, .keep_all = TRUE) %>%
  group_by(first_author, definition_contact_me) %>%
  mutate(n_canonical    = n(),
         is_combined    = n_canonical > 1,
         covered_levels = list(exposure)) %>%
  ungroup()

# Validate: every canonical level must be in the defined hierarchy
unrecognised <- mapping_table %>%
  filter(!exposure %in% canonical_levels) %>%
  distinct(exposure)

# Check that all canonical levels in the mapping table exist
if (nrow(unrecognised) > 0) {
  stop(
    "Unrecognised canonical levels in mapping table — check spelling:\n",
    paste(unrecognised$exposure, collapse = "\n")
  )
} else {
  cat("\nMapping table validated: all canonical levels recognised.\n")
}


# 4. Map exposures onto canonical levels ----------------------------------

dat_mapped <- dat %>%
  left_join(
    mapping_table %>%
      dplyr::select(first_author, definition_contact_me, exposure,
                    is_combined, n_canonical, covered_levels),
    by           = c("first_author", "definition_contact_me"),
    relationship = "many-to-many"
  )

# Flag and remove unmapped rows
unmapped <- dat_mapped %>%
  filter(is.na(exposure)) %>%
  distinct(first_author, definition_contact_me)

if (nrow(unmapped) > 0) {
  warning(
    "No mapping found for the following rows — excluded from analysis:\n",
    paste(unmapped$first_author, unmapped$definition_contact_me,
          sep = " | ", collapse = "\n"))
  dat_mapped <- dat_mapped %>% filter(!is.na(exposure))
}


# 5. Perform fractional encoding ------------------------------------------
# Strategy:
#   (a) Identify fully disaggregated studies (no combined rows anywhere).
#   (b) Within each such study compute the proportion of contacts in each
#       canonical level — proportions sum to 1 over observed levels only.
#   (c) Average across studies (weighted by total contacts) separately per
#       level, using only studies that observed that level. This avoids
#       suppressing rare levels (e.g. Handled corpse) due to structural
#       absence in other studies.
#   (d) Renormalise the global reference distribution to sum to 1.
#   (e) For combined rows, assign each covered level its reference proportion
#       renormalised over only the levels that definition covers.

# --- 5a. Identify fully disaggregated studies --------------------------------

fully_disaggregated_studies <- dat_mapped %>%
  group_by(study_id) %>%
  summarise(any_combined = any(is_combined), .groups = "drop") %>%
  filter(!any_combined) %>%
  pull(study_id)

# --- 5b. Study-level proportions (over observed levels only) -----------------

study_level_props <- dat_mapped %>%
  filter(study_id %in% fully_disaggregated_studies) %>%
  group_by(study_id, exposure) %>%
  summarise(n_exposure = sum(denominator), .groups = "drop") %>%
  group_by(study_id) %>%
  mutate(
    total_n           = sum(n_exposure),
    study_prop        = n_exposure / total_n,
    n_levels_observed = n()
  ) %>%
  ungroup()

cat("\nStudy-level proportions by exposure (disaggregated studies only):\n")
print(study_level_props)

# --- 5c-d. Global reference distribution -------------------------------------
# weighted.mean() per level uses only studies that observed that level,
# because absent levels simply have no rows. Final renormalisation ensures
# the distribution sums to 1 despite unequal study coverage per level.

ref_distribution <- study_level_props %>%
  group_by(exposure) %>%
  summarise(
    ref_prop  = weighted.mean(study_prop, total_n),
    n_studies = n(),
    .groups   = "drop"
  ) %>%
  complete(
    exposure  = canonical_levels,
    fill      = list(ref_prop = 0, n_studies = 0)
  ) %>%
  mutate(
    ref_prop = ref_prop / sum(ref_prop),   # renormalise to sum to 1
    exposure = factor(exposure, levels = canonical_levels)
  ) %>%
  arrange(exposure)

print(ref_distribution)

# Warn if any level has no reference data
zero_ref <- ref_distribution %>% filter(ref_prop == 0)
if (nrow(zero_ref) > 0) {
  warning(
    "The following canonical levels have no reference data:\n",
    paste(zero_ref$exposure, collapse = "\n"),
    "\nFractional encoding for definitions covering these levels will ",
    "renormalise over remaining covered levels."
  )
}

# --- 5e. Apply fractional encoding -------------------------------------------

dat_fractions <- dat_mapped %>%
  left_join(
    ref_distribution %>% dplyr::select(exposure, ref_prop),
    by = "exposure"
  ) %>%
  mutate(
    fraction_raw = if_else(!is_combined, 1, ref_prop)
  ) %>%
  group_by(study_id, definition_contact_me) %>%
  mutate(fraction_sum     = sum(fraction_raw),
         fraction_contact = if_else(
           fraction_sum > 0,
           fraction_raw / fraction_sum,
           1 / n()  )) %>%     # fallback: equal split if all ref_props are zero
  ungroup() %>%
  #   TODO: figure out a more elegant way to handle the fact that the numerators
  #   and denominators are non integers for passing into the Binomial logistic
  mutate(numerator_adj = round(numerator * fraction_contact, digits = 0),
         denominator_adj = round(denominator * fraction_contact, digits = 0),
         uninfected_contacts_adj = round(denominator_adj - numerator_adj, digits = 0))

# Diagnostics
dat_fractions %>%
  filter(is_combined) %>%
  dplyr::select(study_id, definition_contact_me, exposure,
                ref_prop, fraction_contact, denominator, denominator_adj) %>%
  print()

fraction_check <- dat_fractions %>%
  group_by(study_id, definition_contact_me) %>%
  summarise(fraction_sum = sum(fraction_contact), .groups = "drop") %>%
  filter(abs(fraction_sum - 1) > 1e-6)

if (nrow(fraction_check) > 0) {
  warning(
    "Fractions do not sum to 1 for:\n",
    paste(fraction_check$study_id, fraction_check$definition_contact_me,
          sep = " | ", collapse = "\n")
  )
} else {
  cat("\nFraction check passed: all fractions sum to 1 within study.\n")
}


# 6. Build the model data -------------------------------------------------

model_data <- dat_fractions %>%
  group_by(study_id, exposure) %>%
  summarise(
    numerator_adj = sum(numerator_adj),
    denominator_adj  = sum(denominator_adj),
    uninfected_contacts_adj = sum(uninfected_contacts_adj),
    num_index = first(num_index),
    .groups = "drop"
  ) %>%
  mutate(exposure  = factor(exposure, levels = canonical_levels),
         exposure  = relevel(exposure, ref = reference_level)) %>%
  filter(denominator_adj > 0)

# Wide design matrix for inspection only (not passed to model)
design_matrix_wide <- dat_fractions %>%
  group_by(study_id, exposure) %>%
  summarise(numerator_adj = sum(numerator_adj),
            denominator_adj = sum(denominator_adj),
            uninfected_contacts_adj = sum(uninfected_contacts_adj),
            fraction_contact = min(sum(fraction_contact), 1),
            .groups = "drop") %>%
  pivot_wider(names_from   = exposure,
              values_from  = fraction_contact,
              values_fill  = 0,
              id_cols      = c(study_id, numerator_adj, denominator_adj,
                               uninfected_contacts_adj)) %>%
  dplyr::select(-any_of(reference_level)) %>%
  dplyr::select(study_id, numerator_adj, denominator_adj, uninfected_contacts_adj,
                any_of(non_reference_levels))

print(design_matrix_wide)

missing_cols <- setdiff(non_reference_levels, names(design_matrix_wide))
if (length(missing_cols) > 0) {
  warning(
    "The following canonical levels have no observations:\n",
    paste(missing_cols, collapse = "\n")
  )
}


# 7. Logistic regression: Secondary attack rate by exposure level ------------------------------

fit_fixed <- glm(
  cbind(numerator_adj, uninfected_contacts_adj) ~ exposure,
  data   = model_data,
  family = binomial(link = "logit")
)

cat("\n=== Fixed-effects logistic regression (SAR) ===\n")
print(summary(fit_fixed))

# Add in an OR calculation because it's more easily interpretable
# Note: the reference category is "No direct contact"
or_fixed <- broom::tidy(fit_fixed, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(odds_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high)

cat("\nOdds ratios (reference:", reference_level, "):\n")
print(or_fixed)

# --- Mixed-effects model (uncomment when >1 study available) ----------------
# fit_mixed <- glmer(
#   cbind(numerator_adj, uninfected_contacts_adj) ~
#     exposure + household_encoded + (1 | study_id),
#   data   = model_data,
#   family = binomial(link = "logit")
# )
# summary(fit_mixed)


# 8. Predict secondary attack rate by exposure level ----------------------
pred_grid <- tibble(
  exposure          = factor(canonical_levels, levels = canonical_levels),
  household_encoded = 0
)

pred_link <- predict(fit_fixed, newdata = pred_grid, type = "link", se.fit = TRUE)

pred_grid <- pred_grid %>%
  mutate(predicted_sar = plogis(pred_link$fit),
         ci_lower      = plogis(pred_link$fit - 1.96 * pred_link$se.fit),
         ci_upper      = plogis(pred_link$fit + 1.96 * pred_link$se.fit))

cat("\nPredicted SAR by exposure level:\n")
print(pred_grid %>% dplyr::select(exposure, predicted_sar, ci_lower, ci_upper))

# visualise this compared to the data
ggplot() +
  geom_errorbar(data = pred_grid,
                aes(x = exposure, ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_point(data = pred_grid,aes(x = exposure, y = predicted_sar),
             size  = 3,shape = 18   # diamond to distinguish from observed
  ) + geom_point(data = model_data,
                 aes(x = exposure, y = numerator_adj / denominator_adj, col = study_id),
                 size  = 2, shape = 1) +    # open circle
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    x     = NULL,
    y     = "Secondary attack rate",
    col   = "Study",
    title = "Observed and predicted SAR by exposure level"
  ) +
  theme_bw() +
  theme(legend.position  = "bottom", legend.direction = "horizontal")
ggsave("plots/sar_by_exposure.png", dpi = 500, width = 30, height = 20, units = "cm")


# 9. Contacts per exposure ------------------------------------------------

# Uses raw (unadjusted) denominators from dat_mapped 
# fractional adjustment would distort the count structure since num_index is study-level.
# Model: negative binomial with log(num_index) offset 

contact_data <- dat_mapped %>%
  group_by(study_id, exposure) %>%
  summarise(n_contacts = sum(denominator),
            num_index  = first(num_index),
            .groups    = "drop") %>%
  filter(!is.na(num_index), num_index > 0) %>%
  mutate(exposure  = factor(exposure, levels = canonical_levels),
         exposure  = relevel(exposure, ref = reference_level),
         log_index = log(num_index))

# --- 9a. Overall contacts per index case ------------------------------------

contact_data_overall <- contact_data %>%
  group_by(study_id, num_index, log_index) %>%
  summarise(n_contacts = sum(n_contacts), .groups = "drop")

fit_nb_overall <- MASS::glm.nb(n_contacts ~ offset(log_index),
                               data = contact_data_overall)

cat("\n=== Overall contacts per index case (negative binomial) ===\n")
print(summary(fit_nb_overall))

nb_overall_tidy <- broom::tidy(fit_nb_overall, conf.int = TRUE, exponentiate = TRUE)

# Extract and report overdispersion
cat(sprintf("Dispersion parameter theta: %.2f\n", fit_nb_overall$theta))
cat(sprintf("95%% CI: %.2f - %.2f\n",
            fit_nb_overall$theta - 1.96 * fit_nb_overall$SE.theta,
            fit_nb_overall$theta + 1.96 * fit_nb_overall$SE.theta))

overall_rate <- tibble(
  estimate = nb_overall_tidy$estimate[nb_overall_tidy$term == "(Intercept)"],
  ci_lower = nb_overall_tidy$conf.low[nb_overall_tidy$term == "(Intercept)"],
  ci_upper = nb_overall_tidy$conf.high[nb_overall_tidy$term == "(Intercept)"],
  theta    = fit_nb_overall$theta
)
print(overall_rate)

# Descriptive cross-check
contact_data_overall %>%
  mutate(contacts_per_index = n_contacts / num_index) %>%
  summarise(mean_contacts   = weighted.mean(contacts_per_index, num_index),
            median_contacts = median(contacts_per_index),
            min_contacts    = min(contacts_per_index),
            max_contacts    = max(contacts_per_index)) %>%
  print()

# --- 9b. Contacts per index case by exposure level --------------------------

# Poisson dispersion check
fit_pois_check <- glm(
  n_contacts ~ exposure + offset(log_index),
  data   = contact_data,
  family = poisson
)

dispersion_ratio <- fit_pois_check$deviance / fit_pois_check$df.residual

fit_nb_exposure <- MASS::glm.nb(
  n_contacts ~ exposure + offset(log_index),
  data = contact_data
)

print(summary(fit_nb_exposure))

# Rate ratios
rr_contacts <- broom::tidy(fit_nb_exposure, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(rate_ratio = estimate, ci_lower = conf.low, ci_upper = conf.high) %>%
  mutate(term = str_remove(term, "^exposure"))

print(rr_contacts)

# Predicted contacts per index case at each exposure level
pred_contacts <- tibble(exposure  = factor(canonical_levels, levels = canonical_levels),
                        log_index = 0)   # exp(0) = 1 index case

pred_link_contacts <- predict(fit_nb_exposure,
                              newdata = pred_contacts,
                              type    = "link",
                              se.fit  = TRUE)

pred_contacts <- pred_contacts %>%
  mutate(contacts_per_index = exp(pred_link_contacts$fit),
         ci_lower           = exp(pred_link_contacts$fit - 1.96 * pred_link_contacts$se.fit),
         ci_upper           = exp(pred_link_contacts$fit + 1.96 * pred_link_contacts$se.fit))

cat("\nPredicted contacts per index case by exposure level:\n")
print(pred_contacts %>% dplyr::select(exposure, contacts_per_index, ci_lower, ci_upper))

