# sar-exposure.R
# Author: Gina Cuomo-Dannenburg
# Date: 2026-09-30
# Purpose: SAR by exposure from the pseudo-contacts generated in format-data.R
# Input: Data extracted from review of literature with information on SAR across all contacts formatted and cleaned already
# Output: Figures and tables of the analysis of SAR by exposure for use in subsequent analyses/outputs
# Note: used code from Chris to help create initial draft of code base

library(ggplot2)
library(tidyverse)

# read data
dat <- readRDS("data/data-derived/data_logistic.rds")
dat_fractional <- readRDS("data/data-derived/sar_exposure.rds")
ref_distribution <- readRDS("data/data-derived/reference_distribution.rds")


# Exposure names used for modelling and simulation

all_exposures <- c(
  "minimal",
  "indirect",
  "fomites",
  "direct",
  "nursing",
  "body_fluid",
  "corpse"
)

canonical_labels <- c(
  minimal = "No/minimal contact",
  indirect = "Indirect contact only",
  fomites = "Fomite exposure - no direct physical contact",
  direct = "Direct physical contact - no fluids and no nursing",
  nursing = "Nursing care - no body fluids",
  body_fluid = "Body fluid contact",
  corpse = "Handled corpse"
)


# 1. Define logistic ------------------------------------------------------

# logit(p) = log(p / (1 - p)), probabilties to log odds
logit <- function(...) stats::qlogis(...)
# logistic(x) = 1 / (1 + exp(-x)), log odds to probabilties
logistic <- function(...) stats::plogis(...)


# 2. Fit model ------------------------------------------------------------
# currently excluding Swanson et al. because the SAR for direct contact is wildly low 
model <- glm(data = dplyr::filter(dat, label != "Swanson, 2018"), 
             outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse, 
             family = "binomial")
summary(model)

predict.glm(model,
            type = "response", # predictions on the response (outcome) scale, not the predictors scale
            se.fit = TRUE)

exp(coef(model))
# (Intercept)    indirect     fomites      direct     nursing  body_fluid      corpse 
# 0.02067247  2.28573514  3.20986786  4.76380895 11.90567623 25.60688693 24.76833720 

# # TODO: debug why direct contact has a lower OR than the other exposures -- makes no sense
# dat %>% dplyr::filter(direct > 0) %>% group_by(label) %>% reframe(cases = sum(outcome), contacts = denominator) %>% distinct() %>% mutate(sar = (cases/contacts)*100) %>%
#   pull(sar) %>% mean()
# dat %>% dplyr::filter(direct > 0) %>% group_by(label) %>% reframe(cases = sum(outcome), contacts = denominator) %>% 
#   distinct() %>% 
#   summarise(num = sum(cases), denom = sum(contacts)) %>%
#   mutate(sar = num/denom)

# Sensitivity model using only rows with a single known exposure

pseudo_contacts_sensitivity <- dat %>%
  dplyr::filter(label != "Swanson, 2018") %>% 
  filter(
    rowSums(across(all_of(all_exposures), ~ . == 1)) == 1
  )

model_sensitivity <- glm(
  outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse,
  data = pseudo_contacts_sensitivity,
  family = "binomial"
)

# -------------------------------------------------------------------------
# Predicted SAR by exposure: primary and sensitivity analyses
# -------------------------------------------------------------------------

pred_grid <- tibble(
  exposure = all_exposures
)

# Create the design matrix for prediction

for (lvl in all_exposures) {
  pred_grid[[lvl]] <- as.integer(pred_grid$exposure == lvl)
}


# Primary model predictions

pred_link <- predict(
  model,
  newdata = pred_grid,
  type = "link",
  se.fit = TRUE
)

pred_grid <- pred_grid %>%
  mutate(
    predicted_sar = plogis(pred_link$fit),
    ci_lower = plogis(
      pred_link$fit - 1.96 * pred_link$se.fit
    ),
    ci_upper = plogis(
      pred_link$fit + 1.96 * pred_link$se.fit
    )
  )


# Sensitivity model predictions

pred_link_sensitivity <- predict(
  model_sensitivity,
  newdata = pred_grid,
  type = "link",
  se.fit = TRUE
)

pred_grid <- pred_grid %>%
  mutate(
    sensitivity_sar = plogis(pred_link_sensitivity$fit),
    sensitivity_ci_lower = plogis(
      pred_link_sensitivity$fit -
        1.96 * pred_link_sensitivity$se.fit
    ),
    sensitivity_ci_upper = plogis(
      pred_link_sensitivity$fit +
        1.96 * pred_link_sensitivity$se.fit
    )
  )


# Print predictions

cat("\nPrimary analysis: predicted SAR by exposure level:\n")

print(
  pred_grid %>%
    dplyr::select(
      exposure,
      predicted_sar,
      ci_lower,
      ci_upper
    )
)

cat("\nSensitivity analysis: predicted SAR by exposure level:\n")

print(
  pred_grid %>%
    dplyr::select(
      exposure,
      sensitivity_sar,
      sensitivity_ci_lower,
      sensitivity_ci_upper
    )
)


# -------------------------------------------------------------------------
# Primary figure
# -------------------------------------------------------------------------

pred_plot <- pred_grid %>%
  mutate(
    exposure = factor(
      exposure,
      levels = all_exposures,
      labels = unname(canonical_labels[all_exposures])
    )
  )

ggplot(
  pred_plot,
  aes(
    x = exposure,
    y = predicted_sar
  )
) +
  theme_bw() +
  geom_errorbar(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    width = 0.2
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        predicted_sar * 100
      )
    ),
    vjust = -5
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    x = "Maximum exposure category",
    y = "Predicted secondary attack rate"
  )


# -------------------------------------------------------------------------
# Sensitivity figure
# -------------------------------------------------------------------------

ggplot(
  pred_plot,
  aes(
    x = exposure,
    y = sensitivity_sar
  )
) +
  theme_bw() +
  geom_errorbar(
    aes(
      ymin = sensitivity_ci_lower,
      ymax = sensitivity_ci_upper
    ),
    width = 0.2
  ) +
  geom_point(size = 3) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    x = "Maximum exposure category",
    y = "Predicted secondary attack rate"
  )


# -------------------------------------------------------------------------
# Comparison figure
# -------------------------------------------------------------------------

comparison_plot <- pred_grid %>%
  mutate(
    exposure = factor(
      exposure,
      levels = all_exposures,
      labels = canonical_labels
    )
  ) %>%
  dplyr::select(
    exposure,
    predicted_sar,
    ci_lower,
    ci_upper,
    sensitivity_sar,
    sensitivity_ci_lower,
    sensitivity_ci_upper
  ) %>%
  pivot_longer(
    cols = c(
      predicted_sar,
      sensitivity_sar
    ),
    names_to = "analysis",
    values_to = "sar"
  ) %>%
  mutate(
    lower = if_else(
      analysis == "predicted_sar",
      ci_lower,
      sensitivity_ci_lower
    ),
    upper = if_else(
      analysis == "predicted_sar",
      ci_upper,
      sensitivity_ci_upper
    ),
    analysis = recode(
      analysis,
      predicted_sar = "Primary",
      sensitivity_sar = "Sensitivity"
    )
  )

ggplot(
  comparison_plot,
  aes(
    x = exposure,
    y = sar,
    colour = analysis
  )
) +
  theme_bw() +
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper
    ),
    position = position_dodge(width = 0.3),
    width = 0.2
  ) +
  geom_point(
    position = position_dodge(width = 0.3),
    size = 3
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  labs(
    x = "Maximum exposure category",
    y = "Predicted secondary attack rate",
    colour = "Analysis"
  )


# -------------------------------------------------------------------------
# Simulation parameters
# -------------------------------------------------------------------------

set.seed(42)
n_index <- 100

# from Noé's data
mean_contacts <-  13.85
theta <-  0.47

# Exposure proportions from reference distribution (must sum to 1)

exposure_props <- ref_distribution %>%
  arrange(factor(exposure, levels = all_exposures)) %>%
  pull(proportion)
names(exposure_props) <- all_exposures

# Predicted SAR per exposure level

sar_by_exposure <- pred_grid %>%
  arrange(factor(exposure, levels = all_exposures)) %>%
  pull(sensitivity_sar)
names(sar_by_exposure) <- all_exposures

# checks
stopifnot(
  identical(names(exposure_props), all_exposures),
  abs(sum(exposure_props) - 1) < 1e-8,
  identical(names(sar_by_exposure), all_exposures)
)

# --- Simulate ----------------------------------------------------------------

line_list <- map_dfr(1:n_index, function(i) {
  # Step 1: draw total contacts for this index case from negative binomial
  total_contacts <- rnbinom(n = 1, mu = mean_contacts, size = theta)
  # Edge case: if total_contacts = 0, return a single row with zeros
  if (total_contacts == 0) {
    return(tibble(
      index_case               = sprintf("%03d", i),
      total_contacts           = 0L,
      exposure                 = all_exposures,
      total_contacts_exposure  = 0L,
      infected_contacts        = 0L,
      uninfected_contacts      = 0L
    ))
  }
  
  # Step 2: allocate contacts across exposure levels
  # Draw from multinomial to get integer counts that sum to total_contacts
  # Work from highest to lowest risk (all_exposures is already ordered
  # lowest to highest, so reverse for allocation then reorder at the end)
  contacts_per_exposure <- as.integer(
    rmultinom(n = 1, size = total_contacts, prob = exposure_props)
  )
  names(contacts_per_exposure) <- all_exposures
  
  # Step 3: for each exposure level, draw infected contacts from binomial
  # using the predicted SAR for that level
  
  infected_per_exposure <- map_int(all_exposures, function(lvl) {
    n   <- contacts_per_exposure[lvl]
    sar <- sar_by_exposure[lvl]
    if (n == 0) return(0L)
    rbinom(n = 1, size = n, prob = sar)
  })
  
  names(infected_per_exposure) <- all_exposures
  
  # Step 4: assemble into long format, one row per exposure level
  tibble(
    index_case              = sprintf("%03d", i),
    total_contacts          = total_contacts,
    exposure                = all_exposures,
    total_contacts_exposure = as.integer(contacts_per_exposure),
    infected_contacts       = as.integer(infected_per_exposure),
    uninfected_contacts     = as.integer(contacts_per_exposure - infected_per_exposure)
  )
})

# Order by index case and exposure level (highest to lowest risk)
line_list <- line_list %>%
  mutate(exposure = factor(exposure, levels = rev(all_exposures))) %>%
  arrange(index_case, exposure) %>%
  mutate(exposure = as.character(exposure))

# check that things add as they should
# Check 1: total_contacts_exposure sums to total_contacts within each index case
check_totals <- line_list %>%
  group_by(index_case, total_contacts) %>%
  summarise(sum_exposure = sum(total_contacts_exposure), .groups = "drop") %>%
  filter(total_contacts != sum_exposure)

if (nrow(check_totals) > 0) {
  warning("total_contacts_exposure does not sum to total_contacts for:\n",
          paste(check_totals$index_case, collapse = ", "))
} else {
  cat("\nCheck 1 passed: total_contacts_exposure sums to total_contacts for all index cases.\n")
}

# Check 2: infected + uninfected = total_contacts_exposure
check_infected <- line_list %>%
  filter(infected_contacts + uninfected_contacts != total_contacts_exposure)

if (nrow(check_infected) > 0) {
  warning("infected + uninfected != total_contacts_exposure for ",
          nrow(check_infected), " rows.")
} else {
  cat("Check 2 passed: infected + uninfected = total_contacts_exposure for all rows.\n")
}

# check if the linelist proportions of contacts in each category makes sense
line_list$exposure <- factor(line_list$exposure, levels = all_exposures)
line_list %>% dplyr::group_by(exposure) %>%
  dplyr::summarise(total_contacts = sum(total_contacts_exposure)) %>%
  dplyr::mutate(all_contacts = sum(total_contacts), 
                prop = total_contacts/all_contacts)
exposure_props
# passes visual check that this is sensible

# 11. 2x2 summary table for each definition -------------------------------
high_risk_def <- all_exposures[7]

summary_stats <- function(high_risk_def) {
  
  summary_table <- line_list %>%
    mutate(contact_class = if_else(exposure %in% high_risk_def, "high_risk", "non_high_risk")) %>%
    group_by(contact_class) %>%
    summarise(cases = sum(infected_contacts),
              uninfected = sum(uninfected_contacts),
              total_contacts = sum(total_contacts_exposure)) %>%
    bind_rows(summarise(., across(where(is.numeric), sum, na.rm = TRUE),
                        across(where(is.character), ~"total")))
  
  sar <- summary_table$cases[summary_table$contact_class == "high_risk"] /
    summary_table$total_contacts[summary_table$contact_class == "high_risk"]
  
  sensitivity <- summary_table$cases[summary_table$contact_class == "high_risk"] /
    summary_table$cases[summary_table$contact_class == "total"]
  
  # There are no non-high-risk contacts when all exposure categories are high risk
  if (length(summary_table$uninfected[summary_table$contact_class == "non_high_risk"]) == 0) {
    specificity <- NA_real_
    missed_cases <- 0
    prop_missed <- 0
  } else {
    specificity <- summary_table$uninfected[summary_table$contact_class == "non_high_risk"] /
      summary_table$uninfected[summary_table$contact_class == "total"]
    
    missed_cases <- summary_table$cases[summary_table$contact_class == "non_high_risk"]
    prop_missed <- missed_cases/summary_table$cases[summary_table$contact_class == "total"]
  }
  
  found_cases <- summary_table$cases[summary_table$contact_class == "high_risk"]
  
  number_fu <- summary_table$total_contacts[summary_table$contact_class == "high_risk"]
  
  diagnostics <- data.frame(
    sar = sar,
    sensitivity = sensitivity,
    specificity = specificity,
    found_cases = found_cases,
    missed_cases = missed_cases,
    prop_missed = prop_missed,
    number_fu = number_fu
  )
  
  return(list(summary_table, diagnostics))
}

# Define the four cumulative high-risk thresholds (highest to lowest risk)
# These are the possible definitions of a high-risk contact to be followed up
threshold_definitions <- list(
  "class7_handled_corpse" = all_exposures[7],
  "class6to7_body_fluid_and_above" = all_exposures[6:7],
  "class5to7_nursing_and_above" = all_exposures[5:7],
  "class4to7_direct_physical_contact_and_above" = all_exposures[4:7],
  "class3to7_fomites_and_above" = all_exposures[3:7],
  "class2to7_indirect_contact_and_above" = all_exposures[2:7],
  "class1to7_any_possible_contact" = all_exposures[1:7]
)

results <- imap(threshold_definitions, function(high_risk_def, def_name) {
  
  res <- summary_stats(high_risk_def)
  
  list(
    definition = def_name,
    high_risk = paste(high_risk_def, collapse = "; "),
    summary = res[[1]] %>% mutate(definition = def_name, .before = 1),
    diagnostics = res[[2]] %>% mutate(definition = def_name, .before = 1)
  )
})

# combine so they're easy to export
# summary tables stacked
all_summaries <- map_dfr(results, "summary")

# All diagnostics stacked
all_diagnostics <- map_dfr(results, "diagnostics")

# Wide diagnostics table — one row per threshold, easy to read at a glance
diagnostics_wide <- all_diagnostics %>%
  mutate(
    high_risk_classes = map_chr(threshold_definitions[definition],
                                ~ paste(.x, collapse = "; "))
  ) %>%
  dplyr::select(
    definition, high_risk_classes,
    sar, sensitivity, specificity, found_cases, missed_cases, prop_missed, number_fu
  ) %>%
  mutate(across(c(sar, sensitivity, specificity, prop_missed),
                ~ round(.x, 3))) %>%
  mutate(
    definition = factor(
      definition,
      levels = c(
        "class7_handled_corpse",
        "class6to7_body_fluid_and_above",
        "class5to7_nursing_and_above",
        "class4to7_direct_physical_contact_and_above",
        "class3to7_fomites_and_above",
        "class2to7_indirect_contact_and_above",
        "class1to7_any_possible_contact"
      ),
      labels = c(
        "Handled corpse",
        "Body fluid+",
        "Nursing+",
        "Direct contact+",
        "Fomite+",
        "Indirect+",
        "All contacts"
      )
    )
  )

# export the files
write_csv(all_summaries, "contact_risk_summary_tables.csv")
write_csv(all_diagnostics, "contact_risk_diagnostics_long.csv")
write_csv(diagnostics_wide, "contact_risk_diagnostics_wide.csv")


# Make some figures of these results --------------------------------------

ggplot(diagnostics_wide,
       aes(x = definition)) +
  geom_col(
    aes(
      y = number_fu,
      fill = "Contacts followed-up"
    ),
    width = 0.7
  ) +
  geom_col(
    aes(
      y = found_cases,
      fill = "Cases found"
    ),
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "Cases found" = "#d95f02",
      "Contacts followed-up" = "grey80"
    ),
    name = NULL
  ) +
  labs(
    x = "High-risk contact definition",
    y = "Contacts followed up"
  ) +
  theme_bw()
ggsave("plots/contacts_followup.png", dpi = 500, width = 30, height = 20, units = "cm")

plot_diag <- diagnostics_wide %>%
  select(
    definition,
    sensitivity,
    specificity
  ) %>%
  pivot_longer(
    cols = c(
      sensitivity,
      specificity
    ),
    names_to = "diagnostic",
    values_to = "value"
  ) %>%
  # fix specificity for all contacts
  dplyr::mutate(value = if_else(is.na(value), 0, value))

ggplot(
  plot_diag,
  aes(
    x = definition,
    y = value,
    colour = diagnostic,
    group = diagnostic
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  scale_colour_manual(
    values = c(
      sensitivity = "#d95f02",
      specificity = "#1b9e77"
    ),
    labels = c(
      sensitivity = "Sensitivity",
      specificity = "Specificity"
    ),
    name = NULL
  ) +
  labs(
    x = "High-risk contact definition",
    y = "Diagnostic performance"
  ) +
  theme_bw()
ggsave("plots/diagnostics.png", dpi = 500, width = 30, height = 20, units = "cm")

plot_cases <- diagnostics_wide %>%
  select(
    definition,
    found_cases,
    missed_cases
  ) %>%
  pivot_longer(
    cols = c(found_cases, missed_cases),
    names_to = "status",
    values_to = "cases"
  ) %>%
  mutate(
    status = factor(
      status,
      levels = c("missed_cases", "found_cases"),
      labels = c("Cases missed", "Cases found")
    )
  )

ggplot(
  plot_cases,
  aes(
    x = definition,
    y = cases,
    fill = status
  )
) +
  geom_col(
    position = "fill",
    width = 0.7
  ) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  scale_fill_manual(
    values = c(
      "Cases found" = "#d95f02",
      "Cases missed" = "grey80"
    ),
    name = NULL
  ) +
  labs(
    x = "High-risk contact definition",
    y = "Proportion of cases"
  ) +
  theme_bw()
ggsave("plots/prop_cases.png", dpi = 500, width = 30, height = 20, units = "cm")

plot_studies <- dat_fractional %>%
  dplyr::filter(n_exposures_mapped == 1) %>%
  mutate(
    exposure1 = factor(
      exposure1,
      levels = c(
        "minimal",
        "indirect",
        "fomites",
        "direct",
        "nursing",
        "body_fluid",
        "corpse"
      ),
      labels = c(
        "No/minimal contact",
        "Indirect contact only",
        "Fomite exposure",
        "Direct physical contact",
        "Nursing care",
        "Body fluid contact",
        "Handled corpse"
      )
    )
  )

ggplot(
  plot_studies,
  aes(
    x = exposure1,
    y = sar_observed,
    colour = label
  )
) +
  geom_point(
    position = position_jitter(width = 0.1, height = 0),
    size = 2.5
  ) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  labs(
    x = "Maximum exposure category",
    y = "Observed secondary attack rate",
    colour = "Study"
  ) +
  theme_bw()
ggsave("plots/exposure_study.png", dpi = 500, width = 30, height = 20, units = "cm")
