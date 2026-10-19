# logistic-regression.R
# Author: Gina Cuomo-Dannenburg
# Date: 2026-10-16
# Purpose: Perform the logistic regression section of the analysis, including a sensitivity analysis
# Input: Data extracted from review of literature with information on SAR across all contacts, formatted and cleaned already
# Output: Figures and tables of the analysis of SAR by exposure for use in subsequent analyses/outputs
# Note: used code from Chris to help create initial draft of code base

# Stages of the analysis:
# 1. Read in the data, filter columns, etc.
# 2. Define logit/logistic functions
# 3. Standard logistic regression
# 4. Logistic regression with data collection level effects
# 5. Logistic regression with split by data collection level -- 2 separate analyses
# 6. Sensitivity analysis --> no fractional encoding 
# 7. Comparing the outputs 
# 8. Making figures of the outputs requested

# library
library(tidyverse)
library(ggplot2)
library(lme4)


# 1. Read in the data, filter columns, etc. -------------------------------

dat_all <- readRDS("data/data-derived/all_contacts.rds")
pseudo_contacts <- readRDS("data/data-derived/data_logistic.rds")
dat_fractional <- readRDS("data/data-derived/sar_exposure.rds")
reference_distribution <- readRDS("data/data-derived/reference_distribution.rds")

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

# 2. Define logit/logistic functions ------------------------------------

# logit(p) = log(p / (1 - p)), probabilties to log odds
logit <- function(...) stats::qlogis(...)
# logistic(x) = 1 / (1 + exp(-x)), log odds to probabilties
logistic <- function(...) stats::plogis(...)

# 3. Standard logistic regression -----------------------------------------

model <- glm(data = pseudo_contacts, 
             outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse,  
             family = "binomial")
summary(model)

predict.glm(model,
            type = "response", # predictions on the response (outcome) scale, not the predictors scale
            se.fit = TRUE)
exp(coef(model))

# 4. Logistic regression with data collection level effects

model_ct <- glm(data = pseudo_contacts, 
                outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse + ct_type,  
                family = "binomial")
summary(model_ct)

predict.glm(model_ct,
            type = "response", # predictions on the response (outcome) scale, not the predictors scale
            se.fit = TRUE)
exp(coef(model_ct))

# TODO: split the data into two based on ct_level and reanalyse and then compare
# 5. Logistic regression with split by data collection level -- 2 separate analyses

# split data
pseudo_contacts_standard <- pseudo_contacts %>%
  dplyr::filter(ct_type ==  "Real-time Contact Tracing")
pseudo_contacts_enhanced <- pseudo_contacts %>%
  dplyr::filter(ct_type ==  "Research Study")

# standard contact tracing
model_standard <- glm(data = pseudo_contacts_standard, 
                      outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse,  
                      family = "binomial")
summary(model_standard)

predict.glm(model_standard,
            type = "response", # predictions on the response (outcome) scale, not the predictors scale
            se.fit = TRUE)
exp(coef(model_standard))

# enhanced contact tracing
model_enhanced <- glm(data = pseudo_contacts_enhanced, 
                      outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse,  
                      family = "binomial")
summary(model_enhanced)

predict.glm(model_enhanced,
            type = "response", # predictions on the response (outcome) scale, not the predictors scale
            se.fit = TRUE)
exp(coef(model_enhanced))


# 6. Sensitivity analysis --> no fractional encoding  ---------------------
# create a dataset where there is no fractional encoding to use in the below logistic regressions
pseudo_contacts_sensitivity <- pseudo_contacts %>%
  filter(rowSums(across(all_of(all_exposures), ~ . == 1)) == 1)

# 6a. Standard logistic regression -- sensitivity analysis ----------------------------------------

model_sens <- glm(data = pseudo_contacts_sensitivity, 
                  outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse,  
                  family = "binomial")
summary(model_sens)

predict.glm(model_sens,
            type = "response", # predictions on the response (outcome) scale, not the predictors scale
            se.fit = TRUE)
exp(coef(model_sens))

# 6b. Logistic regression with data collection level effects -- sensitivity analysis

model_ct_sens <- glm(data = pseudo_contacts_sensitivity, 
                     outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse + ct_type,  
                     family = "binomial")
summary(model_ct_sens)

predict.glm(model_ct_sens,
            type = "response", # predictions on the response (outcome) scale, not the predictors scale
            se.fit = TRUE)
exp(coef(model_ct_sens))

# 6c. Logistic regression with split data by contact tracing intensity -- sensitivity analysis
pseudo_contacts_sensitivity_standard <- pseudo_contacts %>%
  filter(rowSums(across(all_of(all_exposures), ~ . == 1)) == 1) %>%
  dplyr::filter(ct_type ==  "Real-time Contact Tracing")
pseudo_contacts_sensitivity_enhanced <- pseudo_contacts %>%
  filter(rowSums(across(all_of(all_exposures), ~ . == 1)) == 1) %>%
  dplyr::filter(ct_type ==  "Research Study")

model_standard_sens <- glm(data = pseudo_contacts_sensitivity_standard, 
                           outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse,  
                           family = "binomial")
summary(model_standard_sens)

predict.glm(model_standard_sens,
            type = "response", # predictions on the response (outcome) scale, not the predictors scale
            se.fit = TRUE)
exp(coef(model_standard_sens))

model_enhanced_sens <- glm(data = pseudo_contacts_sensitivity_enhanced, 
                           outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse ,  
                           family = "binomial")
summary(model_enhanced_sens)

predict.glm(model_enhanced_sens,
            type = "response", # predictions on the response (outcome) scale, not the predictors scale
            se.fit = TRUE)
exp(coef(model_enhanced_sens))

# 7. Comparing the outputs ---------------------------------------------------------------------
# we want to compare when all of the data are combined vs split by contact tracing intensity

# List of fitted models
models <- list(
  "All data - Standard" = model,
  "Operational contact tracing - Standard" = model_standard,
  "Enhanced contact tracing - Standard" = model_enhanced,
  "All data - Sensitivity" = model_sens,
  "Operational contact tracing - Sensitivity" = model_standard_sens,
  "Enhanced contact tracing - Sensitivity" = model_enhanced_sens
)

# Create prediction grid
pred_grid <- tibble(exposure = all_exposures)

for (lvl in all_exposures) {
  pred_grid[[lvl]] <- as.integer(pred_grid$exposure == lvl)
}

# Generate predictions for each model
pred_comparison <- purrr::imap_dfr(models, function(fit, model_name) {
  
  # Exclude random effects for mixed models
  pred <- if (inherits(fit, "merMod")) {
    predict(fit, newdata = pred_grid, type = "link", re.form = NA)
  } else {
    predict(fit, newdata = pred_grid, type = "link", se.fit = TRUE)
  }
  
  tibble(
    model = model_name,
    exposure = pred_grid$exposure,
    predicted_sar = plogis(if (is.list(pred)) pred$fit else pred),
    ci_lower = if (is.list(pred)) plogis(pred$fit - 1.96 * pred$se.fit) else NA_real_,
    ci_upper = if (is.list(pred)) plogis(pred$fit + 1.96 * pred$se.fit) else NA_real_
  )
})

# Order and label exposure categories
pred_comparison <- pred_comparison %>%
  # omit the coefficients that are NA because this doesn't reflect an estimated SAR
  mutate(
    across(
      c(predicted_sar, ci_lower, ci_upper),
      ~ if_else(
        model == "Operational contact tracing - Sensitivity" &
          exposure %in% c("indirect", "corpse"),
        NA_real_,
        .x
      )
    )
  ) %>%
  mutate(
    exposure = factor(
      exposure,
      levels = all_exposures,
      labels = unname(canonical_labels[all_exposures])
    )
  ) %>% 
  dplyr::mutate(
    data_type = stringr::str_trim(stringr::str_extract(model, "^[^-]+")),
    analysis_type = case_when(
      str_detect(model, "Standard") ~ "Primary analysis",
      str_detect(model, "Sensitivity") ~ "Sensitivity analysis"
    )
  )

ggplot(pred_comparison,
       aes(x = exposure, y = predicted_sar, colour = model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.5),
    width = 0.2
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Maximum exposure category",
    y = "Predicted secondary attack rate",
    colour = "Model"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      vjust = 1
    )
  )

ggsave("plots/compare-sar.png", units = "cm", width = 30, height = 20, dpi = 300)

ggplot(pred_comparison,
       aes(x = exposure, y = predicted_sar, colour = model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.5),
    width = 0.2
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Maximum exposure category",
    y = "Predicted secondary attack rate",
    colour = "Model"
  ) +
  scale_y_log10()

# 8. Making figures of the outputs requested ---------------------------------

# set theme
theme_set(theme_bw(base_size = 14) +
            theme(
              axis.text.x = element_text(
                angle = 30,
                hjust = 1,
                vjust = 1
              )
            ))

analysis_colours <- c(
  "Primary analysis" = "#3B6FB6",
  "Sensitivity analysis" = "#CC79A7"
)

# Christophe suggested that we present four figures to show this analysis
# plot 1. with the aggregate SAR for the three studies with routine contact tracing,
#   and with all of the studies with research contact tracing.

# use the SAR across all contacts but split by the type of contact tracing and perform these separately
# Operational contact tracing
model_operational <- glm(
  cbind(numerator, denominator - numerator) ~ 1,
  data = dat_all %>% dplyr::filter(ct_type == "Real-time Contact Tracing"),
  family = binomial()
)

# Investigational contact tracing
model_investigational <- glm(
  cbind(numerator, denominator - numerator) ~ 1,
  data = dat_all %>% dplyr::filter(ct_type == "Research Study"),
  family = binomial()
)

# Estimated SARs
plogis(coef(model_operational))
plogis(coef(model_investigational))

# Extract SAR estimates and 95% CIs

models_aggregate <- list(
  "Operational contact tracing" = model_operational,
  "Enhanced contact tracing" = model_investigational
)

sar_aggregate <- purrr::imap_dfr(models_aggregate, function(fit, model_name) {
  
  estimate <- coef(fit)[1]
  se <- sqrt(vcov(fit)[1, 1])
  
  tibble(
    ct_type = model_name,
    sar = plogis(estimate),
    lower = plogis(estimate - 1.96 * se),
    upper = plogis(estimate + 1.96 * se)
  )
})

sar_aggregate$ct_type <- factor(sar_aggregate$ct_type,
                                levels = c("Operational contact tracing",
                                           "Enhanced contact tracing"))

# Plot SAR with 95% CIs

ggplot(sar_aggregate, aes(x = ct_type, y = sar, colour = ct_type)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", sar * 100)),
    vjust = -1,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = stratified_colours) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Contact tracing type",
    y = "Estimated secondary attack rate"
  ) +
  theme(legend.position = "none")

ggsave("plots/plot_1.png", units = "cm", width = 30, height = 20, dpi = 300)


# plot 2. SAR by risk for the research contact tracing studies.
ggplot(filter(pred_comparison, data_type == "Enhanced contact tracing"),
       aes(x = exposure, y = predicted_sar, colour = analysis_type)) +
  geom_point(position = position_dodge(width = 0.5),
             width = 0.2) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.5),
    width = 0.2
  ) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0,1)) +
  scale_colour_manual(values = analysis_colours) +
  labs(
    x = "Maximum exposure category",
    y = "Predicted secondary attack rate",
    colour = "Method",
    subtitle = "Enhanced contact tracing"
  ) 
ggsave("plots/plot_2.png", units = "cm", width = 30, height = 20, dpi = 300)

# plot 3. SAR by risk for the routine contact tracing studies
ggplot(filter(pred_comparison, data_type == "Operational contact tracing"),
       aes(x = exposure, y = predicted_sar, colour = analysis_type)) +
  geom_point(position = position_dodge(width = 0.5),
             width = 0.2) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.5),
    width = 0.2
  ) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0,1)) +
  scale_colour_manual(values = analysis_colours) +
  labs(
    x = "Maximum exposure category",
    y = "Predicted secondary attack rate",
    colour = "Method",
    subtitle = "Operational contact tracing"
  ) 
ggsave("plots/plot_3.png", units = "cm", width = 30, height = 20, dpi = 300)

# plot 4. as above but with odds ratios for a fully stratified analysis.

# Generate odds ratios and 95% CIs for each model

or_comparison <- purrr::imap_dfr(models, function(fit, model_name) {
  
  # Extract fixed effects for mixed models, coefficients for GLMs
  if (inherits(fit, "merMod")) {
    estimates <- lme4::fixef(fit)
    se <- sqrt(diag(vcov(fit)))
  } else {
    estimates <- coef(fit)
    se <- sqrt(diag(vcov(fit)))
  }
  
  tibble(
    model = model_name,
    exposure = names(estimates),
    odds_ratio = exp(estimates),
    ci_lower = exp(estimates - 1.96 * se),
    ci_upper = exp(estimates + 1.96 * se)
  )
})

# Order and label exposure categories

or_comparison <- or_comparison %>%
  # Remove intercept (minimal is the reference)
  dplyr::filter(exposure != "(Intercept)") %>%
  # Keep only exposure coefficients
  dplyr::filter(exposure %in% all_exposures) %>%
  mutate(
    exposure = factor(
      exposure,
      levels = all_exposures,
      labels = unname(canonical_labels[all_exposures])
    )
  ) %>%
  dplyr::mutate(
    data_type = stringr::str_trim(stringr::str_extract(model, "^[^-]+")),
    analysis_type = case_when(
      str_detect(model, "Standard") ~ "Primary analysis",
      str_detect(model, "Sensitivity") ~ "Sensitivity analysis"
    )
  )

# stratified analysis colour palette
stratified_colours <- c(
  "Operational contact tracing" = "#0072B2",
  "Enhanced contact tracing" = "#D55E00"
)

or_comparison$data_type <- factor(or_comparison$data_type,
                                  levels = c(
                                    "Operational contact tracing",
                                    "Enhanced contact tracing"
                                  ))
# make the figure now
ggplot(filter(or_comparison, data_type != "All data" & 
                analysis_type == "Primary analysis"),
       aes(x = exposure, y = odds_ratio, colour = data_type)) +
  geom_point(position = position_dodge(width = 0.5),
             width = 0.2) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.5),
    width = 0.2
  ) +
  scale_colour_manual(values = stratified_colours) +
  labs(
    x = "Maximum exposure category",
    y = "Predicted odds ratio",
    colour = "Stratified analysis",
    subtitle = "Odds ratio - reference: no/minimal contact"
  ) 
ggsave("plots/plot_4.png", units = "cm", width = 30, height = 20, dpi = 300)

ggplot(filter(or_comparison, data_type != "All data" & 
                analysis_type == "Primary analysis"),
       aes(x = exposure, y = odds_ratio, colour = data_type)) +
  geom_point(position = position_dodge(width = 0.5),
             width = 0.2) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.5),
    width = 0.2
  ) +
  scale_colour_manual(values = stratified_colours) +
  labs(
    x = "Maximum exposure category",
    y = "Predicted odds ratio",
    colour = "Stratified analysis",
    subtitle = "Odds ratio - reference: no/minimal contact"
  ) + coord_cartesian(ylim = c(0, 50))
ggsave("plots/plot_4b.png", units = "cm", width = 30, height = 20, dpi = 300)
