# petal-analysis.R
# Author: Gina Cuomo-Dannenburg
# Date: 2026-10-20
# Purpose: Perform the analysis that Charlie's group requires for their analysis
# Input: Data extracted from review of literature with information on SAR across all contacts, formatted and cleaned already
# Output: Figures and tables of the analysis

# Requests from Charlie's group:
# Odds ratios using their exposure mappings rather than ours
# Household SAR based on the enhanced contact tracing studies

# Stages of the analysis
# 1. Read in the data
# 2. Reformat the data with the new exposure mappings
# 3. Estimate the overall household SAR 
# 4. Estimate the ORs for the PETAL exposure mappings

# library
library(tidyverse)
library(ggplot2)
library(lme4)

theme_set(theme_bw(base_size = 14))

# 1. Read in the data -------------------------------

dat_all <- readRDS("data/data-derived/all_contacts.rds")
pseudo_contacts <- readRDS("data/data-derived/data_logistic.rds")
dat_fractional <- readRDS("data/data-derived/sar_exposure.rds")
reference_distribution <- readRDS("data/data-derived/reference_distribution.rds")

# 2. Reformat the data with the new exposure mappings ---------------------
# For PETAL, minimal + indirect + fomites are combined
all_exposures_petal <- c(
  "not_direct",
  "direct",
  "nursing",
  "body_fluid",
  "corpse"
)

canonical_labels_petal <- c(
  not_direct = "No direct physical contact",
  direct = "Direct physical contact - no fluids and no nursing",
  nursing = "Nursing care - no body fluids",
  body_fluid = "Body fluid contact",
  corpse = "Handled corpse"
)

pseudo_contacts_petal <- pseudo_contacts %>%
  mutate(
    not_direct = minimal + indirect + fomites
  ) %>%
  dplyr::select(-minimal, -indirect, -fomites) %>%
  dplyr::filter(ct_type == "Research Study")

# 3. Estimate the overall household SAR 
# Filter for the studies with enhanced contact tracing and only household contacts

dat_household <- dat_all %>%
  dplyr::filter(ct_type == "Research Study") %>%
  dplyr::filter(household == 1) %>%
  # omit the study in infants because this is not representative
  dplyr::filter(label != "Bower, 2016a") %>%
  dplyr::bind_cols(
    binom::binom.confint(
      x = .$numerator,
      n = .$denominator,
      methods = "wilson"
    ) %>%
      dplyr::select(lower, upper)
  ) %>%
  # mean household size
  dplyr::mutate(household_size = if_else(first_author == "Dowell", 7.4, 10.0))

household_model <- glm(
  cbind(numerator, denominator - numerator) ~ 1,
  data = dat_household,
  family = binomial()
)

estimate <- coef(household_model)[1]
se <- sqrt(vcov(household_model)[1, 1])

household_sar <- tibble(
  description = "Household SAR",
  sar = plogis(estimate),
  lower = plogis(estimate - 1.96 * se),
  upper = plogis(estimate + 1.96 * se)
)

# plot this against the data

ggplot(dat_household, aes(x = label, y = sar_observed)) +
  geom_rect(
    data = household_sar,
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = lower, ymax = upper
    ),
    inherit.aes = FALSE,
    fill = "#0072B2",
    alpha = 0.15
  ) +
  geom_hline(
    yintercept = household_sar$sar,
    linetype = 2,
    colour = "#0072B2"
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2
  ) +
  geom_point(size = 2.5) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(
    x = "Study",
    y = "Observed household SAR",
    subtitle = "Inferred mean household size: Bower - 10.0; Dowell - 7.4"
  ) + theme_bw()
ggsave("plots/household_sar-petal.png", dpi = 300, units = "cm", height = 20, width = 30)

# 3b. Household SAR with household size as a covariate

household_model_size <- glm(
  cbind(numerator, denominator - numerator) ~ household_size,
  data = dat_household,
  family = binomial()
)

summary(household_model_size)

# Odds ratio per additional household member
exp(coef(household_model_size)["household_size"])

# Predict household SAR across a range of household sizes

pred_household <- tibble(
  household_size = 4:20
)

pred <- predict(
  household_model_size,
  newdata = pred_household,
  type = "link",
  se.fit = TRUE
)

pred_household <- pred_household %>%
  mutate(
    sar = plogis(pred$fit),
    lower = plogis(pred$fit - 1.96 * pred$se.fit),
    upper = plogis(pred$fit + 1.96 * pred$se.fit)
  )

# Plot predicted SAR and 95% CI

ggplot(pred_household, aes(x = household_size, y = sar)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "#0072B2",
    alpha = 0.2
  ) +
  geom_line(colour = "#0072B2", linewidth = 1) +
  scale_x_continuous(breaks = seq(4, 20, by = 2)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Household size",
    y = "Predicted household SAR"
  ) +
  theme_bw(base_size = 14)

# 4. Estimate the ORs for the PETAL exposure mappings ------------------------
# logistic regression with direct contact as the reference category
# Note: we only have the enhanced contact tracing data included here

# Logistic regression with direct contact as reference

model_or <- glm(
  outcome ~ not_direct + nursing + body_fluid + corpse,
  data = pseudo_contacts_petal,
  family = binomial()
)

summary(model_or)

# Extract odds ratios and 95% confidence intervals

or_estimates <- broom::tidy(
  model_or,
  conf.int = TRUE,
  conf.method = "Wald",
  exponentiate = TRUE
) %>%
  dplyr::filter(term != "(Intercept)") %>%
  dplyr::rename(
    exposure = term,
    odds_ratio = estimate,
    ci_lower = conf.low,
    ci_upper = conf.high
  ) %>%
  bind_rows(
    tibble(
      exposure = "direct",
      odds_ratio = 1,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
  ) %>%
  mutate(
    exposure = factor(
      exposure,
      levels = all_exposures_petal,
      labels = unname(canonical_labels_petal[all_exposures_petal])
    )
  ) %>%
  arrange(exposure)

print(or_estimates)

ggplot(or_estimates,
    aes(x = exposure,
        y = odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.5),
    width = 0.2
  ) +
  scale_y_log10() +
  labs(
    x = "Maximum exposure category",
    y = "Odds ratio",
    colour = "Contact tracing type",
    subtitle = "Reference: Direct physical contact; \nOnly enhanced contact tracing studies"
  )  
ggsave("plots/or-petal.png", dpi = 300, units = "cm", height = 20, width = 30)

# Logistic regression with direct contact as reference with household size as covariate

pseudo_contacts_petal_hh <- pseudo_contacts_petal %>%
  dplyr::filter(household == 1) %>%
  # omit the study in infants because this is not representative
  dplyr::filter(label != "Bower, 2016a") %>%
  dplyr::mutate(household_size = if_else(label == "Dowell, 1999", 7.4, 10.0))

model_or_hh <- glm(
  outcome ~ not_direct + nursing + body_fluid + corpse + household_size,
  data = pseudo_contacts_petal_hh,
  family = binomial()
)

summary(model_or_hh)

# Extract odds ratios and 95% confidence intervals for exposure categories

or_estimates_hh <- broom::tidy(
  model_or_hh,
  conf.int = TRUE,
  conf.method = "Wald",
  exponentiate = TRUE
) %>%
  dplyr::filter(term != "(Intercept)", term != "household_size") %>%
  dplyr::rename(
    exposure = term,
    odds_ratio = estimate,
    ci_lower = conf.low,
    ci_upper = conf.high
  ) %>%
  bind_rows(
    tibble(
      exposure = "direct",
      odds_ratio = 1,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
  ) %>%
  mutate(
    exposure = factor(
      exposure,
      levels = all_exposures_petal,
      labels = unname(canonical_labels_petal[all_exposures_petal])
    )
  ) %>%
  arrange(exposure)

print(or_estimates_hh)

# Extract odds ratio and 95% confidence interval for household size separately

household_size_estimate <- broom::tidy(
  model_or_hh,
  conf.int = TRUE,
  conf.method = "Wald",
  exponentiate = TRUE
) %>%
  dplyr::filter(term == "household_size") %>%
  dplyr::rename(
    odds_ratio = estimate,
    ci_lower = conf.low,
    ci_upper = conf.high
  )

print(household_size_estimate)

# save the outputs so it's easy for Charlie
saveRDS(household_sar, "data/data-out/household-sar.rds")
saveRDS(or_estimates, "data/data-out/or_estimates.rds")