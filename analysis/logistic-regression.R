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
# 5. Logistic regression with study level effects (only standard logistic regression - omitted from sensitivity analysis)
# 6. Sensitivity analysis --> no fractional encoding 
# 7. Comparing the outputs 

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

# 5. Logistic regression with study level effects
# use lme4::glmer in order to have a study-level intercept
model_study <- lme4::glmer(data = pseudo_contacts, 
                outcome ~ indirect + fomites + direct + nursing + body_fluid + corpse + (1 | label),  
                family = "binomial")
summary(model_study)

# Predicted probabilities (SAR), including study-level random effects
predict(model_study,
        type = "response",
        re.form = NULL)

# Predicted probabilities for an average study (random effect = 0)
predict(model_study,
        type = "response",
        re.form = NA)

# Odds ratios for fixed effects
exp(lme4::fixef(model_study))

# 6. Sensitivity analysis --> no fractional encoding  ---------------------
# create a dataset where there is no fractional encoding to use in the below logistic regressions
pseudo_contacts_sensitivity <- pseudo_contacts %>%
  filter(
    rowSums(across(all_of(all_exposures), ~ . == 1)) == 1
  )

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

# 7. Comparing the outputs ---------------------------------------------------------------------


