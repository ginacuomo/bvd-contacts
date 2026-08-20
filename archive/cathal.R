# cathal.R

# Author: Gina Cuomo-Dannenburg
# Date: 2026-06-09
# Purpose: following on from our meeting, I worked on inferring the current epidemic size given
# a) seeding on 15 March across some estimates of doubling time
# b) using current cumulative deaths, inferring when they were exposed, and using case ascertainment
# to estimate both the seeding date and current Nt

# Note: this is very much a "back of the envelope" calculation and isn't designed to be a robust model

# libraries
library(tidyverse)
library(ggplot2)
library(epireview)
library(lubridate)

today <- ymd("2026-06-17")
t <- 94 # 94 days since 15 March 2026 (as of 17 June)
Nt <- function(t, t_d) {
  N <- 2^(t/t_d)
}

dt_inferred <- data.frame(doubling_time = seq(5, 10, by = 0.1))
dt_inferred$Nt <- Nt(t = t, t_d = dt_inferred$doubling_time)

ggplot(dt_inferred, aes(x = doubling_time, y = Nt)) + theme_bw() + 
  geom_point() + labs(x = "Doubling time (days)",
                      y = "N(t)",
                      subtitle = "Inferred epidemic size assuming seeding on 15 March 2026")
ggsave("archive/epidemic_size.png", dpi = 500)
ggplot(dt_inferred, aes(x = doubling_time, y = Nt)) + theme_bw() + 
  geom_point() + labs(x = "Doubling time (days)",
                      y = "N(t)",
                      subtitle = "Inferred epidemic size assuming seeding on 15 March 2026") +
  scale_y_log10()
ggsave("archive/epidemic_size_log.png", dpi = 500)
write.csv(dt_inferred, "archive/inferred_size.csv")

# now working from current deaths

# load epireview and pull out the estimates we want
articles <- epireview::load_epidata_raw("ebola", "article")
models   <- epireview::load_epidata_raw("ebola", "model")
params   <- epireview::load_epidata_raw("ebola", "parameter")

# Assign QA scores and attach to parameters
articles  <- epireview::assign_qa_score(articles = articles)$articles
qa_scores <- articles %>% dplyr::select(covidence_id, qa_score)

params <- params %>% left_join(qa_scores)

# Overview of species in the dataset
# unique(params$ebola_species)

# Join parameters with article-level metadata
df <- left_join(
  params,
  articles[, c("covidence_id", "first_author_surname", "year_publication",
               "article_label", "doi", "notes")],
  by = "covidence_id"
) %>%
  arrange(article_label, -year_publication)

# Identify numeric uncertainty/value columns
unc_cols <- names(df)[
  grepl("uncertainty|upper|lower|value", names(df)) &
    sapply(df, is.numeric)
]

df <- df %>%
  
  # Back-transform inverse parameters
  mutate(
    across(
      all_of(unc_cols),
      ~ if_else(
        inverse_param == TRUE & !is.na(.x) & .x != 0,
        1 / .x,
        .x
      )
    )
  ) %>%
  
  # Apply exponent correction
  mutate(
    across(
      all_of(unc_cols),
      ~ if_else(
        !is.na(exponent) & exponent != 0,
        .x * (10 ^ exponent),
        .x
      )
    )
  ) %>%
  
  # Round values and construct combined uncertainty/range strings
  mutate(
    parameter_value = round(parameter_value, 3),
    parameter_uncertainty_single_value = round(parameter_uncertainty_single_value, 3),
    comb_uncertainty = if_else(
      inverse_param == TRUE,
      paste0(
        round(parameter_uncertainty_upper_value, 3),
        " - ",
        round(parameter_uncertainty_lower_value, 3)
      ),
      paste0(
        round(parameter_uncertainty_lower_value, 3),
        " - ",
        round(parameter_uncertainty_upper_value, 3)
      )
    ),
    comb_range = if_else(
      inverse_param == TRUE,
      paste0(
        round(parameter_upper_bound, 3),
        " - ",
        round(parameter_lower_bound, 3)
      ),
      paste0(
        round(parameter_lower_bound, 3),
        " - ",
        round(parameter_upper_bound, 3)
      )
    )
  )

delays <- df %>% filter(parameter_type %in% c("Human delay - Exposure/Infection to Death",
                                    "Human delay - Exposure/Infection to Death/Recovery",
                                    "Human delay - Exposure/Infection to Symptom Onset/Fever",
                                    "Human delay - Symptom Onset/Fever to Death",
                                    "Human delay - Symptom Onset/Fever to Death in community",                                                     
                                    "Human delay - Symptom Onset/Fever to Death in hospital"))
# look at the delays briefly
epireview::forest_plot_delay_int(delays, ulim = 20, reorder_studies = TRUE) +
  facet_wrap(parameter_type ~ ., scales = "free")

onset_death <- delays %>%
  dplyr::filter(parameter_type == "Human delay - Symptom Onset/Fever to Death") %>%
  dplyr::filter(!(is.na(parameter_value))) %>%
  dplyr::filter(!(is.na(population_sample_size))) %>%
  dplyr::select(parameter_value, population_sample_size) %>%
  dplyr::mutate(weight = population_sample_size/sum(population_sample_size)) %>%
  dplyr::mutate(weighted_val = weight * parameter_value) %>%
  dplyr::summarise(average = weighted.mean(parameter_value, w = weight))
# > onset_death
# A tibble: 1 × 1
# average
# <dbl>
#   1    8.45

exposure_onset <- 11 # estimate from Lindblade et al 2015 for patients who reached ETU

exposure_death <- exposure_onset + onset_death

# data as of today
cum_deaths <- 196 # based on current numbers as of 17 June 2026 - confirmed in DRC

# make some assumptions about the proportions of deaths detected
# TODO: adjust based on linelist and what you expect to be true here
prop_confirmed <- 0.5 # proportion of deaths in an ETU who receive lab confirmation
prop_etu <- 0.7 # proportion of ebola patients who will die who attend an ETU -- again, adjust as you think best

# true deaths as of today - exposure_death days i.e. deaths today were infected at least X days ago
t_deaths <- t - exposure_death
num_deaths <- cum_deaths * 1/prop_etu * 1/prop_confirmed
# [1] 560

# assume some cfr for all ebola cases
cfr <- 0.33
num_cases <- num_deaths * 1/cfr # case estimates from 3 weeks ago
# [1] 1696.97

# Nt = 2^t/td => seeding date depends on doubling time
seeding <- data.frame(doubling_time = seq(5, 10, by = 0.1))
seeding$t <- seeding$doubling_time * (log(num_cases)/log(2))
seeding$days <- round(seeding$t, digits = 0)
seeding$date <- ymd(today) - as.numeric(round(exposure_death, digits = 0)) - days(seeding$days)

ggplot(seeding, aes(x = date, y = doubling_time)) + geom_point() + theme_bw() +
  labs(x = "Seeding date", y = "Doubling time (days)", 
       subtitle = "Date of index case given assumptions about doubling time, death ascertainment, cfr etc.")
ggsave("archive/seeding.png", dpi = 500)
