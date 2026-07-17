# epireview.R
# 
# Author: Gina Cuomo-Dannenburg
# Date: 2026-06-08
# Purpose: to pull out the parameters from epireview for the work to support the
# analysis and understanding of contact tracing, secondary attack rates etc. for
# Bundibugyo

# Uncomment to install the package from GitHub:
# remotes::install_github("mrc-ide/epireview")

library(tidyverse)
library(epireview)
library(ggplot2)
library(DT)

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

attack <- df %>%
  dplyr::filter(parameter_type == "Attack rate") %>%
  dplyr::mutate(attack_median = if_else(parameter_unit == "Percentage (%)",
                                        parameter_value, parameter_value * 100),
                attack_lower = if_else(parameter_unit == "Percentage (%)",
                                       parameter_lower_bound, parameter_lower_bound * 100),
                attack_upper = if_else(parameter_unit == "Percentage (%)",
                                       parameter_upper_bound, parameter_upper_bound * 100)) %>%
  dplyr::select(c(article_label, ebola_species, 
                  parameter_value, parameter_unit, exponent,
                  parameter_lower_bound, parameter_upper_bound, parameter_uncertainty_type,
                  attack_median, attack_lower, attack_upper,
                  survey_start_date, survey_end_date, 
                  population_country, population_location,
                  population_sample_size, population_group, population_sample_type)) %>%
  group_by(article_label) %>%
  mutate(row_num = row_number(),
         n = n(),
         plot_label = if_else(n > 1, 
                              paste0(article_label, ".", row_num), 
                              article_label)) %>%
  ungroup() %>%
  arrange(plot_label)

label_map <- attack %>%
  select(plot_label, article_label) %>%
  deframe()  # named vector: names = plot_label, values = article_label

ggplot(attack, aes(x = attack_median, y = plot_label, col = ebola_species)) + 
  theme_bw() + 
  geom_point() + 
  geom_errorbar(aes(xmin = attack_lower, xmax = attack_upper)) +
  # fix the issues with y axis labels
  scale_y_discrete(breaks = names(label_map), labels = label_map) +
  labs(x = "Attack rate (%)", y = "Source") 

# make a table of risk factors 

risk <- df %>%
  dplyr::filter(parameter_type == "Risk factors") %>%
  dplyr::select(article_label, ebola_species,
                method_disaggregated_by:population_study_end_year, notes, qa_score, doi) %>%
  # now make risk long and split by different riskfactors
  tidyr::separate_rows(riskfactor_name, sep = ";") %>%
  tidyr::separate_rows(riskfactor_occupation, sep = ";") %>%
  dplyr::mutate(riskfactor_occupation = if_else(riskfactor_name == "Occupation",
                                                riskfactor_occupation, NA))


# want to prioritise studies which were:
# - WAE outbreak or later
# - ss and adjusted risk factor
# - high quality according to the quality assessment
# - close contacts or occupations

risk_prioritised <- risk %>%
  dplyr::filter(population_study_end_year > 2012) %>%
  # dplyr::filter(riskfactor_significant == "Significant") %>%
  dplyr::filter(riskfactor_adjusted == "Adjusted") %>% 
  dplyr::filter(qa_score >= 0.5) %>%
  dplyr::filter(riskfactor_name %in% c("Occupation","Close contact",
                                       "Household contact","Non-household contact",
                                       "Funeral","Social gathering"))





