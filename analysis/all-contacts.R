# all-contacts.R
# Author: Gina Cuomo-Dannenburg
# Date: 2026-09-30
# Purpose: Reformat the dataset in order to be able to perform the secondary attack
# rate inference
# Input: Data extracted from review of literature with information on SAR across all contacts
# Output: Figures and tables of the analysis of all contacts from eligible studies
# Data download date: 2026-08-25

# library
library(tidyverse)
require(readxl)
library(ggplot2)

# read data
data_all <- readRDS("data/data-derived/all_contacts.rds")

summary_all <- data_all %>%
  dplyr::arrange(label) %>%
  dplyr::select(label, location, country, year, species, outbreak, 
                study_design:study_population, household, numerator, denominator, sar_observed) %>%
  dplyr::mutate(binom::binom.confint(x = numerator,
                                     n = denominator,
                                     methods = "wilson")) %>%
# map the study designs onto something useful
  dplyr::mutate(study_design = case_when(
    study_design == "Prospective" ~ "Prospective",
    study_design == "Surveillance report" ~ "Prospective",
    study_design == "Concurrent and retrospective" ~ "Mixed",
    study_design %in% c("Retrospective","Cross-sectional","Serosurvey") ~ "Retrospective",
    study_design == "Vaccine ring trial (delayed vaccination arm)" ~ "Vaccine trial",
    TRUE ~ "Other"
  )) %>%
  dplyr::mutate( outbreak_clean = if_else(outbreak == "West Africa", "West Africa", "Other"))

summary_all$study_design_clean = factor(
  summary_all$study_design,
  levels = c(
    "Retrospective",
    "Mixed",
    "Prospective",
    "Vaccine trial",
    "Other"
  )
)

ggplot(summary_all, aes(x = label, y = mean*100, col = study_design_clean)) +
  theme_bw() + facet_grid(. ~ outbreak_clean, scales = "free_x") + 
  geom_point() +
  geom_errorbar(aes(ymin = lower*100, ymax = upper*100)) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) +
  labs(x = "Study", y = "SAR (%)", subtitle = "Secondary attack rate for all contacts") +
  guides(colour = guide_legend(title = "Contact ascertainment")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/sar_all_outbreak.png", dpi = 500, width = 30, height = 20, units = "cm")


# Make a summary table of key information for each study ------------------

table <- summary_all %>% 
  dplyr::mutate(sar = sar_observed * 100) %>%
  dplyr::arrange(year, country, label) %>%
  dplyr::mutate(household = if_else(household == 1, "Household contacts only", "All contacts")) %>%
  dplyr::select(label, country, year, outbreak, study_design, household, numerator, denominator, sar) 

