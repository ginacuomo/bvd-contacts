# # format-data.R
# Author: Gina Cuomo-Dannenburg
# Date: 2026-08-25
# Purpose: Reformat the dataset in order to be able to perform the secondary attack
# rate inference
# Input: Data extracted from review of literature 
# Interim: Unique exposures for inclusion by each study to generate the mapping
# Output: Dataset structured to enable the SAR analysis
# Data download date: 2026-08-25

# library
library(tidyverse)
require(readxl)
library(ggplot2)


# 1. Read and format the data ---------------------------------------------

# data import 
data <- readxl::read_excel("data/data-raw/sar_extraction.xlsx", sheet = "data-fomite")
labels <- readxl::read_excel("data/data-raw/sar_extraction.xlsx", sheet = "study-labels")

# now listing all of the unique exposures for inclusion
exposures_to_map <- data %>%
  dplyr::filter(include == TRUE) %>%
  dplyr::select(first_author:year, household, definition_contact, definition_contact_me:notes) %>%
  filter(definition_contact_me != "All")

length(unique(data$doi))
# Gayedu-Dennis I think should be excluded and some studies also with only all contacts available
length(unique(exposures_to_map$doi)) 

# now output this interim file to enable the mapping onto categories
write.csv(exposures_to_map, "data/data-raw/exposures-to-map.csv", row.names = FALSE)

# This was then manually edited and saved to the data/data-derived folder to read back in and merge with the 
# standard data set in order to proceed with the analysis
# TODO: figure out how to map the Dowell study. 
# we have % of cases with each risk factor but need to convert this into the actual numbers with 
# numerators and denominators

# define the exposure levels
# Classes 0 and 1 merged into a single reference "No direct physical contact"
# because both are sparsely reported and was propagating bias through the analysis

canonical_levels <- c(
  "No/minimal contact",
  "Indirect contact only",
  "Fomite exposure - no direct physical contact ",
  "Direct physical contact - no fluids and no nursing",
  "Nursing care - no body fluids",
  "Body fluid contact",
  "Handled corpse"
)

reference_level      <- canonical_levels[1]
non_reference_levels <- canonical_levels[-1]

exposures <- read.csv("data/data-derived/exposures-mapped.csv")

# reformat the data into the classes that we want/need
dat <- data %>%
  mutate(numerator = as.integer(numerator),
         denominator = as.integer(denominator),
         num_index = as.integer(num_index),
         year = as.character(year),
         study_id = paste(first_author, year_publication, location, sep = "_"),
         sar_observed = numerator / denominator,
         uninfected_contacts = denominator - numerator,
         household = as.integer(household)) %>%
  filter(!is.na(numerator),
         !is.na(denominator),
         denominator > 0,
         numerator <= denominator)

# now reformat the exposures table
# TODO: decide if we need to remap levels 0 and 1 -- separate for the moment
exclude_vars <- c("include", "imputed", "notes", 
                  "definition_contact", "location", "country")

exposures_long <- exposures %>%
  pivot_longer(cols = starts_with("exposure"),
               names_to  = "exposurenum",
               values_to = "exposure") %>%
  dplyr::select(-exposurenum) %>%
  dplyr::filter(!(is.na(exposure))) %>%
  dplyr::select(-any_of(exclude_vars))

# check if there are any duplicate mappings across the different studies
exposures_long %>%
  dplyr::group_by(first_author, doi, exposure) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n > 1) 
# Jezek has duplicate between household and non-household contact exposure for non-caregiving
exposures_long %>%
  dplyr::group_by(first_author, doi, exposure, household) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n > 1)

# reduce number of overlapping columns between the merging
intersect(names(dat), names(exposures_long))

# join the dataframes together
dat_joined <- left_join(dat, labels) %>%
  left_join(exposures_long) %>%
  dplyr::mutate(exposure = if_else(definition_contact_me == "All", "All", exposure))


# 2. Analyse SAR across all contacts --------------------------------------

# studies that didn't have SDBs i.e. included exposures from handling a corpse
# TODO: go to each individual study and see if they mention whether there were unsafe burials included
# non_sdb <- dat_joined %>%
#   dplyr::filter(exposure == "Handled corpse") %>%
#   pull(first_author)

all_contacts <- dat_joined %>%
  dplyr::filter(exposure == "All") %>%
  dplyr::filter(include == TRUE) %>% 
  group_by(label) %>% 
  dplyr::mutate(ci_lower = binom::binom.confint(x = numerator, n = denominator, methods = "wilson")$lower,
                ci_upper = binom::binom.confint(x = numerator, n = denominator, methods = "wilson")$upper) %>%
  dplyr::arrange(desc(sar_observed))

ggplot(all_contacts, aes(x = label, y = sar_observed*100, col = study_design)) +
  theme_bw() + geom_point() + ylim(c(0, 50)) +
  geom_errorbar(aes(ymin = ci_lower*100, ymax = ci_upper*100)) +
  labs(x = "First Author", y = "Observed SAR (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(subtitle = "All contacts") + 
  guides(col = guide_legend(title = "Study design")) + 
  theme(legend.direction = "horizontal", legend.position = "bottom")
ggsave("plots/all_contacts.png", dpi = 500, width = 25, height = 15, units = "cm")
  

# 3. Exposure disaggregation ----------------------------------------------

mapping_table <- dat_joined  %>%
  dplyr::filter(exposure != "All") %>%
  distinct(first_author, definition_contact_me, exposure, household, .keep_all = TRUE) %>%
  group_by(first_author, definition_contact_me) %>%
  mutate(n_canonical = n(),
         disaggregated  = n_canonical == 1,
         covered_levels = list(exposure)) %>%
  mutate(exposure = factor(exposure, levels = canonical_levels)) %>%
  ungroup()

ggplot(mapping_table, aes(y = label, x = exposure, fill = disaggregated)) + 
  theme_bw() + geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("plots/data-source.png", dpi = 500, width = 20, height = 15, units = "cm")

# looking closely at the Bower study because we actually had to aggregate this ourselves
bower <- dat_joined %>%
  dplyr::filter(first_author == "Bower") %>%
  group_by(definition_contact_me) %>%
  dplyr::mutate(ci_lower = binom::binom.confint(x = numerator, n = denominator, methods = "wilson")$lower,
                ci_upper = binom::binom.confint(x = numerator, n = denominator, methods = "wilson")$upper) %>%
  distinct(numerator, denominator, definition_contact_me, .keep_all = TRUE)
bower$definition_contact_me <- factor(bower$definition_contact_me,
                                      levels = c("All", 
                                                 "Handled corpse",
                                                 "Handled fluids", 
                                                 "Direct wet contact",        
                                                 "Direct dry contact",
                                                 "Direct contact combined",
                                                 "Indirect wet contact",        
                                                 "Indirect dry contact",
                                                 "Indirect contact combined",
                                                 "Minimal/no contact"))

ggplot(bower, aes(x = definition_contact_me, y = (numerator/denominator)*100, col = imputed)) + 
  theme_bw() + geom_point() + 
  geom_errorbar(aes(ymin = ci_lower*100, ymax = ci_upper*100)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "Exposure", y = "SAR (%)", subtitle = "Note: imputed exposures were to enable comparison with other papers") +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) #+ ylim(c(0,100))
ggsave("plots/bower.png", dpi = 500, width = 20, height = 15, units = "cm")
