# # format-data.R
# Author: Gina Cuomo-Dannenburg
# Date: 2026-09-28
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
data <- readxl::read_excel("data/data-raw/data-extracted.xlsx", sheet = "complete-data")
labels <- readxl::read_excel("data/data-raw/data-extracted.xlsx", sheet = "study-labels")

# check how all contacts are encoded so that we exclude these from the mapping
sort(unique(data$definition_contact_me))[1:6]

# now listing all of the unique exposures for inclusion
dat <- data %>%
  dplyr::mutate(numerator = as.integer(numerator),
                denominator = as.integer(denominator),
                num_index = as.integer(num_index),
                year = as.character(year),
                sar_observed = numerator / denominator,
                uninfected_contacts = denominator - numerator) %>%
  # additional checks but this has manually been checked already
  dplyr::filter(!is.na(numerator),
                !is.na(denominator),
                denominator > 0,
                numerator <= denominator) %>%
  dplyr::filter(include == TRUE)

# check that we have all rows for inclusion here 
nrow(dat) == sum(data$include)

exposures_to_map <- dat %>%
  dplyr::filter(include == TRUE) %>%
  dplyr::select(first_author:year, household, definition_contact, definition_contact_me:notes) %>%
  filter(!(definition_contact_me  %in% c("All", "All contacts")))

length(unique(data$doi)) == length(unique(exposures_to_map$doi)) 

# now output this interim file to enable the mapping onto categories
write.csv(exposures_to_map, "data/data-raw/exposures-to-map.csv", row.names = FALSE)

# This was then manually edited and saved to the data/data-derived folder to read back in and merge with the 
# standard data set in order to proceed with the analysis
exposures <- read.csv("data/data-derived/exposures-mapped.csv")

# check we aren't missing any exposures
dat$definition_contact_me[which(!(dat$definition_contact_me %in% exposures$definition_contact_me))]
# only missing the All / Al contacts contacts which is correct => something is going wrong in the joining further down

# join the datasets together
dat_joined <- dat %>%
  left_join(
    exposures %>% select(first_author, year, doi, definition_contact_me, exposure1:exposure5),
    by = c("first_author", "year", "doi", "definition_contact_me"))

# now reformat the data to be long
dat_joined <- dat_joined %>% 
  pivot_longer(exposure1:exposure5, names_to = "exposurenum", values_to = "exposure") %>%
  dplyr::filter(exposure != "") %>% 
  dplyr::filter(!(is.na(exposure))) %>%
  dplyr::select(-exposurenum)


# 2. Map exposures onto the canonical levels correctly --------------------

# define the canonical levels
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

# reformat it so that it matches the canonical levels
# make the lookup table
unique(dat_joined$exposure)

lookup <- tibble(exposure = c("handled corpse", "body fluid exposure", "nursing care", "direct contact", "fomites", "indirect", "minimal"),
                 levels = rev(canonical_levels))

dat_joined <- dat_joined %>%
  mutate(exposure = lookup$levels[match(exposure, lookup$exposure)])

# add the study labels to the dataset
dat_joined <- left_join(dat_joined, labels) %>%
  relocate(label) %>%
  arrange(label)


# 3. Save data output -----------------------------------------------------

saveRDS(dat_joined, "data/data-derived/dat_joined.rds")

