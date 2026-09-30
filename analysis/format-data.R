# format-data.R
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
  dplyr::select(first_author:year, outbreak, household, definition_contact, definition_contact_me:notes) %>%
  filter(!(definition_contact_me  %in% c("All", "All contacts")))

length(unique(data$doi)) == length(unique(exposures_to_map$doi)) 

# now output this interim file to enable the mapping onto categories
write.csv(exposures_to_map, "data/data-raw/exposures-to-map.csv", row.names = FALSE)

# This was then manually edited and saved to the data/data-derived folder to read back in and merge with the 
# standard data set in order to proceed with the analysis
exposures <- read.csv("data/data-derived/exposures-mapped.csv")

# adding in some checks
# Check that the exposure mapping is unique for each original row
mapping_check <- exposures %>%
  count(first_author, year, doi, household, definition_contact_me) %>%
  filter(n > 1)

stopifnot(nrow(mapping_check) == 0)

# check we aren't missing any exposures
dat$definition_contact_me[which(!(dat$definition_contact_me %in% exposures$definition_contact_me))]
# only missing the All / Al contacts contacts which is correct => something is going wrong in the joining further down

# join the datasets together
dat_joined <- dat %>%
  left_join(
    exposures %>% select(first_author, year, doi, definition_contact_me, exposure1:exposure5),
    by = c("first_author", "year", "doi", "definition_contact_me")) %>%
  left_join(labels) %>%
  relocate(label)

# Check that all non-"All" rows have a mapped exposure
missing_mapping <- dat_joined %>%
  filter(!definition_contact_me %in% c("All", "All contacts")) %>%
  filter(if_all(exposure1:exposure5, ~ is.na(.) | . == ""))

stopifnot(nrow(missing_mapping) == 0)

# Map exposure CSV labels onto syntactic analysis variable names

# syntactic labels to make the coding easier
exposure_mapping <- c(
  "minimal" = "minimal",
  "indirect" = "indirect",
  "fomites" = "fomites",
  "direct contact" = "direct",
  "nursing care" = "nursing",
  "body fluid exposure" = "body_fluid",
  "handled corpse" = "corpse"
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

dat_joined <- dat_joined %>%
  mutate(across(exposure1:exposure5,
                ~ unname(exposure_mapping[.x]) ) )

sort(unique(unlist(dat_joined %>% select(exposure1:exposure5))))

# Identify how many canonical exposure categories are represented in each row
# (ignoring blanks)

dat_joined <- dat_joined %>%
  mutate(n_exposures_mapped = rowSums(!is.na(dplyr::across(exposure1:exposure5)) & dplyr::across(exposure1:exposure5) != ""))

# check that all with no mapping are for "all contacts
dat_joined %>%
  dplyr::filter(n_exposures_mapped == 0) %>%
  dplyr::pull(definition_contact_me)
# passes check

# Identify studies with a completely disaggregated exposure distribution:
# - every exposure row maps to exactly one exposure category
# - the exposure-specific denominators sum to the total number of contacts
# - rows labelled "All" are excluded from this calculation

disaggregated_studies <- dat_joined %>%
  filter(!definition_contact_me %in% c("All", "All contacts")) %>%
  group_by(label, total_contacts) %>%
  summarise(all_single_exposure = all(n_exposures_mapped == 1),
            sum_contacts = sum(denominator),
            .groups = "drop") %>%
  filter(all_single_exposure,
         sum_contacts == total_contacts) %>%
  pull(label)

disaggregated_studies

# [1] "Bower, 2016a" "Bower, 2016b"

# Note however, that one of these Bower studies only includes children < 3 and therefore is not
# representative => only use "Bower, 2016b"
disaggregated_studies <- "Bower, 2016b"

# Due to only one study, we now use the empirical information
reference_distribution <- dat_joined %>%
  filter(label %in% disaggregated_studies, !definition_contact_me %in% c("All", "All contacts")) %>%
  pivot_longer(cols = exposure1:exposure5, values_to = "exposure",values_drop_na = TRUE) %>%
  filter(exposure != "") %>%
  group_by(exposure) %>%
  summarise(n_contacts = sum(denominator), .groups = "drop") %>%
  mutate(proportion = n_contacts / sum(n_contacts))

print(reference_distribution)

# Check that proportions sum to 1

stopifnot(abs(sum(reference_distribution$proportion) - 1) < 1e-8)

# 3. Create fractional exposure vectors -------------------------------

# Put the reference proportions into a named vector for easy lookup

reference_props <- reference_distribution$proportion
names(reference_props) <- reference_distribution$exposure

# For each row, identify the canonical exposure categories it covers and

# assign the reference proportions within those categories.

dat_fractional <- dat_joined %>%
  dplyr::filter(!(definition_contact_me %in% c("All", "All contacts"))) %>%
  rowwise() %>%
  mutate(mapped_exposures = list(c(exposure1, exposure2, exposure3, exposure4, exposure5) %>%
      na.omit() )) %>%
  mutate(mapped_exposures = list(mapped_exposures[mapped_exposures != ""] )) %>%
  mutate(mapped_props = list(reference_props[match(mapped_exposures, names(reference_props))] )) %>%
  mutate(mapped_props = list(mapped_props / sum(mapped_props) )) %>%
  ungroup()


# 4. Create wide fractional exposure design matrix ---------------------

# Create one column for each exposure category

all_exposures <- names(canonical_labels)

# Create one column per exposure

for (lvl in all_exposures) {
  dat_fractional[[lvl]] <- purrr::map2_dbl(dat_fractional$mapped_exposures, 
                                           dat_fractional$mapped_props, 
                                           ~{ idx <- match(lvl, .x) 
    if (is.na(idx)) {
      0
    } else {
      .y[idx]
    } }
  )
}

# Check that the rows sum to 1 -- check passed
dat_fractional %>%
  dplyr::mutate(fraction_sum = rowSums(across(all_of(all_exposures)))) %>%
  pull(fraction_sum)

# 5. Expand to pseudo-individual contacts -------------------------------
# In order to perform the logistic regression, we need to reformat the data into the form
# where each row represents an 'individual'
# Give each original row a unique identifier before expanding

dat_fractional <- dat_fractional %>%
  mutate(source_row = row_number())

# Expand each row to one pseudo-individual per contact

pseudo_contacts <- dat_fractional %>%
  dplyr::select(source_row, label, household, outbreak,  
                numerator, denominator, total_contacts, dplyr::all_of(all_exposures)) %>%
  uncount(denominator, .id = "contact", .remove = FALSE) %>%
  group_by(source_row) %>%
  mutate(outcome = as.integer(contact <= first(numerator))) %>%
  ungroup()


# 6. All contacts ---------------------------------------------------------

# Can use the original dataset to filter for the information we are interested in for "all contacts"
unique(dat$definition_contact_me) %>% sort()

dat_all <- dat %>%
  dplyr::filter(definition_contact_me %in% c("All", "All contacts")) %>%
  left_join(labels) %>%
  relocate(label)


# 7. Save data ------------------------------------------------------------

saveRDS(dat_all, "data/data-derived/all_contacts.rds")
saveRDS(pseudo_contacts, "data/data-derived/data_logistic.rds")
saveRDS(dat_fractional, "data/data-derived/sar_exposure.rds")
