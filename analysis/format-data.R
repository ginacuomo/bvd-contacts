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

# data import 
data <- readxl::read_excel("data/data-raw/sar_extraction.xlsx", sheet = "data-fomite")

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
