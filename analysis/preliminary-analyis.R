# preliminary-analyis.R
# Author: Gina Cuomo-Dannenburg
# Date: 2026-09-02
# Purpose: Perform the preliminary analysis of SARs across studies and exposures based on the current dataset (see date)
# Input: Data extracted and then cleaned and mapped in format-data.R
# Output: Preliminary results of the analysis
# Data download date: 2026-08-25

# libraries
library(ggplot2)
library(tidyverse)

dat <- readRDS("data/data-derived/dat_joined.rds")

