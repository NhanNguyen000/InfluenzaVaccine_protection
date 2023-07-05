rm(list = ls())

library(tidyverse)
# load data =======================================================================
proPred_dat <- read.csv2("data/Olink_trait_validation_results_with_OMICSPRED_ID.csv", sep = "\t")
load("selected_DAPs.RData")

selected_proPred <- proPred_dat %>% 
  filter(Gene %in% selected_DAs) %>% 
  select(OMICSPRED.ID, Gene, X.SNP, Internal_R2, NSPHS_R2, ORCADES_R2)
