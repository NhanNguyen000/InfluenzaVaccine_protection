# clean the environment
rm(list = ls())

library(tidyverse)

# load data from 2 cohorts -------------------------------
load("../ZirFlu_NhanNguyen/ZirFlu.RData")
load("../iMED_NhanNguyen/iMED.RData")

# load processed and cleanned data -------------------------------
load("cohorts_dat.RData")
load("protein_normDat.RData")
load("proteinDat_impute.RData")
load("HAIreclassify.RData")
load("HAIreclassify2.RData")
# check data -----------------------
a <- cohorts_dat$donorInfo_all %>% full_join(HAIreclassify2$all_cohorts)
a2 <- HAIreclassify2$iMED %>% left_join(cohorts_dat$donorInfo_all)
a3 <- HAIreclassify2$ZirFlu %>% left_join(cohorts_dat$donorInfo_all)

