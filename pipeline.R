# clean the environment
rm(list = ls())

library(tidyverse)
# # load data from 2 cohorts -------------------------------
# load("../ZirFlu_NhanNguyen/ZirFlu.RData")
# load("../iMED_NhanNguyen/iMED.RData")

# load processed and cleanned data -------------------------------
load("cohorts_dat.RData") # after run the code: 00_combine_patientInfo_v1.R
load("protein_normDat.RData") # after run the code: 00_normalize_proteinDat_2cohorts_v1.R
load("proteinDat_impute.RData") # after run the code: 01_imputeProteinNA_datNorm2cohort_v1.R
load("HAIreclassify.RData") # 01_reclasify_HAIgroups_v1.R
load("HAIreclassify2.RData") # 01_reclasify_HAIgroups_changeBaseline_v1.R
# check data -----------------------
a <- cohorts_dat$donorInfo_all %>% full_join(HAIreclassify2$all_cohorts)
a2 <- HAIreclassify2$iMED %>% left_join(cohorts_dat$donorInfo_all)
a3 <- HAIreclassify2$ZirFlu %>% left_join(cohorts_dat$donorInfo_all)

