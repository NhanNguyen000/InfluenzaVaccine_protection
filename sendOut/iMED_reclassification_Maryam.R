rm(list = ls())
library(tidyverse)

load("/vol/projects/BIIM/Influenza/ZirrFlu/InfluenzaCohorts_NhanNguyen/cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
iMED_2015 <- cohorts$HAI_all %>% 
  filter(cohort == "iMED", season == "2015") %>% 
  select(-matches("victoria|yamagata|vaccine")) %>%
  left_join(cohorts$donorInfo_all %>% select(probandID, sex, age)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>% 
  relocate(name, sex, age, 
           H1N1_reclassify, H3N2_reclassify, B_reclassify,
           .after = probandID)

write.table(iMED_2015, file = "sendOut/iMED_2015_cohort_reclassify.txt", quote = FALSE, row.names = FALSE)


# Column name explanation ---------------
# ProbandID: participant ID
# name: sample ID from the participant
# H1N1, H3N2, B are 3 strain
# H1N1_reclassify  = the 4 class for H1N1 strain, etc
# H1N1_T1, H1N1_T4, H1N1_abFC = the antibody titer at T1, T4 and antibody foldchange for H1N1 strain
# H1N1_T1_log2, H1N1_T4_log2, H1N1_abFC_log2 = the log2(value) for H1N1 strain
# season: the year of the cohort
# cohort: the name of the cohort
# responder: the groups based on responder definition across 3 strains