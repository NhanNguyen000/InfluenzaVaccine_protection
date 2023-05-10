library(tidyverse)

# load data from 2 cohorts
load("../ZirFlu_NhanNguyen/ZirFlu.RData")
load("../iMED_NhanNguyen/iMED.RData")

cohorts_dat <- list()
# donor info ---------------------------------------
cohorts_dat$donorInfo_all <- ZirFlu$donorInfo %>% mutate(cohort = "ZirFlu") %>%
  full_join(iMED$donorInfo %>% 
              select(patientID, gender, age) %>%
              rename("sex" = "gender") %>%
              mutate(season = NA, 
                     disease = "healthy", 
                     condition = "healthy", 
                     cohort = "iMED")) %>%
  mutate(sex = substring(sex, 1, 1)) %>%
  mutate(age_group = ifelse(age >= 60, "old", "young"))

# check groups
cohorts_dat$donorInfo_all %>% 
  count(cohort, season, disease, age_group)

# donor sample ---------------------------------------
cohorts_dat$donorSample_all <- ZirFlu$donorSamples %>% 
  mutate(cohort = "ZirFlu", name = probenID, probenID = as.numeric(probenID)) %>%
  full_join(iMED$donorSample %>% 
              dplyr::select(patientID, probenID, name, time) %>% 
              mutate(cohort = "iMED"))

# HAI titer ---------------------------------------
cohorts_dat$HAItiter_all <- cohorts_dat$donorInfo_all %>%
  left_join(ZirFlu$HAItiter %>% select(patientID, season, category, matches("ab|T1|T2"))) %>%
  left_join(iMED$HAItiter %>% select(patientID, responder, matches("ab|T1|T4")) %>%
              rename("iMED_H1N1_T1" = "H1N1_T1", "iMED_H1N1_T4" = "H1N1_T4",
                     "iMED_H3N2_T1" = "H3N2_T1", "iMED_H3N2_T4" = "H3N2_T4",
                     "iMED_B_T1" = "B_T1", "iMED_B_T4" = "B_T4")) %>%
  mutate(category = ifelse(is.na(category), responder, category)) %>% 
  select(-responder) %>%
  mutate(category = ifelse(category == "other", "Other", category)) %>%
  mutate(H1N1_abFC_combine = ifelse(is.na(H1N1_abFC), ab_H1N1, H1N1_abFC)) %>%
  mutate(H3N2_abFC_combine = ifelse(is.na(H3N2_abFC), ab_H3N2, H3N2_abFC)) %>%
  mutate(Bvictoria_abFC_combine = ifelse(is.na(Bvictoria_abFC), ab_B, Bvictoria_abFC)) %>%
  mutate(Byamagata_abFC_combine = ifelse(is.na(Byamagata_abFC), ab_B, Byamagata_abFC)) %>%
  mutate(H1N1_T1_combine = ifelse(is.na(H1N1_T1), iMED_H1N1_T1, H1N1_T1),
         H3N2_T1_combine = ifelse(is.na(H3N2_T1), iMED_H3N2_T1, H3N2_T1))

# check sample size per group
cohorts_dat$HAItiter_all %>% 
  count(cohort, season, disease, age_group, category)

cohorts_dat$HAItiter_all %>% 
  count(cohort, disease, age_group, category)

cohorts_dat$HAItiter_all %>% filter(disease == "healthy") %>%
  count(cohort, age_group, category)

# protein annotation ---------------------------------------
cohorts_dat$proteinAnnot_all <- ZirFlu$protein_annot %>% # ZirFlu has 313 proteins, iMED has 311 proteins
  full_join(iMED$protein_annot) # 2 cohort has total 311 proteins

# save data ------------------------------------------------
# save(cohorts_dat, file = "cohorts_dat.RData")

