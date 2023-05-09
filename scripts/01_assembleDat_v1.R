rm(list = ls())

library(tidyverse)
library(readxl)

# load data ---------------------------------------
load("metaDat_antibody.RData")
load("protein_normDat.RData")
load("metaboliteDat.RData")

# donor info ---------------------------------------
cohorts <- list()

cohorts$donorInfo_all <- iMED$meta2014 %>% 
  # iMED pilot cohort (n = 34)
  select(probandID, gender, age) %>% distinct() %>% mutate(season = "2014") %>%
  
  # iMED discovery cohort (n = 200)
  full_join(iMED$meta2015 %>% 
      select(probandID, gender, age) %>% distinct() %>% mutate(season = "2015")) %>%
  mutate(condition = "Healthy", cohort = "iMED") %>%
  
  # ZirFlu cohort, 2 seasons
  full_join(
    ZirFlu$HAI_2019 %>% mutate(season = "2019") %>%
      full_join(ZirFlu$HAI_2020 %>% mutate(season = "2020")) %>% 
      select(patientID, condition, season) %>%
      left_join(
        read_excel("/vol/projects/CIIM/Influenza/ZirrFlu/proteomic/raw_data/ZirrFlu plates final.xlsx",
                   sheet = "Phenotypes") %>% # get the participant's age and sex information
          select(patientID, Age, Sex, Condition, Season) %>% distinct() %>% 
          rename("season" = "Season", "age" = "Age", "gender" = "Sex")) %>%
      mutate(cohort = "ZirFlu") %>% rename("probandID" = "patientID")
  ) %>%
  # clean the data
  select(-Condition) %>%
  mutate(sex = substring(gender, 1, 1),
         age_group = ifelse(age >= 60, "old", "young"),
         disease = ifelse(condition == "Healthy", "healthy", "cirrhosis"))

# check groups
cohorts$donorInfo_all %>% 
  count(cohort, season, disease, age_group)

# donor sample ---------------------------------------
cohorts$donorSample_all <- iMED$meta2014 %>% 
  # iMED pilot cohort (n = 34)
  select(probandID, ProbenID, name, time) %>% mutate(season = "2014") %>%
  
  # iMED discovery cohort (n = 200)
  full_join(iMED$meta2015 %>% 
      select(probandID, ProbenID, name, time) %>% mutate(season = "2015")) %>% 
  mutate(cohort = "iMED") %>%
  
  # ZirFlu cohort, 2 seasons
  full_join(
    read_excel("/vol/projects/CIIM/Influenza/ZirrFlu/proteomic/raw_data/ZirrFlu plates final.xlsx",
               sheet = "Phenotypes") %>%
      select(patientID, probenID, Season, Time) %>% filter(Time != "Bridge") %>%
      rename("probandID" = "patientID", "season" = "Season", "time" = "Time") %>%
      mutate(name = probenID, cohort = "ZirFu",
             time = ifelse(time == "Baseline", "T1", 
                           ifelse(time == "T1", "T2", 
                                  ifelse(time == "T2", "T3", time)))) # check time match with visit 1 and 2, and between cohorts
  ) %>%
  # clean the data
  select(name, time, season, probandID, cohort) %>% relocate(probandID)

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

cohorts$HAI_all <- iMED$HAI_2014 %>% 
  # iMED pilot cohort (n = 34)
  rename("probandID" = "pID") %>% mutate(season = "2014") %>%
  full_join(iMED$meta2014 %>% select(probandID, responder) %>% distinct()) %>%
 
  # iMED discovery cohort (n = 200)
  full_join(
    iMED$HAI_2015 %>% rename("probandID" = "X...sample") %>% 
     # full_join(iMED$meta2015 %>% select(probandID, responder, matches("ab")) %>% drop_na()) %>%, the abFC are strange, need to check
      mutate(season = "2015")
  ) %>% mutate(cohort = "iMED") %>%
  full_join(
    ZirFlu$HAI_2019 %>% mutate(season = "2019") %>%
      full_join(ZirFlu$HAI_2020 %>% mutate(season = "2020")) %>%
      select(patientID, season, matches("d0|d21|abFC"), vaccine_response) %>%
      setNames(names(.) %>% 
                 stringr::str_replace("_d0", "_T1") %>% 
                 stringr::str_replace("_d21-35", "_T4")) %>%
      mutate(responder = ifelse(vaccine_response == "Non", "NR", 
                                ifelse(vaccine_response == "Single" | vaccine_response == "Double", 
                                       "other", "TR"))) %>%
      mutate(cohort = "ZirFlu") %>% rename("probandID" = "patientID")
  ) %>%
  # clean the data
  mutate(responder = ifelse(responder == "other", "Other", responder))

# protein ---------------------------------------
length(intersect(names(protein_normDat$iMED), names(protein_normDat$ZirFlu)))
# iMED has 311 proteins, ZirFlu has 308 proteins -> overlapped proteins: 306 proteins

# metabolites ---------------------------------------
intersect(iMED_meboAnnot$Formula, ZirFlu_meboAnnot$Formula))
# iMED has 1345 metabolites, ZirFlu has 786 metabolites -> overlapped metabolites: 549 metabolites

a <- iMED_meboAnnot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_")) %>%
  group_by(ionIdx, Formula) %>%
  mutate(Formula_v2 = Formula + row_number() -1)
#     id value onemore

double_Formula <- a %>% count(Formula)

mebo_Dat <- iMED_mebo %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("ionIdx") %>% mutate(ionIdx = as.numeric(ionIdx)) %>%
  full_join(
    iMED_meboAnnot %>% select(ionIdx, Formula) %>% distinct() %>%
      group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))
  ) %>% select(-ionIdx) %>% column_to_rownames("Formula")

b <- ZirFlu_meboAnnot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

## iMED metabolites -------------------------
iMED_metabolites <- iMED$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

iMED_metaboliteDat <- iMED$metabolite
if (identical(iMED_metabolites$ionIdx, as.numeric(colnames(iMED$metabolite))) == TRUE) {
  colnames(iMED_metaboliteDat) <- iMED_metabolites$Formula
}


## ZirFlu metabolites -------------------------
ZirFlu_metabolites <- ZirFlu$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

ZirFlu_metaboliteDat <- ZirFlu$metabolite_dat
if (identical(ZirFlu_metabolites$ionIdx, as.numeric(colnames(ZirFlu$metabolite_dat))) == TRUE) {
  colnames(ZirFlu_metaboliteDat) <- ZirFlu_metabolites$Formula
}

## overlap metabolites between iMED and ZirFlu -------------------------
overlapped_metabolites <- intersect(iMED_metabolites$Formula, ZirFlu_metabolites$Formula)

metadat_iMED <- metadat_healthy %>% filter(cohort == "iMED")
inputDat_iMED <- iMED_metaboliteDat[, overlapped_metabolites] %>% 
  t() %>% as.data.frame() %>% select(metadat_iMED$name)

metadat_ZirFlu_healthy <- metadat_healthy %>% filter(cohort == "ZirFlu")
inputDat_ZirFlu_healthy <- ZirFlu_metaboliteDat[, overlapped_metabolites] %>% 
  t() %>% as.data.frame() %>% select(metadat_ZirFlu_healthy$name)

iMED_mebo, iMED_meboAnnot, ZirFlu_mebo, ZirFlu_meboAnnot
# save data ------------------------------------------------
# save(cohorts, file = "cohorts.RData")

