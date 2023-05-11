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
      mutate(name = probenID, cohort = "ZirFlu",
             time = ifelse(time == "Baseline", "T1", 
                           ifelse(time == "T1", "T4", 
                                  ifelse(time == "T2", "T5", time)))) # check time match with visit 1 and 2, and between cohorts
  ) %>%
  # clean the data
  select(name, time, season, probandID, cohort) %>% relocate(probandID)

# HAI titer ---------------------------------------
cohorts$HAI_all <- iMED$HAI_2014 %>% 
  # iMED pilot cohort (n = 34)
  rename("probandID" = "pID") %>% mutate(season = "2014") %>%
  full_join(iMED$meta2014 %>% select(probandID, responder) %>% distinct()) %>%
  mutate(ab_H1N1 = ifelse(H1N1_T1 == 0, H1N1_T4, H1N1_T4/H1N1_T1),
         ab_H3N2 = ifelse(H3N2_T1 == 0, H3N2_T4, H3N2_T4/H3N2_T1),
         ab_B = ifelse(B_T1 == 0, B_T4, B_T4/B_T1)) %>%
 
  # iMED discovery cohort (n = 200)
  full_join(
    iMED$HAI_2015 %>% rename("probandID" = "X...sample") %>% 
     full_join(iMED$meta2015 %>% select(probandID, matches("ab")) %>% drop_na()) %>%
      mutate(season = "2015")
  ) %>% 
  mutate(cohort = "iMED") %>% rename("H1N1_abFC" = "ab_H1N1", "H3N2_abFC" = "ab_H3N2", "B_abFC" = "ab_B") %>%
  
  # ZirFlu cohorts, 2 season
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
overlapped_proteins <- intersect(names(protein_normDat$iMED), names(protein_normDat$ZirFlu))
length(overlapped_proteins)
# iMED has 311 proteins, ZirFlu has 308 proteins -> overlapped proteins: 306 proteins

protein_Dat <- list()

protein_Dat$iMED_2014 <- protein_normDat$iMED[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2014"], overlapped_proteins]

protein_Dat$iMED_2015 <- protein_normDat$iMED[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2015"], overlapped_proteins]

protein_Dat$ZirFlu_2019 <- protein_normDat$ZirFlu[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2019"], overlapped_proteins]

protein_Dat$ZirFlu_2020 <- protein_normDat$ZirFlu[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2020"], overlapped_proteins]

# metabolites ---------------------------------------
overlapped_metabolites <- intersect(iMED_meboAnnot$Formula_v2, ZirFlu_meboAnnot$Formula)
length(overlapped_metabolites)
# iMED has 1345 metabolites, ZirFlu has 786 metabolites -> overlapped metabolites: 508 metabolites

mebo_Dat <- list()
# iMED have case: 1 formula - multiple ionIdxs
overlapped_metabolites_iMED <- iMED_meboAnnot %>% filter(Formula %in% overlapped_metabolites)

mebo_Dat$iMED_2014 <- iMED_mebo[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2014"], overlapped_metabolites_iMED$Formula_v2]

mebo_Dat$iMED_2015 <- iMED_mebo[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2015"], overlapped_metabolites_iMED$Formula_v2]

mebo_Dat$ZirFlu_2019 <- ZirFlu_mebo[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2019"], overlapped_metabolites]

mebo_Dat$ZirFlu_2020 <- ZirFlu_mebo[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2020"], overlapped_metabolites]

# save data ------------------------------------------------
# save(cohorts, protein_Dat, mebo_Dat, file = "cohorts_dat.RData")

