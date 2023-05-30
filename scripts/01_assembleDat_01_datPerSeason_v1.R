# NOTE: I log2 the antibody titer and antibodi foldchange for all strains and log2 the metabolite data to normalize the data distributioin
# protein data is already in the normal distribution, so I do not log2 protein data
rm(list = ls())

library(tidyverse)
library(readxl)

get.log2 <- function(dat) {
  library(tidyverse)
  library(rstatix)
  outcome <- dat %>% mutate(across(where(is.numeric), ~log2(.x))) %>%
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))
  return(outcome)
}

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
  mutate(responder = ifelse(responder == "other", "Other", responder)) %>%
  mutate(across(ends_with("T1") | ends_with("T4") | ends_with("abFC"), # convert antibody titer and abFC to normalize the value distribution
                ~log2(.x), .names = "{.col}_log2")) %>%
  mutate(across(ends_with("log2"), ~ifelse(is.infinite(.x), 0, .x)))


# summary(cohorts$HAI_all$H1N1_T1) # check the value range
# summary(cohorts$HAI_all$H1N1_T1_log2)

# hist(cohorts$HAI_all$H1N1_T1) # show skew distribution
# hist(cohorts$HAI_all$H1N1_T1_log2) # show less skew distribution
# hist(cohorts$HAI_all$H1N1_abFC)
# hist(cohorts$HAI_all$H1N1_abFC_log2)

# MN titer ---------------------------------------
cohorts$MN_all <- iMED$MNbothCohorts %>% rename("probandID" = "X..." )%>%
  # iMED pilot and discovery cohort (n = 34 + n = 200, total n = 234)
  full_join(iMED$meta2014 %>% select(probandID) %>% 
              distinct() %>% mutate(season = "2014") %>%
              full_join(iMED$meta2015 %>% select(probandID) %>% 
                          distinct() %>% mutate(season = "2015"))) %>%
  mutate(cohort = "iMED", 
         "H1N1_abFC" = ifelse(H1N1_T1 == 0, H1N1_T4, H1N1_T4/H1N1_T1), 
         "H3N2_abFC" = ifelse(H3N2_T1 == 0, H3N2_T4, H3N2_T4/H3N2_T1), 
         "B_abFC" = ifelse(B_T1 == 0, B_T4, B_T4/B_T1)) %>%
  # ZirFlu cohorts, season 2019
  full_join(
    ZirFlu$MN_2019 %>% select(-condition, -matches("T3")) %>%
      rename("probandID" = "patientID") %>% 
      setNames(names(.) %>% stringr::str_replace("_T2", "_T4")) %>%
      mutate(cohort = "ZirFlu")
  ) %>%
  # clean the data
  mutate(across(ends_with("T1") | ends_with("T4") | ends_with("abFC"), # convert antibody titer and abFC to normalize the value distribution
                ~log2(.x), .names = "{.col}_log2")) %>%
  mutate(across(ends_with("log2"), ~ifelse(is.infinite(.x), 0, .x)))


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

# summary(unlist(protein_Dat$iMED_2014)) # check the value range
# summary(protein_Dat$iMED_2014$ADA) # check the value range
# 
# hist(unlist(protein_Dat$iMED_2014))
# hist(protein_Dat$iMED_2014$ADA) # show quite normal distribution

# metabolites ---------------------------------------
overlapped_metabolites <- intersect(unique(iMED_meboAnnot$Formula), ZirFlu_meboAnnot$Formula)
length(overlapped_metabolites)
# iMED has 1326 metabolites (with 1345 metabolite Idx), ZirFlu has 786 metabolites -> overlapped metabolites: 508 metabolites

mebo_Dat <- list()
# iMED have case: 1 formula - multiple ionIdxs
overlapped_metabolites_iMED <- iMED_meboAnnot %>% filter(Formula %in% overlapped_metabolites)

mebo_Dat$iMED_2014 <- iMED_mebo[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2014"], unique(overlapped_metabolites_iMED$Formula)] %>%
  get.log2()

mebo_Dat$iMED_2015 <- iMED_mebo[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2015"], unique(overlapped_metabolites_iMED$Formula)] %>%
  get.log2()

mebo_Dat$ZirFlu_2019 <- ZirFlu_mebo[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2019"], overlapped_metabolites] %>%
  get.log2()

mebo_Dat$ZirFlu_2020 <- ZirFlu_mebo[
  cohorts$donorSample_all$name[cohorts$donorSample_all$season == "2020"], overlapped_metabolites] %>%
  get.log2()

summary(unlist(mebo_Dat$iMED_2014)) # check the value range
summary(mebo_Dat$iMED_2014$C3H4O) # check the value range

hist(unlist(mebo_Dat$iMED_2014)) # show less scale distribution compared to the raw data without log2
hist(mebo_Dat$iMED_2014$C3H4O) # show quite normal distribution


# save data ------------------------------------------------
save(cohorts, protein_Dat, mebo_Dat, file = "cohorts_dat.RData")

