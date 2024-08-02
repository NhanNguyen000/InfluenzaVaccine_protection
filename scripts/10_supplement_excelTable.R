rm(list = ls())

library(openxlsx)
library(tidyverse)

# sample information ===============================================================
load("processedDat/cohorts_dat.RData")

convert_protectees <- function(input) {
  # Aim: convert the vecotr of 4 groups reclassification (LL, LH, HL, HH) into 2 groups: LL and protectee (LH, HL HH)
  outcome = ifelse(input %in% c("LH", "HL", "HH"), "protectee", input) 
  return(outcome)
}

sampleDat <- cohorts$HAI_all %>% ## metadata for all healthy subjects (we only use healthy subject)
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~convert_protectees(.x)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "protectee"))) %>%
  select(-name, -time)

# check numbers
sampleDat %>% count(season)

sampleDat %>% count(season, H1N1_reclassify)
sampleDat %>% count(season, H3N2_reclassify)
sampleDat %>% count(season, B_reclassify)
sampleDat %>% count(season, Bvictoria_reclassify)
sampleDat %>% count(season, Byamagata_reclassify)

rm(convert_protectees, cohorts, mebo_Dat, protein_Dat)

# sample information ===============================================================
load("processedDat/selectedProteins_statOutcome.RData")

selectedProteins_statOutcome <- tstat_longDat_consistVars
rm(tstat_longDat_consistVars)



# save as excel files ====================================================================================
write.xlsx(list(sampleDat = sampleDat,
                selectedProteins_statOutcome = selectedProteins_statOutcome, 
           "output/suppTables.xlsx", rowNames = TRUE)
