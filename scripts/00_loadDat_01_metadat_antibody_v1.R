library(tidyverse)
library(readxl)

get.log2 <- function(dat) {
  library(tidyverse)
  library(rstatix)
  outcome <- dat %>% mutate(across(where(is.numeric), ~log2(.x))) %>%
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))
  return(outcome)
}

# iMED ================================================================================
iMED <- list()
## metadata  --------------------------------------------------------
iMED$metaCohort1 <- read.csv('/vol/projects/CIIM/Influenza/iMED/metadata/meta_cohort1.csv', row.names = 1)
iMED$metaCohort2 <- read.csv('/vol/projects/CIIM/Influenza/iMED/metadata/meta_cohort2.csv', row.names = 1)

## HAI antibody titers -----------------------------------------------
iMED$HAIcohort1 <- read.table('/vol/projects/CIIM/Influenza/iMED/metadata/hai_titers_pilot.tsv', header = TRUE)
iMED$HAIcohort2 <- read.csv2('/vol/projects/CIIM/Influenza/iMED/metadata/hai_titers.csv')

## MN antibody titers ---------------------------------------------
iMED$MNbothCohorts <- read.csv2('/vol/projects/CIIM/Influenza/iMED/metadata/MNtiters_raw.csv')

# ZirFlu ================================================================================
ZirFlu <- list()
## metadata  ------------------------------------------------------------------
ZirFlu$meta2019 <- read.table(file = '/vol/projects/CIIM/Influenza/ZirrFlu/metadata/zirrflu_meta.csv', sep = '\t', header=T)

## HAI antibody titers ------------------------------------------------------------
### season 2019 --------------------------------------
HAItiter_2019 <- read_excel("/vol/projects/CIIM/Influenza/ZirrFlu/metadata/20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx",
  sheet = "2019-2020 HAI Titer") %>% 
  select(-c(15:19)) %>%
  drop_na("HAI Titer against Influenza A/H1N1") %>% 
  as.data.frame()

names(HAItiter_2019) <- c("condition", "patientID", 
                          paste0(rep(c("H1N1_", "H3N2_", "Bvictoria_", "Byamagata_"), each = 3),
                                 HAItiter_2019[1,][3:14]))

HAIabFC_2019 <- read_excel("/vol/projects/CIIM/Influenza/ZirrFlu/metadata/20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx",
                            sheet = "2019-2020 Fold-increase") %>% 
  drop_na("Fold-increase") %>%
  as.data.frame()

names(HAIabFC_2019) <- c("condition", "patientID", 
                            "H1N1_abFC", "H3N2_abFC", "Bvictoria_abFC", "Byamagata_abFC",
                            "vaccine_response")
ZirFlu$HAI_2019 <- HAItiter_2019 %>% slice(-1) %>% fill(condition) %>%
  full_join(HAIabFC_2019 %>% slice(-1) %>% fill(condition))

rm(HAIabFC_2019, HAItiter_2019)

### season 2020 ----------------------------------
HAItiter_2020 <- read_excel("/vol/projects/CIIM/Influenza/ZirrFlu/metadata/20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx",
                          sheet = "2020-2021 HAI Titer") %>% 
  select(-c(15:17)) %>%
  drop_na("HAI Titer against Influenza A/H1N1") %>% 
  as.data.frame()

names(HAItiter_2020) <- c("condition", "patientID", 
                          paste0(rep(c("H1N1_", "H3N2_", "Bvictoria_", "Byamagata_"), each = 3),
                                 HAItiter_2020[1,][3:14]))


HAIabFC_2020 <- read_excel("/vol/projects/CIIM/Influenza/ZirrFlu/metadata/20221101_Final_Info_HAI_Titer_Valerie&Nhan.xlsx",
                            sheet = "2020-2021 Fold-increase") %>% 
  drop_na("Fold-increase") %>% as.data.frame()

names(HAIabFC_2020) <- c("condition", "patientID", 
                            "H1N1_abFC", "H3N2_abFC", "Bvictoria_abFC", "Byamagata_abFC",
                            "vaccine_response")


ZirFlu$HAI_2020 <- HAItiter_2020 %>% slice(-1) %>% fill(condition) %>%
  full_join(HAIabFC_2020 %>% slice(-1) %>% fill(condition))

rm(HAIabFC_2020, HAItiter_2020)
## MN antibody titers ---------------------------------------------------------
### season 2019 --------------------------
MNtiter_2019 <- read_excel("/vol/projects/CIIM/Influenza/ZirrFlu/metadata/20230105_ZirFlu-2019-2020_MN_Titer_Valerie&Nhan.xlsx",
                            sheet = "2019-2020 MN Titer") %>%
  select(1:14) %>% 
  slice(-c(1, 2, 36, 52)) %>%
  fill(1, .direction = "down") %>%
  as.data.frame() %>%
  mutate_at(c(3:14), as.numeric)

names(MNtiter_2019) <- c("condition", "patientID",
                          paste0("H1N1_", c("T1", "T2", "T3")),
                          paste0("H3N2_", c("T1", "T2", "T3")),
                          paste0("Bvictoria_", c("T1", "T2", "T3")),
                          paste0("Byamagata_", c("T1", "T2", "T3")))

MNabFC_2019 <- read_excel("/vol/projects/CIIM/Influenza/ZirrFlu/metadata/metadata/20230105_ZirFlu-2019-2020_MN_Titer_Valerie&Nhan.xlsx",
                                sheet = "2019-2020 MN Fold-increase") %>%
  slice(-c(1, 2, 36, 52)) %>%
  fill(1, .direction = "down") %>%
  as.data.frame() %>%
  mutate_at(c(3:6), as.numeric)

names(MNabFC_2019) <- c("condition", "patientID",
                        paste0(c("H1N1_", "H3N2_", "Bvictoria_", "Byamagata_"), "abFC"))

ZirFlu$MN_2019 <- MN_iter_2019 %>% full_join(MNabFC_2019) %>%
  mutate(season = "2019", condition = str_to_lower(condition)) %>%
  relocate(season)

# save data ================================================================================
cohorts_rawDat <- list("iMED" = iMED, "ZirFlu" = ZirFlu)
# save data
save(cohorts, file = "cohorts.RData")



