rm(list = ls())

library(tidyverse)
library(OlinkAnalyze)
library(openxlsx)
library(readxl)

# load data =========================================================================
# iMED data --------------------------------------------------------------------
rawDat_iMEDprotein <- read_NPX('/vol/projects/CIIM/Influenza/iMED/proteomic/raw_data/20210222_Li_NPX_2021-05-18.csv') 

iMED_meta2015 <- read.csv('/vol/projects/CIIM/Influenza/iMED/metadata/meta_cohort2.csv', row.names = 1)
iMED_HAI_2014 <- read.table('/vol/projects/CIIM/Influenza/iMED/metadata/hai_titers_pilot.tsv', header = TRUE)
iMED_HAI_2015 <- read.csv2('/vol/projects/CIIM/Influenza/iMED/metadata/hai_titers.csv')


# ZirFlu --------------------------------------------------------------------
ZirFlu_proteinRawDat <- read_NPX(filename = "/vol/projects/CIIM/Influenza/ZirrFlu/proteomic/raw_data/20212645_Li_NPX_2022-02-02.csv")
ZirFlu_proteinPlate <- read_excel("/vol/projects/CIIM/Influenza/ZirrFlu/proteomic/raw_data/ZirrFlu plates final.xlsx",
                                  sheet = "Phenotypes")
# Martijn use this file to process the protein Olink sample
ZirFlu_samples <- read_excel("/vol/projects/BIIM/Influenza/ZirrFlu/proteomic/raw_data/Auszug-ZirFlu 1 und 2 CentraXX23112021_YangLi.xlsx")

# check ZirFlu

ZirFLu_allinfo <- ZirFlu_samples %>% 
  select(ProbandenID, ProbenID, Geschlecht,
         DatumProbe, DatumEpisode, DatumBasis, 
         Therapiephase, Zeitpunkt) %>% 
  full_join(ZirFlu_proteinPlate %>% 
              select(patientID, probenID, probeDatum, Sex, pID_time, Season, Condition), 
            by = c("ProbandenID" = "patientID", "ProbenID" = "probenID")) %>% # have 10 more bridging-samples from iMED
  full_join(ZirFlu_proteinRawDat %>% 
              select(SampleID, PlateID) %>% distinct(),
            by = c("ProbenID" = "SampleID")) # have 8 more CONTROL_SAMPLE --> total 1478 samples at the end

ZirFlu_proteinSample <- ZirFLu_allinfo %>% 
  drop_na(PlateID) %>% # 318 sample = 300 (ZirFlu sample) + 10 (briding sample) + 8 CONTROL_SAMPLE
  drop_na(pID_time) %>% # remove 10 (briding sample) and 8 CONTROL_SAMPLE
  separate(pID_time, c("pID_proteinPlate", "time_proteinPlate"), sep = " ") %>%
  mutate(Zeitpunkt_v2 = ifelse(Zeitpunkt == "-", "Baseline", 
                               ifelse(Zeitpunkt == "01", "T1", 
                                      ifelse(Zeitpunkt == "02", "T2", 
                                             ifelse(Zeitpunkt == "03", "T3", Zeitpunkt))))) %>% # 300 samples, but have some double sample
  distinct()

all.equal(ZirFlu_proteinSample$ProbandenID, ZirFlu_proteinSample$pID_proteinPlate) # TRUE, the same 277 samples

missmatch_timePoint <- ZirFlu_proteinSample %>% 
  slice(which(ZirFlu_proteinSample$Zeitpunkt_v2 != ZirFlu_proteinSample$time_proteinPlate)) # 4 miss matched in the Zeitpunk, (3 visit 3 become T1, 1 maybe baseline become baseline)

# Plot the time point
which(ZirFlu_proteinSample$DatumProbe != ZirFlu_proteinSample$DatumEpisode)
which(ZirFlu_proteinSample$DatumProbe != ZirFlu_proteinSample$DatumBasis)
all.equal(ZirFlu_proteinSample$DatumProbe, ZirFlu_proteinSample$probeDatum) # TRUE

ZirFlu_timePoint <- ZirFlu_proteinSample %>% 
  group_by(ProbandenID, Season) %>% arrange(Zeitpunkt_v2) %>% 
  mutate(Zeitpunkt_timeDiff = difftime(DatumProbe, lag(DatumProbe), units = "days")) %>% # have 1 negative value
  group_by(ProbandenID, Season) %>% arrange(time_proteinPlate) %>% 
  mutate(proteinPlate_timeDiff = difftime(DatumProbe, lag(DatumProbe), units = "days")) # all seems to be correct

missmatch_timeRange <- ZirFlu_timePoint %>%
  slice(which(Zeitpunkt_timeDiff != proteinPlate_timeDiff)) # du to the set baseline timepoint for the Z-01-99-085 donor.

