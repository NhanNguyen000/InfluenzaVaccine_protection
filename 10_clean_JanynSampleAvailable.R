library(readxl)
library(tidyverse)

# load ZirFlu samples ------------------------------
rawDat <- read_xlsx("../metadata/Auszug ZirFlu CentraXX 08082022.xlsx")
avai_samples <- rawDat %>%
  filter(ProbeOEName == "ZirFlu (P-2024-GAS)",
         ProbenArt == "Zellen, PBMC, vital [ZZZ(pbm)]",
         Projekt == "ZirFlu")  %>% # based on Janyn's filter to selected only ZirFlu participants
  mutate(season = ifelse(DatumProbe < "2020-06-01 UTC", "2019",
                         ifelse(DatumProbe < "2021-06-01 UTC", "2020", "2021"))) # the collection date of the sample
# The participants from the first season (2019-2020) donated the samples from end of 2019 until the beginning of 2020. 
# While the participants from the second season (2020-2021) have samples from end of 2020 until beginning of 2021. 
# The third season we didn’t analyze, because the cohort was too small.

avai_samples %>% count(season)

# matches with samples from previous ZirFlu data  ------------------------------
View(ZirFlu$donorSamples)

avaiSample_matchPrevious <- avai_samples %>% 
  filter(ProbandenID %in% ZirFlu$donorSamples$patientID) %>%
  mutate(time = ifelse(Zeitpunkt == "-", "T1",
                       ifelse(Zeitpunkt == "01", "T2",
                              ifelse(Zeitpunkt == "02", "T3",
                                     ifelse(Zeitpunkt == "03", "T4", NA))))) %>%
  filter(season %in% c("2019", "2020"))

avaiSample_matchPrevious %>% count(season)
avaiSample_matchPrevious %>% count(season, time)

avaiSample_matchPrevious_addCondition <- ZirFlu$donorInfo %>% distinct() %>%
  full_join(avaiSample_matchPrevious %>% select(ProbandenID, ProbenID, season, time), 
            by = c("patientID" = "ProbandenID", "season"))

# note: some samples appear in 1 season but not in the other  ------------------------------
sample_NAs <- avaiSample_matchPrevious_addCondition[which(is.na(avaiSample_matchPrevious_addCondition$condition)), ] 
matchedParticipant_avaiSamples <- avaiSample_matchPrevious_addCondition %>% drop_na()

participant_avaiSamples <- matchedParticipant_avaiSamples %>% select(-ProbenID) %>% distinct()

participant_avaiSamples %>% filter(season == "2019") %>% count(condition, time)
participant_avaiSamples %>% filter(season == "2020") %>% count(condition, time)

# check longitude available samples ----------------------
season_2019 <- participant_avaiSamples %>% 
  filter(season == "2019") %>%
  group_by(patientID, condition) %>% count()
a <- season_2019 %>% filter(condition == "decompensated cirrhosis")

# 2019
decomCirrhosis_2019_T2 <- participant_avaiSamples %>% 
  filter(season == "2019", time == "T2", condition == "decompensated cirrhosis") 
decomCirrhosis_2019_longitude <- participant_avaiSamples %>% 
  filter(season == "2019", patientID %in% decomCirrhosis_2019_T2$patientID)

decomCirrhosis_2019_T3 <- participant_avaiSamples %>%
  filter(season == "2019", time == "T3", condition == "decompensated cirrhosis")
decomCirrhosis_2019_3timePoint <- participant_avaiSamples %>%
  filter(season == "2019", time %in% c("T1", "T3", "T4"), patientID %in% decomCirrhosis_2019_T3$patientID)
a <- decomCirrhosis_2019_3timePoint %>% group_by(patientID) %>% count()

decomCirrhosis_2019_2timePoint <- participant_avaiSamples %>%
  filter(season == "2019", time %in% c("T1", "T3"), patientID %in% decomCirrhosis_2019_T3$patientID)
a <- decomCirrhosis_2019_2timePoint %>% group_by(patientID) %>% count()

comCirrhosis_2019_T2 <- participant_avaiSamples %>% 
  filter(season == "2019", time == "T2", condition == "compensated cirrhosis") 
comCirrhosis_2019_longitude <- participant_avaiSamples %>% 
  filter(season == "2019", patientID %in% comCirrhosis_2019_T2$patientID)
a <- comCirrhosis_2019_longitude %>% group_by(patientID) %>% count()

comCirrhosis_2019_T3 <- participant_avaiSamples %>%
  filter(season == "2019", time == "T3", condition == "compensated cirrhosis")
comCirrhosis_2019_3timePoint <- participant_avaiSamples %>%
  filter(season == "2019", time %in% c("T1", "T3", "T4"), patientID %in% comCirrhosis_2019_T3$patientID)
a <- comCirrhosis_2019_3timePoint %>% group_by(patientID) %>% count()

comCirrhosis_2019_2timePoint <- participant_avaiSamples %>%
  filter(season == "2019", time %in% c("T1", "T3"), patientID %in% comCirrhosis_2019_T3$patientID)
a <- comCirrhosis_2019_2timePoint %>% group_by(patientID) %>% count()

healthy_2019_T2 <- participant_avaiSamples %>%
  filter(season == "2019", time == "T2", condition == "healthy")
healthy_2019_longitude <- participant_avaiSamples %>%
  filter(season == "2019", patientID %in% healthy_2019_T2$patientID)
a <- healthy_2019_longitude %>% group_by(patientID) %>% count()

healthy_2019_T3 <- participant_avaiSamples %>%
  filter(season == "2019", time == "T3", condition == "healthy")
healthy_2019_3timepoint <- participant_avaiSamples %>%
  filter(season == "2019", time %in% c("T1", "T3", "T4"), patientID %in% healthy_2019_T3$patientID)
a <- healthy_2019_3timepoint %>% group_by(patientID) %>% count()

healthy_2019_2timepoint <- participant_avaiSamples %>%
  filter(season == "2019", time %in% c("T1", "T3"), patientID %in% healthy_2019_T3$patientID)
a <- healthy_2019_2timepoint %>% group_by(patientID) %>% count()

# 2020

decomCirrhosis_2020_T3 <- participant_avaiSamples %>%
  filter(season == "2020", time == "T3", condition == "decompensated cirrhosis")
decomCirrhosis_2020_3timePoint <- participant_avaiSamples %>%
  filter(season == "2020", time %in% c("T1", "T3", "T4"), patientID %in% decomCirrhosis_2020_T3$patientID)
a <- decomCirrhosis_2020_3timePoint %>% group_by(patientID) %>% count()

decomCirrhosis_2020_2timePoint <- participant_avaiSamples %>%
  filter(season == "2020", time %in% c("T1", "T3"), patientID %in% decomCirrhosis_2020_T3$patientID)
a <- decomCirrhosis_2020_2timePoint %>% group_by(patientID) %>% count()


comCirrhosis_2020_T2 <- participant_avaiSamples %>%
  filter(season == "2020", time == "T2", condition == "compensated cirrhosis")
comCirrhosis_2020_longitude <- participant_avaiSamples%>%
  filter(season == "2020", patientID %in% comCirrhosis_2020_T2$patientID)
a <- comCirrhosis_2020_longitude %>% group_by(patientID) %>% count()

comCirrhosis_2020_T3 <- participant_avaiSamples %>%
  filter(season == "2020", time == "T3", condition == "compensated cirrhosis")

comCirrhosis_2020_3timepoint <- participant_avaiSamples%>%
  filter(season == "2020", time %in% c("T1", "T3", "T4"), patientID %in% comCirrhosis_2020_T3$patientID)
a <- comCirrhosis_2020_3timepoint %>% group_by(patientID) %>% count()

comCirrhosis_2020_2timepoint <- participant_avaiSamples%>%
  filter(season == "2020", time %in% c("T1", "T3"), patientID %in% comCirrhosis_2020_T3$patientID)
a <- comCirrhosis_2020_2timepoint %>% group_by(patientID) %>% count()

healthy_2020_T2 <- participant_avaiSamples %>%
  filter(season == "2020", time == "T2", condition == "healthy")
healthy_2020_longitude <- participant_avaiSamples %>%
  filter(season == "2020", patientID %in% healthy_2020_T2$patientID)
a <- healthy_2020_longitude %>% group_by(patientID) %>% count()

healthy_2020_T3 <- participant_avaiSamples %>%
  filter(season == "2020", time == "T3", condition == "healthy")
healthy_2020_3timePoint <- participant_avaiSamples %>%
  filter(season == "2020", time %in% c("T1", "T3", "T4"), patientID %in% healthy_2020_T3$patientID)
a <- healthy_2020_3timePoint %>% group_by(patientID) %>% count()

healthy_2020_2timePoint <- participant_avaiSamples %>%
  filter(season == "2020", time %in% c("T1", "T3"), patientID %in% healthy_2020_T3$patientID)
a <- healthy_2020_2timePoint %>% group_by(patientID) %>% count()

# count the sex ------------
a <- participant_avaiSamples %>% filter(season == "2019") %>% count(condition, sex, time)

# save information as an excel file ------------------------------------------------
library(openxlsx)
dat <- list("participants_information" = participant_avaiSamples,
            "samplesFromUsedParticipants" = matchedParticipant_avaiSamples,
            "availableSamples_inZirFlu" = avai_samples)
write.xlsx(dat, file = "sendOut/ZirFlu_availableSamples_fromAuszugCentraXX08082022_NhanNguyen.xlsx")

interested_participants <- participant_avaiSamples %>% select(-time) %>% distinct()
dat <- list("interested_participants" = interested_participants,
            "participant_probenID" = matchedParticipant_avaiSamples,
            "availableSamples_inZirFlu" = avai_samples)
write.xlsx(dat, file = "sendOut/ZirFlu_interestedParticipant_fromAuszugCentraXX08082022_NhanNguyen.xlsx")


