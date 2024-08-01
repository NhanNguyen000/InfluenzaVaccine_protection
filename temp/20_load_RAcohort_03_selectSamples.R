rm(list = ls())

library(tidyverse)
library(rstatix)
library(openxlsx)

convert.titer <- function(input_titers) {
  titers <- ifelse(input_titers == "<10", 0,
                   ifelse(input_titers == ">1280", "1280", input_titers)) %>% 
    as.numeric()
  return(titers)
}

get.abFC <- function(titer_baseline, titer_postVac) {
  abFC <- ifelse(titer_baseline == 0, titer_postVac, titer_postVac/titer_baseline)
  return(abFC)
}


get.reclassify <- function(dat, baseline_col, abFC_col, baseline_cutoff, abFC_cutoff) {
  # Description: reclassify vaccine response based on the baseline and abFC
  
  # Arguments: 
  # dat: data with patient ID, HAI value per time point, abFC of HAI between time point per column
  # baseline_col, abFC_col: names of the baseline and abFC columns
  # baseline_cutoff, abFC_cutoff: threshold for the classification
  
  # Returns: 
  # 4 reclassified groups: HH - high baseline & high abFC, 
  # HL - high baseline & log abFC, LH - low baseline & high abFC, LL - low baseline & low abFC
  
  reclassify_temp <- paste0(ifelse(dat[, baseline_col] >= baseline_cutoff, "H", "L"),
                            ifelse(dat[, abFC_col] >= abFC_cutoff, "H", "L"))
  reclassify <- ifelse(grepl("NA", reclassify_temp), NA, reclassify_temp)
  return(reclassify)
}

# load data ==============================================================================
load("/vol/projects/skumar/influenza_RA/RA_LateSeasons_Antibody_OLINK_v1.RData" )
# This is 1 dataframe with antibody titres at each time point, and OLINK dataframe filtered after 1st QC as Martijn did. Once I fit other 20 donors somewhere in here (or not), I will update you

RAcohort_info <- updated_LateSeason_MetaData_deduplicatedMerged_dfFlt %>% 
  select(Sample_ID, Donor_Date, bioBank_olink_NSC, NSC, Sample_Number, 
         season, condition, sex..1..F..0..M., age, TimePoint,
         B.Colorado.B.Vic, B.Phuket.B.Yam, A.Michigan.A.H1N1.pdm09, A.Singapore.A.H3N2.) %>%
  distinct() # 721 samples
length(unique(RAcohort_info$NSC)) # 233 donors across season

# adjust the titer and donor condition information
RAcohort_long <- RAcohort_info %>%
  mutate(H1N1_Michigan = convert.titer(A.Michigan.A.H1N1.pdm09), 
         H3N2_Singapore = convert.titer(A.Singapore.A.H3N2.),
         B_ColoradoVic = convert.titer(B.Colorado.B.Vic), 
         B_PhuketYam = convert.titer(B.Phuket.B.Yam)) %>%
  mutate(condition = ifelse(condition == "RA", "Patient", 
                            ifelse(condition == "Healthy", "Control", condition)))

RAcohort_wide <- RAcohort_long %>% 
  select(-Sample_ID, -Donor_Date, -bioBank_olink_NSC, -Sample_Number,
         -A.Michigan.A.H1N1.pdm09, -A.Singapore.A.H3N2.,
         -B.Colorado.B.Vic, -B.Phuket.B.Yam) %>%
  mutate(TimePoint = ifelse(TimePoint == "V1", "d0",
                            ifelse(TimePoint == "V2", "d7",
                                   ifelse(TimePoint == "V3", "d28", TimePoint)))) %>%
  pivot_wider(names_from = TimePoint, 
              values_from = c(H1N1_Michigan, H3N2_Singapore, B_ColoradoVic, B_PhuketYam)) # 251 rows, not 265/269 rows after using only Control vs. Patient for donor's condition
#intersect(RAcohort_wide$NSC, RAcohort_info$NSC) # still 233 identical donors

# remove dublicated donor sample
RAcohort_check <- RAcohort_wide %>% group_by(season) %>% add_count(NSC)
RAcohort_check_v2 <- RAcohort_check %>% filter(n == 1) # 243 donors

# calculate abFC
RAcohort_abFC <- RAcohort_check_v2 %>% select(-n) %>%
  mutate(H1N1_Michigan_abFC = get.abFC(H1N1_Michigan_d0, H1N1_Michigan_d28), 
         H3N2_Singapore_abFC = get.abFC(H3N2_Singapore_d0, H3N2_Singapore_d28),
         B_ColoradoVic_abFC = get.abFC(B_ColoradoVic_d0, B_ColoradoVic_d28), 
         B_PhuketYam_abFC = get.abFC(B_PhuketYam_d0, B_PhuketYam_d28))

# reclassification ------------------------------------------------------------------------
## cutoff 
baseline <- 40
abFC <- 4

strains <- c("H1N1_Michigan", "H3N2_Singapore", "B_ColoradoVic", "B_PhuketYam")
reclassify_temp <- list()
for (strain in strains) {
  reclassify_temp[[strain]] <- get.reclassify(dat = RAcohort_abFC, 
                                              baseline_col = paste0(strain, "_d0"),
                                              abFC_col = paste0(strain, "_abFC"),
                                              baseline_cutoff = baseline,
                                              abFC_cutoff = abFC)
  
}

RAcohort_reclassify <-  RAcohort_abFC %>% 
  cbind(as.data.frame(reclassify_temp) %>% 
          rename_with(~paste0(.x, "_reclassify")))  %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

rm(RAcohort_abFC, RAcohort_check, RAcohort_check_v2, RAcohort_long, RAcohort_wide)

# prepare the data (reclassification, antibody titer, CD83 level) ------------------------------------------------------
cd83_V1 <- updated_LateSeason_MetaData_deduplicatedMerged_dfFlt %>% 
  filter(Assay == "CD83", TimePoint == "V1") %>% rename("CD83" = NPX) %>% 
  select(NSC, season, sex..1..F..0..M., age, CD83) # CD83 before vaccination

RA_dat <- RAcohort_info %>% 
  select(-Sample_ID, -Donor_Date, -bioBank_olink_NSC, -Sample_Number) %>%
  mutate(condition = ifelse(condition == "RA", "Patient", 
                            ifelse(condition == "Healthy", "Control", condition)))  %>%
  pivot_wider(names_from = TimePoint, 
              values_from = c(A.Michigan.A.H1N1.pdm09, A.Singapore.A.H3N2.,
                              B.Colorado.B.Vic, B.Phuket.B.Yam)) %>%
  right_join(RAcohort_reclassify %>% select(-matches("d0|d7|d28"))) %>% 
  left_join(cd83_V1)

rm(cd83_V1)

# donors same reclassified groups across 4 strains ----------------------------------
donors_sameGroup <- RA_dat %>% 
  slice(which(H1N1_Michigan_reclassify == H3N2_Singapore_reclassify)) %>% 
  slice(which(H1N1_Michigan_reclassify == B_ColoradoVic_reclassify)) %>% 
  slice(which(H1N1_Michigan_reclassify == B_PhuketYam_reclassify))

table(donors_sameGroup$H1N1_Michigan_reclassify, donors_sameGroup$condition)
write.xlsx(donors_sameGroup, file = "sendOut/donors_sameGroup.xlsx")

rm(donors_sameGroup)
# LL, LH, HH donors per strain, with antibody titer  ----------------------------------
## LL group per strain ------------------------------------------------------
# LL donors in H1N1, titer d0
LL_H1N1_Michigan <- RA_dat %>% 
  filter(H1N1_Michigan_reclassify == "LL", A.Michigan.A.H1N1.pdm09_V1 %in% c("<10", "10")) 

# LL donors in H3N2, titer d0
LL_H3N2_Singapore <- RA_dat %>% 
  filter(H3N2_Singapore_reclassify == "LL", A.Singapore.A.H3N2._V1 %in% c("<10", "10")) 

# LL donors in B_ColoradoVic, titer d0
LL_B_ColoradoVic <- RA_dat %>% 
  filter(B_ColoradoVic_reclassify == "LL", B.Colorado.B.Vic_V1 %in% c("<10")) 

# LL donors in B_PhuketYam, titer d0
LL_B_PhuketYam <- RA_dat %>% 
  filter(B_PhuketYam_reclassify == "LL", B.Phuket.B.Yam_V1 %in% c("<10")) 

## LH group per strain ------------------------------------------------------
# LH donors in H1N1, tier d0 and titer d28
LH_H1N1_Michigan <- RA_dat %>% 
  filter(H1N1_Michigan_reclassify == "LH", 
         A.Michigan.A.H1N1.pdm09_V3 %in% c("320", "1280", ">1280"))

# LH donors in H3N2, tier d0 and titer d28
LH_H3N2_Singapore <- RA_dat %>% 
  filter(H3N2_Singapore_reclassify == "LH", 
         A.Singapore.A.H3N2._V3 %in% c("320", "1280", ">1280"))

# LH donors in B_ColoradoVic, tier d0 and titer d28
LH_B_ColoradoVic <- RA_dat %>% 
  filter(B_ColoradoVic_reclassify == "LH", 
         B.Colorado.B.Vic_V3 %in% c("160", "320", "1280", ">1280")) # lower antibody titer cutoff to 160, to have enough donors

# LH donors in B_PhuketYam, tier d0 and titer d28
LH_B_PhuketYam <- RA_dat %>% 
  filter(B_PhuketYam_reclassify == "LH", 
         B.Phuket.B.Yam_V3 %in% c("160", "320", "1280", ">1280")) # lower antibody titer cutoff to 160, to have enough donors

## HH group per strain ------------------------------------------------------
# HH donors in H1N1, tier d0 and titer d28
HH_H1N1_Michigan <- RA_dat %>% 
  filter(H1N1_Michigan_reclassify == "HH", 
         A.Michigan.A.H1N1.pdm09_V3 %in% c("320", "1280", ">1280"))

# HH donors in H3N2, tier d0 and titer d28
HH_H3N2_Singapore <- RA_dat %>% 
  filter(H3N2_Singapore_reclassify == "HH", 
         A.Singapore.A.H3N2._V3 %in% c("320", "1280", ">1280"))

# HH donors in B_ColoradoVic, tier d0 and titer d28
HH_B_ColoradoVic <- RA_dat %>% 
  filter(B_ColoradoVic_reclassify == "HH") # only 2 HH donors in this strain

# HH donors in B_PhuketYam, tier d0 and titer d28
HH_B_PhuketYam <- RA_dat %>% 
  filter(B_PhuketYam_reclassify == "HH") # only 2 HH donors in this strain

## save file ------------------------------------------------------
donors_titer_baseline <- list('LL_H1N1_Michigan' = LL_H1N1_Michigan, 'LL_H3N2_Singapore' = LL_H3N2_Singapore, 
              'LL_B_ColoradoVic' = LL_B_ColoradoVic, "LL_B_PhuketYam" = LL_B_PhuketYam,
              'LH_H1N1_Michigan' = LH_H1N1_Michigan, 'LH_H3N2_Singapore' = LH_H3N2_Singapore, 
              'LH_B_ColoradoVic' = LH_B_ColoradoVic, "LH_B_PhuketYam" = LH_B_PhuketYam,
              'HH_H1N1_Michigan' = HH_H1N1_Michigan, 'HH_H3N2_Singapore' = HH_H3N2_Singapore, 
              'HH_B_ColoradoVic' = HH_B_ColoradoVic, "HH_B_PhuketYam" = HH_B_PhuketYam)

write.xlsx(donors_titer_baseline, file = "sendOut/donors_titer_baseline.xlsx")

rm(donors_titer_baseline, 
   LL_B_ColoradoVic, LL_B_PhuketYam, LL_H1N1_Michigan, LL_H3N2_Singapore,
   LH_B_ColoradoVic, LH_B_PhuketYam, LH_H1N1_Michigan, LH_H3N2_Singapore,
   HH_B_ColoradoVic, HH_B_PhuketYam, HH_H1N1_Michigan, HH_H3N2_Singapore)

# LL, LH, HH donors per strain, with CD83 level  ----------------------------------
## LL group per strain ------------------------------------------------------
# LL donors in H1N1, 10 donors with lowest level of CD83 at baseline
LL_H1N1_Michigan <- RA_dat %>% 
  filter(H1N1_Michigan_reclassify == "LL") %>% slice_min(CD83, n = 10)

# LL donors in H3N2, 10 donors with lowest level of CD83 
LL_H3N2_Singapore <- RA_dat %>% 
  filter(H3N2_Singapore_reclassify == "LL") %>% slice_min(CD83, n = 10)

# LL donors in B_ColoradoVic, 10 donors with lowest level of CD83 
LL_B_ColoradoVic <- RA_dat %>% 
  filter(B_ColoradoVic_reclassify == "LL") %>% slice_min(CD83, n = 10)

# LL donors in B_PhuketYam, 10 donors with lowest level of CD83 
LL_B_PhuketYam <- RA_dat %>% 
  filter(B_PhuketYam_reclassify == "LL") %>% slice_min(CD83, n = 10) 

## LH group per strain ------------------------------------------------------
# LH donors in H1N1, 10 donors with highest level of CD83 at baseline
LH_H1N1_Michigan <- RA_dat %>% 
  filter(H1N1_Michigan_reclassify == "LH") %>% slice_max(CD83, n = 10)

# LH donors in H3N2, tier d0 and titer d28
LH_H3N2_Singapore <- RA_dat %>% 
  filter(H3N2_Singapore_reclassify == "LH") %>% slice_max(CD83, n = 10)

# LH donors in B_ColoradoVic, tier d0 and titer d28
LH_B_ColoradoVic <- RA_dat %>% 
  filter(B_ColoradoVic_reclassify == "LH") %>% slice_max(CD83, n = 10)

# LH donors in B_PhuketYam, tier d0 and titer d28
LH_B_PhuketYam <- RA_dat %>% 
  filter(B_PhuketYam_reclassify == "LH") %>% slice_max(CD83, n = 10)

## HH group per strain ------------------------------------------------------
# HH donors in H1N1, 10 donors with highest level of CD83 at baseline
HH_H1N1_Michigan <- RA_dat %>% 
  filter(H1N1_Michigan_reclassify == "HH") %>% slice_max(CD83, n = 10)

# HH donors in H3N2, 10 donors with highest level of CD83 at baseline
HH_H3N2_Singapore <- RA_dat %>% 
  filter(H3N2_Singapore_reclassify == "HH") %>% slice_max(CD83, n = 10)

# HH donors in B_ColoradoVic, 10 donors with highest level of CD83 at baseline
HH_B_ColoradoVic <- RA_dat %>% 
  filter(B_ColoradoVic_reclassify == "HH") # only 2 HH donors in this strain

# HH donors in B_PhuketYam, 10 donors with highest level of CD83 at baseline
HH_B_PhuketYam <- RA_dat %>% 
  filter(B_PhuketYam_reclassify == "HH") # only 2 HH donors in this strain

## save file ------------------------------------------------------
donors_CD83_baseline <- list('LL_H1N1_Michigan' = LL_H1N1_Michigan, 'LL_H3N2_Singapore' = LL_H3N2_Singapore, 
                             'LL_B_ColoradoVic' = LL_B_ColoradoVic, "LL_B_PhuketYam" = LL_B_PhuketYam,
                             'LH_H1N1_Michigan' = LH_H1N1_Michigan, 'LH_H3N2_Singapore' = LH_H3N2_Singapore, 
                             'LH_B_ColoradoVic' = LH_B_ColoradoVic, "LH_B_PhuketYam" = LH_B_PhuketYam,
                             'HH_H1N1_Michigan' = HH_H1N1_Michigan, 'HH_H3N2_Singapore' = HH_H3N2_Singapore, 
                             'HH_B_ColoradoVic' = HH_B_ColoradoVic, "HH_B_PhuketYam" = HH_B_PhuketYam)

write.xlsx(donors_CD83_baseline, file = "sendOut/donors_CD83_baseline.xlsx")

rm(donors_CD83_baseline, 
   LL_B_ColoradoVic, LL_B_PhuketYam, LL_H1N1_Michigan, LL_H3N2_Singapore,
   LH_B_ColoradoVic, LH_B_PhuketYam, LH_H1N1_Michigan, LH_H3N2_Singapore,
   HH_B_ColoradoVic, HH_B_PhuketYam, HH_H1N1_Michigan, HH_H3N2_Singapore)
