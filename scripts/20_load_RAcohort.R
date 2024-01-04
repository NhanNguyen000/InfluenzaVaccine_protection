rm(list = ls())

library(tidyverse)
library(rstatix)

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
View(RAcohort_info %>% group_by(season) %>% count(NSC))

# test duplicated donors
a <- RAcohort_info %>% group_by(season, NSC) %>% count(condition) # 233 identical donors, 32 donors appears in 2 season --> 265 donors per season
length(which(duplicated(a$NSC) == TRUE))
b <- a[(which(duplicated(a$NSC) == TRUE)), ]
table(b$condition)

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
intersect(RAcohort_wide$NSC, RAcohort_info$NSC) # still 233 identical donors

# duplicate donor
RAcohort_check <- RAcohort_wide %>% group_by(season) %>% add_count(NSC)
unique(RAcohort_check$n)

unique((RAcohort_check %>% filter(n > 1))$NSC)
#[1] "300415"  "2352537" "1885881"

duplicated_samples <- RAcohort_long %>% filter(NSC %in% unique((RAcohort_check %>% filter(n > 1))$NSC))
#save(duplicated_samples, file = "sendOut/RAcohort_duplicatedSample_20231415.RData")
unique(duplicated_samples$NSC)
# [1] "300415"  "2352537" "1885881"


# remove dublicated donor sample
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

#save(RAcohort_reclassify, file = "sendOut/RAcohort_reclassify_20231218.RData")
table(RAcohort_reclassify$season, RAcohort_reclassify$condition)

RAcohort_reclassify %>% group_by(condition, season) %>% 
  drop_na(H1N1_Michigan_reclassify) %>%
  count(H1N1_Michigan_reclassify)

RAcohort_reclassify %>% group_by(condition, season) %>% 
  drop_na(H3N2_Singapore_reclassify) %>%
  count(H3N2_Singapore_reclassify)

RAcohort_reclassify %>% group_by(condition, season) %>% 
  drop_na(B_ColoradoVic_reclassify) %>%
  count(B_ColoradoVic_reclassify)

RAcohort_reclassify %>% group_by(condition, season) %>% 
  drop_na(B_PhuketYam_reclassify) %>%
  count(B_PhuketYam_reclassify)

RAcohort_reclassify %>% drop_na(B_PhuketYam_reclassify) %>%
  filter(condition == "Control") %>%
  group_by(season) %>% 
  count(B_PhuketYam_reclassify)

dat_temp <- RAcohort_reclassify %>% filter(condition == "Control")
table(dat_temp$season, dat_temp$H1N1_Michigan_reclassify)
table(dat_temp$season, dat_temp$H3N2_Singapore_reclassify)
table(dat_temp$season, dat_temp$B_ColoradoVic_reclassify)
table(dat_temp$season, dat_temp$B_PhuketYam_reclassify)

dat_temp <- RAcohort_reclassify %>% filter(condition == "Patient")
table(dat_temp$season, dat_temp$H1N1_Michigan_reclassify)
table(dat_temp$season, dat_temp$H3N2_Singapore_reclassify)
table(dat_temp$season, dat_temp$B_ColoradoVic_reclassify)
table(dat_temp$season, dat_temp$B_PhuketYam_reclassify)

# plot the violine plot of antibody titer ------------------------------------------------------
strain <- "H1N1_Michigan"
strain <- "H3N2_Singapore"
strain <- "B_ColoradoVic"
strain <- "B_PhuketYam"
year <- "2018/2019"
year <- "2020/2021"

strain_dat <- RAcohort_reclassify %>% 
  filter(season == year) %>%
  dplyr::select(NSC, condition, season, matches(strain)) %>%
  gather(matches("_d"), key = time, value = HAItiter) %>%
  rename_at(vars(ends_with("_abFC")), ~ "abFC") %>%
  rename_at(vars(ends_with("_reclassify")), ~ "reclassify") %>%
  mutate(category = ifelse(abFC >= 4, "Response", "Non-response"),
         time = gsub(paste0(strain, "_"), "", time),
         HAItiter = ifelse(HAItiter == 0, 0, log2(HAItiter))) %>% 
  mutate_at(vars(contains("time")), ~factor(.x, levels = c("d0", "d7", "d28"))) %>%
  drop_na()

plotDat_NRvsR <- strain_dat  %>%
  ggplot(aes(time, HAItiter)) + 
  geom_violin(width = 0.9) +
  geom_boxplot(width = 0.2, color="gray41") + 
  geom_point() +
  geom_line(aes(group = NSC), color="grey70") + 
  facet_wrap(vars(category), nrow = 1) +
  theme_bw() +
  ylab(paste0("Log2 of HAI titer - ", strain))

plotDat_NRvsR

plotDat_4groups <- strain_dat  %>%
  ggplot(aes(time, HAItiter)) + 
  geom_violin(width = 0.9) +
  geom_boxplot(width = 0.2, color="gray41") + 
  geom_point() +
  geom_line(aes(group = NSC), color="grey70") + 
  facet_wrap(vars(reclassify), nrow = 1) +
  theme_bw() +
  ylab(paste0("Log2 of HAI titer - ", strain))

plotDat_4groups 
