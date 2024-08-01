rm(list = ls())

library(tidyverse)
library(rstatix)
library(ggpubr)

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

# remove duplicated donor sample
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

# CD83 data ------------------------------------------------------
cd83_dat <- updated_LateSeason_MetaData_deduplicatedMerged_dfFlt %>% 
  filter(Assay == "CD83") %>% rename("CD83" = NPX) %>% 
  select(NSC, season, sex..1..F..0..M., age, TimePoint, CD83) %>%
  mutate(TimePoint = ifelse(TimePoint == "V1", "d0",
                            ifelse(TimePoint == "V2", "d7",
                                   ifelse(TimePoint == "V3", "d28", TimePoint)))) 
# CD83 boxtplot at baseline -----------------------------------------------------
inputDat <- RAcohort_reclassify %>% 
  select(NSC, season, sex..1..F..0..M., age, condition, matches("_reclassify")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season)) %>%
  mutate(reclassify = factor(reclassify, levels = c("LL", "LH", "HL", "HH")),
         abFC = ifelse(reclassify == "LL" |reclassify == "HL", "NR", "R"),
         strainSeason = factor(strainSeason, 
                               levels = c("H1N1_Michigan_2018/2019", "H3N2_Singapore_2018/2019",
                                          "B_ColoradoVic_2018/2019", "B_PhuketYam_2018/2019",
                                          "H1N1_Michigan_2020/2021", "H3N2_Singapore_2020/2021", 
                                          "B_ColoradoVic_2020/2021", "B_PhuketYam_2020/2021"))) %>%
  left_join(cd83_dat %>% filter(TimePoint == "d0"))

protein <- "CD83"


## compare_reClass -----------------------------------------------------
compare_reClass <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))

inputDat %>% 
  ggboxplot(x = "reclassify", y = protein,
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 2) +
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

inputDat %>% filter(condition == "Control") %>%
  ggboxplot(x = "reclassify", y = protein,
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 2) +
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

inputDat %>% filter(condition == "Patient") %>%
  ggboxplot(x = "reclassify", y = protein,
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 2) +
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

## compare_abFC -----------------------------------------------------
compare_abFC <- list(c("R", "NR"))

inputDat %>% 
  ggboxplot(x = "abFC", y = protein,
            paletter = "jco", add = "jitter") + 
  facet_grid(season~strain) +
  stat_compare_means(comparisons = compare_abFC , method = "t.test")

inputDat %>% filter(condition == "Control") %>%
  ggboxplot(x = "abFC", y = protein,
            paletter = "jco", add = "jitter") + 
  facet_grid(season~strain) +
  stat_compare_means(comparisons = compare_abFC , method = "t.test")

inputDat %>% filter(condition == "Patient") %>%
  ggboxplot(x = "abFC", y = protein,
            paletter = "jco", add = "jitter") + 
  facet_grid(season~strain) +
  stat_compare_means(comparisons = compare_abFC , method = "t.test")


# CD83 overtime: d0 vs. d28 -----------------------------------------------------
inputDat <- RAcohort_reclassify %>% 
  select(NSC, season, sex..1..F..0..M., age, condition, matches("_reclassify")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season)) %>%
  mutate(reclassify = factor(reclassify, levels = c("LL", "LH", "HL", "HH")),
         strainSeason = factor(strainSeason, 
                               levels = c("H1N1_Michigan_2018/2019", "H3N2_Singapore_2018/2019",
                                          "B_ColoradoVic_2018/2019", "B_PhuketYam_2018/2019",
                                          "H1N1_Michigan_2020/2021", "H3N2_Singapore_2020/2021", 
                                          "B_ColoradoVic_2020/2021", "B_PhuketYam_2020/2021"))) %>%
  left_join(cd83_dat %>% filter(TimePoint %in% c("d0", "d28"))) %>% 
  mutate(TimePoint =  factor(TimePoint, levels = c("d0","d28")))

# plot data
inputDat2 <- inputDat %>% add_count(NSC, strainSeason) %>% filter(n == 2)

stat.test <- inputDat2 %>% 
  group_by(strainSeason, reclassify) %>%  
  t_test(CD83 ~ TimePoint) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(reclassify == "LL", 0.9, 
                       ifelse(reclassify == "LH", 1.9,
                              ifelse(reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(reclassify == "LL", 1.1, 
                       ifelse(reclassify == "LH", 2.1,
                              ifelse(reclassify == "HL", 3.1, 4.1))))

bxp <- inputDat %>% 
  ggboxplot(x = "reclassify", y = protein, color = "TimePoint",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 2)
bxp + stat_pvalue_manual(stat.test, label = "p.signif") + theme(legend.position = "top")

# plot only RA patients
inputDat_RA <- inputDat %>% filter(condition == "Patient")
inputDat_RA2 <- inputDat_RA %>% add_count(NSC, strainSeason) %>% filter(n == 2)

stat.test <- inputDat_RA2 %>%  #need to check
  group_by(strainSeason, reclassify) %>%  
  t_test(CD83 ~ TimePoint) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(reclassify == "LL", 0.9, 
                       ifelse(reclassify == "LH", 1.9,
                              ifelse(reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(reclassify == "LL", 1.1, 
                       ifelse(reclassify == "LH", 2.1,
                              ifelse(reclassify == "HL", 3.1, 4.1))))

bxp <- inputDat_RA %>% 
  ggboxplot(x = "reclassify", y = protein, color = "TimePoint",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 2)
bxp + stat_pvalue_manual(stat.test, label = "p.signif") + theme(legend.position = "top")

# CD83 overtime: d0 vs. d7 -----------------------------------------------------
inputDat <- RAcohort_reclassify %>% 
  select(NSC, season, sex..1..F..0..M., age, condition, matches("_reclassify")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season)) %>%
  mutate(reclassify = factor(reclassify, levels = c("LL", "LH", "HL", "HH")),
         strainSeason = factor(strainSeason, 
                               levels = c("H1N1_Michigan_2018/2019", "H3N2_Singapore_2018/2019",
                                          "B_ColoradoVic_2018/2019", "B_PhuketYam_2018/2019",
                                          "H1N1_Michigan_2020/2021", "H3N2_Singapore_2020/2021", 
                                          "B_ColoradoVic_2020/2021", "B_PhuketYam_2020/2021"))) %>%
  left_join(cd83_dat %>% filter(TimePoint %in% c("d0", "d7"))) %>% 
  mutate(TimePoint =  factor(TimePoint, levels = c("d0","d7"))) %>% drop_na()

# plot data 
inputDat2 <- inputDat %>% add_count(NSC, strainSeason) %>% filter(n == 2)

stat.test <- inputDat2 %>% 
  group_by(strainSeason, reclassify) %>%  
  t_test(CD83 ~ TimePoint) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(reclassify == "LL", 0.9, 
                       ifelse(reclassify == "LH", 1.9,
                              ifelse(reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(reclassify == "LL", 1.1, 
                       ifelse(reclassify == "LH", 2.1,
                              ifelse(reclassify == "HL", 3.1, 4.1))))

bxp <- inputDat %>% 
  ggboxplot(x = "reclassify", y = protein, color = "TimePoint",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 2)
bxp + stat_pvalue_manual(stat.test, label = "p.signif") + theme(legend.position = "top")
