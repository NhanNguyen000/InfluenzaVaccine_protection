library(readxl)
get.log2 <- function(dat) {
  library(tidyverse)
  library(rstatix)
  outcome <- dat %>% mutate(across(where(is.numeric), ~log2(.x))) %>%
    mutate(across(where(is.numeric), ~ifelse(is.infinite(.x), 0, .x)))
  return(outcome)
}


# get MN titer data ---------------
MN_rawDat2019 <- read_excel("../metadata/20230105_ZirFlu-2019-2020_MN_Titer_Valerie&Nhan.xlsx",
                       sheet = "2019-2020 MN Titer")
MN_titer_2019 <- MN_rawDat2019 %>% 
  select(1:14) %>% slice(-c(1, 2, 36, 52)) %>%
  fill(1, .direction = "down") %>%
  as.data.frame() %>%
  mutate_at(c(3:14), as.numeric) %>% get.log2()
names(MN_titer_2019) <- c("condition", "patientID", 
                     paste0("H1N1_", c("T1", "T2", "T3")),
                     paste0("H3N2_", c("T1", "T2", "T3")),
                     paste0("Bvictoria_", c("T1", "T2", "T3")),
                     paste0("Byamagata_", c("T1", "T2", "T3")))

MNabFC_rawDat2019 <- read_excel("../metadata/20230105_ZirFlu-2019-2020_MN_Titer_Valerie&Nhan.xlsx",
                        sheet = "2019-2020 MN Fold-increase")
MN_abFC2019 <- MNabFC_rawDat2019 %>% 
  slice(-c(1, 2, 36, 52)) %>%
  fill(1, .direction = "down") %>%
  as.data.frame() %>%
  mutate_at(c(3:6), as.numeric) %>% get.log2()

names(MN_abFC2019) <- c("condition", "patientID",
                        paste0(c("H1N1_", "H3N2_", "Bvictoria_", "Byamagata_"), "abFC"))

MN_2019 <- MN_titer_2019 %>% full_join(MN_abFC2019) %>% 
  mutate(season = "2019",
         condition = str_to_lower(condition)) %>%
  relocate(season)
names(MN_2019)[4:19] <- paste0("MN_", names(MN_2019)[4:19])

# compare HAI and MN value (abFC)--------------
cols <- c("LL" = "red", "LH" = "blue", "HL" = "darkgreen", "HH" = "orange")

HAIvsMN_ZirFlu_2019 <- MN_2019 %>% left_join(metadat)

HAIvsMN_ZirFlu_2019 %>% count(H1N1_reclassify)

HAIvsMN_ZirFlu_2019 %>%
  ggplot(aes(x = H1N1_abFC, y = MN_H1N1_abFC, col = H1N1_reclassify)) + 
  geom_point(size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top")

cor.test(HAIvsMN_ZirFlu_2019$H1N1_abFC, HAIvsMN_ZirFlu_2019$MN_H1N1abFC)

HAIvsMN_ZirFlu_2019 %>%
  ggplot(aes(x = H1N1_abFC, y = MN_H1N1_abFC)) + 
  geom_point(aes(col = H1N1_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

HAIvsMN_ZirFlu_2019 %>%
  ggplot(aes(x = H3N2_abFC, y = MN_H3N2_abFC)) + 
  geom_point(aes(col = H3N2_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)


HAIvsMN_ZirFlu_2019 %>%
  ggplot(aes(x = Bvictoria_abFC, y = MN_Bvictoria_abFC)) + 
  geom_point(aes(col = Bvictoria_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

HAIvsMN_ZirFlu_2019 %>%
  ggplot(aes(x = Byamagata_abFC, y = MN_Byamagata_abFC)) + 
  geom_point(aes(col = Byamagata_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

# compare HAI and MM titer at T1 -----------------------
HAIvsMN_ZirFlu_2019 %>%
  ggplot(aes(x = H1N1_T1, y = MN_H1N1_T1)) + 
  geom_point(aes(col = H1N1_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

HAIvsMN_ZirFlu_2019 %>%
  ggplot(aes(x = H3N2_T1, y = MN_H3N2_T1)) + 
  geom_point(aes(col = H3N2_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

HAIvsMN_ZirFlu_2019 %>%
  ggplot(aes(x = Bvictoria_T1, y = MN_Bvictoria_T1)) + 
  geom_point(aes(col = Bvictoria_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

HAIvsMN_ZirFlu_2019 %>%
  ggplot(aes(x = Byamagata_T1, y = MN_Byamagata_T1)) + 
  geom_point(aes(col = Byamagata_reclassify),
             size = 2, 
             position = position_jitter()) + theme_bw() +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme(legend.position = "top") +
  scale_color_manual(values = cols)

# boxplot of MN value with HAI reclassifycation -----------------
get.boxplot_reclassify_perStrain <- function(dat, strain_list) {
  # Description: boxplot for 4 reclassified vaccine response groups per strain
  
  # Arguments: 
  # dat: data with patient ID, HAI values, abFC of HAI, vaccine reclassify per column
  # strain_list: names of strain in the season
  # baseline_cutoff, abFC_cutoff: threshold for the classification
  
  # Returns: Boxplot of 4 reclassified groups per strain, with each dot is patients, 
  # and dot connection project the change of HAI level in patient 
  # 4 reclassified groups: HH - high baseline & high abFC, HL - high baseline & log abFC, 
  # LH - low baseline & high abFC, LL - low baseline & low abFC
  
  plot_list <- list()
  
  for (strain in strain_list) {
    plotDat <- dat %>% dplyr::select(patientID, matches(strain), -matches("ab")) %>%
      gather(matches("_T"), key = time, value = HAItiter) %>%
      arrange(patientID) %>%
      rename_at(vars(ends_with("_reclassify")), ~ "reclassify")
    
    plot_list[[strain]] <- plotDat %>%
      ggplot(aes(time, HAItiter)) + geom_boxplot() + geom_point() +
      geom_line(aes(group = patientID)) + theme_bw() +
      facet_wrap(vars(reclassify), nrow = 1)
  }
  
  return(plot_list)
}

plotDat <- HAIvsMN_ZirFlu_2019 %>%
  select(patientID, MN_H1N1_T1, MN_H1N1_T2, H1N1_reclassify) %>%
  gather(matches("_T"), key = time, value = MNtiter) %>% 
  arrange(patientID) %>%
  rename_at(vars(ends_with("_reclassify")), ~"reclassify") %>%
  mutate(time = str_replace(time, "MN_", ""),
         reclassify = factor(reclassify, levels = c("HH", "HL", "LH", "LL")))

plotDat %>% ggplot(aes(time, MNtiter)) + 
  geom_boxplot() + geom_point() +
  geom_line(aes(group = patientID)) + theme_bw() +
  facet_wrap(vars(reclassify), nrow = 1) +
  scale_y_continuous(breaks = seq(2, 15, 2.5), limits = c(2, 15))

