rm(list = ls())

library(tidyverse)
library(ggpubr)
library(rstatix)
# load data =======================================================================
load("cohorts_dat.RData")
load("sigMebo_dynamic_from30selectedMebo.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- mebo_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

# boxplot for all strains, time d0 vs. d7 ================================
mebo <- "C18H32O2"
mebo <- "C20H32O2"
mebo <- "C22H36O2"
mebo <- "C3H7NO2S"
mebo <- "C2H2O3"

metadat_boxplot <- inputDat %>% 
  select(probandID, season, responder, time, c(mebo),
         matches("_abFC|_d0|_d7|_reclassify")) %>% filter(time %in% c("d0", "d7")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season))

## boxplots for each timepoint, per group, per season -------------
metadat_boxplot %>% 
  filter(strain == "H1N1" | strainSeason %in% c("B_2014", "B_2015", "H3N2_2015", "Byamagata_2020")) %>%
  ggboxplot(x = "reclassify", y = mebo, color = "time",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 1)

## add non sig. paired t_test to the graph -------------
metadat_boxplot2 <- metadat_boxplot %>% 
  add_count(probandID, strainSeason) %>% filter(n == 2) %>% # only use donor have preand post data point
  filter(strain == "H1N1" | strainSeason %in% c("B_2014", "B_2015", "H3N2_2015", "Byamagata_2020")) %>% 
  group_by(strainSeason, reclassify) %>%  add_count(reclassify) %>% filter(nn >= 3) # only use comparision have >= 3 donors 

stat.test <- metadat_boxplot2 %>% 
  #t_test(C18H32O2 ~ time) %>% 
  #t_test(C20H32O2 ~ time) %>% 
  #t_test(C22H36O2 ~ time) %>%
  #t_test(C3H7NO2S ~ time) %>%
  t_test(C2H2O3 ~ time) %>%
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(reclassify == "LL", 0.9, 
                       ifelse(reclassify == "LH", 1.9,
                              ifelse(reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(reclassify == "LL", 1.1, 
                       ifelse(reclassify == "LH", 2.1,
                              ifelse(reclassify == "HL", 3.1, 4.1))))

bxp <- metadat_boxplot %>% 
  filter(strain == "H1N1" | strainSeason %in% c("B_2014", "B_2015", "H3N2_2015", "Byamagata_2020")) %>%
  ggboxplot(x = "reclassify", y = mebo, color = "time",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 2)
bxp + stat_pvalue_manual(stat.test, label = "p.signif") + theme(legend.position = "top")

# boxplot for all strains, time d0 vs. d28 ================================
mebo <- "C18H32O2"
mebo <- "C20H32O2"
mebo <- "C22H36O2"
mebo <- "C3H7NO2S"
mebo <- "C2H2O3"

metadat_boxplot <- inputDat %>% 
  select(probandID, season, responder, time, c(mebo),
         matches("_abFC|_d0|_d28|_reclassify")) %>% filter(time %in% c("d0", "d28")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season))

## boxplots for each timepoint, per group, per season -------------
metadat_boxplot %>% 
  filter(strain == "H1N1" | strainSeason %in% c("B_2014", "B_2015", "H3N2_2015", "Byamagata_2020")) %>%
  ggboxplot(x = "reclassify", y = mebo, color = "time",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 1)

## add non sig. paired t_test to the graph -------------
metadat_boxplot2 <- metadat_boxplot %>% 
  add_count(probandID, strainSeason) %>% filter(n == 2) %>% # only use donor have preand post data point
  filter(strain == "H1N1" | strainSeason %in% c("B_2014", "B_2015", "H3N2_2015", "Byamagata_2020")) %>% 
  group_by(strainSeason, reclassify) %>%  add_count(reclassify) %>% filter(nn >= 3) # only use comparision have >= 3 donors 

stat.test <- metadat_boxplot2 %>% 
  #t_test(C18H32O2 ~ time) %>% 
  #t_test(C20H32O2 ~ time) %>% 
  #t_test(C22H36O2 ~ time) %>%
  #t_test(C3H7NO2S ~ time) %>%
  t_test(C2H2O3 ~ time) %>%
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(reclassify == "LL", 0.9, 
                       ifelse(reclassify == "LH", 1.9,
                              ifelse(reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(reclassify == "LL", 1.1, 
                       ifelse(reclassify == "LH", 2.1,
                              ifelse(reclassify == "HL", 3.1, 4.1))))

bxp <- metadat_boxplot %>% 
  filter(strain == "H1N1" | strainSeason %in% c("B_2014", "B_2015", "H3N2_2015", "Byamagata_2020")) %>%
  ggboxplot(x = "reclassify", y = mebo, color = "time",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 2)
bxp + stat_pvalue_manual(stat.test, label = "p.signif") + theme(legend.position = "top")

