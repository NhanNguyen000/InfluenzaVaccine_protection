rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data =======================================================================
load("processedDat/cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy) %>% 
  full_join(mebo_Dat %>% 
              lapply(function(x) x %>% rownames_to_column("name")) %>%
              purrr::reduce(full_join) %>% right_join(metadata_healthy))


# scater plot protein / antibody and metabolite---------------------------
mebo <- "C3H7NO2S"
mebo <- "C18H32O2"
mebo <- "C20H32O2"
mebo <- "C22H36O2"
protein <- "IL6"
protein <- "IL10"
protein <- "CD83"

inputDat %>% 
  ggplot(aes_string(x = protein, y = mebo)) +
  geom_point() + facet_wrap(~season) + 
  stat_cor(label.y = 13) + 
  geom_smooth( method = "lm", se = FALSE) + 
  theme_classic() 

inputDat %>% filter(season == "2015") %>%
  ggplot(aes_string(x = protein, y = mebo)) +
  geom_point() + facet_wrap(~H1N1_reclassify) + 
  stat_cor(label.y = 13) + 
  geom_smooth( method = "lm", se = FALSE) + 
  theme_classic() 


inputDat %>% 
  ggplot(aes_string(x = protein, y = mebo)) +
  geom_point() + facet_wrap(~season) + 
  stat_cor(method = "spearman", label.y = 13) + 
  geom_smooth( method = "lm", se = FALSE) + 
  theme_classic() 

ab <- "H1N1_d0_log2"
ab <- "H1N1_abFC_log2"
inputDat %>% 
  ggplot(aes_string(x = ab, y = mebo)) +
  geom_point() + facet_wrap(~season) + 
  stat_cor() + 
  geom_smooth( method = "lm", se = FALSE) + 
  theme_classic() 
