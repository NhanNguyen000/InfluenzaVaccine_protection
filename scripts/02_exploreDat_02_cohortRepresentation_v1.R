rm(list = ls())
library(tidyverse)
library(ggpubr)
library(webr)

load("cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# donut chart -------------------------
PieDonut(test, aes(cohort, count = n))

PieDonut(metadata_healthy, aes(season, sex), 
         r0 = 0.5, r2 = 1.1, start = -120,
         title = "Distribution of gender per season")
