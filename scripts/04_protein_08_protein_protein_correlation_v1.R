rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data =======================================================================
load("cohorts_dat.RData")

## metadata for all healthy subjects 
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(age_v2 = ifelse(age >=64, "old", "young"))

inputDat <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

# CHECK PROTEINS =======================================================================
protein <- c("CD83", "IL15")
plotDat <- inputDat %>% 
  select(season, responder, age_v2, c(protein), matches("_abFC|_T1|_T4|_reclassify")) %>%
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_T1 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high")))

# all donors and all seasons in 1 plot -------------------------
plotDat %>% 
  ggplot() +
  geom_jitter(aes(x = CD83, y = IL15), width = 0.2, size = 1.2) +
  geom_smooth(aes(x = CD83, y = IL15), method = "lm", se = FALSE) +
  stat_cor(aes(x = CD83, y = IL15)) +
  ylim(0, 3) + theme_classic()

# donors per season  -------------------------
plotDat %>% 
  ggplot() +
  geom_jitter(aes(x = CD83, y = IL15), width = 0.2, size = 1.2) + 
  facet_wrap(~season, nrow = 1) +
  geom_smooth(aes(x = CD83, y = IL15), method = "lm", se = FALSE) +
  stat_cor(aes(x = CD83, y = IL15)) +
  theme_classic()

# donors per season and group  -------------------------
plotDat %>% 
  ggplot(aes(x = CD83, y = IL15, color = H1N1_reclassify)) +
  geom_jitter(width = 0.2, size = 1.2) + 
  facet_wrap(~season, nrow = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor() +
  theme_classic()

plotDat %>% 
  ggplot(aes(x = CD83, y = IL15, color = H1N1_abFC)) +
  geom_jitter(width = 0.2, size = 1.2) + 
  facet_wrap(~season, nrow = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor() +
  theme_classic()

plotDat %>% 
  ggplot(aes(x = CD83, y = IL15, color = H1N1_abBaseline)) +
  geom_jitter(width = 0.2, size = 1.2) + 
  facet_wrap(~season, nrow = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor() +
  theme_classic()

plotDat %>% 
  ggplot(aes(x = CD83, y = IL15, color = age_v2)) +
  geom_jitter(width = 0.2, size = 1.2) + 
  facet_wrap(~season, nrow = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor() +
  theme_classic()
