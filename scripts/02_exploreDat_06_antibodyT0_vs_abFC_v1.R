rm(list = ls())
library(tidyverse)
library(ggpubr)
library(webr)
library(patchwork)

load("cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

metadata_healthy %>%
  ggplot() +
  geom_bar(aes(x = name, y = H1N1_T1_log2), stat = "identity", fill="cyan",colour="cyan") +
  geom_line(aes(x = name, y = H1N1_abFC_log2, group = 1), stat = "identity", color="red",size=1) 

a <- metadata_healthy %>%
  ggplot() +
  geom_bar(aes(x = name, y = 13, fill = H1N1_T1_log2), stat = "identity") +
  scale_fill_gradient2(low = "white", high = "red") +
  geom_line(aes(x = name, y = H1N1_abFC_log2, group = 1), stat = "identity", color="black",size=0.7) +
  ylim(c(-10, 15))
a
a + coord_polar("x") + 
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
  ) 

