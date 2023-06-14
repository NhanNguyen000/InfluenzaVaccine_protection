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

# antibody T1 and abFC heterogeneity ------------------
H1N1_T1_abFC_diverPlot <- metadata_healthy %>% filter(season == "2015") %>%
  arrange(H1N1_T1_log2) %>% mutate(name = fct_reorder(name, H1N1_T1_log2)) %>%
  ggplot() +
  geom_bar(aes(x = name, y = 13, fill = H1N1_T1_log2), stat = "identity") +
  scale_fill_gradient2(low = "white", high = "red") +
  geom_line(aes(x = name, y = H1N1_abFC_log2, group = 1), stat = "identity", color="black",size=0.7) +
  ylim(c(-10, 15)) +
  geom_line(aes(x = name, y = log2(4), group = 1), color = "blue")

H1N1_T1_abFC_diverPlot 
H1N1_T1_abFC_diverPlot  + coord_polar("x") + 
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
  ) 
# antibody T1 and abFC negative correlation ------------------
H1N1_T1_abFC_corr <- metadata_healthy %>% filter(season == "2015") %>%
  mutate(group = ifelse(H1N1_T1 > 10, as.character(H1N1_T1), "10")) %>%
  mutate(group = factor(group, levels = c("10", "20", "40", "80", "160", "320", "640", "1280", "2560")))

H1N1_T1_abFC_corr %>% 
  ggplot() +
  geom_boxplot(aes(x = group, y = H1N1_abFC_log2)) + 
  geom_jitter(aes(x = group, y = H1N1_abFC_log2), width = 0.2, size = 1.2) +
  geom_smooth(aes(x = H1N1_T1_log2, y = H1N1_abFC_log2), method = "lm", se = FALSE) +
  stat_cor(aes(x = H1N1_T1_log2, y = H1N1_abFC_log2), label.x = 5, label.y = 10) +
  ylim(0, 12.5) + theme_classic()


H1N1_T1_abFC_corr %>% 
  ggplot(aes(x = group, y = H1N1_abFC_log2)) +
  geom_boxplot() + 
  geom_jitter(aes(color = H1N1_reclassify), width = 0.2, size = 1.2) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE) +
  stat_cor(aes(group = 1), label.x = 5, label.y = 10) +
  ylim(0, 12.5) + theme_classic()

# add color based on group reclassification
H1N1_T1_abFC_corr %>% 
  ggplot() +
  geom_boxplot(aes(x = group, y = H1N1_abFC_log2), outlier.color = NA) + 
  geom_jitter(aes(x = group, y = H1N1_abFC_log2, color = H1N1_reclassify), 
              width = 0.2, size = 2) +
  geom_smooth(aes(x = H1N1_T1_log2, y = H1N1_abFC_log2), method = "lm", se = FALSE) +
  stat_cor(aes(x = H1N1_T1_log2, y = H1N1_abFC_log2), label.x = 5, label.y = 10) +
  ylim(0, 12.5) + theme_classic() + 
  theme(legend.position = "top")
