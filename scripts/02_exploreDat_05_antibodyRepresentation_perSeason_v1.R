rm(list = ls())
library(tidyverse)
library(ggpubr)
library(webr)

load("cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# for season 2015 -----------------------
abFC_2015 <- metadata_healthy %>% 
  select(name, cohort, season, group, sex, matches("abFC_log2")) %>%
  filter(season == "2015")

abFC_2015_longDat <- abFC_2015 %>% 
  pivot_longer(cols = ends_with("abFC_log2"), names_to = "strain", values_to = "abFC") %>%
  slice(-which(is.na(abFC)))

abFC_2015_longDat %>%
  ggboxplot(x = "sex", y = "abFC", color = "sex",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strain, nrow = 1) + 
  stat_compare_means(comparisons = list(c("f", "m")), method = "t.test")

abFC_2015_longDat %>% 
  ggplot(aes(x = strain, y = name, fill = abFC)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "white", high = "red") +
  geom_col(width = 0.7, just = -0.5) +
  coord_polar("y") +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        ) +
  labs(title = "Season 2015 - log2(abFC) for H3N2, H1N1, B")
