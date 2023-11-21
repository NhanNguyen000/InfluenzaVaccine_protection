rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
library(reshape2)

load("cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  arrange(season,  H1N1_abFC_log2, H1N1_d0_log2)

# plot --------------------------------------------------
H1N1_d0_abFC_diverPlot <- metadata_healthy %>% filter(season == "2015") %>%
  arrange(H1N1_d0_log2) %>% mutate(name = fct_reorder(name, H1N1_d0_log2)) %>%
  ggplot() +
  geom_bar(aes(x = name, y = 13, fill = H1N1_d0_log2), stat = "identity") +
  scale_fill_gradient2(low = "white", high = "red") +
  geom_line(aes(x = name, y = H1N1_abFC_log2, group = 1), stat = "identity", color="black",size=0.7) +
  ylim(c(-10, 15)) +
  geom_line(aes(x = name, y = log2(4), group = 1), color = "blue")

H1N1_d0_abFC_diverPlot 

a <- metadata_healthy %>%
  select(season, name, matches("d0_log2")) %>%
  pivot_longer(cols = ends_with("d0_log2"), names_to = "strain", values_to = "d0_log2") %>%
  mutate(strain = gsub("_d0_log2", "", strain)) %>% drop_na(d0_log2) %>%
  full_join(
    metadata_healthy %>%
      select(season, name, matches("abFC_log2")) %>%
      pivot_longer(cols = ends_with("abFC_log2"), names_to = "strain", values_to = "abFC_log2") %>%
      mutate(strain = gsub("_abFC_log2", "", strain)) %>% drop_na(abFC_log2)
  )

a %>% filter(season == "2015") %>%
  ggplot(aes(x = name, y = strain, fill = d0_log2)) + 
  geom_tile() +
#  geom_text(aes(label = ifelse(is.na(p.value), NA, "*"))) +
  scale_fill_gradient(low = "white", high = "red") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


metadata_healthy %>% filter(season == "2015") %>%
  ggplot(aes(x = name, y = 13, fill = H1N1_d0_log2)) +
  geom_bar(stat = "identity") + 
  scale_fill_gradient(low = "white", high = "red") + 
  theme(axis.text.x = element_text(angle = 90))

metadata_healthy %>%  arrange(H1N1_d0_log2) %>% mutate(name = fct_reorder(name, H1N1_d0_log2)) %>%
  ggplot() +
  geom_bar(aes(x = name, y = 13, fill = H1N1_d0_log2), stat = "identity") +
  scale_fill_gradient2(low = "white", high = "red") +
  geom_line(aes(x = name, y = H1N1_abFC_log2, group = 1), stat = "identity", color="black",size=0.7) +
  ylim(c(-10, 15)) +
  geom_line(aes(x = name, y = log2(4), group = 1), color = "blue")
