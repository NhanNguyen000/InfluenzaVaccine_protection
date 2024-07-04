rm(list = ls())
library(tidyverse)
library(ggpubr)
library(webr)
library(patchwork)

load("processedDat/cohorts_dat.RData")

# metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate(group = paste0(cohort, "_", season)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# stack bar chart -------------
metadata_NRvsR <- metadata_healthy %>%
  mutate_at(vars(contains("abFC")), ~ifelse(.x >= 4, "responder", "non-responder")) %>%
  select(name, cohort, season, group, sex, responder, ends_with("abFC")) %>%
  pivot_longer(ends_with("abFC"), names_to = "strain", values_to = "strain_response") %>%
  drop_na("strain_response") %>%
  mutate(strain = gsub("_abFC", "", strain),
         strain = factor(strain, levels = c("Bvictoria", "Byamagata", "B", "H3N2", "H1N1")))

bars_NRvsR <- metadata_NRvsR %>%
  ggplot(aes(y = strain, fill = strain_response)) +
  geom_bar() + 
  #facet_grid(season ~., scales = "free", space = "free") + # the facet doesn have equal size with this code
  facet_grid(season ~., scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) +
  theme_bw() + 
  labs(x = "Number of vaccinees", y = "Strain", fill = "Per strain")

# pie chart --------------------------------
metadata_NRotherTR <- metadata_healthy %>% 
  select(name, cohort, season, group, sex, responder) %>%
  add_count(season) %>%  
  group_by(season, responder, n) %>% summarise(count = n())  %>%
  mutate(n2 = ifelse(n > 100, n, 100)) # for later use to increase the size of pie chart for 3 season


pies_NRotherTR <- metadata_NRotherTR  %>% 
  #ggplot(aes(x = n/2, y = count, fill = responder, width = n)) + # using raw data, the size of pie chart for 3 season is too small
  ggplot(aes(x = n2/2, y = count, fill = responder, width = n2)) + # increase the size of pie chart for 3 season
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  facet_grid(season ~.) +
  coord_polar("y") + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()) + 
  labs(fill = "Category")

bars_NRvsR + pies_NRotherTR +
  plot_layout(widths= c(4, -1, 1), guides = "collect") & theme(legend.position = "top")
