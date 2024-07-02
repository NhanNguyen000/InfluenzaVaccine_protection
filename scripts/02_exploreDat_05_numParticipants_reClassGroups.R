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
metadata_4groups <- metadata_healthy %>%
  select(name, cohort, season, group, sex, responder, ends_with("reclassify")) %>%
  pivot_longer(ends_with("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na("reclassify") %>%
  mutate(strain = gsub("_reclassify", "", strain),
         strain = factor(strain, levels = c("Bvictoria", "Byamagata", "B", "H3N2", "H1N1")))

plot_reClassGroups <- metadata_4groups %>%
  ggplot(aes(y = strain, fill = reclassify)) +
  geom_bar() + 
  #facet_grid(season ~., scales = "free", space = "free") + # the facet doesn have equal size with this code
  facet_grid(season ~., scales = "free") + 
  scale_fill_manual(values=c( c("#792C74", "#65771E", "#B65008", "#036879"))) +
  theme(strip.text.y = element_text(angle = 0)) +
  theme_classic() + 
  labs(x = "Number of vaccinees", y = "Strain", fill = "Per strain") + 
  theme(legend.position = "top", text = element_text(size = 24))

# violin plot to demonstrate the relcassification concept
png("output/numParticipants_per_reClassGroups.png", width = 720)
plot_reClassGroups 
dev.off()
