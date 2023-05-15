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
PieDonut(metadata_healthy, aes(season, sex), r0 = 0.5, r2 = 1.1, start = -120,
         title = "Distribution of gender by season")

# check more: https://stackoverflow.com/questions/17069436/hierarchical-multilevel-pie-chart
https://bookdown.org/content/b298e479-b1ab-49fa-b83d-a57c2b034d49/part.html#circular-packing
https://stackoverflow.com/questions/50004058/multiple-dependent-level-sunburst-doughnut-chart-using-ggplot2


metadata_healthy %>% 
  ggplot() + 
  geom_bar(aes(x = season, y = cohort), stat = "identity") + coord_polar("y")

metadata_healthy %>% 
  ggplot() +  
  geom_col(aes(x = season, y = sex)) + coord_polar("y")

metadata_healthy %>% 
  ggplot(aes(x = name, y = cohort, fill = cohort)) +  
  geom_bar(stat = "identity") + coord_polar(theta = "y")

ggplot() + 
  geom_col(data = metadata_healthy, aes(x = 2, y = season), color = "white")

a<- metadata_healthy %>% rownames_to_column("rowId") 
a %>%
  ggplot() +
  geom_col(aes(x = 2, y = rowId, fill = H1N1_T1_log2)) 

    ggplot() +
    geom_col(data = metadata_healthy, aes(x = 1, y = name, fill = H1N1_T1_log2)) +
    geom_text(data = metadata_healthy, )
    geom_col(data = metadata_healthy, aes(x = 2, y = name, fill = H3N2_T1_log2))
  