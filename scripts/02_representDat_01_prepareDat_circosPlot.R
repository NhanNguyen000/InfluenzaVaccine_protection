rm(list = ls())
try(dev.off())

library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
library(reshape2)

load("processedDat/cohorts_dat.RData")

# metadata for all healthy subjects --------------------------------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  arrange(season, H1N1_reclassify, H1N1_abFC)

# prepare data for circo plot --------------------------------------------------
## assign the positions (start and end) for each participant per seasons -----------------------------
season_cohorts <- list()
for (year in unique(metadata_healthy$season)) {
  season_cohorts[[year]] <- metadata_healthy %>% 
    filter(season == year) %>% 
    mutate(start = as.numeric(rownames(.)) - 1, end = as.numeric(rownames(.)))
}
metadata_healthy_v2 <- season_cohorts %>% purrr::reduce(rbind)

## step 1: the overview file for cohort participant and their sex -----------------------------
metadata_healthy_v2  %>% count(cohort, season)

part1 <- metadata_healthy_v2 %>% 
  count(cohort, season) %>%
  mutate(col_1 = "chr", col_2 = "-", col_3 = season, 
         valName = season, 
         start = 0, end = n, color = "white") %>%
  select(-cohort, -season, -n)

part2 <- metadata_healthy_v2 %>% 
  select(probandID, season, sex, start, end) %>%
  mutate(col_1 = "band", col_2 = season, 
         col_3 = paste0(season, "_", probandID),
         valName = probandID, 
         color = ifelse(sex == "f", "pink", "blue")) %>% 
  select(-probandID, -season, -sex) %>%
  relocate(start, end, .before = color)

cohorts_sex <- rbind(part1, part2)

write.table(cohorts_sex, file = 'processedDat/circos_input/cohorts_sex.txt', row.names = F, col.names = F, quote = F)

# step 2: H1N1 strain files ----------------------------------------------------------
# H1N1 abFC
H1N1_abFC <- metadata_healthy_v2 %>% select(season, start, end, H1N1_abFC)
write.table(H1N1_abFC, file = 'processedDat/circos_input/H1N1_abFC.txt', row.names = F, col.names = F, quote = F)

H1N1_abFC_log2 <- metadata_healthy_v2 %>% select(season, start, end, H1N1_abFC_log2)
write.table(H1N1_abFC_log2, file = 'processedDat/circos_input/H1N1_abFC_log2.txt', row.names = F, col.names = F, quote = F)

# H1N1 d0
H1N1_d0 <- metadata_healthy_v2 %>% select(season, start, end, H1N1_d0)
write.table(H1N1_d0, file = 'processedDat/circos_input/H1N1_d0.txt', row.names = F, col.names = F, quote = F)

H1N1_d0_log2 <- metadata_healthy_v2 %>% select(season, start, end, H1N1_d0_log2)
write.table(H1N1_d0_log2, file = 'processedDat/circos_input/H1N1_d0_log2.txt', row.names = F, col.names = F, quote = F)

# step 2: H3N2 strain files ----------------------------------------------------------
# H3N2 abFC
# H3N2_abFC <- metadata_healthy_v2 %>% select(season, start, end, H3N2_abFC)
# write.table(H3N2_abFC, file = 'processedDat/circos_input/H3N2_abFC.txt', row.names = F, col.names = F, quote = F)

H3N2_abFC_log2 <- metadata_healthy_v2 %>% select(season, start, end, H3N2_abFC_log2)
write.table(H3N2_abFC_log2, file = 'processedDat/circos_input/H3N2_abFC_log2.txt', row.names = F, col.names = F, quote = F)

# H3N2 d0
H3N2_d0 <- metadata_healthy_v2 %>% select(season, start, end, H3N2_d0)
write.table(H3N2_d0, file = 'processedDat/circos_input/H3N2_d0.txt', row.names = F, col.names = F, quote = F)

H3N2_d0_log2 <- metadata_healthy_v2 %>% select(season, start, end, H3N2_d0_log2)
write.table(H3N2_d0_log2, file = 'processedDat/circos_input/H3N2_d0_log2.txt', row.names = F, col.names = F, quote = F)

# step 3: B strains files ----------------------------------------------------------
# B abFC for season 2014, 2015
B_abFC_log2 <- metadata_healthy_v2 %>% select(season, start, end, B_abFC_log2) %>% drop_na()
write.table(B_abFC_log2, file = 'processedDat/circos_input/B_abFC_log2.txt', row.names = F, col.names = F, quote = F)

# B d0 for season 2014, 2015
B_d0 <- metadata_healthy_v2 %>% select(season, start, end, B_d0) %>% drop_na()
write.table(B_d0, file = 'processedDat/circos_input/B_d0.txt', row.names = F, col.names = F, quote = F)

B_d0_log2 <- metadata_healthy_v2 %>% select(season, start, end, B_d0_log2) %>% drop_na()
write.table(B_d0_log2, file = 'processedDat/circos_input/B_d0_log2.txt', row.names = F, col.names = F, quote = F)

# Bvictoria abFC for season 2019, 2020
Bvictoria_abFC_log2 <- metadata_healthy_v2 %>% select(season, start, end, Bvictoria_abFC_log2) %>% drop_na()
write.table(Bvictoria_abFC_log2, file = 'processedDat/circos_input/Bvictoria_abFC_log2.txt', row.names = F, col.names = F, quote = F)

# Bvictoria d0 for season 2019, 2020
Bvictoria_d0 <- metadata_healthy_v2 %>% select(season, start, end, Bvictoria_d0) %>% drop_na()
write.table(Bvictoria_d0, file = 'processedDat/circos_input/Bvictoria_d0.txt', row.names = F, col.names = F, quote = F)

Bvictoria_d0_log2 <- metadata_healthy_v2 %>% select(season, start, end, Bvictoria_d0_log2) %>% drop_na()
write.table(Bvictoria_d0_log2, file = 'processedDat/circos_input/Bvictoria_d0_log2.txt', row.names = F, col.names = F, quote = F)

# Byamagata abFC for season 2019, 2020
Byamagata_abFC_log2 <- metadata_healthy_v2 %>% select(season, start, end, Byamagata_abFC_log2) %>% drop_na()
write.table(Byamagata_abFC_log2, file = 'processedDat/circos_input/Byamagata_abFC_log2.txt', row.names = F, col.names = F, quote = F)

# B d0 for season 2019, 2020
Byamagata_d0 <- metadata_healthy_v2 %>% select(season, start, end, Byamagata_d0) %>% drop_na()
write.table(Byamagata_d0, file = 'processedDat/circos_input/Byamagata_d0.txt', row.names = F, col.names = F, quote = F)

Byamagata_d0_log2 <- metadata_healthy_v2 %>% select(season, start, end, Byamagata_d0_log2) %>% drop_na()
write.table(Byamagata_d0_log2, file = 'processedDat/circos_input/Byamagata_d0_log2.txt', row.names = F, col.names = F, quote = F)
