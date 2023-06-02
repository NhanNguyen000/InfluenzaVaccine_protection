rm(list = ls())
try(dev.off())

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
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy")

# assign the positions (start and end) for each participant per seasons
season_cohorts <- list()
for (year in unique(metadata_healthy$season)) {
  season_cohorts[[year]] <- metadata_healthy %>% 
    filter(season == year) %>% 
    mutate(start = as.numeric(rownames(.)) - 1, end = as.numeric(rownames(.)))
}
metadata_healthy_v2 <- season_cohorts %>% purrr::reduce(rbind)

# step 1: the soverview file for cohort participant and their sex
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

write.table(cohorts_sex, file = 'processedDat/cohorts_sex.txt', row.names = F, col.names = F, quote = F)

# step 2: H1N1 strain files
# H1N1 abFC
H1N1_abFC <- metadata_healthy_v2 %>% select(season, start, end, H1N1_abFC)
write.table(H1N1_abFC, file = 'processedDat/H1N1_abFC.txt', row.names = F, col.names = F, quote = F)

H1N1_abFC_log2 <- metadata_healthy_v2 %>% select(season, start, end, H1N1_abFC_log2)
write.table(H1N1_abFC_log2, file = 'processedDat/H1N1_abFC_log2.txt', row.names = F, col.names = F, quote = F)

# H1N1 T1
H1N1_T1 <- metadata_healthy_v2 %>% select(season, start, end, H1N1_T1)
write.table(H1N1_T1, file = 'processedDat/H1N1_T1.txt', row.names = F, col.names = F, quote = F)

# step 2: H3N2 strain files
# H3N2 abFC
# H3N2_abFC <- metadata_healthy_v2 %>% select(season, start, end, H3N2_abFC)
# write.table(H3N2_abFC, file = 'processedDat/H3N2_abFC.txt', row.names = F, col.names = F, quote = F)

H3N2_abFC_log2 <- metadata_healthy_v2 %>% select(season, start, end, H3N2_abFC_log2)
write.table(H3N2_abFC_log2, file = 'processedDat/H3N2_abFC_log2.txt', row.names = F, col.names = F, quote = F)

# H3N2 T1
H3N2_T1 <- metadata_healthy_v2 %>% select(season, start, end, H3N2_T1)
write.table(H3N2_T1, file = 'processedDat/H3N2_T1.txt', row.names = F, col.names = F, quote = F)

# step 3: B strains files
# B abFC for season 2014, 2015
B_abFC_log2 <- metadata_healthy_v2 %>% select(season, start, end, B_abFC_log2) %>% drop_na()
write.table(B_abFC_log2, file = 'processedDat/B_abFC_log2.txt', row.names = F, col.names = F, quote = F)

# B T1 for season 2014, 2015
B_T1 <- metadata_healthy_v2 %>% select(season, start, end, B_T1) %>% drop_na()
write.table(B_T1, file = 'processedDat/B_T1.txt', row.names = F, col.names = F, quote = F)

# Bvictoria abFC for season 2019, 2020
Bvictoria_abFC_log2 <- metadata_healthy_v2 %>% select(season, start, end, Bvictoria_abFC_log2) %>% drop_na()
write.table(Bvictoria_abFC_log2, file = 'processedDat/Bvictoria_abFC_log2.txt', row.names = F, col.names = F, quote = F)

# Bvictoria T1 for season 2019, 2020
Bvictoria_T1 <- metadata_healthy_v2 %>% select(season, start, end, Bvictoria_T1) %>% drop_na()
write.table(Bvictoria_T1, file = 'processedDat/Bvictoria_T1.txt', row.names = F, col.names = F, quote = F)

# Byamagata abFC for season 2019, 2020
Byamagata_abFC_log2 <- metadata_healthy_v2 %>% select(season, start, end, Byamagata_abFC_log2) %>% drop_na()
write.table(Byamagata_abFC_log2, file = 'processedDat/Byamagata_abFC_log2.txt', row.names = F, col.names = F, quote = F)

# B T1 for season 2019, 2020
Byamagata_T1 <- metadata_healthy_v2 %>% select(season, start, end, Byamagata_T1) %>% drop_na()
write.table(Byamagata_T1, file = 'processedDat/Byamagata_T1.txt', row.names = F, col.names = F, quote = F)


# martijn code -----------------
setwd('/Users/martijnzoodsma/Documents/PhD/influenza/iMED/circos/Ab titers 3/')
meta1 <- read.csv('data/meta_cohort1.csv') %>% filter(time == 'T1') %>% mutate(cohort = 'cohort1')
meta2 <- read.csv('data/meta_cohort2.csv') %>% filter(time == 'T1') %>% mutate(cohort = 'cohort2')
​
# HAI titers
hai1 <- read.table('data/hai_titers_pilot.tsv', sep = '\t', header=T)
hai1$pID <- lapply(strsplit(hai1$pID, split='-', fixed=T), function(x) {x[[2]]}) %>% unlist() %>% stringr::str_remove(., "^0+")
​
hai2 <- read.table('data/hai_titers.csv', sep = ';', header=T)
hai2$sample <- lapply(strsplit(hai2$sample, split='-', fixed=T), function(x) {x[[2]]}) %>% unlist() %>% stringr::str_remove(., "^0+")
​
meta1 <- merge(meta1, hai1, by.x = 'ProbandID', by.y = 'pID', sort=F, all.x = T, all.y = T) 
meta2 <- merge(meta2, hai2 %>% select(-responder), by.x = 'ProbandID', by.y = 'sample', sort=F, all.x = T, all.y = T) 
​
​
# Responder categories
resp1 <- read.csv('data/cohort1_responder_classes.tsv', sep = ' ', header=T) %>% 
  distinct(ProbandID, .keep_all = T)
​
# resp2 <- read.csv('data/resp_classes_cohort2.csv', sep = '\t')
# resp2$Pseudonym <- lapply(strsplit(resp2$Pseudonym, split='-', fixed=T), function(x) {x[[2]]}) %>% unlist() %>% stringr::str_remove(., "^0+")
# resp2$H1N1_response <- ifelse(resp2$H1N1_response == 'low', F, T)
# resp2$H3N2_response <- ifelse(resp2$H3N2_response == 'low', F, T)
# resp2$B_response <- ifelse(resp2$B_response == 'low', F, T)
# 
# meta1 <- merge(meta1, resp1, by = 'ProbandID', all.x=T, all.y=T)
# meta2 <- merge(meta2, resp2, by.x = 'ProbandID', by.y = 'Pseudonym', all.x=T, all.y=T)
​
df <- rbind(
  meta1 %>% select(cohort, ProbandID, age, gender, responder, H1N1_T1, H1N1_T4, H3N2_T1, H3N2_T4, B_T1, B_T4) %>% mutate(gender = ifelse(gender == 'female', 1, 0)), 
  meta2 %>% select(cohort, ProbandID, age, gender, responder, H1N1_T1, H1N1_T4, H3N2_T1, H3N2_T4, B_T1, B_T4) %>% mutate(gender = ifelse(gender == 'female', 1, 0))
)
​
df$cohort <- ifelse(df$cohort == 'cohort1', 'Replication', 'Discovery')
​
df$H1N1_T1 <- ifelse(df$H1N1_T1 == 0, 1, df$H1N1_T1)
df$H1N1_T4 <- ifelse(df$H1N1_T4 == 0, 1, df$H1N1_T4)
df$H3N2_T1 <- ifelse(df$H3N2_T1 == 0, 1, df$H3N2_T1)
df$H3N2_T4 <- ifelse(df$H3N2_T4 == 0, 1, df$H3N2_T4)
df$B_T1 <- ifelse(df$B_T1 == 0, 1, df$B_T1)
df$B_T4 <- ifelse(df$B_T4 == 0, 1, df$B_T4)
​
df$abH1N1 <- df$H1N1_T4 / df$H1N1_T1
df$abH3N2 <- df$H3N2_T4 / df$H3N2_T1
df$abB <- df$B_T4 / df$B_T1
​
df$H1N1_T1 <- log2(df$H1N1_T1)
df$H1N1_T4 <- log2(df$H1N1_T4)
df$H3N2_T1 <- log2(df$H3N2_T1)
df$H3N2_T4 <- log2(df$H3N2_T4)
df$B_T1 <- log2(df$B_T1)
df$B_T4 <- log2(df$B_T4)
​
df$abH1N1 <- log2(df$abH1N1)
df$abH3N2 <- log2(df$abH3N2)
df$abB <- log2(df$abB)
​
head(df)
df$B_response <- ifelse(df$abB >= 4, T, F)
df$H3N2_response <- ifelse(df$abH3N2 >= 4, T, F)
df$H1N1_response <- ifelse(df$abH1N1 >= 4, T, F)
​
# Define the classes
df$resp_class <- NA
df$resp_class <- ifelse(df$H1N1_response & df$H3N2_response & df$B_response, 'TR', df$resp_class)
df$resp_class <- ifelse(!df$H1N1_response & !df$H3N2_response & !df$B_response, 'NR', df$resp_class)
df$resp_class <- ifelse(df$H1N1_response & !df$H3N2_response & !df$B_response, 'SR (H1N1)', df$resp_class)
df$resp_class <- ifelse(!df$H1N1_response & df$H3N2_response & !df$B_response, 'SR (H3N2)', df$resp_class)
df$resp_class <- ifelse(!df$H1N1_response & !df$H3N2_response & df$B_response, 'SR (B)', df$resp_class)
​
df$resp_class <- ifelse(df$H1N1_response & !df$H3N2_response & df$B_response, 'DR (B & H1N1)', df$resp_class)
df$resp_class <- ifelse(!df$H1N1_response & df$H3N2_response & df$B_response, 'DR (B & H3N2)', df$resp_class)
df$resp_class <- ifelse(df$H1N1_response & df$H3N2_response & !df$B_response, 'DR (H1N1 & H3N2)', df$resp_class)
​
# Sorting
# df$cohort <- factor(df$cohort, levels = c('cohort1', 'cohort2'))
df$cohort <- factor(df$cohort, levels = c('Replication', 'Discovery'))
df$resp_class <- factor(df$resp_class, levels = c('NR',"SR (H1N1)", "SR (H3N2)", "SR (B)", "DR (H1N1 & H3N2)", "DR (B & H3N2)",  "DR (B & H1N1)", 'TR'))
df$responder_broad <- df$responder
df$responder <- df$resp_class
​
# df %<>% group_by(cohort) %>% arrange(abH1N1) %>% ungroup() %>% 
df$med <- apply(df[, c('abB', 'abH3N2', 'abH1N1')], 1, median)
# df %<>% group_by(cohort) %>% arrange(abH1N1, abH3N2, abB, .by_group = TRUE)  %>% ungroup() %>% mutate(start = as.numeric(rownames(.)) - 1, end = as.numeric(rownames(.)))
df %<>% group_by(cohort, responder_broad) %>% arrange(med, .by_group = T) %>% ungroup() %>% mutate(start = as.numeric(rownames(.)) - 1, end = as.numeric(rownames(.)))
​
# View(df)
# df %<>% arrange(cohort, responder, gender, age) %>%
#   mutate(start = as.numeric(rownames(.)) - 1, end = as.numeric(rownames(.)))
​
# Export all
# Gender, age
write.table(df %>% select(cohort, start, end, gender) %>% mutate(gender = ifelse(gender == 1, 'm', 'f')), file = 'data/genders2.txt', row.names = F, col.names = F, quote = F)
# write.table(df %>% select(cohort, start, end, responder) %>% mutate(responder = as.numeric(responder)), file = 'data/responder.txt', row.names = F, col.names = F, quote = F)
​
write.table(df %>% select(cohort, start, end, gender), file = 'data/genders.txt', row.names = F, col.names = F, quote = F)
write.table(df %>% select(cohort, start, end, age), file = 'data/ages.txt', row.names = F, col.names = F, quote = F)
​
head(df)
​
​
cats <- df %>% 
  select(cohort, responder_broad, start, end) %>% 
  mutate(var = ifelse(responder_broad == 'NR', 1, ifelse(responder_broad == 'other', 2, 3))) %>% 
  select(-responder_broad)
write.table(cats, file = 'data/responder.txt', row.names = F, col.names = F, quote = F)
​
# Strains
write.table(df %>% select(cohort, start, end, H1N1_T1), file = 'data/H1N1_T1.txt', row.names = F, col.names = F, quote = F)
write.table(df %>% select(cohort, start, end, H1N1_T4), file = 'data/H1N1_T4.txt', row.names = F, col.names = F, quote = F)
write.table(df %>% select(cohort, start, end, H3N2_T1), file = 'data/H3N2_T1.txt', row.names = F, col.names = F, quote = F)
write.table(df %>% select(cohort, start, end, H3N2_T4), file = 'data/H3N2_T4.txt', row.names = F, col.names = F, quote = F)
write.table(df %>% select(cohort, start, end, B_T1), file = 'data/B_T1.txt', row.names = F, col.names = F, quote = F)
write.table(df %>% select(cohort, start, end, B_T4), file = 'data/B_T4.txt', row.names = F, col.names = F, quote = F)
​
# Fold changes
write.table(df %>% select(cohort, start, end, abB), file = 'data/B_FC.txt', row.names = F, col.names = F, quote = F)
write.table(df %>% select(cohort, start, end, abH1N1), file = 'data/H1N1_FC.txt', row.names = F, col.names = F, quote = F)
write.table(df %>% select(cohort, start, end, abH3N2), file = 'data/H3N2_FC.txt', row.names = F, col.names = F, quote = F)
​
​
# Line plot
df$lim = log2(4)
write.table(df %>% select(cohort, start, end, lim), file = 'data/lim_FC.txt', row.names = F, col.names = F, quote = F)
​
# What are the limits we need to set for the Ab tracks?
# logFC
summary(df$abB)
summary(df$abH1N1)
summary(df$abH3N2)
​
# Pre-vacc titers
summary(df$B_T1)
summary(df$H1N1_T1)
summary(df$H3N2_T1)
​
# Responder categories
write.table(df %>% select(cohort, start, end, H1N1_response) %>% mutate(ifelse(H1N1_response, 1, 0)) %>% select(-H1N1_response), file = 'data/H1N1_response.txt', row.names = F, col.names = F, quote = F)
write.table(df %>% select(cohort, start, end, H3N2_response) %>% mutate(ifelse(H3N2_response, 1, 0))%>% select(-H3N2_response), file = 'data/H3N2_response.txt', row.names = F, col.names = F, quote = F)
write.table(df %>% select(cohort, start, end, B_response) %>% mutate(ifelse(B_response, 1, 0))%>% select(-B_response), file = 'data/B_response.txt', row.names = F, col.names = F, quote = F)
​
​
​
​
​
###### Color legends
library(cowplot)
​
# Genders
genders <- ifelse(df$gender == 0, 'Male', 'Female')
p1 <- ggplot() + geom_point(aes(x=1, y=genders, color = genders, fill = genders), shape=22) + 
  scale_color_manual(values = c('Male' = '#c6dbef', 'Female' = '#fee5d9')) + 
  scale_fill_manual(values = c('Male' = '#c6dbef', 'Female' = '#fee5d9')) + 
  theme_bw() + 
  labs(color = 'Sex', fill = 'Sex') + 
  theme(legend.title = element_text(size=8), legend.text = element_text(size=6)) +
  guides(color = guide_legend(override.aes = list(size = 8)))
p1
leg <- get_legend(p1)
​
​
pdf('legend_sex.pdf', width = 1, height = 1)
plot(leg)
dev.off()
​
# Responders
p1 <- ggplot(df) + geom_point(aes(x=1, y=responder_broad, color = responder_broad, fill = responder_broad), shape=22) + 
  scale_color_manual(values = c('NR' = 'dodgerblue', 'other' = 'goldenrod1', 'TR' = 'indianred2')) + 
  scale_fill_manual(values = c('NR' = 'dodgerblue', 'other' = 'goldenrod1', 'TR' = 'indianred2')) + 
  theme_bw() + 
  labs(color = 'Category', fill = 'Category')+ 
  guides(color = guide_legend(override.aes = list(size = 8))) +
  theme(legend.title = element_text(size=8), legend.text = element_text(size=6))
p1
leg <- get_legend(p1)
​
pdf('legend_resp.pdf', width = 1, height = 2)
plot(leg)
dev.off()
​
# Ab titers
library(RColorBrewer)
df %>% select(abH1N1, abH3N2, abB) %>% melt() -> x
​
​
# ggplot() +
#   geom_point(data = filter(x, variable == 'abB'), aes(variable, value, color=value)) +
#   scale_color_distiller(palette = 'BuGn', direction = 1, breaks = c(0, 12), labels = c('Low', 'High')) + 
#   
#   ggnewscale::new_scale_colour() +
#   geom_point(data = filter(x, variable == 'abH3N2'), aes(variable, value, color=value)) +
#   scale_color_distiller(palette = 'OrRd', direction = 1, breaks = c(-3, 14), labels = c('Low', 'High')) +
#   
#   ggnewscale::new_scale_colour() +
#   geom_point(data = filter(x, variable == 'abH1N1'), aes(variable, value, color=value)) +
#   scale_color_distiller(palette = 'BuPu', direction = 1, breaks = c(-7, 12), labels = c('Low', 'High')) +
#   theme(legend.position = 'top') +
#   guides(color=guide_legend(title = 'test'))
​
ggplot(filter(x, variable == 'abB')) + geom_point(aes(variable, value, color=value)) + 
  scale_color_distiller(palette = 'BuGn', direction = 1, breaks = c(0, 12), #labels = c('Low', 'High'), 
                        guide = guide_colourbar(direction = "horizontal", label = F,
                                                ticks = F, title.vjust = .70, title = 'B',
                                                title.position = "left", frame.colour = "black")) -> p1
cowplot::get_legend(p1) -> l1
​
ggplot(filter(x, variable == 'abH1N1')) + geom_point(aes(variable, value, color=value)) + 
  scale_color_distiller(palette = 'BuPu', direction = 1, breaks = c(-7, 12), labels = c('Low', 'High'), 
                        guide = guide_colourbar(direction = "horizontal", label = F,
                                                ticks = F, title.vjust = .70, title = 'H1N1',
                                                title.position = "left", frame.colour = "black")) -> p2
cowplot::get_legend(p2) -> l2
​
ggplot(filter(x, variable == 'abH3N2')) + geom_point(aes(variable, value, color=value)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.margin = margin(t=-.1,b=-.1,unit="cm")) +
  scale_color_distiller(palette = 'OrRd', direction = 1, breaks = c(-3, 14), labels = c('Low', 'High'), 
                        guide = guide_colourbar(direction = "horizontal", label = T,
                                                ticks = F, title.vjust = .70, title = 'H3N2',
                                                title.position = "left", frame.colour = "black")) -> p3
cowplot::get_legend(p3) -> l3
​
library(patchwork)
plot(l3)
plot(l1) / plot(l2) / plot(l3)
​
pdf('colorlegends.pdf', width = 2, height = 3)
ggpubr::ggarrange(plotlist = list(l1, l2, l3), nrow = 3, ncol=1)
dev.off()
​