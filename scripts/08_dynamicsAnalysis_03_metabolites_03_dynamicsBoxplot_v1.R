rm(list = ls())

library(tidyverse)
library(ggpubr)
library(rstatix)
# load data =======================================================================
load("cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- mebo_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

# boxplot for T1, T4 ================================
mebo <- "C3H6O4"
metadat_boxplot <- inputDat %>% 
  select(probandID, season, responder, time, c(mebo),
         matches("_abFC|_T1|_T4|_reclassify")) %>% filter(time %in% c("T1", "T4"))

# boxplots for each timepoint, per group, per season -------------
ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo, color = "time",
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)

# add non sig. paired t_test to the graph -------------
metadat_boxplot2 <- metadat_boxplot %>% add_count(probandID) %>% filter(n == 2)
stat.test <- metadat_boxplot2 %>% 
  group_by(season, H1N1_reclassify) %>%  
  t_test(C3H6O4 ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(H1N1_reclassify == "LL", 0.9, 
                       ifelse(H1N1_reclassify == "LH", 1.9,
                              ifelse(H1N1_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(H1N1_reclassify == "LL", 1.1, 
                       ifelse(H1N1_reclassify == "LH", 2.1,
                              ifelse(H1N1_reclassify == "HL", 3.1, 4.1))))

bxp <- ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo, color = "time",
                 paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)
bxp + stat_pvalue_manual(stat.test, label = "p.signif")

# compare 2 independent samples -------------------
compare_means(C3H6O4 ~ time, data = metadat_boxplot, 
              group.by = c("season", "H1N1_reclassify"), method = "wilcox.test") # no sig. (same result as using t.test, kruskal.tests) 

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo, color = "time",
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +
  stat_compare_means(label = "p.signif")

# boxplot for all time point ================================
mebo <- "C3H6O4"
mebo <- "C8H8O5"
mebo <- "C3H7NO2"
mebo <- "C18H30O2"
mebo <- "C18H32O2"

metadat_boxplot <- inputDat %>% 
  select(probandID, season, responder, time, c(mebo),
         matches("_abFC|_T1|_T4|_reclassify"))

# boxplots for each timepoint, per group, per season -------------
ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo, color = "time",
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)

# add non sig. paired t_test to the graph -------------
metadat_boxplot2 <- metadat_boxplot %>% add_count(probandID) %>% filter(n == 3)
stat.test <- metadat_boxplot2 %>% 
  group_by(season, H1N1_reclassify) %>%  
  #t_test(C3H6O4 ~ time, ref.group = "T1") %>% 
  #t_test(C8H8O5 ~ time, ref.group = "T1") %>% # no sig. p value
 # t_test(C3H7NO2 ~ time, ref.group = "T1") %>% 
 # t_test(C18H30O2 ~ time, ref.group = "T1") %>% 
  t_test(C18H32O2 ~ time, ref.group = "T1") %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(H1N1_reclassify == "LL", 0.9, 
                       ifelse(H1N1_reclassify == "LH", 1.9,
                              ifelse(H1N1_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(H1N1_reclassify == "LL", 1.1, 
                       ifelse(H1N1_reclassify == "LH", 2.1,
                              ifelse(H1N1_reclassify == "HL", 3.1, 4.1)))) 

bxp <- ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo, color = "time",
                 paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)
bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif")

#  modify plot
stat.test_manual  <- stat.test %>% filter(p < 0.05) %>%
  mutate(xmin = c(0.8, 0.8, 0.8, 0.8), 
         xmax = c(1, 1, 1.2, 1.1),
         y.position = c(14.2, 14.5, 14.7, 12.5))

stat.test_manual  <- stat.test %>% filter(p < 0.05) %>%
  mutate(xmin = c(2.7, 2.7, 0.7, 0.7, 1.7, 1.7, 2.7, 2.7, 3.7, 3.7, 0.7, 0.7), 
         xmax = c(3, 3.2, 1, 1.2, 2, 2.2, 3, 3.2, 4, 4.2, 1, 1.2))

stat.test_manual  <- stat.test %>% filter(p < 0.05) %>%
  mutate(xmin = c(2.7, 2.7, 0.7, 0.7, 1.7, 1.7, 2.7, 2.7, 3.7, 3.7, 1.7), 
         xmax = c(3, 3.2, 1, 1.2, 2, 2.2, 3, 3.2, 4, 4.2, 2.2))

stat.test_manual  <- stat.test %>% filter(p < 0.05) %>%
  mutate(xmin = c(1.7, 2.7, 2.7, 0.7, 0.7, 1.7, 1.7, 2.7, 2.7, 3.7, 3.7), 
         xmax = c(2, 3, 3.2, 1, 1.2, 2, 2.2, 3, 3.2, 4, 4.2))

bxp + stat_pvalue_manual(stat.test_manual, label = "p.adj.signif")


