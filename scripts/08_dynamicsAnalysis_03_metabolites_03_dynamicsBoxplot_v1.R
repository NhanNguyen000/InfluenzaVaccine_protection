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
mebo <- "C18H30O2"
mebo <- "C18H32O2"
mebo <- "C18H36O2"

mebo <- "C6H14N2O2" # lysine
mebo <- "C3H7NO2S" # cysteine

metadat_boxplot <- inputDat %>% 
  select(probandID, season, responder, time, c(mebo),
         matches("_abFC|_d0|_d28|_reclassify")) %>% filter(time %in% c("d0", "d28"))

# boxplots for each timepoint, per group, per season -------------
ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo, color = "time",
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)

# add non sig. paired t_test to the graph -------------
metadat_boxplot2 <- metadat_boxplot %>% add_count(probandID) %>% filter(n == 2)
stat.test <- metadat_boxplot2 %>% 
  group_by(season, H1N1_reclassify) %>%  
 # t_test(C3H6O4 ~ time) %>% 
 # t_test(C18H30O2 ~ time) %>% 
 # t_test(C18H32O2 ~ time) %>%
 # t_test(C18H36O2 ~ time) %>% 
 # t_test(C6H14N2O2 ~ time) %>%
  t_test(C3H7NO2S ~ time) %>%
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

compare_means(C18H30O2 ~ time, data = metadat_boxplot, 
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
  t_test(C3H6O4 ~ time, ref.group = "T1") %>% 
  #t_test(C8H8O5 ~ time, ref.group = "T1") %>% # no sig. p value
 # t_test(C3H7NO2 ~ time, ref.group = "T1") %>% 
 # t_test(C18H30O2 ~ time, ref.group = "T1") %>% 
  #t_test(C18H32O2 ~ time, ref.group = "T1") %>% 
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


# Calculate the ratio dynamics ---------------------
mebos <- c("C20H30O2", "C20H32O2", "C20H36O2", "C20H38O2", "C20H40O2",
           "C18H28O2", "C18H30O2", "C18H32O2", "C18H36O2",
           "C14H22O2", "C14H24O2", "C14H26O2", "C14H28O2")
mebos <- "C3H6O4"
metadat_boxplot <- inputDat %>% 
  select(probandID, season, responder, time, c(mebos), matches("_abFC|_T1|_T4|_reclassify")) %>%
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_T1 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high"))) %>% 
  filter(time %in% c("T1", "T4"))

# boxplots for each timepoint, per group, per season -------------
ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo, color = "time",
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)

# add non sig. paired t_test to the graph -------------
metadat_boxplot2 <- metadat_boxplot %>% add_count(probandID) %>% filter(n == 2)
stat.test <- metadat_boxplot2 %>% 
  group_by(season, H1N1_reclassify) %>%  
  # t_test(C3H6O4 ~ time) %>% 
  # t_test(C18H30O2 ~ time) %>% 
  #t_test(C18H36O2 ~ time) %>% 
  #t_test(C18H32O2 ~ time) %>% 
  #t_test(C20H40O2 ~ time) %>% 
 # t_test(C20H38O2 ~ time) %>% 
#  t_test(C20H36O2 ~ time) %>% 
 # t_test(C14H28O2 ~ time) %>%
 # t_test(C14H24O2 ~ time) %>%
 # t_test(C14H22O2 ~ time) %>%
  t_test(C3H6O4 ~ time) %>%
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(H1N1_reclassify == "LL", 0.9, 
                       ifelse(H1N1_reclassify == "LH", 1.9,
                              ifelse(H1N1_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(H1N1_reclassify == "LL", 1.1, 
                       ifelse(H1N1_reclassify == "LH", 2.1,
                              ifelse(H1N1_reclassify == "HL", 3.1, 4.1))))
mebo <- "C18H36O2"
mebo <- "C18H32O2"
mebo <- "C20H40O2"
mebo <- "C20H38O2"
mebo <- "C20H36O2"
mebo <- "C14H28O2"
mebo <- "C14H24O2"
mebo <- "C14H22O2"
mebo <- "C3H6O4"
bxp <- ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo, color = "time",
                 paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)
bxp + stat_pvalue_manual(stat.test, label = "p.signif")

a <- metadat_boxplot %>% mutate(ratio = C18H32O2 - C18H36O2)

metadat_boxplot3 <- metadat_boxplot %>% add_count(probandID) %>% filter(n == 2) %>%
#  mutate(ratio = C18H32O2 - C18H36O2) # sig.
#  mutate(ratio = C18H30O2 - C18H36O2) # no sig.
#  mutate(ratio = C18H28O2 - C18H36O2) # no sig.
#  mutate(ratio = C20H38O2 - C20H40O2) # sig. in LH 2015, HL 2014
#  mutate(ratio = C20H36O2 - C20H40O2) # sig. in LL and LH 2015, HL 2014
#  mutate(ratio = C20H32O2 - C20H40O2) # no sig.
#  mutate(ratio = C20H30O2 - C20H40O2) # no sig.
#  mutate(ratio = C14H22O2 - C14H28O2) # sig.
  mutate(ratio = C14H24O2 - C14H28O2) # sig. in LL and LH 2015
#  mutate(ratio = C14H26O2 - C14H28O2) # no sig.

stat.test <- metadat_boxplot3 %>% 
  group_by(season, H1N1_reclassify) %>%  
  t_test(ratio ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(H1N1_reclassify == "LL", 0.9, 
                       ifelse(H1N1_reclassify == "LH", 1.9,
                              ifelse(H1N1_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(H1N1_reclassify == "LL", 1.1, 
                       ifelse(H1N1_reclassify == "LH", 2.1,
                              ifelse(H1N1_reclassify == "HL", 3.1, 4.1))))

bxp <- ggboxplot(metadat_boxplot3, x = "H1N1_reclassify", y = "ratio", color = "time",
                 paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)
bxp + stat_pvalue_manual(stat.test, label = "p.signif")
