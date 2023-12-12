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

# boxplot ================================
mebo <- "C3H7NO2S"

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
  t_test(C3H7NO2S ~ time) %>% 
  #t_test(IL15 ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(H1N1_reclassify == "LL", 0.9, 
                       ifelse(H1N1_reclassify == "LH", 1.9,
                              ifelse(H1N1_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(H1N1_reclassify == "LL", 1.1, 
                       ifelse(H1N1_reclassify == "LH", 2.1,
                              ifelse(H1N1_reclassify == "HL", 3.1, 4.1))))

bxp <- ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo, color = "time",
                 paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)
bxp + stat_pvalue_manual(stat.test, label = "p.signif") + theme(legend.position = "top")

# compare 2 independent samples -------------------
compare_means(C3H7NO2S ~ time, data = metadat_boxplot, 
              group.by = c("season", "H1N1_reclassify"), method = "wilcox.test") # no sig. (same result as using t.test, kruskal.tests) 

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = mebo, color = "time",
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +
  stat_compare_means(label = "p.signif")

# compare 2 dependent samples ----------------------
metadat_boxplot2 <- metadat_boxplot %>% add_count(probandID) %>% filter(n == 2)
compare_means(C3H7NO2S ~ time, data = metadat_boxplot2, 
              group.by = c("season", "H1N1_reclassify"), method = "wilcox.test", paired = TRUE)

t.test(C3H7NO2S ~time, data = metadat_boxplot2, paired = TRUE)

metadat_boxplot2 %>% 
  group_by(season, H1N1_reclassify) %>% summarise(p_paired.test = t.test(C3H7NO2S ~ time)$p.value)
stat.test <- metadat_boxplot2 %>% group_by(season, H1N1_reclassify) %>% t_test(C3H7NO2S ~ time) # show the same p.value result as the above command


# boxplot with other strains -------------------


metadat_boxplot <- inputDat %>% 
  select(probandID, season, responder, time, c(mebo),
         matches("_abFC|_d0|_d28|_reclassify")) %>% filter(time %in% c("d0", "d28")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season))

# boxplots for each timepoint, per group, per season -------------
metadat_boxplot %>% 
  filter(strain == "H1N1" | strainSeason %in% c("B_2014", "B_2015", "H3N2_2015", "Byamagata_2020")) %>%
  ggboxplot(x = "reclassify", y = mebo, color = "time",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 1)

# add non sig. paired t_test to the graph -------------
metadat_boxplot2 <- metadat_boxplot %>% add_count(probandID, strainSeason) %>% filter(n == 2)

stat.test <- metadat_boxplot2 %>% 
  filter(strain == "H1N1" | strainSeason %in% c("B_2014", "B_2015", "H3N2_2015", "Byamagata_2020")) %>%
  group_by(strainSeason, reclassify) %>%  
  t_test(C3H7NO2S ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(reclassify == "LL", 0.9, 
                       ifelse(reclassify == "LH", 1.9,
                              ifelse(reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(reclassify == "LL", 1.1, 
                       ifelse(reclassify == "LH", 2.1,
                              ifelse(reclassify == "HL", 3.1, 4.1))))

bxp <- metadat_boxplot %>% 
  filter(strain == "H1N1" | strainSeason %in% c("B_2014", "B_2015", "H3N2_2015", "Byamagata_2020")) %>%
  ggboxplot(x = "reclassify", y = mebo, color = "time",
            paletter = "jco", add = "jitter") + 
  facet_wrap(~strainSeason, nrow = 2)
bxp + stat_pvalue_manual(stat.test, label = "p.signif") + theme(legend.position = "top")


