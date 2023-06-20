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

inputDat <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

# boxplot ================================
protein <- "CD83"
metadat_boxplot <- inputDat %>% 
  select(probandID, season, responder, time, c(protein), matches("_abFC|_T1|_T4|_reclassify")) %>%
  filter(time %in% c("T1", "T4"))

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = protein, color = "time",
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = protein, 
          color = "time", 
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)
# compare 2 independent samples -------------------
compare_means(CD83 ~ time, data = metadat_boxplot, 
              group.by = c("season", "H1N1_reclassify"), method = "wilcox.test") # no sig. (same result as using t.test, kruskal.tests) 

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = protein, color = "time",
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +
  stat_compare_means(label = "p.signif")

# compare 2 dependent samples ----------------------
metadat_boxplot2 <- metadat_boxplot %>% add_count(probandID) %>% filter(n == 2)
compare_means(CD83 ~ time, data = metadat_boxplot2, 
              group.by = c("season", "H1N1_reclassify"), method = "wilcox.test", paired = TRUE)

t.test(CD83 ~time, data = metadat_boxplot2, paired = TRUE)

metadat_boxplot2 %>% 
  group_by(season, H1N1_reclassify) %>% summarise(p_paired.test = t.test(CD83 ~ time)$p.value)
stat.test <- metadat_boxplot2 %>% group_by(season, H1N1_reclassify) %>% t_test(CD83 ~ time) # show the same p.value result as the above command

# add non sig. paired t_test to the graph
stat.test <- metadat_boxplot2 %>% 
  group_by(season, H1N1_reclassify) %>%  t_test(CD83 ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(H1N1_reclassify == "LL", 0.9, 
                       ifelse(H1N1_reclassify == "LH", 1.9,
                              ifelse(H1N1_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(H1N1_reclassify == "LL", 1.1, 
                       ifelse(H1N1_reclassify == "LH", 2.1,
                              ifelse(H1N1_reclassify == "HL", 3.1, 4.1))))

bxp <- ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = protein, color = "time",
                 paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1)
bxp + stat_pvalue_manual(stat.test, label = "p.signif")s
