rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data =======================================================================
load("cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

# boxplot ================================
my_comparisons <- list( c("R", "NR") )
my_comparisons_v1 <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))

protein <- "CD83"
metadat_boxplot <- inputDat %>% 
  select(season, responder, c(protein), matches("_abFC|_T1|_T4|_reclassify")) %>%
  mutate(H1N1_class = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NRs"))

ggboxplot(metadat_boxplot, x = "H1N1_class", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggboxplot(metadat_boxplot, x = "responder", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = my_comparisons_v1s, method = "t.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")


