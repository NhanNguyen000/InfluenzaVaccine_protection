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
compare_responder <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
compare_reClass <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
compare_abFC <- list(c("R", "NR"))
compare_abT1 <- list(c("low", "high"))

protein <- "PGF"
protein <- "CXCL9"
protein <- "CD83"

metadat_boxplot <- inputDat %>% 
  select(season, responder, c(protein), matches("_abFC|_T1|_T4|_reclassify")) %>%
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_T1 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high")))

# based on reclassification 
ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

# based on abFC: NR vs R
ggboxplot(metadat_boxplot, x = "H1N1_abFC", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_abFC, method = "t.test")

# based on abT1: low Ab vs high Ab
ggboxplot(metadat_boxplot, x = "H1N1_abBaseline", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_abT1, method = "t.test")

# based on responders group: NR, other, TR
ggboxplot(metadat_boxplot, x = "responder", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_responder, method = "t.test")
