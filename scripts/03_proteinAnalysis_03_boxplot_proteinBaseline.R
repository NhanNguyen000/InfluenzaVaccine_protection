rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data =======================================================================
load("processedDat/cohorts_dat.RData")

## metadata for all healthy subjects --------------------------------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

# boxplot ================================================================
compare_responder <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
compare_reClass <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
compare_abFC <- list(c("R", "NR"))
compare_abD0 <- list(c("low", "high"))

## select season and proteins --------------------------------------------------
strainSeasons <- c("H1N1_2014", "H1N1_2015", "H1N1_2019", "H1N1_2020",
                   "B_2014", "B_2015", "H3N2_2015", "Byamagata_2020") # season have at least >= 3 participant per reclassification group
protein <- "CD83"

## prepare the ploting data --------------------------------------------------
metadat_boxplot <- inputDat %>% 
  select(season, responder, c(protein), matches("_abFC|_d0|_d28|_reclassify")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season)) %>%
  mutate(abFC = ifelse(reclassify == "LH" | reclassify == "HH", "R", "NR")) %>%
  mutate(abBaseline = ifelse(reclassify == "LL" | reclassify == "LH", "low", "high"))


## make boxplot --------------------------------------------------
# based on reclassification 
boxplot_reClass <- metadat_boxplot %>% 
  filter(strainSeason %in% strainSeasons) %>%
  ggboxplot(x = "reclassify", y = protein,
            color = "#898366", add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  facet_wrap(~strainSeason, nrow = 2) +
  stat_compare_means(comparisons = compare_reClass, size = 5, method = "t.test")+
  theme(text = element_text(size = 18))

boxplot_reClass 

# based on abFC: NR vs R
boxplot_abFC <- metadat_boxplot %>% 
  filter(strainSeason %in% strainSeasons) %>%
  ggboxplot(x = "abFC", y = protein,
            color = "#898366", add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  facet_wrap(~strainSeason, nrow = 2) + 
  stat_compare_means(comparisons = compare_abFC, size = 5, method = "t.test")+
  theme(text = element_text(size = 18))

boxplot_abFC

# based on abD0: low Ab vs high Ab
boxplot_abD0 <- metadat_boxplot %>% 
  filter(strainSeason %in% strainSeasons) %>%
  ggboxplot(x = "abBaseline", y = protein,
            color = "#898366", add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  facet_wrap(~strainSeason, nrow = 2) + 
  stat_compare_means(comparisons = compare_abD0, size = 5, method = "t.test")+
  theme(text = element_text(size = 18))

boxplot_abD0

## save the plot --------------------------------------------------
png(paste0("output/boxplotProtein_reClass_", protein, ".png"), width = 720, height = 696)
boxplot_reClass 
dev.off()

png(paste0("output/boxplotProtein_abFC_", protein, ".png"), width = 576, height = 696)
boxplot_abFC
dev.off()

png(paste0("output/boxplotProtein_abD0_", protein, ".png"), width = 576, height = 696)
boxplot_abD0 
dev.off()
