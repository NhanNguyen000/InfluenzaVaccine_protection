rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data ====================================================================
load("processedDat/cohorts_dat.RData")

## metadata for all healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- mebo_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

# boxplot =======================================================================
compare_responder <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
compare_reClass <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
compare_abFC <- list(c("R", "NR"))
compare_abD0 <- list(c("low", "high"))

## select season and metabolites --------------------------------------------------
strainSeasons <- c("H1N1_2014", "H1N1_2015", "H1N1_2019", "H1N1_2020",
                   "B_2014", "B_2015", "H3N2_2015", "Byamagata_2020") # season have at least >= 3 participant per reclassification group

mebo <- "C18H32O2"
mebo <- "C3H7NO2S"

# selected_vars
mebo <- "C9H14O4"
mebo <- "C19H40NO7P"
mebo <- "C22H36O2"
mebo <- "C20H32O2"
mebo <- "C12H22O2"
mebo <- "C23H32O2"

# consistent trend across 8 strain in lm() models
mebo <- "C3H6O2"     
mebo <- "C7H10"      # maybe
mebo <-"C5H10O3" 
mebo <- "C7H7NO"
mebo <-"C4H10O4" # maybe
mebo <- "C2H6O4S" 
mebo <- "C7H13NO2"  # a good one
mebo <- "C9H20NO2" # a good one
mebo <- "C10H16O3" # maybe
mebo <- "C3H8NO6P"  
mebo <- "C13H18O" # maybe
mebo <- "C10H12O5S"  
mebo <- "C11H21N3O5"
mebo <- "C16H14O6" # maybe
mebo <- "C19H16O10"
mebo <- "C26H44O5" # maybe
mebo <- "C27H46O5"

# after check pathway analysis
mebo <- "C18H32O2" # belong to linoleic acid metabolism
mebo <- "C18H32O3"
mebo <- "C5H10N2O3" # belong to the glutamate pathway
mebo <- "C5H9NO4"
mebo <- "C2H2O3" # belong to glycine, serine and threoline metabolism
mebo <- "C3H7NO3"
mebo <- "C3H7NO2S"
mebo <- "C4H9NO3"
mebo <- "C3H6O4"
mebo <- "C2H5NO2"
mebo <- "C4H9NO3"
mebo <- "C3H7NO2"

mebo <- "C20H32O2" # arachidonic acid metabolism
mebo <- "C20H34O4"
mebo <- "C2H2O3" # glyoxylate and dicarboxylate metabolism
mebo <- "C2H4O2"
mebo <- "C5H10N2O3"

mebo <- "C7H14N2O3" # theanine
mebo <- "C10H17N3O6S" # gluathione - could not measure
mebo <- "C5H9NO4" # glutamate

## prepare the ploting data --------------------------------------------------
metadat_boxplot <- inputDat %>% 
  select(season, responder, c(mebo), matches("_abFC|_d0|_d28|_reclassify")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season)) %>%
  mutate(abFC = ifelse(reclassify == "LH" | reclassify == "HH", "R", "NR")) %>%
  mutate(abBaseline = ifelse(reclassify == "LL" | reclassify == "LH", "low", "high"))

## make boxplot --------------------------------------------------
# based on reclassification 
boxplot_reClass  <- metadat_boxplot %>% 
  filter(strainSeason %in% strainSeasons) %>%
  ggboxplot(x = "reclassify", y = mebo,
            color = "#898366", add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  facet_wrap(~strainSeason, nrow = 2) +
  stat_compare_means(comparisons = compare_reClass, size = 5, method = "t.test")+
  theme(text = element_text(size = 18))

boxplot_reClass 

# based on abFC: NR vs R
boxplot_abFC <- metadat_boxplot %>% 
  filter(strainSeason %in% strainSeasons) %>%
  ggboxplot(x = "abFC", y = mebo,
            color = "#898366", add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  facet_wrap(~strainSeason, nrow = 2) + 
  stat_compare_means(comparisons = compare_abFC, size = 5, method = "t.test")+
  theme(text = element_text(size = 18))

boxplot_abFC

# based on abT1: low Ab vs high Ab
boxplot_abD0 <- metadat_boxplot %>% 
  filter(strainSeason %in% strainSeasons) %>%
  ggboxplot(x = "abBaseline", y = mebo,
            color = "#898366", add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  facet_wrap(~strainSeason, nrow = 2) + 
  stat_compare_means(comparisons = compare_abD0, size = 5, method = "t.test")+
  theme(text = element_text(size = 18))

boxplot_abD0 


## save the plot --------------------------------------------------
png(paste0("output/boxplotMetabolite_reClass_", mebo, ".png"), width = 720, height = 696)
boxplot_reClass 
dev.off()

png(paste0("output/boxplotMetabolite_abFC_", mebo, ".png"), width = 576, height = 696)
boxplot_abFC
dev.off()

png(paste0("output/boxplotMetabolite_abD0_", mebo, ".png"), width = 576, height = 696)
boxplot_abD0 
dev.off()


