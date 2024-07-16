rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data =======================================================================
load("processedDat/cohorts_dat.RData")
load("processedDat/omicsPred_output.RData")

## predicted value from OmicsPred -------------------
pred_dat <- pred_Tab %>% 
  t() %>% as.data.frame %>% rownames_to_column("OMICSPRED.ID") %>% 
  full_join(selected_models %>% dplyr::select(OMICSPRED.ID, model, type, Gene, Ensembl.ID))

## metadata for healthy subjects -------------------------
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              dplyr::select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
  filter(condition == "Healthy", probandID %in% names(pred_dat)) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

measureDat <- protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% right_join(metadata_healthy)

measureDat_allSubjects <-  protein_Dat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>% 
  right_join(cohorts$HAI_all %>% 
               full_join(cohorts$donorInfo_all %>% 
                           dplyr::select(probandID, season, cohort, sex, age, condition)) %>%
               left_join(cohorts$donorSample_all %>% filter(time == "d0")) %>%
               filter(condition == "Healthy") %>%
               mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))))


# make boxplot ------------------------------------------------------------------
compare_responder <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
compare_reClass <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
compare_abFC <- list(c("R", "NR"))
compare_abD0 <- list(c("low", "high"))

protein <- "CCL7"

## measured value - all subjects ------------------------------------------------------------------
strainSeasons <- c("H1N1_2014", "H1N1_2015", "H1N1_2019", "H1N1_2020",
                   "B_2014", "B_2015", "H3N2_2015", "Byamagata_2020") # season have at least >= 3 participant per reclassification group

metadat_boxplot <- measureDat_allSubjects %>% 
  dplyr::select(season, responder, c(protein), matches("_abFC|_d0|_d28|_reclassify")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season)) %>%
  mutate(abFC = ifelse(reclassify == "LH" | reclassify == "HH", "R", "NR")) %>%
  mutate(abBaseline = ifelse(reclassify == "LL" | reclassify == "LH", "low", "high"))

boxplot_reClass <- metadat_boxplot %>% 
  filter(strainSeason %in% strainSeasons) %>%
  ggboxplot(x = "reclassify", y = protein,
            color = "#898366", add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  facet_wrap(~strainSeason, nrow = 2) +
  stat_compare_means(comparisons = compare_reClass, size = 5, method = "t.test")+
  theme(text = element_text(size = 18))

## save the plot
png(paste0("output/boxplotProtein_reClass_", protein, ".png"), width = 720, height = 696)
boxplot_reClass 
dev.off()


## measured value and predicted value - 159 subjects have genetic data ------------------------------------------------------------------
model_name <- "Somalogic"

metadat_boxplot <- measureDat %>% 
  # add prediction value from omicsPred
  inner_join(
    as.data.frame(pred_dat %>% 
      filter(Gene %in% protein, model %in% model_name) %>% 
      mutate(varName =  paste0(Gene, "_", model)) %>% 
      dplyr::select(-Gene, -model, -OMICSPRED.ID, -type, -Ensembl.ID) %>% 
        column_to_rownames("varName") %>% t()) %>% 
  tibble::rownames_to_column(var = "probandID")) %>%
  # arrange the data format
  dplyr::select( 
    season, responder, c(protein, paste0(protein, "_", model_name)), 
    matches("_abFC|_d0|_d28|_reclassify")) %>%
  pivot_longer(matches("reclassify"), names_to = "strain", values_to = "reclassify") %>%
  drop_na(reclassify) %>% 
  mutate(strain = gsub("_reclassify", "", strain),
         strainSeason = paste0(strain, "_", season)) %>%
  mutate(abFC = ifelse(reclassify == "LH" | reclassify == "HH", "R", "NR")) %>%
  mutate(abBaseline = ifelse(reclassify == "LL" | reclassify == "LH", "low", "high"))

# based on reclassification 
boxplot_reClass_measure <- metadat_boxplot %>% 
  ggboxplot(x = "reclassify", y = protein,
            color = "#898366", add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  facet_wrap(~strainSeason, nrow = 1) +
  stat_compare_means(comparisons = compare_reClass, size = 5, method = "t.test")+
  theme(text = element_text(size = 18))

boxplot_reClass_predict <- metadat_boxplot %>% 
  ggboxplot(x = "reclassify", y = paste0(protein, "_", model_name),
            color = "#898366", add = "jitter", add.params = list(size = 3, alpha = 0.5)) + 
  facet_wrap(~strainSeason, nrow = 1) +
  stat_compare_means(comparisons = compare_reClass, size = 5, method = "t.test")+
  theme(text = element_text(size = 18))

## save the plot
png(paste0("output/boxplotProtein_reClass_", protein, "_159subjects_season2015.png"), width = 720, height = 432)
boxplot_reClass_measure
dev.off()

png(paste0("output/omicPred_boxplotProtein_reClass_", protein, "_159subjects_season2015.png"), width = 720, height = 432)
boxplot_reClass_predict
dev.off()

