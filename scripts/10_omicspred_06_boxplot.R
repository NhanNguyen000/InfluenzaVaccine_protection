rm(list = ls())

library(tidyverse)
library(ggpubr)
# load data =======================================================================
load("cohorts_dat.RData")
load("omicsPred_predTab.RData")
load("cohorts_dat.RData")

## predicted value from OmicsPred -------------------
pred_dat <- pred_Tab %>% 
  t() %>% as.data.frame %>% rownames_to_column("OMICSPRED.ID") %>% 
  full_join(selected_models %>% dplyr::select(OMICSPRED.ID, model, type, Gene, Ensembl.ID))

# for FCRL3 (OPGS003959 - RNAseq, OPGS000406 - protein), TNFSF12(OPGS002615 - protein, OPGS006022 - RNAseq) , PTX3 (OPGS002455 - protein)
predDat_v2 <- pred_dat %>% as.data.frame() %>%
  filter(OMICSPRED.ID %in% c("OPGS003959", "OPGS000406", "OPGS002615", "OPGS006022", "OPGS002455")) %>%
  mutate(valName = ifelse(OMICSPRED.ID == "OPGS003959", "FCRL3_RNAseq",
                          ifelse(OMICSPRED.ID == "OPGS000406", "FCRL3_protein", 
                                 ifelse(OMICSPRED.ID == "OPGS002455", "PTX3_protein", 
                                        ifelse(OMICSPRED.ID == "OPGS006022", "TNFSF12_RNAseq", "TNFSF12_protein"))))) %>%
  dplyr::select(-OMICSPRED.ID, -model, -type, -Ensembl.ID, -Gene) %>% 
  column_to_rownames("valName") %>% t() %>% as.data.frame()


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


# boxplot ================================
compare_responder <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
compare_reClass <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
compare_abFC <- list(c("R", "NR"))
compare_abD0 <- list(c("low", "high"))

protein <- "FCRL3"
protein <- "PTX3"

# measured value - all subjects 
metadat_boxplot <- measureDat_allSubjects %>% 
  dplyr::select(season, responder, c(protein), matches("_abFC|_d0|_d28|_reclassify")) %>%
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_d0 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high")))

# based on reclassification 
ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

ggboxplot(metadat_boxplot, x = "B_reclassify", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +  
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

ggboxplot(metadat_boxplot, x = "Bvictoria_reclassify", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +  
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

ggboxplot(metadat_boxplot, x = "Byamagata_reclassify", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) +  
  stat_compare_means(comparisons = compare_reClass, method = "t.test")


# measured value - 159 subjects have genetic data
metadat_boxplot <- measureDat %>% 
  dplyr::select(season, responder, c(protein), matches("_abFC|_d0|_d28|_reclassify")) %>%
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_d0 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high")))

# based on reclassification 
ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = protein,
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

ggboxplot(metadat_boxplot, x = "B_reclassify", y = protein,
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = compare_reClass, method = "t.test")


# predicted value
metadat_boxplot_v2 <- predDat_v2 %>% rownames_to_column("probandID") %>% 
  right_join(metadata_healthy) %>% 
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_d0 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high")))

# based on reclassification 
ggboxplot(metadat_boxplot_v2, x = "H1N1_reclassify", y = paste0(protein, "_protein"),
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

ggboxplot(metadat_boxplot_v2, x = "B_reclassify", y = paste0(protein, "_protein"),
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = compare_reClass, method = "t.test")

## LL vs. other -----------------------------------------------
compare_reClass_v2 <- list( c("LL", "Other"))

# measured value - 159 subjects have genetic data
metadat_boxplot_v3 <- metadat_boxplot %>% 
  mutate(H1N1_reclassify_v2 = ifelse(H1N1_reclassify == "LL", "LL", "Other"),
         H3N2_reclassify_v2 = ifelse(H3N2_reclassify == "LL", "LL", "Other"),
         B_reclassify_v2 = ifelse(B_reclassify == "LL", "LL", "Other")) %>%
  mutate(H1N1_reclassify_v2 = factor(H1N1_reclassify_v2, levels = c("LL", "Other")),
         H3N2_reclassify_v2 = factor(H3N2_reclassify_v2, levels = c("LL", "Other")),
         B_reclassify_v2 = factor(B_reclassify_v2, levels = c("LL", "Other")))

# based on reclassification 
ggboxplot(metadat_boxplot_v3, x = "H1N1_reclassify_v2", y = protein,
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = compare_reClass_v2 , method = "t.test")

ggboxplot(metadat_boxplot_v3, x = "H3N2_reclassify_v2", y = protein,
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = compare_reClass_v2 , method = "t.test")

ggboxplot(metadat_boxplot_v3, x = "B_reclassify_v2", y = protein,
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = compare_reClass_v2 , method = "t.test")


# predicted value
metadat_boxplot_v4 <- metadat_boxplot_v2 %>% 
  mutate(H1N1_reclassify_v2 = ifelse(H1N1_reclassify == "LL", "LL", "Other"),
         H3N2_reclassify_v2 = ifelse(H3N2_reclassify == "LL", "LL", "Other"),
         B_reclassify_v2 = ifelse(B_reclassify == "LL", "LL", "Other"))

# based on reclassification 
ggboxplot(metadat_boxplot_v4, x = "H1N1_reclassify_v2", y = paste0(protein, "_protein"),
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = compare_reClass_v2, method = "t.test")

ggboxplot(metadat_boxplot_v4, x = "H3N2_reclassify_v2", y = paste0(protein, "_protein"),
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = compare_reClass_v2, method = "t.test")

ggboxplot(metadat_boxplot_v4, x = "B_reclassify_v2", y = paste0(protein, "_protein"),
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = compare_reClass_v2, method = "t.test")
