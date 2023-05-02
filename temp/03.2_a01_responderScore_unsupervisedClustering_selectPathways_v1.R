library(ggpubr)
# prepare data ----------------------------------------------------------------------------
load("inflamScore_NESfgsae.RData")

inputDat <- fgseaResTab %>% 
  lapply(function(x) x[["hallmarkSets"]] %>% rownames_to_column("pathway")) %>%
  purrr::reduce(full_join) %>%
  column_to_rownames("pathway") %>% select(metadat_healthy$name)

hallmarkSubsets_v0 <- c("HALLMARK_HYPOXIA", 
                        "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                        hallmarkSubsets_v1)
picked_pathways <- c("HALLMARK_XENOBIOTIC_METABOLISM", 
                     "HALLMARK_SPERMATOGENESIS",
                     "HALLMARK_GLYCOLYSIS", 
                     "HALLMARK_ANGIOGENESIS",
                     "HALLMARK_APICAL_SURFACE", 
                     "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                     "HALLMARK_KRAS_SIGNALING_UP",
                     "HALLMARK_ALLOGRAFT_REJECTION",
                     "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                     "HALLMARK_E2F_TARGETS",
                     "HALLMARK_PROTEIN_SECRETION", 
                     "HALLMARK_INFLAMMATORY_RESPONSE",
                     "HALLMARK_IL2_STAT5_SIGNALING")

inputDat_selectedPathway <- inputDat[picked_pathways,]
temp.Score <- inputDat_selectedPathway %>% colMeans(na.rm = TRUE)

metadat2 <- metadat_healthy %>% 
  full_join(as.data.frame(temp.Score) %>% rownames_to_column("name")) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  left_join(HAIreclassify$all_cohorts) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("HL", "HH", "LL", "LH")))

# boxplot ----------------------------------------------------------
my_comparisons <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
ggboxplot(metadat2, x = "category", y = "temp.Score",
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

my_comparisons_v2 <- list( c("HH", "HL"), c("LL", "LH") )

ggboxplot(metadat2, x = "H1N1_reclassify", y = "temp.Score",
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

ggboxplot(metadat2, x = "H3N2_reclassify", y = "temp.Score",
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")
