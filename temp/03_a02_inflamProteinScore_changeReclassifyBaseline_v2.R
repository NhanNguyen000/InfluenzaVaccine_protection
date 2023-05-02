# libraries & functions
# load libraries & functions --------------------------
library(ggpubr)

## focus on healthy subject -------------------------------------------------------------
# metadata for all healthy subjects
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

metadat_healthy %>% count(age_group, category)

# input data
inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>%
  column_to_rownames("name")

# metadata for all subjects
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>%
  filter(name %in% metadat$name) %>%
  column_to_rownames("name")
# inflam.Score for 3 categories (NR, Other, TR) -----------------------------------


# calculate inflammation score - using average (NES values of inflammation hallmark sets from GSEA) at baseline  -----------------------
load("gsea_hallmarktSet.RData") # load hallmark gene sets (prepare from from MSigDB) & convert protein - Entrez gene ID

# inflamScore use differnt hallmark pathways -------------------------------
# hallmarkSubsets_v2
picked_hallmarks <- hallmarkSubsets_v2 # pan-vaccine paper --> 58 proteins
picked_hallmarks <- c(hallmarkSubsets_v1, 
                      "HALLMARK_HYPOXIA",
                      "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
picked_proteins <- protein_Entrez %>% 
  filter(ENTREZID %in% unique(unlist(hallmarkSets[picked_hallmarks]))) # 84 proteins

inflamPro.Score <- inputDat[, picked_proteins$OlinkID] %>% rowMeans(na.rm = TRUE)

inputDat_Zscore <-  inputDat %>% scale() 
inflamPro.Zscore <- inputDat_Zscore[, picked_proteins$OlinkID] %>% rowMeans(na.rm = TRUE)
inflamPro.Zscore_abs <- inputDat_Zscore[, picked_proteins$OlinkID] %>% abs() %>% rowMeans(na.rm = TRUE)

library(psych)
inflamPro.GeomScore <- inputDat[, picked_proteins$OlinkID] %>% t() %>% geometric.mean(na.rm = TRUE)
# plots --------------------------------------------------------------
# boxplot
my_comparisons <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("HL", "HH"), c("LL", "HH"))

# boxplot all_together (healthy & cirrhosis)
metadat_boxplot <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("old", "young", "cirrhosis"))) %>%
  left_join(as.data.frame(inflamPro.Score) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamPro.Zscore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamPro.Zscore_abs) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamPro.GeomScore) %>% rownames_to_column("name"))

score <- "inflamPro.Score"
score <- "inflamPro.Zscore"
score <- "inflamPro.Zscore_abs"
score <- "inflamPro.GeomScore"

ggboxplot(metadat_boxplot, x = "category", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  stat_compare_means(method = "anova", label.y = 2.3) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label.y = c(1.3, 1.5, 1.7))
  # stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  stat_compare_means(method = "anova", label.y = 2.3) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test", label.y = c(1.3, 1.5, 1.7, 1.9))
  # stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H3N2_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  stat_compare_means(method = "anova", label.y = 2.3) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test", label.y = c(1.3, 1.5, 1.7, 1.9))
  # stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")
# save(fgseaResTab, fgseaRes_table, file = "20230227_inflamScore.RData")
