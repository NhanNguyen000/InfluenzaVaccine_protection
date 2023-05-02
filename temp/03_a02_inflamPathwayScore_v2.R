# libraries & functions
# load libraries & functions --------------------------
library(fgsea)
#library('org.Hs.eg.db')
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

get.fgsea_singleInput <- function(dat, col_input, col_geneID) {
  # Description: prepare the input for fgsea
  
  # Arguments: 
  # dat - data with gene ID in one column and gene value/expression in other columns
  # col_input - name of the input column
  # col_geneID - name of the column present gene ID
  
  # Returns: fgsaInput - input vector for the fgsea analysis
  
  fgseaInput <- dat[, col_input]
  names(fgseaInput) <- dat[, col_geneID]
  
  return(fgseaInput)
}

get.fgseaRes <- function(dat, convertID, hallmarkSet) {
  # Description: run fgsea (GSEA) for all subjects/patients
  
  # Arguments: 
  # dat - data with patients in row, proteins in column
  # convertID - table to convert protein OlinkID to gene Entrez ID
  # hallmarkSet - the hallmark / reference gene set to run GSEA
  
  # Returns: fgsaRes - a list of fgsea result for each patient
  
  fgseaRes <- list()
  for (patient in rownames(dat)) {
    
    ranks_temp <- get.fgsea_singleInput(
      dat = dat %>% t() %>% as.data.frame() %>% 
        rownames_to_column("OlinkID") %>% left_join(convertID), 
      col_input = patient, col_geneID = "ENTREZID")
    
    fgseaRes[[patient]] <- fgsea(hallmarkSet, ranks_temp) %>% 
      as_tibble() %>% dplyr::select(pathway, size, leadingEdge, NES) %>%
      mutate(name = patient)
  }
  
  return(fgseaRes)
}

## focus on healthy subject -------------------------------------------------------------
# metadata for all healthy subjects
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

metadat_healthy %>% count(age_group, category)

# metadata forcirrhosis subjects
metadat_cirrhosis <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify$all_cohorts) %>%
  filter(disease == "cirrhosis") %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

metadat_cirrhosis %>% count(age_group, category)
# inflam.Score for 3 categories (NR, Other, TR) -----------------------------------


# calculate inflammation score - using average (NES values of inflammation hallmark sets from GSEA) at baseline  -----------------------
load("proteinDat_impute.RData") # load imputed protein data using median impute
load("gsea_hallmarktSet.RData") # load hallmark gene sets (prepare from from MSigDB) & convert protein - Entrez gene ID

# run code for all cohort
fgseaResTab <- list()
for (cohort in names(proteinDat_impute)) {
  inputDat <- proteinDat_impute[[cohort]]
  
  fgseaRes <- get.fgseaRes(dat = inputDat,
                           convertID = protein_Entrez, 
                           hallmarkSet = hallmarkSets) 
  
  fgseaResTab[[cohort]] <- fgseaRes %>% 
    purrr::reduce(full_join) %>% # long format
    select(-size, -leadingEdge) %>%
    pivot_wider(names_from = name, values_from = NES) %>% # convert to wide format
    column_to_rownames("pathway")
  
}

fgseaRes_table <- fgseaResTab %>% 
  lapply(function(x) x %>% rownames_to_column("pathway")) %>%
  purrr::reduce(full_join) %>% column_to_rownames("pathway")

# load calculated data 
load("20230227_inflamScore.RData")

# calculate Z score
fgseaRes_table_Zscore <-  fgseaRes_table %>% 
  t() %>% scale() %>% as.data.frame() %>% t()
# inflamScore use differnt hallmark pathways -------------------------------
# hallmarkSubsets_v2
picked_hallmarks <- hallmarkSubsets_v2 # pan-vaccine paper - figure 2
picked_hallmarks <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING",
                      "HALLMARK_INFLAMMATORY_RESPONSE",
                      "HALLMARK_TNFA_SIGNALING_VIA_NFKB") # pan-vaccine paper - figure 4
picked_hallmarks <- c(hallmarkSubsets_v1, 
                      "HALLMARK_HYPOXIA",
                      "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")

inflamScore <- fgseaRes_table[picked_hallmarks, ] %>% colMeans(na.rm = TRUE)
inflamScore_Zscore <- fgseaRes_table_Zscore[picked_hallmarks, ] %>%
  colMeans(na.rm = TRUE)
inflamScore_Zscore_abs <- fgseaRes_table_Zscore[picked_hallmarks, ] %>% abs() %>%
  colMeans(na.rm = TRUE)

inflamScore_norm_allTab <- fgseaRes_table_Zscore[picked_hallmarks, ] -  mean(fgseaRes_table %>% as.matrix(), na.rm = TRUE)
inflamScore_Zscore_allTab <- inflamScore_norm_allTab %>% colMeans(na.rm = TRUE)

library(psych)
inflamGeomScore <- fgseaRes_table[picked_hallmarks, ] %>% geometric.mean(na.rm = TRUE)
# plots --------------------------------------------------------------
# heatmap for average inflamScore
heatmap_dat <- metadat_healthy %>% select(name, age_group, category) %>%
  left_join(fgseaRes_table %>% t() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  # left_join(fgseaRes_table_Zscore %>% t() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  group_by(category, age_group) %>%
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  mutate(name2 = paste0(category, "_", age_group)) %>% 
  column_to_rownames("name2") %>%
  select(picked_hallmarks)

heatmap_dat %>% t() %>% Heatmap() # scale column
fgseaRes_table %>% Heatmap()
fgseaRes_table_Zscore %>% Heatmap()

# boxplot
my_comparisons <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"), c("HL", "HH") )

metadat_boxplot <- metadat_healthy %>% 
  left_join(as.data.frame(inflamScore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamScore_Zscore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamScore_Zscore_abs) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamScore_Zscore_allTab) %>% rownames_to_column("name"))

score <- "inflamScore"
score <- "inflamScore_Zscore"
score <- "inflamScore_Zscore_abs"
score <- "inflamScore_Zscore_allTab"

ggboxplot(metadat_boxplot, x = "category", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(method = "anova", label.y = 2.1) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggboxplot(metadat_boxplot, x = "category", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = score,
          paletter = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H3N2_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H3N2_reclassify", y = score,
          paletter = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

metadat_boxplot %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "B_reclassify", y = score,
          paletter = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "Bvictoria_reclassify", y = score,
            paletter = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "Byamagata_reclassify", y = score,
            paletter = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

# boxplot cirrhosis
metadat_boxplot <- metadat_cirrhosis %>% 
  left_join(as.data.frame(inflamScore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamScore_Zscore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamScore_Zscore_abs) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamScore_Zscore_allTab) %>% rownames_to_column("name"))

ggboxplot(metadat_boxplot, x = "category", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "category", y = score,
          paletter = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = score,
          paletter = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H3N2_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

# boxplot all_together (healthy & cirrhosis)
metadat_boxplot <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  # full_join(HAIreclassify$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  # mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("young", "old", "cirrhosis"))) %>%
  left_join(as.data.frame(inflamScore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamScore_Zscore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamScore_Zscore_abs) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamScore_Zscore_allTab) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(inflamGeomScore) %>% rownames_to_column("name"))

score <- "inflamScore"
score <- "inflamScore_Zscore"
score <- "inflamScore_Zscore_abs"
score <- "inflamScore_Zscore_allTab"
score <- "inflamGeomScore"

compare_means(inflamScore ~ group,  data = metadat_boxplot, ref.group = "young",
              method = "t.test")
ggboxplot(metadat_boxplot, x = "group", y = score,
          paletter = "jco", add = "jitter") + 
  stat_compare_means(paired = TRUE, method = "t.test")


# save(fgseaResTab, fgseaRes_table, file = "20230227_inflamScore.RData")
