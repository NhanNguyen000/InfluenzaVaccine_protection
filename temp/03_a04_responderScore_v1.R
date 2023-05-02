library(limma)
library(psych)

get.proteinName <- function(OlinkIDs, OlinkDat) {
  # Aim: get protein names from thier Olink ID
  # input: OlinkIDs - vector of interested Olink ID protein, OlinkDat - Olink ID to protein name data 
  # output: protein names
  
  outcome <- OlinkDat %>% 
    filter(OlinkID %in% OlinkIDs)
  return(outcome$Assay)
}

# responderScore ------------------------------
# all subjects (healhty + cirrhosis)
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

# calculate score based on DEPs ---------------
pickProteins <- intersect(res.venn$H1N1, res.venn$B) # LL vs. Protector, 10 proteins
pickProteins <- unique(c(res.venn$H1N1, res.venn$B)) # LL vs. Protector, 57 proteins

picked_proteinsDat <- cohorts_dat$proteinAnnot_all %>% 
  filter(Assay %in% pickProteins)

picked_proteins <- picked_proteinsDat$OlinkID

responsePro.Score <- inputDat[, picked_proteins] %>% rowMeans(na.rm = TRUE)
inputDat_Zscore <-  inputDat %>% scale() 
responsePro.Zscore <- inputDat_Zscore[, picked_proteins] %>% rowMeans(na.rm = TRUE)
responsePro.Zscore_abs <- inputDat_Zscore[, picked_proteins] %>% abs() %>% rowMeans(na.rm = TRUE)


responsePro.GeomScore <- inputDat[, picked_proteins] %>% t() %>% geometric.mean(na.rm = TRUE)

# calculate score based on elasitic model ---------------
pickProteins <- a3 # elastic model,  LL vs.LH, HL, HH, 94 proteins
pickProteins <- a3 # elastic model,  LL vs. Protector, 126 proteins

picked_proteinsDat <- cohorts_dat$proteinAnnot_all %>% 
  right_join(pickProteins %>% rename("Assay" = "proteins"))

picked_proteins <- picked_proteinsDat$OlinkID

dat_v1 <- inputDat[, picked_proteins] %>% t()
identical(rownames(dat_v1), picked_proteinsDat$OlinkID)

dat_v2 <- dat_v1 * picked_proteinsDat$avg_coef

responsePro.Score <- dat_v2 %>% t() %>% rowMeans(na.rm = TRUE)

dat_v2_Zscore <-  dat_v2 %>% t() %>% scale() 
responsePro.Zscore <- dat_v2_Zscore %>% rowMeans(na.rm = TRUE)
responsePro.Zscore_abs <- dat_v2_Zscore %>% abs() %>% rowMeans(na.rm = TRUE)

responsePro.GeomScore <- dat_v2 %>% geometric.mean(na.rm = TRUE)

# plots --------------------------------------------------------------
# boxplot
metadat_boxplot <- metadat_healthy %>% 
  left_join(as.data.frame(responsePro.Score) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.Zscore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.Zscore_abs) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.GeomScore) %>% rownames_to_column("name"))

#my_comparisons <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("HL", "HH"), c("LL", "HH"))
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HH"))

score <- "responsePro.Score"
score <- "responsePro.Zscore"
score <- "responsePro.Zscore_abs"
score <- "responsePro.GeomScore"

# ggboxplot(metadat_boxplot, x = "category", y = score,
#           paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
#   # stat_compare_means(comparisons = my_comparisons, method = "t.test")
#   stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

ggboxplot(metadat_boxplot, x = "H3N2_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

metadat_boxplot %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "B_reclassify", y = score,
            paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

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

# boxplot all_together (healthy & cirrhosis)
metadat_boxplot <- metadat %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("old", "young", "cirrhosis"))) %>%
  left_join(as.data.frame(responsePro.Score) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.Zscore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.Zscore_abs) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.GeomScore) %>% rownames_to_column("name"))

score <- "responsePro.Score"
score <- "responsePro.Zscore"
score <- "responsePro.Zscore_abs"
score <- "responsePro.GeomScore"

ggboxplot(metadat_boxplot, x = "category", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~group) + 
  stat_compare_means(method = "anova", label.y = 2.4) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label.y = c(1.1, 1.7, 2))

# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  # stat_compare_means(method = "anova", label.y = 2.8) + # y = 3.3 or 2.8
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

ggboxplot(metadat_boxplot, x = "H3N2_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  # stat_compare_means(method = "anova", label.y = 2.8) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")
# stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

metadat_boxplot %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "B_reclassify", y = score,
            paletter = "jco", add = "jitter") + facet_wrap(~group) +
  # stat_compare_means(method = "anova", label.y = 3.3) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "Bvictoria_reclassify", y = score,
            paletter = "jco", add = "jitter") + facet_wrap(~group) +
  stat_compare_means(method = "anova", label.y = 3.3) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "Byamagata_reclassify", y = score,
            paletter = "jco", add = "jitter") + facet_wrap(~group) +
  stat_compare_means(method = "anova", label.y = 3.3) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

# heatmap of DEPs ------------------------
library(ComplexHeatmap)
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("young", "old","cirrhosis")))

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>%
  filter(name %in% metadat$name) %>%
  column_to_rownames("name")


heatmap_H1N1 <- metadat %>% select(name, group, H1N1_reclassify) %>%
  left_join(inputDat %>% rownames_to_column("name")) %>%
  group_by(group, H1N1_reclassify) %>%
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  mutate(name2 = paste0(group, "_", H1N1_reclassify)) %>% 
  column_to_rownames("name2") %>%
  select(unique(unlist(resSigPro$H1N1))) %>% t()

heatmap_H1N1_proName <- heatmap_H1N1 %>% as.data.frame() %>% 
  rownames_to_column("OlinkID") %>% 
  left_join(cohorts_dat$proteinAnnot_all %>% select(OlinkID, Assay)) %>%
  column_to_rownames("Assay") %>% select(-OlinkID)

heatmap_B <- metadat %>% select(name, group, B_reclassify) %>% drop_na(B_reclassify) %>%
  left_join(inputDat %>% rownames_to_column("name")) %>%
  group_by(group, B_reclassify) %>%
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  mutate(name2 = paste0(group, "_", B_reclassify)) %>% 
  column_to_rownames("name2") %>%
  select(unique(unlist(resSigPro$B))) %>% t()

heatmap_B_proName <- heatmap_B %>% as.data.frame() %>% 
  rownames_to_column("OlinkID") %>% 
  left_join(cohorts_dat$proteinAnnot_all %>% select(OlinkID, Assay)) %>%
  column_to_rownames("Assay") %>% select(-OlinkID)

heatmap_dat <- heatmap_H1N1_proName
heatmap_dat <- heatmap_B_proName
heatmap_dat %>% Heatmap()
heatmap_dat %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))
heatmap_dat %>% t() %>% scale() %>% t() %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))



heatmap_dat %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"),
          right_annotation = row_ha)

H1N1_annot <- resSig$H1N1 %>% 
  lapply(function(x) x %>% rownames_to_column("OlinkID")) %>% 
  purrr::reduce(full_join) %>%
  arrange(match(OlinkID, rownames(heatmap_H1N1))) %>% # TRUE, identical(rownames(H1N1_annot), rownames(heatmap_H1N1))
  column_to_rownames("OlinkID") %>% 
  lapply(function(x) x %>% as.logical)

row_ha <- rowAnnotation(
  DEP_LLvsLH = H1N1_annot$H1N1_reclassifyLH,
  DEP_LLvsHL = H1N1_annot$H1N1_reclassifyHL,
  DEP_LLvsHH = H1N1_annot$H1N1_reclassifyHH,
  col = list(
    DEP_LLvsLH = c("TRUE" = "black"),
    DEP_LLvsHL = c("TRUE" = "black"),
    DEP_LLvsHH = c("TRUE" = "black"))
)

heatmap_H1N1_proName %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"),
          right_annotation = row_ha)

heatmap_H1N1_proName %>% t() %>% scale() %>% t() %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"),
          right_annotation = row_ha)

B_annot <- resSig$B %>% 
  lapply(function(x) x %>% rownames_to_column("OlinkID")) %>% 
  purrr::reduce(full_join) %>%
  arrange(match(OlinkID, rownames(heatmap_B))) %>% # TRUE, identical(rownames(B_annot), rownames(heatmap_B))
  column_to_rownames("OlinkID") %>% 
  lapply(function(x) x %>% as.logical)

row_ha <- rowAnnotation(
  DEP_LLvsLH = B_annot$B_reclassifyLH,
  DEP_LLvsHL = B_annot$B_reclassifyHL,
  DEP_LLvsHH = B_annot$B_reclassifyHH,
  col = list(
    DEP_LLvsLH = c("TRUE" = "black"),
    DEP_LLvsHL = c("TRUE" = "black"),
    DEP_LLvsHH = c("TRUE" = "black"))
)

heatmap_B_proName %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"),
          right_annotation = row_ha)

heatmap_B_proName %>% t() %>% scale() %>% t() %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"),
          right_annotation = row_ha)


