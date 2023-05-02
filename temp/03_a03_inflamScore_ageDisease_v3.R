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

#save(fgseaResTab, fgseaRes_table, file = "20230227_inflamScore.RData")
# load calculated data
load("20230227_inflamScore.RData")

# calculate Z score
fgseaRes_table_Zscore <-  fgseaRes_table %>% 
  t() %>% scale() %>% as.data.frame() %>% t()

# inflamScore use differnt hallmark pathways -------------------------------
# hallmarkSubsets_v2
picked_hallmarks <- hallmarkSubsets_v2 # pan-vaccine paper - figure 2
# picked_hallmarks <- c("HALLMARK_IL6_JAK_STAT3_SIGNALING",
#                       "HALLMARK_INFLAMMATORY_RESPONSE",
#                       "HALLMARK_TNFA_SIGNALING_VIA_NFKB") # pan-vaccine paper - figure 4
# picked_hallmarks <- c(hallmarkSubsets_v1, 
#                       "HALLMARK_HYPOXIA",
#                       "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")

inflamScore <- fgseaRes_table[picked_hallmarks, ] %>% colMeans(na.rm = TRUE)
inflamScore_Zscore <- fgseaRes_table_Zscore[picked_hallmarks, ] %>%
  colMeans(na.rm = TRUE)
inflamScore_Zscore_abs <- fgseaRes_table_Zscore[picked_hallmarks, ] %>% abs() %>%
  colMeans(na.rm = TRUE)

inflamScore_norm_allTab <- fgseaRes_table_Zscore[picked_hallmarks, ] -  mean(fgseaRes_table %>% as.matrix(), na.rm = TRUE)
inflamScore_Zscore_allTab <- inflamScore_norm_allTab %>% colMeans(na.rm = TRUE)

library(psych)
inflamGeomScore <- fgseaRes_table[picked_hallmarks, ] %>% geometric.mean(na.rm = TRUE)

# plot ------------------------
# boxplot all_together (healthy & cirrhosis)
metadat_boxplot <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH"))) %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("young","old",  "cirrhosis"))) %>%
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

my_comparisons_group <- list(c("young", "old"), c("old", "cirrhosis"), c("young", "cirrhosis"))

ggboxplot(metadat_boxplot, x = "group", y = score,
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons_group ,
                     method = "t.test")

# compare_means(inflamScore ~ group, data = metadat_boxplot, ref.group = ".all.", method = "t.test")
stat.test <- compare_means(inflamScore ~ group, data = metadat_boxplot, 
                           comparisons = my_comparisons_group, 
                           method = "t.test")
stat.test <- compare_means(inflamScore_Zscore ~ group, data = metadat_boxplot, 
                           comparisons = my_comparisons_group, 
                           method = "t.test")
stat.test <- compare_means(inflamGeomScore ~ group, data = metadat_boxplot, 
                           comparisons = my_comparisons_group, 
                           method = "t.test")
ggboxplot(metadat_boxplot, x = "group", y = score,
          paletter = "jco", add = "jitter") + 
  stat_pvalue_manual(stat.test, label = "p.adj", 
                     # y.position = c(1.3, 1.7, 1.5))
                     y.position = c(1.6, 1.8, 1.7))

# compare abFC 
stat.test <- compare_means(H1N1_abFC_combine ~ group, data = metadat_boxplot, 
                           comparisons = my_comparisons_group, 
                           method = "t.test")
ggboxplot(metadat_boxplot, x = "group", y = "H1N1_abFC_combine",
          paletter = "jco", add = "jitter") + 
  stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(13, 17, 15))

stat.test <- compare_means(H3N2_abFC_combine ~ group, data = metadat_boxplot, 
                           comparisons = my_comparisons_group, 
                           method = "t.test")
ggboxplot(metadat_boxplot, x = "group", y = "H3N2_abFC_combine",
          paletter = "jco", add = "jitter") + 
  stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(15, 19, 17))

# compare baseline
stat.test <- compare_means(H1N1_T1_combine ~group, data = metadat_boxplot,
                           comparisons = my_comparisons_group,
                           method = "t.test")
ggboxplot(metadat_boxplot, x = "group", y = "H1N1_T1_combine",
          palette = "jco", add = "jitter") +
  stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(12, 14, 13))

stat.test <- compare_means(H3N2_T1_combine ~group, data = metadat_boxplot,
                           comparisons = my_comparisons_group,
                           method = "t.test")
ggboxplot(metadat_boxplot, x = "group", y = "H3N2_T1_combine",
          palette = "jco", add = "jitter") +
  stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(11, 13, 12))

# heatmap for average inflamScore -----------
heatmap_dat <- metadat_boxplot %>% select(name, group) %>%
  left_join(fgseaRes_table %>% t() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  # left_join(fgseaRes_table_Zscore %>% t() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  group_by(group) %>%
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  column_to_rownames("group") %>%
  select(picked_hallmarks)

fgseaRes_table_Zscore_v2 <-  fgseaRes_table %>%  scale() %>% as.data.frame()

heatmap_dat %>% t() %>% Heatmap()
heatmap_dat %>% t() %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))
heatmap_dat %>% scale() %>% t() %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))
fgseaRes_table %>% Heatmap()
fgseaRes_table_Zscore %>% Heatmap()

# inflam.Score to abFC ------------------------------
metadat_boxplot %>%
  ggplot(aes(x = H1N1_abFC_combine, y = inflamScore)) + 
  geom_point(aes(col = group),
    position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.25)) +
  stat_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme_bw()

metadat_boxplot %>%
  ggplot(aes(x = H1N1_abFC_combine, y = inflamScore, col = group)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.25)) +
  stat_smooth(method = "lm", se = FALSE) + stat_cor(label.x = 7.5, label.y = c(-0.2, -0.4, -0.6)) +
  theme_bw() + theme(legend.position = "top")

metadat_boxplot %>%
  ggplot(aes(x = H3N2_abFC_combine, y = inflamScore, col = group)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25, jitter.height = 0.25)) +
  stat_smooth(method = "lm", se = FALSE) + stat_cor(label.x = 7.5, label.y = c(-0.2, -0.4, -0.6)) +
  theme_bw() + theme(legend.position = "top")
# 