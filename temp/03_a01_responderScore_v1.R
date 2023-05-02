# Aim: pick proteins/pathways differnt between 3 categories (NR, other, TR) or 4 reclassification (HH, HL, LH, LL)
# to calculate the inflammation score (based on average protein expresison, average pathway NES, ....)

library(ComplexHeatmap)
library(circlize)

# global function & setup variable -------------------------------------------------------------------------
## functions -------------------------------------------------------------
get.datExpr <- function(inputDat, metadat, OlinkIDs, proteinAnnot) {
  # Aim: extract the protein data with selected samples and protein from a list of protein data per cohort
  
  # input: inputDat - a list of different dataframes, each dataframe is for one cohort.
  # metadat - metadata of interested sample, with "name" column overlap with the sample name in protein data
  # OlinkIDs - list of interested portein IDs 
  # proteinAnnot - a dataframe with OlinkID and protein names (Assay)
  
  # output: outcome - a protein expression table with protein name in rowname, sample in colnames
  
  inputTable <- inputDat %>% 
    lapply(function(x) x %>% rownames_to_column("name")) %>% 
    purrr::reduce(full_join) # convert a list of different table per cohort to one table
  
  table_selectSamples <- inputTable %>%
    filter(name %in% metadat$name) %>% arrange(match(name, metadat$name)) %>%
    column_to_rownames("name") %>% t()# select only samples on the metadata
  
  table_selectSamples_selectProteins <-  table_selectSamples[OlinkIDs, ]  # select only interested proteins
  
  outcome <- table_selectSamples_selectProteins %>%
    as.data.frame() %>% rownames_to_column("OlinkID") %>%
    left_join(proteinAnnot) %>% dplyr::select(-OlinkID, -UniProt) %>%
    column_to_rownames("Assay") # convert OlinkID to protein names
  
  return(outcome)
}

## annotation color -------------------------------------------------------------
col_cohort <- c("iMED" = "red", "ZirFlu" = "blue")
col_age <- c("young" = "blue", "old" = "red")
col_condition <- c("healthy" = "orange", 
                   "compensated cirrhosis" = "darkblue", 
                   "decompensated cirrhosis"  = "darkred")
col_reclasify <- c("LL" = "black", "LH" = "green", "HL" = "blue", "HH" = "red")
col_category <- c("NR" = "grey90", "Other" = "khaki3", "TR" = "peru")
col_inflam <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))

## focus on healthy subject -------------------------------------------------------------
# metadata for all healthy subjects
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# options------------------------------------------------------------------------------
## option 1: --------------------------
# inflammation hallmarks (from literature) + added more pathways --> proteins


## option 2: --------------------------
# unsupervised clustering --> pick pathways --> proteins 


## option 3: --------------------------
# limma: protein ~ 3 categories --> proteins


## option 4: --------------------------
# limma: protein ~ 4 reclassification --> protein per strains


## option 5: --------------------------
# limma: protein ~ low vs. high baseline --> protein per strains
# limma: protein ~ abFC  --> protein per strains

# calculate respond.Score = average expression of selected proteins -----------
avgProtein.score <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>% 
  filter(name %in% metadat_healthy$name) %>%
  column_to_rownames("name") %>% 
  select(unique(as.vector(unlist(res.venn)))) %>%
  #select(picked_proteinsDat$OlinkID) %>%
  #select(rownames(H1N1abFC_pro)) %>%
  #select(rownames(H1N1_T1pro)) %>%
  rowMeans(na.rm = TRUE)

metadat2 <- metadat_healthy %>% 
  full_join(as.data.frame(avgProtein.score) %>% rownames_to_column("name")) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("HL", "HH", "LL", "LH")))


library(ggpubr)
my_comparisons <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
my_comparisons_v2 <- list( c("HH", "HL"), c("LL", "LH") )
ggboxplot(metadat2, x = "H1N1_reclassify", y = "avgProtein.score",
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

ggboxplot(metadat2, x = "category", y = "avgProtein.score",
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

ggboxplot(metadat2, x = "H3N2_reclassify", y = "avgProtein.score",
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

# check specific proteins -----------------------------------------
selectPro <- cohorts_dat$proteinAnnot_all %>% 
  # filter(Assay %in% c("CXCL8", "IL6", "NBN", "CCL26", "TIMP3")) 
  # filter(Assay %in% c("CD83", "TNFRSF4", "SIT1", "ADAM23", "CXCL8", "SULT2A1"))
  filter(Assay %in% c("LSP1", "HSD11B1", "TNFRSF11A", "FGF19"))

proteinDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>% 
  filter(name %in% metadat_healthy$name) %>%
  column_to_rownames("name") %>% 
  select(selectPro$OlinkID)

if (identical(names(proteinDat), selectPro$OlinkID)) {
  names(proteinDat) <- selectPro$Assay
}

metadat2 <- metadat_healthy %>% 
  full_join(as.data.frame(proteinDat) %>% rownames_to_column("name")) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("HL", "HH", "LL", "LH")))


proId <- "NBN" 
proId <- "CCL26"
proId <- "IL6"
proId <- "CXCL8"
proId <- "TIMP3"

cowplot::plot_grid(
  ggboxplot(metadat2, x = "H1N1_reclassify", y = proId,
            paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  ggboxplot(metadat2, x = "category", y = proId,
            paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test"),
  nrow = 2
)

proId <- "CD83"
proId <- "TNFRSF4"
proId <- "SIT1"
proId <- "ADAM23"
proId <- "CXCL8"
proId <- "SULT2A1"
ggboxplot(metadat2, x = "H1N1_reclassify", y = proId,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")

proId <- "LSP1"
proId <- "HSD11B1"
proId <- "TNFRSF11A"
proId <- "FGF19"
cowplot::plot_grid(
  ggboxplot(metadat2, x = "H1N1_reclassify", y = proId,
            paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  ggboxplot(metadat2, x = "H3N2_reclassify", y = proId,
            paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
    stat_compare_means(comparisons = my_comparisons_v2, method = "t.test"),
  nrow = 2
)
