library(ComplexHeatmap)
library(circlize)

# prepare data ----------------------------------------------------------------------------
# select hallmarkt pathways
load("gsea_hallmarktSet.RData")
hallmarkSubsets_v0 <- c("HALLMARK_HYPOXIA", 
                        "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                        hallmarkSubsets_v1)

hallmarkSubsets <- hallmarkSets[hallmarkSubsets_v0]

hallmarkSubsets_existProteins <- hallmarkSubsets %>% 
  lapply(function(x) x[which(x %in% protein_Entrez$ENTREZID)])

hallmarkSubsets_existProteinsID <- hallmarkSubsets %>% 
  lapply(function(x) x <- protein_Entrez %>% filter(ENTREZID %in% x) %>% select(OlinkID))

hallmarkSubsets_existProteinsNames <- hallmarkSubsets %>% 
  lapply(function(x) x <- protein_Entrez %>% filter(ENTREZID %in% x) %>% select(Assay))

# input protein data (normalized data 2 cohorts)
inputDat <- get.datExpr(inputDat = proteinDat_impute,
                        metadat = metadat_healthy,
                        OlinkIDs = unique(unlist(hallmarkSubsets_existProteinsID)),
                        proteinAnnot = cohorts_dat$proteinAnnot_all)


# heatmap -----------------------------------
identical(metadat_healthy$name, colnames(inputDat)) # check the sampleorder = TRUE
metadat <- metadat_healthy

inputDat %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)))

column_ha = HeatmapAnnotation(
  age_group = metadat$age_group,
  category = metadat$category,
  col = list(
    age_group = col_age,
    category = col_category))

inputDat %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          show_column_dend = FALSE,
          row_dend_width = unit(6, "cm"),
          show_column_names = FALSE,
          top_annotation = column_ha)

column_ha = HeatmapAnnotation(
  cohort = metadat$cohort,
  age_group = metadat$age_group,
  condition = metadat$condition,
  H1N1 = metadat$H1N1_reclassify,
  H3N2 = metadat$H3N2_reclassify,
  Bvic = metadat$Bvictoria_reclassify2,
  Byam = metadat$Byamagata_reclassify2,
  category = metadat$category,
  col = list(
    cohort = col_cohort,
    age_group = col_age,
    condition = col_condition,
    H1N1 = col_reclasify,
    H3N2 = col_reclasify,
    Bvic = col_reclasify,
    Byam = col_reclasify,
    category = col_category))

inputDat %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)

# only NR vs TR
metadat_NRvsTR <- metadat %>% filter(category != "Other")
inputDat_NRvsTR <- inputDat %>% select(metadat_NRvsTR$name)
identical(metadat_NRvsTR$name, colnames(inputDat_NRvsTR)) # check the sampleorder = TRUE

column_ha = HeatmapAnnotation(
  cohort = metadat_NRvsTR$cohort,
  age_group = metadat_NRvsTR$age_group,
  category = metadat_NRvsTR$category,
  col = list(
    cohort = col_cohort,
    age_group = col_age,
    category = col_category))

inputDat_NRvsTR %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)

## Heatmap for average proteins -----------------------------------
picked_proteins <- c("CXCL3", "BANK1", "GMPR", "PRDX5", "PDGFB", "CXCL6", "CXCL1",
                     "IL6", "SLAMF7", "CXCL8", "CXCL9", "CXCL10",
                     "CCL4", "KYNU", "HSPA1A", "IL18", "IL17RB", "LAMP3", "TLR3", "SMPDL3A",
                     "GZMB", "TNFSF11", "PTH1R", "CCL20", "CCL24", "IL10")
picked_proteinsDat <- cohorts_dat$proteinAnnot_all %>% filter(Assay %in% picked_proteins)

proInput <- get.datExpr(inputDat = protein_normDat,
                        metadat = metadat_healthy,
                        OlinkIDs = picked_proteinsDat$OlinkID,
                        proteinAnnot = cohorts_dat$proteinAnnot_all)

heatmap_dat <- metadat_healthy %>% select(name, age_group, H1N1_reclassify) %>%
  full_join(proInput %>% t() %>% 
              as.data.frame() %>% rownames_to_column("name")) %>%
  group_by(H1N1_reclassify, age_group) %>% 
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  mutate(name2 = paste0(H1N1_reclassify, "_", age_group)) %>% 
  column_to_rownames("name2")

heatmap_dat[, -c(1,2)] %>% t() %>% Heatmap()
