library(ComplexHeatmap)
library(circlize)

# prepare data ----------------------------------------------------------------------------
load("inflamScore_NESfgsae.RData")

# annotation color -------------------------------------------------------------
col_cohort <- c("iMED" = "red", "ZirFlu" = "blue")
col_age <- c("young" = "blue", "old" = "red")
col_condition <- c("healthy" = "orange", 
                   "compensated cirrhosis" = "darkblue", 
                   "decompensated cirrhosis"  = "darkred")
col_reclasify <- c("LL" = "black", "LH" = "green", "HL" = "blue", "HH" = "red")
col_category <- c("NR" = "grey90", "Other" = "khaki3", "TR" = "peru")
col_inflam <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))

# input data (no NA) & metadata -------------------------------------------------
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  filter(disease == "healthy")

inputDat <- fgseaResTab %>% 
  lapply(function(x) x[["hallmarkSets"]] %>% rownames_to_column("pathway")) %>%
  purrr::reduce(full_join) %>%
  column_to_rownames("pathway") %>% select(metadat$name)

hallmarkSubsets_v0 <- c("HALLMARK_HYPOXIA", 
                        "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                        hallmarkSubsets_v1)
inputDat_selectedPathway <- inputDat[hallmarkSubsets_v0,]

metadat_NRvsTR <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  filter(disease == "healthy") %>%
  filter(category %in% c("NR", "TR"))

inputDat_NRvsTR <- fgseaResTab %>% 
  lapply(function(x) x[["hallmarkSets"]] %>% rownames_to_column("pathway")) %>%
  purrr::reduce(full_join) %>%
  column_to_rownames("pathway") %>% select(metadat_NRvsTR$name)
# heatmap -----------------------------------
identical(metadat$name, colnames(inputDat)) # check the sampleorder = TRUE
metadat <- metadat_NRvsTR
column_ha = HeatmapAnnotation(
  cohort = metadat$cohort,
  age_group = metadat$age_group,
  # condition = metadat$condition,
  # H1N1 = metadat$H1N1$reclassify,
  # H3N2 = metadat$H3N2$reclassify,
  # Bvic = metadat$Bvictoria$reclassify,
  # Byam = metadat$Byamagata$reclassify,
  category = metadat$category,
  col = list(
    cohort = col_cohort,
    age_group = col_age,
    # condition = col_condition,
    # H1N1 = col_reclasify, 
    # H3N2 = col_reclasify, 
    # Bvic = col_reclasify, 
    # Byam = col_reclasify,
    category = col_category))

inputDat %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)))

inputDat %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)

inputDat_selectedPathway %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)))

inputDat_selectedPathway %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)

inputDat_NRvsTR %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)

inputDat_NRvsTR[hallmarkSubsets_v0,] %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)

picked_pathways <- c("HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_SPERMATOGENESIS",
                     "HALLMARK_GLYCOLYSIS", "HALLMARK_ANGIOGENESIS",
                     "HALLMARK_APICAL_SURFACE", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"  ,
                     "HALLMARK_KRAS_SIGNALING_UP" ,"HALLMARK_ALLOGRAFT_REJECTION"    ,
                     "HALLMARK_CHOLESTEROL_HOMEOSTASIS"   ,"HALLMARK_E2F_TARGETS" ,
                     "HALLMARK_PROTEIN_SECRETION", "HALLMARK_INFLAMMATORY_RESPONSE" ,
                     "HALLMARK_IL2_STAT5_SIGNALING"
                     )
inputDat_NRvsTR[picked_pathways,] %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)

inputDat_NRvsTR[unique(c(picked_pathways, hallmarkSubsets_v0)), ] %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)

inputDat[unique(c(picked_pathways, hallmarkSubsets_v0)), ] %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)


lt <- apply(inputDat_selectedPathway, 2, 
            function(x) data.frame(density(x, na.rm = TRUE)[c("x", "y")]))
ha = rowAnnotation(foo = anno_joyplot(lt, width = unit(4, "cm"), 
                                      gp = gpar(fill = 1:10), transparency = 0.75))
plot(ha)
# can use NR vs. TR --> see clear patterns --> pick pathways


inflamScore_v2 <- list()
for (cohort in names(fgseaResTab)) {
  inflamScore_v2[[cohort]] <- fgseaResTab[[cohort]] %>% 
    lapply(function(x) x[["hallmarkSets"]] %>% 
             rownames_to_column("pathway") %>% 
             filter(pathway %in% c(picked_pathways, hallmarkSubsets_v0)) %>% # inflamScore use picked & literature pathways
             t() %>% as.data.frame() %>%
             dplyr::select(all_of(hallmarkSubsets_v2)) %>%  # inflamScore use picked & literature pathways
             rowSums() ) %>%
    as.data.frame() %>% rename_with(~paste0("inflamScore.", .x))
}

inflamScore_v2 <- list()
for (cohort in names(fgseaResTab)) {
  inflamScore_v2[[cohort]] <- fgseaResTab[[cohort]] %>% 
    lapply(function(x) x[["hallmarkSets"]] %>% rownames_to_column("pathway"))
}

a <- fgseaResTab %>% 
  lapply(function(x) x[["hallmarkSets"]] %>% 
           rownames_to_column("pathway")%>% 
           filter(pathway %in% c(picked_pathways, hallmarkSubsets_v0))) %>%
  purrr::reduce(full_join) %>%
  column_to_rownames("pathway")

a2 <- a %>% select(inflamScore_allCohort$name) %>% t() %>% rowMeans()
cor.test(as.numeric(a2), inflamScore_allCohort$inflamScore.hallmarkSubsets_v2)
