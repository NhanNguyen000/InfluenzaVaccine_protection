library(limma)

get.proteinName <- function(OlinkIDs, OlinkDat) {
  # Aim: get protein names from thier Olink ID
  # input: OlinkIDs - vector of interested Olink ID protein, OlinkDat - Olink ID to protein name data 
  # output: protein names
  
  outcome <- OlinkDat %>% 
    filter(OlinkID %in% OlinkIDs)
  return(outcome$Assay)
}

# input data --------------------------
inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t()

# run limmar model --------------------------
identical(metadat_healthy$name, colnames(inputDat))

designTable <- model.matrix(~ + sex + age + H3N2_reclassify, metadat_healthy)
res <- lmFit(inputDat, design = designTable) %>% eBayes() 

# check result
View(res$p.value)
H3N2_reclassify_LLvsLH <-res$p.value %>% as.data.frame() %>% 
  select(H3N2_reclassifyLH) %>% filter(H3N2_reclassifyLH < 0.05)

H3N2_reclassify_LLvsHL <-res$p.value %>% as.data.frame() %>% 
  select(H3N2_reclassifyHL) %>% filter(H3N2_reclassifyHL < 0.05)

H3N2_reclassify_LLvsHH <-res$p.value %>% as.data.frame() %>% 
  select(H3N2_reclassifyHH) %>% filter(H3N2_reclassifyHH < 0.05)

library(ggvenn)
res.venn <- list("LLvsLH" = rownames(H3N2_reclassify_LLvsLH),
                 "LLvsHL" = rownames(H3N2_reclassify_LLvsHL),
                 "LLvsHH" = rownames(H3N2_reclassify_LLvsHH))
ggvenn(res.venn,  fill_color = c("#999999", "#E69F00", "#56B4E9"))


library(gplots)
res.venn2 <- list(
  "LLvsLH" = get.proteinName(rownames(H3N2_reclassify_LLvsLH), cohorts_dat$proteinAnnot_all),
  "LLvsHL" = get.proteinName(rownames(H3N2_reclassify_LLvsHL), cohorts_dat$proteinAnnot_all),
  "LLvsHH" = get.proteinName(rownames(H3N2_reclassify_LLvsHH), cohorts_dat$proteinAnnot_all))
v.table <- venn(res.venn2)
print(v.table)

## Heatmap for average proteins -----------------------------------
proInput <- get.datExpr(inputDat = protein_normDat,
                        metadat = metadat_healthy,
                        OlinkIDs = unique(unlist(res.venn)),
                        proteinAnnot = cohorts_dat$proteinAnnot_all)

heatmap_dat <- metadat_healthy %>% select(name, age_group, H1N1_reclassify) %>%
  full_join(proInput %>% t() %>% 
              as.data.frame() %>% rownames_to_column("name")) %>%
  group_by(H1N1_reclassify, age_group) %>% 
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  mutate(name2 = paste0(H1N1_reclassify, "_", age_group)) %>% 
  column_to_rownames("name2")

heatmap_dat[, -c(1,2)] %>% t() %>% Heatmap()

# heatmap per samples-----------------------------------
inputDat <- get.datExpr(inputDat = proteinDat_impute,
                        metadat = metadat_healthy,
                        OlinkIDs = unique(unlist(res.venn)),
                        proteinAnnot = cohorts_dat$proteinAnnot_all)


identical(metadat$name, colnames(inputDat)) # check the sampleorder = TRUE

inputDat %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)))

column_ha = HeatmapAnnotation(
  cohort = metadat$cohort,
  age_group = metadat$age_group,
  condition = metadat$condition,
  H1N1 = metadat$H1N1_reclassify,
  category = metadat$category,
  col = list(
    cohort = col_cohort,
    age_group = col_age,
    condition = col_condition,
    H1N1 = col_reclasify,
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
  H1N1 = metadat_NRvsTR$H1N1_reclassify,
  col = list(
    cohort = col_cohort,
    age_group = col_age,
    category = col_category,
    H1N1 = col_reclasify))

inputDat_NRvsTR %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)

# only LL vs LH
metadat_LLvsLH <- metadat %>% filter(H1N1_reclassify %in% c("LL", "LH"))
inputDat_LLvsLH <- inputDat %>% select(metadat_LLvsLH$name)
identical(metadat_LLvsLH$name, colnames(inputDat_LLvsLH)) # check the sampleorder = TRUE

column_ha = HeatmapAnnotation(
  cohort = metadat_LLvsLH$cohort,
  age_group = metadat_LLvsLH$age_group,
  category = metadat_LLvsLH$category,
  H1N1 = metadat_LLvsLH$H1N1_reclassify,
  col = list(
    cohort = col_cohort,
    age_group = col_age,
    category = col_category,
    H1N1 = col_reclasify))

inputDat_LLvsLH %>%
  Heatmap(heatmap_legend_param = list(at = c(-5, 0, 5)),
          top_annotation = column_ha)
