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
# metadata for all healthy subjects
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t()

# run limma model for LL vs. (LH, HH) --------------------------
metadat_temp <- metadat_healthy %>% 
  filter(H1N1_reclassify %in% c("LL", "LH", "HH")) %>%
  mutate(H1N1_reclassify_LLvsResponder = ifelse(H1N1_reclassify == "LL", "LL", "responder"))

inputDat_temp <- inputDat[, metadat_temp$name]
identical(metadat_temp$name, colnames(inputDat_temp))

res <- lmFit(inputDat_temp, 
             design =  model.matrix(~ + sex + age + H1N1_reclassify_LLvsResponder, metadat_temp)) %>% 
  eBayes() 

resSig <- res$p.value %>% as.data.frame() %>%  select(4) %>% filter(. <0.05) 
#top.table <- topTable(res, coef= "H1N1_reclassify_LLvsResponderresponder", number = 17) # check again
# res %>% decideTests() %>% summary() # apply for padj, we do not have sig. padj
resSigPro <- rownames(resSig)
get.proteinName(OlinkIDs = resSigPro, OlinkDat = cohorts_dat$proteinAnnot_all)

get.proteinName(OlinkIDs = "OID20624", OlinkDat = cohorts_dat$proteinAnnot_all)

volcanoplot(res, coef= "H1N1_reclassify_LLvsResponderresponder", names = rownames(res), 
            highlight = 5)

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

metadat_heatmap <- metadat %>% filter(name %in% metadat_temp$name)

inputDat_heatmap <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_heatmap$name) %>%
  column_to_rownames("name")

heatmap_temp <- metadat_heatmap %>% select(name, group, H1N1_reclassify) %>%
  left_join(inputDat_heatmap %>% rownames_to_column("name")) %>%
  group_by(group, H1N1_reclassify) %>%
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  mutate(name2 = paste0(group, "_", H1N1_reclassify)) %>% 
  column_to_rownames("name2") %>%
  select(resSigPro) %>% t()

heatmap_temp <- metadat_heatmap %>% select(name, group, H1N1_reclassify) %>%
  left_join(inputDat_heatmap %>% rownames_to_column("name")) %>%
  mutate(H1N1_reclassify = ifelse(H1N1_reclassify == "LL", "LL", "LL_andHH")) %>%
  group_by(group, H1N1_reclassify) %>%
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  mutate(name2 = paste0(group, "_", H1N1_reclassify)) %>% 
  column_to_rownames("name2") %>%
  select(resSigPro) %>% t()

heatmap_proName <- heatmap_temp %>% as.data.frame() %>% 
  rownames_to_column("OlinkID") %>% 
  left_join(cohorts_dat$proteinAnnot_all %>% select(OlinkID, Assay)) %>%
  column_to_rownames("Assay") %>% select(-OlinkID)

heatmap_proName %>% Heatmap()
heatmap_proName %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))
heatmap_proName %>% t() %>% scale() %>% t() %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))


