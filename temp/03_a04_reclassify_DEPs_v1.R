library(limma)
library(gplots)
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

# run limma model for LL vs. protector (LH, HL, HH), and LL vs. LH --------------------------
strain_groups <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")

res_LLvsLH <- list()
resSig_LLvsLH <- list()

res_LLvsProtector <- list()
resSig_LLvsProtector <- list()

for (strain_group in strain_groups) {
  # LL vs LH
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    drop_na(reclassify) %>% 
    select(-matches("H1N1|H3N2|B")) %>% 
    filter(reclassify %in% c("LL", "LH")) %>% 
    mutate(reclassify = factor(reclassify, levels = c("LL", "LH")))
  
  inputDat_temp <- inputDat[, metadat_temp$name]
  identical(metadat_temp$name, colnames(inputDat_temp))
  
  res_LLvsLH[[strain_group]] <- lmFit(inputDat_temp, 
               design =  model.matrix(~ + sex + age + reclassify, metadat_temp)) %>% 
    eBayes() 
  resSig_LLvsLH[[strain_group]] <- res_LLvsLH[[strain_group]]$p.value %>% 
    as.data.frame() %>%  select(4) %>% filter(. <0.05) 
  rm(metadat_temp, inputDat_temp)
  
  # LL vs Protector
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    drop_na(reclassify) %>% 
    select(-matches("H1N1|H3N2|B")) %>% 
    mutate(reclassify_LLvsProtector = ifelse(reclassify == "LL", "LL", "Protector"))
  
  inputDat_temp <- inputDat[, metadat_temp$name]
  identical(metadat_temp$name, colnames(inputDat_temp))
  
  res_LLvsProtector[[strain_group]] <- lmFit(inputDat_temp, 
               design =  model.matrix(~ + sex + age + reclassify_LLvsProtector, metadat_temp)) %>% 
    eBayes()
  resSig_LLvsProtector[[strain_group]] <- res_LLvsProtector[[strain_group]]$p.value %>% 
    as.data.frame() %>%  select(4) %>% filter(. <0.05)
  rm(metadat_temp, inputDat_temp)
}

# save(res_LLvsLH, res_LLvsProtector, resSig_LLvsLH, resSig_LLvsProtector,
#      file = "DEPs_combineCohorts.RData")

load("DEPs_combineCohorts.RData")

# top.table <- topTable(res_LLvsLH$H1N1_reclassify, coef= "reclassifyLH", number = 39) # check again
# res_LLvsProtector$B_reclassify %>% decideTests() %>% summary() # apply for padj, we do not have sig. padj for any strain

resSigPro <- rownames(resSig_LLvsProtector$H1N1_reclassify)
get.proteinName(OlinkIDs = resSigPro, OlinkDat = cohorts_dat$proteinAnnot_all)

volcanoplot(res_LLvsLH$H1N1_reclassify, coef= "reclassifyLH", 
            names = rownames(res_LLvsLH$H1N1_reclassify), 
            highlight = 5)
# venn diagram -------------------------------
# LL vs. Protector
res.venn <- list(
  "H1N1" = get.proteinName(OlinkIDs = rownames(resSig_LLvsProtector$H1N1_reclassify), 
                           OlinkDat = cohorts_dat$proteinAnnot_all),
  "H3N2" = get.proteinName(OlinkIDs = rownames(resSig_LLvsProtector$H3N2_reclassify), 
                           OlinkDat = cohorts_dat$proteinAnnot_all),
  "B" = get.proteinName(OlinkIDs = rownames(resSig_LLvsProtector$B_reclassify), 
                        OlinkDat = cohorts_dat$proteinAnnot_all))
v.table <- venn(res.venn)

intersect(res.venn$H1N1, res.venn$H3N2)
intersect(res.venn$H1N1, res.venn$B)
intersect(res.venn$H3N2, res.venn$B)

# LL vs. LH
res.venn2 <- list(
  "H1N1" = get.proteinName(OlinkIDs = rownames(resSig_LLvsLH$H1N1_reclassify), 
                           OlinkDat = cohorts_dat$proteinAnnot_all),
  "H3N2" = get.proteinName(OlinkIDs = rownames(resSig_LLvsLH$H3N2_reclassify), 
                           OlinkDat = cohorts_dat$proteinAnnot_all),
  "B" = get.proteinName(OlinkIDs = rownames(resSig_LLvsLH$B_reclassify), 
                        OlinkDat = cohorts_dat$proteinAnnot_all))
v.table <- venn(res.venn2)

intersect(res.venn2$H1N1, res.venn2$H3N2)
intersect(res.venn2$H1N1, res.venn2$B)
intersect(res.venn2$H3N2, res.venn2$B)

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

# LL vs. Other
strain_group <- "H1N1_reclassify"
strain_group <- "H3N2_reclassify"
strain_group <- "B_reclassify"

metadat_heatmap <- metadat %>% filter(disease == "healthy") %>% 
  rename("reclassify" := strain_group) %>% drop_na(reclassify) %>%
  mutate(class = ifelse(reclassify == "LL", "LL", "Others")) %>% 
  select(name, class)

heatmap_temp <- metadat_heatmap %>% 
  left_join(protein_normDat %>% 
              lapply(function(x) x %>% rownames_to_column("name")) %>%
              purrr::reduce(full_join)) %>%
  group_by(class) %>%
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  column_to_rownames("class") %>%
  select(rownames(resSig_LLvsProtector[[strain_group]])) %>% t()

heatmap_proName <- heatmap_temp %>% as.data.frame() %>% 
  rownames_to_column("OlinkID") %>% 
  left_join(cohorts_dat$proteinAnnot_all %>% select(OlinkID, Assay)) %>%
  column_to_rownames("Assay") %>% select(-OlinkID)

heatmap_proName %>% Heatmap()
heatmap_proName %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))
heatmap_proName %>% t() %>% scale() %>% t() %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))


# LL vs. LH
strain_group <- "H1N1_reclassify"
strain_group <- "H3N2_reclassify"
strain_group <- "B_reclassify"

metadat_heatmap <- metadat %>% filter(disease == "healthy") %>% 
  rename("class" := strain_group) %>% drop_na(class) %>%
  filter(class %in% c("LL", "LH")) %>%
  select(name, class)

heatmap_temp <- metadat_heatmap %>% 
  left_join(protein_normDat %>% 
              lapply(function(x) x %>% rownames_to_column("name")) %>%
              purrr::reduce(full_join)) %>%
  group_by(class) %>%
  summarise_at(vars(-name), funs(mean(., na.rm = TRUE))) %>%
  column_to_rownames("class") %>%
  select(rownames(resSig_LLvsLH[[strain_group]])) %>% t()

heatmap_proName <- heatmap_temp %>% as.data.frame() %>% 
  rownames_to_column("OlinkID") %>% 
  left_join(cohorts_dat$proteinAnnot_all %>% select(OlinkID, Assay)) %>%
  column_to_rownames("Assay") %>% select(-OlinkID)

heatmap_proName %>% Heatmap()
heatmap_proName %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))
heatmap_proName %>% t() %>% scale() %>% t() %>% 
  Heatmap(cluster_columns=FALSE, cluster_rows = FALSE, row_names_max_width = unit(10, "cm"))
