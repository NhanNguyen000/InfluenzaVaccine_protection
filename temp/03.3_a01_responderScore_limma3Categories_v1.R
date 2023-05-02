library(limma)
library(ComplexHeatmap)

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

# run limma model --------------------------
identical(metadat_healthy$name, colnames(inputDat))

designTable <- model.matrix(~ + sex + age + category, metadat_healthy)
res <- lmFit(inputDat, design = designTable) %>% eBayes() 

# check result
# View(res$p.value)
Other_vsNR <-res$p.value %>% as.data.frame() %>% 
  select(categoryOther) %>% filter(categoryOther < 0.05)

TR_vsNR <-res$p.value %>% as.data.frame() %>% 
  select(categoryTR) %>% filter(categoryTR < 0.05)

library(ggvenn)
res.venn <- list("Other_vsNR" = rownames(Other_vsNR),
                 "TR_vsNR" = rownames(TR_vsNR))
ggvenn(res.venn,  fill_color = c("#999999", "#E69F00"))

res.venn2 <- list(
  "Other_vsNR" = get.proteinName(rownames(Other_vsNR), cohorts_dat$proteinAnnot_all),
  "TR_vsNR" = get.proteinName(rownames(TR_vsNR), cohorts_dat$proteinAnnot_all))

## Heatmap -----------------------------------
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
