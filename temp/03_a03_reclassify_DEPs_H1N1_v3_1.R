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
  mutate(H1N1_reclassify_lowVsHighBaseline = ifelse(H1N1_reclassify == "LL" | H1N1_reclassify == "LH", 
                                                    "lowBaseline", "highBaseline"))

metadat_temp <- metadat_healthy %>% 
  filter(H1N1_reclassify %in% c("LL", "HL")) %>%
  mutate(H1N1_reclassify_lowVsHighBaseline = ifelse(H1N1_reclassify == "LL", 
                                                    "lowBaseline", "highBaseline"))

inputDat_temp <- inputDat[, metadat_temp$name]
identical(metadat_temp$name, colnames(inputDat_temp))

res <- lmFit(inputDat_temp, 
             design =  model.matrix(~ + sex + age + H1N1_reclassify_lowVsHighBaseline, metadat_temp)) %>% 
  eBayes() 

resSig <- res$p.value %>% as.data.frame() %>%  select(4) %>% filter(. <0.05) 
# top.table <- topTable(res, coef= "H1N1_reclassify_lowVsHighBaselinelowBaseline", number = 41) # check again
# res %>% decideTests() %>% summary() # apply for padj, we do not have sig. padj
resSigPro <- rownames(resSig)
get.proteinName(OlinkIDs = resSigPro, OlinkDat = cohorts_dat$proteinAnnot_all)

# check abFC-------------

metadat_temp <- metadat_healthy %>% 
  filter(H1N1_reclassify %in% c("LL", "LH")) %>%
  mutate(H1N1_reclassify_lowVsHighAbFC = ifelse(H1N1_reclassify == "LL", 
                                                    "lowAbFC", "highAbFC"))

metadat_temp <- metadat_healthy %>% 
  filter(H1N1_reclassify %in% c("HL", "HH")) %>%
  mutate(H1N1_reclassify_lowVsHighAbFC = ifelse(H1N1_reclassify == "HL", 
                                                "lowAbFC", "highAbFC"))

inputDat_temp <- inputDat[, metadat_temp$name]
identical(metadat_temp$name, colnames(inputDat_temp))

res <- lmFit(inputDat_temp, 
             design =  model.matrix(~ + sex + age + H1N1_reclassify_lowVsHighAbFC, metadat_temp)) %>% 
  eBayes() 

resSig <- res$p.value %>% as.data.frame() %>%  select(4) %>% filter(. <0.05) 
resSigPro <- rownames(resSig)
get.proteinName(OlinkIDs = resSigPro, OlinkDat = cohorts_dat$proteinAnnot_all)


