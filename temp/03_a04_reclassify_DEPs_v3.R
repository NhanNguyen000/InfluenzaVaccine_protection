library(limma)
library(gplots)

get.limmaRes <- function(metaDat, inputDat) {
  # Aim: identify DE proteins/metabolites correct with sex, age, and reclassify (the interested vaccine response reclassification groups) 
  #by running the linear model (limma package)
  
  # input: metadata - interested participant (in "name" column) with sex, age, and reclassify group, 
  #inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  
  # output: res - limma output
  
  inputDat_temp <- inputDat[, metaDat$name]
  
  if (identical(metaDat$name, colnames(inputDat_temp)) == TRUE) {
    res <- lmFit(inputDat_temp,
                 design = model.matrix(~ + sex + age + reclassify, metaDat)) %>% 
      eBayes()
  } else res <- "Error: check input"
  
  return(res)
}

# input data --------------------------
# metadata for all healthy subjects
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

metadat_healthy <- metadat_healthy %>% filter(cohort == "iMED")
metadat_healthy <- metadat_healthy %>% filter(cohort == "ZirFlu")

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t()

# run limma model per cohort (iMED or ZirFlu) data --------------------------
strain_groups <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")
strain_groups <- "H1N1_reclassify"
res <- list()
for (strain_group in strain_groups) {
  # LL vs LH, HL, HH together
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na()
  
  res[[strain_group]]$LL_vsOthers <- get.limmaRes(metadat_temp, inputDat)
  
  # LL vs (LH, HL, HH) as Protector
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na() %>%
    mutate(reclassify = ifelse(reclassify == "LL", "LL", "Protector")) %>%
    mutate(reclassify = factor(reclassify, levels = c("LL", "Protector")))
  
  res[[strain_group]]$LL_vsProtectors <- get.limmaRes(metadat_temp, inputDat)
  
  # LL vs (LH an HH) together
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na() %>%
    filter(reclassify %in% c("LL", "LH", "HH")) %>%
    mutate(reclassify = ifelse(reclassify == "LL", "LL", "LH_HH")) %>%
    mutate(reclassify = factor(reclassify, levels = c("LL", "LH_HH")))
  
  res[[strain_group]]$LL_vsLH_HH <- get.limmaRes(metadat_temp, inputDat)
  
  # HL vs (LH an HH) together
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na() %>%
    filter(reclassify %in% c("HL", "LH", "HH")) %>%
    mutate(reclassify = ifelse(reclassify == "HL", "HL", "LH_HH")) %>%
    mutate(reclassify = factor(reclassify, levels = c("LL", "HL", "LH_HH")))
  
  res[[strain_group]]$HL_vsLH_HH <- get.limmaRes(metadat_temp, inputDat)
  
  # HL vs HH 
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na() %>%
    filter(reclassify %in% c("HL", "HH")) %>%
    mutate(reclassify = ifelse(reclassify == "HL", "HL", "HH"))  %>%
    mutate(reclassify = factor(reclassify, levels = c("LL", "HL", "HH")))
  
  res[[strain_group]]$HL_vsHH <- get.limmaRes(metadat_temp, inputDat)
  
}
# sigProteins ---------------------------------------------------------
resSig_LLvsOther <- list()
for (i in c("reclassifyLH", "reclassifyHL", "reclassifyHH")) {
  resSig_LLvsOther[[i]] <- res %>%
    lapply(function(x) x$LL_vsOther$p.value %>% 
             as.data.frame() %>% select(i) %>% filter(. <0.05)%>%
             rownames_to_column("OlinkID") %>% left_join(cohorts_dat$proteinAnnot_all) %>%
             select(-OlinkID, -UniProt) %>% relocate(Assay))
}

resSig_LLvsOther_perGroup <- res_4groups %>%
  lapply(function(x) x$LL_vsOther %>% 
           lapply(function(y) y$p.value %>% 
                    as.data.frame() %>% select(4) %>% filter(. < 0.05) %>%
                    rownames_to_column("OlinkID") %>% left_join(cohorts_dat$proteinAnnot_all) %>%
                    select(-OlinkID, -UniProt) %>% relocate(Assay)))

resSig_Protector <- res %>% 
  lapply(function(x) x$LL_vsProtectors$p.value %>% 
           as.data.frame() %>%  select(4) %>% filter(. <0.05)%>%
           rownames_to_column("OlinkID") %>% left_join(cohorts_dat$proteinAnnot_all) %>%
           select(-OlinkID, -UniProt) %>% relocate(Assay))

resSig_LLvsLH_HH <- res %>%
  lapply(function(x) x$LL_vsLH_HH$p.value %>%
           as.data.frame() %>% select(4) %>% filter(. < 0.05)%>%
           rownames_to_column("OlinkID") %>% left_join(cohorts_dat$proteinAnnot_all) %>%
           select(-OlinkID, -UniProt) %>% relocate(Assay))

resSig_HLvsLH_HH <- res %>%
  lapply(function(x) x$HL_vsLH_HH $p.value %>%
           as.data.frame() %>% select(4) %>% filter(. < 0.05)%>%
           rownames_to_column("OlinkID") %>% left_join(cohorts_dat$proteinAnnot_all) %>%
           select(-OlinkID, -UniProt) %>% relocate(Assay))

resSig_HLvsHH <- res %>%
  lapply(function(x) x$HL_vsHH$p.value %>%
           as.data.frame() %>% select(4) %>% filter(. < 0.05)%>%
           rownames_to_column("OlinkID") %>% left_join(cohorts_dat$proteinAnnot_all) %>%
           select(-OlinkID, -UniProt) %>% relocate(Assay))
# Check venn diagram --------------------
# LL vs. Protector
res.venn <- list(
  "H1N1" = resSig_Protector$H1N1_reclassify$Assay,
  "H3N2" = resSig_Protector$H3N2_reclassify$Assay,
  "B" = resSig_Protector$B_reclassify$Assay)
v.table <- venn(res.venn)

intersect(res.venn$H1N1, res.venn$H3N2)
intersect(res.venn$H1N1, res.venn$B)
intersect(res.venn$H3N2, res.venn$B)
