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

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t() %>%
  as.data.frame() %>% rownames_to_column("OlinkID") %>% 
  left_join(cohorts_dat$proteinAnnot_all) %>% select(-OlinkID, -UniProt) %>%
  column_to_rownames("Assay")

# run limma model for LL vs. protector (LH, HL, HH), and LL vs. LH --------------------------
strain_groups <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")

res <- list()
res_4groups <- list()
for (strain_group in strain_groups) {
  # LL vs LH, HL, HH together
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na()
  
  res[[strain_group]]$LL_vsOthers <- get.limmaRes(metadat_temp, inputDat)

  # LL vs LH, HL, HH per group
  reclass <- c("LH", "HL", "HH")
  for (i in reclass) {
    metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      filter(reclassify %in% c("LL", i)) %>%
      mutate(reclassify = factor(reclassify, levels = c("LL", i)))
    
    res_4groups[[strain_group]]$LL_vsOthers[[i]] <- get.limmaRes(metadat_temp, inputDat)
  }

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
             as.data.frame() %>% select(i) %>% filter(. <0.05))
}

resSig_LLvsOther_perGroup <- res_4groups %>%
  lapply(function(x) x$LL_vsOther %>% 
           lapply(function(y) y$p.value %>% 
                    as.data.frame() %>% select(4) %>% filter(. < 0.05)))

resSig_Protector <- res %>% 
  lapply(function(x) x$LL_vsProtectors$p.value %>% 
           as.data.frame() %>%  select(4) %>% filter(. <0.05))

resSig_LLvsLH_HH <- res %>%
  lapply(function(x) x$LL_vsLH_HH$p.value %>%
           as.data.frame() %>% select(4) %>% filter(. < 0.05))

resSig_HLvsLH_HH <- res %>%
  lapply(function(x) x$HL_vsLH_HH $p.value %>%
           as.data.frame() %>% select(4) %>% filter(. < 0.05))

resSig_HLvsHH <- res %>%
  lapply(function(x) x$HL_vsHH$p.value %>%
           as.data.frame() %>% select(4) %>% filter(. < 0.05))
# vocanol plot ----------------------

volcanoplot(res$H1N1_reclassify$LL_vsProtectors, coef= "reclassifyProtector", 
            names = rownames(res$H1N1_reclassify$LL_vsProtectors), 
            highlight = 25, pch = 10, cex = 0.25)

volcanoplot(res$H3N2_reclassify$LL_vsProtectors, coef= "reclassifyProtector", 
            names = rownames(res$H1N1_reclassify$LL_vsProtectors), 
            highlight = 19, pch = 10, cex = 0.25)

volcanoplot(res$H3N2_reclassify$LL_vsProtectors, coef= "reclassifyProtector", 
            names = rownames(res$H1N1_reclassify$LL_vsProtectors), 
            highlight = 42, pch = 10, cex = 0.25)

# other volcano plot
library(EnhancedVolcano)
resTable <- topTable(res$H1N1_reclassify$LL_vsProtectors, coef = "reclassifyProtector", n = Inf)
resTable <- topTable(res$B_reclassify$LL_vsProtectors, coef = "reclassifyProtector", n = Inf)
# test: 
resTable2 <- resTable %>% as.data.frame %>% rownames_to_column("Protein") %>%
  full_join(res$H1N1_reclassify$LL_vsProtectors$coefficients %>% 
              as.data.frame %>% rownames_to_column("Protein")) %>%
  full_join(res$H1N1_reclassify$LL_vsProtectors$p.value %>% 
              as.data.frame %>% rownames_to_column("Protein") %>% 
              select(Protein, reclassifyProtector) %>% rename("pvalue" = 2))
identical(resTable2$logFC, resTable2$reclassifyProtector)
identical(resTable2$P.Value, resTable2$pvalue)

EnhancedVolcano(resTable,
                lab = rownames(resTable),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05, FCcutoff = 0.1,
                xlim = c(-0.5, 0.5), ylim = c(0, 3.5))

EnhancedVolcano(resTable,
                lab = ifelse(rownames(resTable) %in% H1N1_B_proteins, rownames(resTable), NA),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05, FCcutoff = 0,
                xlim = c(-0.5, 0.5), ylim = c(0, 3.5))


# Check venn diagram --------------------
# LL vs. Protector
res.venn <- list(
  "H1N1" = rownames(resSig_Protector$H1N1_reclassify),
  "H3N2" = rownames(resSig_Protector$H3N2_reclassify),
  "B" = rownames(resSig_Protector$B_reclassify))
v.table <- venn(res.venn)

intersect(res.venn$H1N1, res.venn$H3N2)
H1N1_B_proteins <- intersect(res.venn$H1N1, res.venn$B)
intersect(res.venn$H3N2, res.venn$B)

res.venn <- list(
  "H1N1" = rownames(resSig_Protector$H1N1_reclassify),
  "B" = rownames(resSig_Protector$B_reclassify))
v.table <- venn(res.venn)

# check the t-statistic ----------------
resSig_LLvsOther_tstatistic <- res %>% 
  lapply(function(x) x$LL_vsOthers$t %>% 
           as.data.frame() %>% select(matches("reclassify")))

resSig_LLvsOther_perGroup_tstatistic <- res_4groups %>%
  lapply(function(x) x$LL_vsOthers %>% 
           lapply(function(y) y$t %>% as.data.frame() %>% 
                    select(4) %>% rownames_to_column("Assay")) %>% 
           purrr::reduce(full_join) %>% relocate(Assay)
  )

resSig_Protector_tstatistic <- res %>%
  lapply(function(x) x$LL_vsProtectors$t %>% as.data.frame() %>% select(4)) %>%
  imap(~.x %>% rename_with(function(x) paste(x, .y, sep = "_"))) %>%
  lapply(function(x) x %>% rownames_to_column("Assay")) %>% 
  purrr::reduce(full_join)

resSig_LLvsLH_HH_tstatistic <- res %>%
  lapply(function(x) x$LL_vsLH_HH$t %>% as.data.frame() %>% select(4))%>%
  imap(~.x %>% rename_with(function(x) paste(x, .y, sep = "_"))) %>%
  lapply(function(x) x %>% rownames_to_column("Assay")) %>% 
  purrr::reduce(full_join)

resSig_HLvsLH_HH_tstatistic <- res %>%
  lapply(function(x) x$HL_vsLH_HH $t %>% as.data.frame() %>% select(4))%>%
  imap(~.x %>% rename_with(function(x) paste(x, .y, sep = "_"))) %>%
  lapply(function(x) x %>% rownames_to_column("Assay")) %>% 
  purrr::reduce(full_join)

resSig_HLvsHH_tstatistic <- res %>%
  lapply(function(x) x$HL_vsHH$t %>% as.data.frame() %>% select(4))%>%
  imap(~.x %>% rename_with(function(x) paste(x, .y, sep = "_"))) %>%
  lapply(function(x) x %>% rownames_to_column("Assay")) %>% 
  purrr::reduce(full_join)

# heatmap: LL vs. LH, HL, HH ------------- ------------- ------------- -------------
# t-statistic data
tstat_longDat <- resSig_LLvsOther_tstatistic %>% 
  lapply(function(x) x %>% as.data.frame %>% rownames_to_column("Assay") %>%
           gather(key = "compare", value = "tstatistic", -1))

# pick strain
strain <- "H1N1_reclassify"
strain <- "H3N2_reclassify"
strain <- "B_reclassify"

# prepare data
LLvsOther_wide <- resSig_LLvsOther %>%
  lapply(function(x) x[[strain]] %>% rownames_to_column("Assay")) %>% purrr::reduce(full_join)
dim(LLvsOther_wide)

LLvsOther_long <- LLvsOther_wide %>% 
  gather(key = "compare", value = "sig.p.value", -1)

plotDat <- tstat_longDat[[strain]] %>% 
  right_join(LLvsOther_long) %>% 
  mutate(compare = paste0("LLvs", substr(compare, 11, 12))) %>% 
  mutate(compare = factor(compare, levels = c("LLvsLH", "LLvsHL", "LLvsHH")))

# plot heatmap
plotDat %>%
  ggplot(aes(x = compare, y = Assay, fill = tstatistic)) + 
  geom_tile() +
  geom_text(aes(label = ifelse(is.na(sig.p.value), NA, "*"))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw()

# heatmap: LL vs. Protector (LH, HL, HH)--------------
View(resSig_Protector)

# sum up the data --------------------
pval_dat <- list()
for (compareType in names(res$H1N1_reclassify)) {
  pval_dat[[compareType]] <- res %>% 
    lapply(function(x) x[[compareType]]$p.value %>% 
             as.data.frame %>% select(matches("reclassify")) %>% 
             rownames_to_column("Assay") %>%
             gather(key = "group", value = "p.value", -1) %>%
             mutate(compare = compareType)) %>% 
    imap(~mutate(.x, strain = .y)) %>%
    purrr::reduce(full_join)
}

tstat_dat <- list()
for (compareType in names(res$H1N1_reclassify)) {
  tstat_dat[[compareType]] <- res %>% 
    lapply(function(x) x[[compareType]]$t %>% 
             as.data.frame %>% select(matches("reclassify")) %>% 
             rownames_to_column("Assay") %>%
             gather(key = "group", value = "t", -1) %>%
             mutate(compare = compareType)) %>% 
    imap(~mutate(.x, strain = .y)) %>%
    purrr::reduce(full_join)
}

pval_Tab <- pval_dat %>% purrr::reduce(full_join) %>%
  full_join(tstat_dat %>% purrr::reduce(full_join)) %>%
  relocate(p.value, .after = t) %>%
  mutate(group = gsub("reclassify", "", group),
         strain = gsub("_reclassify", "", strain)) %>% 
  mutate(lab_temp = ifelse(compare == "LL_vsOthers", 
                           paste(strain, compare, group, sep = "_"),
                           paste(strain, compare, sep = "_")))

# all strains
DEPs <- pval_Tab %>% filter(p.value < 0.05) %>% 
  select(Assay) %>% unlist() %>% unique()

pval_Tab %>% filter(Assay %in% DEPs)  %>%
  ggplot(aes(x = lab_temp, y = Assay, fill = t)) +
  geom_tile() +
  geom_text(aes(label = ifelse(p.value < 0.05, "*", NA))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# per strain 
strainType <- "H1N1"
strainType <- "H3N2"
strainType <- "B"

DEPs_temp <- pval_Tab %>% filter(p.value < 0.05, strain == strainType) %>%
  select(Assay) %>% unlist() %>% unique()

pval_Tab %>% filter(Assay %in% DEPs_temp, strain == strainType)  %>%
  ggplot(aes(x = lab_temp, y = Assay, fill = t)) +
  geom_tile() +
  geom_text(aes(label = ifelse(p.value < 0.05, "*", NA))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# selected protein and info
DEPs_temp <- intersect(res.venn$H1N1, res.venn$B)
compareType <- c("LL_vsOthers", "LL_vsProtectors")

pval_Tab %>% filter(Assay %in% DEPs_temp, 
                    compare %in% compareType)  %>%
  ggplot(aes(x = lab_temp, y = Assay, fill = t)) +
  geom_tile() +
  geom_text(aes(label = ifelse(p.value < 0.05, "*", NA))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

