library(tidyverse)
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

get.limmaRes_multiComparesion <- function(metadat, inputDat, strain_groups) {
  # Aim: run the linear model (using get.limaRes function) with multiple comparison 
  
  # input: metadata - interested participant (in "name" column) with sex, age, and reclassify group, 
  # inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  # strain_groups - strains which the test apply to
  
  # output: res - list of outcome from limmar model compare 2 groups,
  # res_4groups - list of outcome from limmar model compare 4 groups
  
  res <- list()
  res_4groups <- list()
  for (strain_group in strain_groups) {
    # LL vs LH, HL, HH together
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na()
    
    res[[strain_group]]$LL_vsOthers <- get.limmaRes(metadat_temp, inputDat)
    
    # LL vs LH, HL, HH per group
    reclass <- c("LH", "HL", "HH")
    for (i in reclass) {
      metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
        select(name, sex, age, reclassify) %>% drop_na() %>%
        filter(reclassify %in% c("LL", i)) %>%
        mutate(reclassify = factor(reclassify, levels = c("LL", i)))
      
      res_4groups[[strain_group]]$LL_vsOthers[[i]] <- get.limmaRes(metadat_temp, inputDat)
    }
    
    # LL vs (LH, HL, HH) as Protector
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      mutate(reclassify = ifelse(reclassify == "LL", "LL", "Protector")) %>%
      mutate(reclassify = factor(reclassify, levels = c("LL", "Protector")))
    
    res[[strain_group]]$LL_vsProtectors <- get.limmaRes(metadat_temp, inputDat)
    
    # LL vs (LH an HH) together
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      filter(reclassify %in% c("LL", "LH", "HH")) %>%
      mutate(reclassify = ifelse(reclassify == "LL", "LL", "LH_HH")) %>%
      mutate(reclassify = factor(reclassify, levels = c("LL", "LH_HH")))
    
    res[[strain_group]]$LL_vsLH_HH <- get.limmaRes(metadat_temp, inputDat)
    
    # HL vs LH, HL, HH together
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>% 
      mutate(reclassify = factor(reclassify, levels = c("HL", "LL","LH", "HH")))
    
    res[[strain_group]]$HL_vsOthers <- get.limmaRes(metadat_temp, inputDat)
    
    # HL vs (LH an HH) together
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      filter(reclassify %in% c("HL", "LH", "HH")) %>%
      mutate(reclassify = ifelse(reclassify == "HL", "HL", "LH_HH")) %>%
      mutate(reclassify = factor(reclassify, levels = c("HL", "LH_HH")))
    
    res[[strain_group]]$HL_vsLH_HH <- get.limmaRes(metadat_temp, inputDat)
    
    # HL vs HH 
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      filter(reclassify %in% c("HL", "HH")) %>%
      mutate(reclassify = ifelse(reclassify == "HL", "HL", "HH"))  %>%
      mutate(reclassify = factor(reclassify, levels = c("HL", "HH")))
    
    res[[strain_group]]$HL_vsHH <- get.limmaRes(metadat_temp, inputDat)
    
    # (LL and HL as low reponse group) vs. (LH and HH as high response group)
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      mutate(reclassify = substr(reclassify, 2, 2)) %>%
      mutate(reclassify = factor(reclassify, levels = c("L", "H"))) 
    res[[strain_group]]$LvsH_response <- get.limmaRes(metadat_temp, inputDat)
    
    # (LL and LH as low baseline group) vs. (HL and HH as high baseline group)
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      mutate(reclassify = substr(reclassify, 1, 1)) %>%
      mutate(reclassify = factor(reclassify, levels = c("L", "H"))) 
    res[[strain_group]]$LvsH_baseline <- get.limmaRes(metadat_temp, inputDat)
    
  }
  return(list("res" = res, "res_4groups" = res_4groups))
}

# note ---------------------------
# use iMED is a discovery cohort, and ZirFlu as a replicate

# input data -------------------------
# metadata for all healthy subjects
metadat_healthy <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

# iMED metabolites
iMED_metabolites <- iMED$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

identical(iMED_metabolites$ionIdx, as.numeric(colnames(iMED$metabolite)))

iMED_metaboliteDat <- iMED$metabolite
colnames(iMED_metaboliteDat) <- iMED_metabolites$Formula

# ZirFlu metabolites
ZirFlu_metabolites <- ZirFlu$metabolite_annot %>% 
  select(ionIdx, Formula) %>% distinct() %>%
  group_by(ionIdx) %>% summarise(Formula = paste(Formula, collapse = "_"))

identical(ZirFlu_metabolites$ionIdx, as.numeric(colnames(ZirFlu$metabolite_dat)))

ZirFlu_metaboliteDat <- ZirFlu$metabolite_dat
colnames(ZirFlu_metaboliteDat) <- ZirFlu_metabolites$Formula

# input data: overlap metabolites in iMED + ZirFlu --------------------------
# overlap metabolites
overlapped_metabolites <- intersect(iMED_metabolites$Formula, ZirFlu_metabolites$Formula)

metadat_healthy_iMED <- metadat_healthy %>% filter(cohort == "iMED")
metadat_healthy <- metadat_healthy_iMED
inputDat <- iMED_metaboliteDat[, overlapped_metabolites] %>% 
  t() %>% as.data.frame() %>% select(metadat_healthy_iMED$name)

metadat_healthy_ZirFlu <- metadat_healthy %>% filter(cohort == "ZirFlu")

metadat_healthy <- metadat_healthy_ZirFlu
inputDat <- ZirFlu_metaboliteDat[, overlapped_metabolites] %>% 
  t() %>% as.data.frame() %>% select(metadat_healthy_ZirFlu$name)

# run limma model for LL vs. protector (LH, HL, HH), and LL vs. LH --------------------------
strain_groups <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")
strain_groups <- "H1N1_reclassify"

resMeta <- list()
resMeta_4groups <- list()
for (strain_group in strain_groups) {
  # LL vs LH, HL, HH together
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na()
  
  resMeta[[strain_group]]$LL_vsOthers <- get.limmaRes(metadat_temp, inputDat)
  
  # LL vs LH, HL, HH per group
  reclass <- c("LH", "HL", "HH")
  for (i in reclass) {
    metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
      select(name, sex, age, reclassify) %>% drop_na() %>%
      filter(reclassify %in% c("LL", i)) %>%
      mutate(reclassify = factor(reclassify, levels = c("LL", i)))
    
    resMeta_4groups[[strain_group]]$LL_vsOthers[[i]] <- get.limmaRes(metadat_temp, inputDat)
  }
  
  # LL vs (LH, HL, HH) as Protector
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na() %>%
    mutate(reclassify = ifelse(reclassify == "LL", "LL", "Protector")) %>%
    mutate(reclassify = factor(reclassify, levels = c("LL", "Protector")))
  
  resMeta[[strain_group]]$LL_vsProtectors <- get.limmaRes(metadat_temp, inputDat)
  
  # LL vs (LH an HH) together
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na() %>%
    filter(reclassify %in% c("LL", "LH", "HH")) %>%
    mutate(reclassify = ifelse(reclassify == "LL", "LL", "LH_HH")) %>%
    mutate(reclassify = factor(reclassify, levels = c("LL", "LH_HH")))
  
  resMeta[[strain_group]]$LL_vsLH_HH <- get.limmaRes(metadat_temp, inputDat)
  
  # HL vs (LH an HH) together
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na() %>%
    filter(reclassify %in% c("HL", "LH", "HH")) %>%
    mutate(reclassify = ifelse(reclassify == "HL", "HL", "LH_HH")) %>%
    mutate(reclassify = factor(reclassify, levels = c("LL", "HL", "LH_HH")))
  
  resMeta[[strain_group]]$HL_vsLH_HH <- get.limmaRes(metadat_temp, inputDat)
  
  # HL vs HH 
  metadat_temp <- metadat_healthy %>% rename("reclassify" := strain_group) %>%
    select(name, sex, age, reclassify) %>% drop_na() %>%
    filter(reclassify %in% c("HL", "HH")) %>%
    mutate(reclassify = ifelse(reclassify == "HL", "HL", "HH"))  %>%
    mutate(reclassify = factor(reclassify, levels = c("LL", "HL", "HH")))
  
  resMeta[[strain_group]]$HL_vsHH <- get.limmaRes(metadat_temp, inputDat)
  
}

# sum up the data --------------------
pval_dat <- list()
for (compareType in names(resMeta$H1N1_reclassify)) {
  pval_dat[[compareType]] <- resMeta %>% 
    lapply(function(x) x[[compareType]]$p.value %>% 
             as.data.frame %>% select(matches("reclassify")) %>% 
             rownames_to_column("Formula") %>%
             gather(key = "group", value = "p.value", -1) %>%
             mutate(compare = compareType)) %>% 
    imap(~mutate(.x, strain = .y)) %>%
    purrr::reduce(full_join)
}

tstat_dat <- list()
for (compareType in names(resMeta$H1N1_reclassify)) {
  tstat_dat[[compareType]] <- resMeta %>% 
    lapply(function(x) x[[compareType]]$t %>% 
             as.data.frame %>% select(matches("reclassify")) %>% 
             rownames_to_column("Formula") %>%
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
DEMs <- pval_Tab %>% filter(p.value < 0.05) %>% 
  select(Formula) %>% unlist() %>% unique()

pval_Tab %>% filter(Formula %in% DEMs)  %>%
  ggplot(aes(x = lab_temp, y = Formula, fill = t)) +
  geom_tile() +
  geom_text(aes(label = ifelse(p.value < 0.05, "*", NA))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# per strain 
strainType <- "H1N1"
strainType <- "H3N2"
strainType <- "B"

DEMs_temp <- pval_Tab %>% filter(p.value < 0.05, strain == strainType) %>%
  select(Formula) %>% unlist() %>% unique()

pval_Tab %>% filter(Formula %in% DEMs_temp, strain == strainType)  %>%
  ggplot(aes(x = lab_temp, y = Formula, fill = t)) +
  geom_tile() +
  geom_text(aes(label = ifelse(p.value < 0.05, "*", NA))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
# sigMetabolites ---------------------------------------------------------
resMetaSig_Protector <- resMeta %>% 
  lapply(function(x) x$LL_vsProtectors$p.value %>% 
           as.data.frame() %>%  select(4) %>% filter(. <0.05))

resMetaSig_LLvsOther_perGroup <- resMeta_4groups %>%
  lapply(function(x) x$LL_vsOther %>% 
           lapply(function(y) y$p.value %>% 
                    as.data.frame() %>% select(4) %>% filter(. < 0.05)))

resMetaSig_Protector <- resMeta %>% 
  lapply(function(x) x$LL_vsProtectors$p.value %>% 
           as.data.frame() %>%  select(4) %>% filter(. <0.05))

resMetaSig_LLvsLH_HH <- resMeta %>%
  lapply(function(x) x$LL_vsLH_HH$p.value %>%
           as.data.frame() %>% select(4) %>% filter(. < 0.05))

resMetaSig_HLvsLH_HH <- resMeta %>%
  lapply(function(x) x$HL_vsLH_HH $p.value %>%
           as.data.frame() %>% select(4) %>% filter(. < 0.05))

resMetaSig_HLvsHH <- resMeta %>%
  lapply(function(x) x$HL_vsHH$p.value %>%
           as.data.frame() %>% select(4) %>% filter(. < 0.05))

# Check venn diagram --------------------
# LL vs. Protector
res.venn <- list(
  "H1N1" = rownames(resMetaSig_Protector$H1N1_reclassify),
  "H3N2" = rownames(resMetaSig_Protector$H3N2_reclassify),
  "B" = rownames(resMetaSig_Protector$B_reclassify))
v.table <- venn(res.venn)

intersect(res.venn$H1N1, res.venn$H3N2)
H1N1_B_metabolites <- intersect(res.venn$H1N1, res.venn$B)
intersect(res.venn$H3N2, res.venn$B)

H1N1_B_H3N2_metabolites <- intersect(H1N1_B_metabolites, res.venn$H3N2)

res.venn <- list(
  "H1N1" = rownames(resSig_Protector$H1N1_reclassify),
  "B" = rownames(resSig_Protector$B_reclassify))
v.table <- venn(res.venn)

# vocano plot -------------------
volcanoplot(resMeta$H3N2_reclassify$LL_vsProtectors, coef= "reclassifyProtector", 
            names = rownames(resMeta$H3N2_reclassify$LL_vsProtectors), 
            highlight = 20, pch = 10, cex = 0.25)

library(EnhancedVolcano)
resTable <- topTable(resMeta$H1N1_reclassify$LL_vsProtectors, coef = "reclassifyProtector", n = Inf)
resTable <- topTable(resMeta$B_reclassify$LL_vsProtectors, coef = "reclassifyProtector", n = Inf)
resTable <- topTable(resMeta$H3N2_reclassify$LL_vsProtectors, coef = "reclassifyProtector", n = Inf)
# test: 
resTable2 <- resTable %>% as.data.frame %>% rownames_to_column("Formula") %>%
  full_join(resMeta$H1N1_reclassify$LL_vsProtectors$coefficients %>% 
              as.data.frame %>% rownames_to_column("Formula")) %>%
  full_join(resMeta$H1N1_reclassify$LL_vsProtectors$p.value %>% 
              as.data.frame %>% rownames_to_column("Formula") %>% 
              select(Formula, reclassifyProtector) %>% rename("pvalue" = 2))
identical(resTable2$logFC, resTable2$reclassifyProtector)
identical(resTable2$P.Value, resTable2$pvalue)

EnhancedVolcano(resTable,
                lab = rownames(resTable),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05, FCcutoff = 0.1,
                xlim = c(-0.5, 0.5), ylim = c(0, 3.5))

EnhancedVolcano(resTable,
                lab = ifelse(rownames(resTable) %in% H1N1_B_H3N2_metabolites, rownames(resTable), NA),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05, FCcutoff = 0,
                #xlim = c(-0.5, 0.5), ylim = c(0, 3.5)
                xlim = c(-2, 1), ylim = c(0, 5)
                )
EnhancedVolcano(resTable,
                lab = rownames(resTable),
                selectLab = H1N1_B_H3N2_metabolites,
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05, FCcutoff = 0,
                #xlim = c(-0.5, 0.5), ylim = c(0, 3.5)
                xlim = c(-2, 1), ylim = c(0, 5)
)

g <- resTable %>% rownames_to_column("Formula")
# HEATMAP ---------
# sum up the data 
pval_dat <- list()
for (compareType in names(resMeta$H1N1_reclassify)) {
  pval_dat[[compareType]] <- resMeta %>% 
    lapply(function(x) x[[compareType]]$p.value %>% 
             as.data.frame %>% select(matches("reclassify")) %>% 
             rownames_to_column("Formula") %>%
             gather(key = "group", value = "p.value", -1) %>%
             mutate(compare = compareType)) %>% 
    imap(~mutate(.x, strain = .y)) %>%
    purrr::reduce(full_join)
}

tstat_dat <- list()
for (compareType in names(resMeta$H1N1_reclassify)) {
  tstat_dat[[compareType]] <- resMeta %>% 
    lapply(function(x) x[[compareType]]$t %>% 
             as.data.frame %>% select(matches("reclassify")) %>% 
             rownames_to_column("Formula") %>%
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

# selected protein and info
DEMs_temp <- intersect(res.venn$H1N1, res.venn$B)
DEMs_temp <- c("C9H18O", "C7H12N2O3", "C7H10O6", "C4H8O3", "C3H7NO2S", "C12H22O3")
compareType <- c("LL_vsOthers", "LL_vsProtectors")

pval_Tab %>% filter(Formula %in% DEMs_temp, 
                    compare %in% compareType)  %>%
  ggplot(aes(x = lab_temp, y = Formula, fill = t)) +
  geom_tile() +
  geom_text(aes(label = ifelse(p.value < 0.05, "*", NA))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# correlation ------------------------------
inputDat_protein <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy_iMED$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t()

inputDat_protein <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy_ZirFlu$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t()


DEmetabolites_dat <- inputDat[H1N1_B_metabolites,]

DEPs_dat <- inputDat_protein %>% as.data.frame() %>%
  rownames_to_column("OlinkID") %>% left_join(cohorts_dat$proteinAnnot_all) %>%
  select(-OlinkID, -UniProt) %>% filter(Assay %in% H1N1_B_proteins) %>% 
  column_to_rownames("Assay")


# check
identical(colnames(DEmetabolites_dat), colnames(DEPs_dat))


corRes <- list()

corResTab <- as.data.frame(matrix(NA, ncol = 4))
colnames(corResTab) <- c("protein", "metabolites", "p.value", "cor")

for(protein in rownames(DEPs_dat)) {
  for (metabolite in rownames(DEmetabolites_dat)) {
    res_temp <- cor.test(as.numeric(DEPs_dat[protein,]), 
                         as.numeric(DEmetabolites_dat[metabolite, ]))
    corResTab <- rbind(corResTab, 
                       c(protein, metabolite, res_temp$p.value, res_temp$estimate))
    s
    corRes[[protein]][[metabolite]] <- res_temp
  }
}
corResTab <- corResTab[-1, ] %>% 
  mutate(cor = as.numeric(cor), 
         p.value = as.numeric(p.value))

corResTab %>% 
  ggplot(aes(x = cor, y = -log10(p.value))) + 
  geom_jitter() + theme_bw()

# check protein vs. metabolites --------------------
plot_dat <- DEPs_dat %>% t() %>% as.data.frame() %>% rownames_to_column("name") %>%
  full_join(DEmetabolites_dat %>% t() %>% as.data.frame() %>% rownames_to_column("name"))

plot_dat %>%
  ggplot(aes(x = C4H8O3, y = CD83)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme_bw()

plot_dat %>%
  ggplot(aes(x = C12H22O3, y = CD83)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme_bw() + xlim(12.5, 14)

library(ggpubr)
# CHECK CIRRHOSIS ----------
# metadata for all healthy subjects
metadat_cirrhosis <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify2$all_cohorts) %>%
  filter(disease == "cirrhosis") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))


inputDat_protein_cirrhosis <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_cirrhosis$name) %>% 
  arrange(factor(name, levels = metadat_cirrhosis$name)) %>%
  column_to_rownames("name") %>% t()

inputDat_Meta_cirrhosis <- ZirFlu_metaboliteDat %>% 
  t() %>% as.data.frame() %>% select(metadat_cirrhosis$name)

cirrhosis_proDat <- inputDat_protein_cirrhosis %>% as.data.frame() %>%
  rownames_to_column("OlinkID") %>% left_join(cohorts_dat$proteinAnnot_all) %>%
  select(-OlinkID, -UniProt) %>% column_to_rownames("Assay") %>% t()

plot_dat <- cirrhosis_proDat %>% as.data.frame() %>% rownames_to_column("name") %>%
  full_join(inputDat_Meta_cirrhosis  %>% t() %>% as.data.frame() %>% rownames_to_column("name"))

plot_dat %>%
  ggplot(aes(x = C4H8O3, y = CD83)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() +
  theme_bw()

