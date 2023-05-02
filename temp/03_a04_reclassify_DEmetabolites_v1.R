library(tidyverse)
library(limma)

# input data --------------------------
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

# check if the data in the same range & normalization
range(iMED_metaboliteDat)
range(ZirFlu_metaboliteDat)

hist(unlist(iMED_metaboliteDat), xlim = range(9, 24))
hist(unlist(ZirFlu_metaboliteDat), xlim = range(9, 24))

# combine iMED + ZirFlu metabolites --------------------------
inputDat <- iMED_metaboliteDat %>% rownames_to_column("name") %>%
  full_join(ZirFlu_metaboliteDat %>% rownames_to_column("name")) %>%
  column_to_rownames("name") %>% t() %>% as.data.frame() %>%
  select(metadat_healthy$name)

# overlap metabolites
overlapped_metabolites <- intersect(iMED_metabolites$Formula, ZirFlu_metabolites$Formula)

inputDat <- inputDat[overlapped_metabolites, ]

# combine iMED + ZirFlu metabolites (normalization Z-score)--------------------------
inputDat_scaled <- iMED_metaboliteDat %>% scale() %>% as.data.frame() %>% rownames_to_column("name") %>%
  full_join(ZirFlu_metaboliteDat %>% scale() %>% as.data.frame() %>% rownames_to_column("name")) %>%
  column_to_rownames("name") %>% t() %>% as.data.frame() %>%
  select(metadat_healthy$name)
inputDat_scaled <- inputDat_scaled[overlapped_metabolites, ]

# run limma model for LL vs. (LH, HH) --------------------------
strain_groups <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")

res_LLvsLH <- list()
resSig_LLvsLH <- list()

res_LLvsProtector <- list()
resSig_LLvsProtector <- list()

inputDat <- inputDat_scaled
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
  resSig_LLvsLH <- res_LLvsLH[[strain_group]]$p.value %>% 
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

save(res_LLvsLH, res_LLvsProtector, resSig_LLvsLH, resSig_LLvsProtector,
     file = "DEmetabolites_combineCohorts.RData")

load("DEmetabolites_combineCohorts.RData")

# top.table <- topTable(res_LLvsLH$H1N1_reclassify, coef= "reclassifyLH", number = 39) # check again

res_LLvsProtector$B_reclassify %>% decideTests() %>% summary() # apply for padj
res_LLvsLH$B_reclassify %>% decideTests() %>% summary()

resSig_metabolites <- rownames(resSig_LLvsProtector$H1N1_reclassify)

# venn diagram -------------------------------
# LL vs. Protector
res.venn <- list(
  "H1N1" = rownames(resSig_LLvsProtector$H1N1_reclassify), 
  "H3N2" = rownames(resSig_LLvsProtector$H3N2_reclassify), 
  "B" = rownames(resSig_LLvsProtector$B_reclassify))
v.table <- venn(res.venn)

intersect(res.venn$H1N1, res.venn$H3N2)
intersect(res.venn$H1N1, res.venn$B)
intersect(res.venn$H3N2, res.venn$B)

# LL vs. LH
res.venn2 <- list(
  "H1N1" = rownames(resSig_LLvsLH$H1N1_reclassify), 
  "H3N2" = rownames(resSig_LLvsLH$H3N2_reclassify), 
  "B" = rownames(resSig_LLvsLH$B_reclassify))
v.table <- venn(res.venn2)

intersect(res.venn2$H1N1, res.venn2$H3N2)
intersect(res.venn2$H1N1, res.venn2$B)
intersect(res.venn2$H3N2, res.venn2$B)

# check metabolites with protein expression -------------
# metabolites
DEmetabolites <- unique(c(rownames(resSig_LLvsLH$H1N1_reclassify),
                   rownames(resSig_LLvsLH$H3N2_reclassify),
                   rownames(resSig_LLvsLH$B_reclassify)))
DEmetabolites_dat <- inputDat[DEmetabolites,]

# proteins
load("DEPs_combineCohorts.RData")
DEPs <- unique(c(rownames(resSig_LLvsLH$H1N1_reclassify),
                 rownames(resSig_LLvsLH$H3N2_reclassify),
                 rownames(resSig_LLvsLH$B_reclassify)))
inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t()

DEPs_dat <- inputDat[DEPs, ] %>% as.data.frame() %>%
  rownames_to_column("OlinkID") %>% left_join(cohorts_dat$proteinAnnot_all) %>%
  select(-OlinkID, -UniProt) %>% column_to_rownames("Assay")

# check
identical(colnames(DEmetabolites_dat), colnames((DEPs_dat)))


corRes <- list()

corResTab <- as.data.frame(matrix(NA, ncol = 4))
colnames(corResTab) <- c("protein", "metabolites", "p.value", "cor")

for(protein in rownames(DEPs_dat)) {
  for (metabolite in rownames(DEmetabolites_dat)) {
    res_temp <- cor.test(as.numeric(DEPs_dat[protein,]), 
                         as.numeric(DEmetabolites_dat[metabolite, ]))
    corResTab <- rbind(corResTab, 
                       c(protein, metabolite, res_temp$p.value, res_temp$estimate))
    
    corRes[[protein]][[metabolite]] <- res_temp
  }
}
corResTab <- corResTab[-1, ] %>% 
  mutate(cor = as.numeric(cor), 
         p.value = as.numeric(p.value))

corResTab %>% 
  ggplot(aes(x = cor, y = -log10(p.value))) + 
  geom_jitter() + theme_bw()

a <- corResTab %>% filter(abs(cor) >= 0.5)
a1 <- corResTab %>% filter(protein %in% c("CD83", "CXCL8"))

corResTab %>% 
  # filter(protein == "CD83") %>% 
  filter(protein == "CXCL8") %>%
  ggplot(aes(x = cor, y = -log10(p.value))) + 
  geom_jitter() + theme_bw()

b1 <- corResTab %>% filter(protein == "CD83") %>% filter(abs(cor) >= 0.2)
b2 <- corResTab %>% filter(protein == "CXCL8") %>% filter(abs(cor) >= 0.2)

g1 <- iMED_metabolites %>% filter(Formula %in% b1$metabolites) %>% 
  left_join(iMED$metabolite_annot) %>%
  slice(-grep("CHEBI|HMDB", CompoundID))

g2 <- iMED_metabolites %>% filter(Formula %in% b2$metabolites) %>% 
  left_join(iMED$metabolite_annot) %>%
  slice(-grep("CHEBI|HMDB", CompoundID))

metabolites_combineCohorts <- iMED_metabolites %>% 
  filter(Formula %in% overlapped_metabolites) %>% 
  left_join(iMED$metabolite_annot) %>%
  select(Formula, CompoundID) %>%
  slice(-grep("CHEBI|HMDB", CompoundID))

write.table(metabolites_combineCohorts$CompoundID, "metabolites_combineCohorts_KeggID.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(g1$CompoundID, "CD83_cor02_KeggID.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(g2$CompoundID, "CXCL8_cor02_KeggID.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# boxplot check the metabolites level ---------------

# boxplot all_together (healthy & cirrhosis)
metadat_boxplot <- metadat_healthy  %>%
  left_join(inputDat %>% t() %>% as.data.frame() %>% rownames_to_column("name"))

metadat_boxplot_scaled <- metadat_healthy  %>%
  left_join(inputDat_scaled %>% t() %>% as.data.frame() %>% rownames_to_column("name"))

score <- "inflamScore"

my_comparisons_group <- list(c("LL", "LH"), c("HL", "HH"))

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = "C4H8O3",
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons_group ,
                     method = "t.test")

metadat_boxplot %>%
  ggplot(aes(x = C4H8O3, y = H1N1_abFC_combine)) +
  geom_jitter() +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() + theme_bw()

metadat_boxplot %>%
  ggplot(aes(x = C4H8O3, y = H1N1_abFC_combine, col = cohort)) +
  geom_jitter() +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() + 
  theme_bw() + theme(legend.position = "top")

metadat_boxplot_scaled %>%
  ggplot(aes(x = C4H8O3, y = H1N1_abFC_combine)) +
  geom_jitter() +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() + theme_bw()

metadat_boxplot_scaled %>%
  ggplot(aes(x = C4H8O3, y = H1N1_abFC_combine, col = cohort)) +
  geom_jitter() +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() + theme_bw()

# add protein data
inputProDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name))

metadat_boxplot2 <- metadat_boxplot %>% left_join(inputProDat)
metadat_boxplot2 %>%
  ggplot(aes(x = C4H8O3, y = OID20565)) +
  geom_jitter() +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() + theme_bw()

metadat_boxplot2 %>%
  ggplot(aes(x = C4H8O3, y = OID20565, col = cohort)) +
  geom_jitter() +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() + 
  theme_bw() + theme(legend.position = "top")


metadat_boxplot2 %>%
  ggplot(aes(x = OID20565, y = H1N1_abFC_combine)) +
  geom_jitter() +
  geom_smooth(method = "lm", se = FALSE) + stat_cor() + theme_bw()

ggboxplot(metadat_boxplot2, x = "H1N1_reclassify", y = "OID20565",
          paletter = "jco", add = "jitter") + 
  stat_compare_means(comparisons = my_comparisons_group ,
                     method = "t.test")
