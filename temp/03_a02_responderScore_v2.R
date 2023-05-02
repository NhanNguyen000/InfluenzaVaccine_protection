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
  full_join(HAIreclassify$all_cohorts) %>%
  filter(disease == "healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") %>% t()

# run limma model --------------------------
identical(metadat_healthy$name, colnames(inputDat))

res <- list()

res$H1N1 <- lmFit(inputDat, 
                  design =  model.matrix(~ + sex + age + H1N1_reclassify, metadat_healthy)) %>% eBayes() 
res$H3N2 <- lmFit(inputDat, 
                  design =  model.matrix(~ + sex + age + H3N2_reclassify, metadat_healthy)) %>% eBayes() 

metadat_healthy_iMED <- metadat_healthy %>% filter(cohort == "iMED")
inputDat_iMED <- inputDat %>% as.data.frame %>% select(metadat_healthy_iMED$name)
res$B <- lmFit(inputDat_iMED,
               design = model.matrix(~ + sex + age + B_reclassify, metadat_healthy_iMED)) %>% eBayes()

metadat_healthy_ZirFlu <- metadat_healthy %>% filter(cohort == "ZirFlu")
inputDat_ZirFlu <- inputDat %>% as.data.frame %>% select(metadat_healthy_ZirFlu$name)
res$Bvictoria <- lmFit(inputDat_ZirFlu,
                       design = model.matrix(~ + sex + age + Bvictoria_reclassify, 
                                             metadat_healthy_ZirFlu)) %>% eBayes()
res$Byamagata <- lmFit(inputDat_ZirFlu,
                       design = model.matrix(~ + sex + age + Byamagata_reclassify,
                                             metadat_healthy_ZirFlu)) %>% eBayes()

# get results
resSig <- list()
resSigPro <- list()
for (strain in names(res)) {
  comparisions <- c("LH", "HL", "HH")
  for(type in comparisions) {
    resSig[[strain]][[paste0("LLvs", type)]] <- res[[strain]]$p.value %>% as.data.frame() %>%
      select(matches(type)) %>% filter(. <0.05)
    resSigPro[[strain]][[paste0("LLvs", type)]] <- rownames(resSig[[strain]][[paste0("LLvs", type)]])
  }
}

a <- unique(unlist(resSigPro$H1N1)) # 33 proteins
b <- unique(unlist(resSigPro$B)) # 68 proteins
proPverlap_H1N1_B <- intersect(a, b) # 20 proteins
get.proteinName(OlinkIDs = proPverlap_H1N1_B, OlinkDat = cohorts_dat$proteinAnnot_all)

# responderScore -------------------------------
# only healthy subjects
inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>% 
  purrr::reduce(full_join) %>%
  filter(name %in% metadat_healthy$name) %>% 
  arrange(factor(name, levels = metadat_healthy$name)) %>%
  column_to_rownames("name") 

# all subjects (healhty + cirrhosis)
metadat <- cohorts_dat$donorInfo_all %>% 
  left_join(cohorts_dat$donorSample_all %>% filter(time == "T1")) %>%
  full_join(cohorts_dat$HAItiter_all) %>%
  full_join(HAIreclassify$all_cohorts) %>%
  mutate(category = factor(category, levels = c("NR", "Other", "TR"))) %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

inputDat <- protein_normDat %>% 
  lapply(function(x) x %>% rownames_to_column("name")) %>%
  purrr::reduce(full_join) %>%
  filter(name %in% metadat$name) %>%
  column_to_rownames("name")

# calculate score
picked_proteins <- proPverlap_H1N1_B

responsePro.Score <- inputDat[, picked_proteins] %>% rowMeans(na.rm = TRUE)

inputDat_Zscore <-  inputDat %>% scale() 
responsePro.Zscore <- inputDat_Zscore[, picked_proteins] %>% rowMeans(na.rm = TRUE)
responsePro.Zscore_abs <- inputDat_Zscore[, picked_proteins] %>% abs() %>% rowMeans(na.rm = TRUE)

library(psych)
responsePro.GeomScore <- inputDat[, picked_proteins] %>% t() %>% geometric.mean(na.rm = TRUE)

# plots --------------------------------------------------------------
# boxplot
metadat_boxplot <- metadat_healthy %>% 
  left_join(as.data.frame(responsePro.Score) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.Zscore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.Zscore_abs) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.GeomScore) %>% rownames_to_column("name"))

my_comparisons <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
my_comparisons_v2 <- list( c("LL", "LH"), c("LL", "HL"), c("HL", "HH"), c("LL", "HH"))

score <- "responsePro.Score"
score <- "responsePro.Zscore"
score <- "responsePro.Zscore_abs"
score <- "responsePro.GeomScore"

# ggboxplot(metadat_boxplot, x = "category", y = score,
#           paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
#   # stat_compare_means(comparisons = my_comparisons, method = "t.test")
#   stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H3N2_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~age_group) +
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

metadat_boxplot %>% filter(cohort == "iMED") %>%
  ggboxplot(x = "B_reclassify", y = score,
            paletter = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "Bvictoria_reclassify", y = score,
            paletter = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

metadat_boxplot %>% filter(cohort == "ZirFlu") %>%
  ggboxplot(x = "Byamagata_reclassify", y = score,
            paletter = "jco", add = "jitter") + 
  # stat_compare_means(comparisons = my_comparisons, method = "t.test")
  stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

# boxplot all_together (healthy & cirrhosis)
metadat_boxplot <- metadat %>%
  mutate(group = ifelse(disease == "healthy", age_group, disease)) %>%
  mutate(group = factor(group, levels = c("old", "young", "cirrhosis"))) %>%
  left_join(as.data.frame(responsePro.Score) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.Zscore) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.Zscore_abs) %>% rownames_to_column("name")) %>%
  left_join(as.data.frame(responsePro.GeomScore) %>% rownames_to_column("name"))

score <- "responsePro.Score"
score <- "responsePro.Zscore"
score <- "responsePro.Zscore_abs"
score <- "responsePro.GeomScore"

ggboxplot(metadat_boxplot, x = "category", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")
# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H1N1_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")
# stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")

ggboxplot(metadat_boxplot, x = "H3N2_reclassify", y = score,
          paletter = "jco", add = "jitter") + facet_wrap(~group) +
  stat_compare_means(comparisons = my_comparisons_v2, method = "t.test")
# stat_compare_means(comparisons = my_comparisons_v2, method = "wilcox.test")
