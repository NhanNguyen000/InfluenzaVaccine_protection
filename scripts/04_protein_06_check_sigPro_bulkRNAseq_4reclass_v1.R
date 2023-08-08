rm(list = ls())
library(ggpubr)

get.limmaRes <- function(metaDat, inputDat) {
  # Aim: identify DE proteins/metabolites correct with sex, age, and reclassify (the interested vaccine response reclassification groups) 
  #by running the linear model (limma package)
  
  # input: metadata - interested participant (in "SampleName" column) with sex, age, and reclassify group, 
  #inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  
  # output: res - limma output
  
  inputDat_temp <- inputDat[, metaDat$SampleName]
  
  if (identical(metaDat$SampleName, colnames(inputDat_temp)) == TRUE) {
    res <- lmFit(inputDat_temp,
                 design = model.matrix(~ + sex + age + reclassify, metaDat)) %>% 
      eBayes()
  } else res <- "Error: check input"
  
  return(res)
}

get.limmaRes_perStrain <- function(metadat, inputDat, strain_groups) {
  # Aim: run the linear model (using get.limaRes function) with for multiple strain
  
  # input: metadata - interested participant (in "name" column) with sex, age, and reclassify group, 
  # inputDat - data with participant in colnames and protein/metabolites in rownames (need to select only the intested participants)
  # strain_groups - strains which the test apply to
  
  # output: res - list of outcome from limmar model per strains
  
  resList <- list()
  for (strain_group in strain_groups) {
    metadat_temp <- metadat %>% rename("reclassify" := strain_group) %>%
      select(SampleName, sex, age, reclassify) %>% drop_na()
    
    resList[[strain_group]] <- get.limmaRes(metadat_temp, inputDat)
  }
  return(resList)
}

# load data =======================================================================

## load data from iMED transcriptome ---------------
load("/vol/projects/CIIM/Influenza/iMED/transcriptomic/transcriptomeDF.RData")
# View(transcriptomeList[[1]])
# View(transcriptomeList[[2]])
iMED_transcrip_T1 <- transcriptomeList[[2]] %>% filter(SampleTime == "T1")
iMED_transcripDat <- transcriptomeList[[1]] %>% as.data.frame %>% 
  select(iMED_transcrip_T1$SampleName)


## metadata for all healthy subjects -------------------------
load("cohorts_dat.RData")
metadata_healthy <- cohorts$HAI_all %>% 
  full_join(cohorts$donorInfo_all %>% 
              select(probandID, season, cohort, sex, age, condition)) %>%
  left_join(cohorts$donorSample_all %>% filter(time == "T1")) %>%
  filter(condition == "Healthy") %>%
  mutate_at(vars(contains("reclassify")), ~factor(.x, levels = c("LL", "LH", "HL", "HH")))

## run the limma model -------------------------
iMED_strains <- c("H1N1_reclassify", "H3N2_reclassify", "B_reclassify")

# iMED cohort 2015
inputDat_iMED_2015 <- iMED_transcripDat %>% as.data.frame()

metadat_iMED_2015 <- metadata_healthy %>% filter(season == "2015", cohort == "iMED") %>%
  right_join(iMED_transcrip_T1, by = c("probandID" = "patientID")) %>% 
  arrange(factor(SampleName, levels = colnames(inputDat_iMED_2015)))

resPro_2015_bulkRNAseq <- get.limmaRes_perStrain(metadat =  metadat_iMED_2015,
                                                 inputDat = inputDat_iMED_2015, 
                                                 strain_groups = iMED_strains)

# save data ------------------------------------------------
save(resPro_2015_bulkRNAseq, file = "resPro_2015_bulkRNAseq_4groups.RData")

# heatmap ------------------------------------------------
load("resPro_2015_bulkRNAseq_4groups.RData")
load("selected_DAPs_padj2015.RData")

pval_dat <- resPro_2015_bulkRNAseq %>% 
  lapply(function(x) x$p.value %>% as.data.frame %>% 
           select(matches("reclassify")) %>% 
           rownames_to_column("valName")) %>% 
  bind_rows(.id = "strain") %>% mutate(strain = gsub("_reclassify", "", strain)) %>%
  pivot_longer(matches("reclassify"), names_to = "reclassify", values_to = "p.value")

tstat_dat <- resPro_2015_bulkRNAseq %>% 
  lapply(function(x) x$t %>% as.data.frame %>% 
           select(matches("reclassify")) %>% 
           rownames_to_column("valName")) %>% 
  bind_rows(.id = "strain") %>% mutate(strain = gsub("_reclassify", "", strain)) %>%
  pivot_longer(matches("reclassify"), names_to = "reclassify", values_to = "relative_diff")

plotDat_DAPs_bulkRNAseq <- tstat_dat %>% 
  left_join(pval_dat) %>% 
  mutate(group = paste0(strain, "_", reclassify), 
         group = gsub("reclassify", "LLvs", group)) %>% 
  filter(valName %in% selected_DAs)

plotDat_DAPs_bulkRNAseq %>%
  ggplot(aes(x = group, y = valName, fill = relative_diff)) + 
  geom_tile() +
 # geom_text(aes(label = ifelse(is.na(p.value), NA, "*"))) +
  scale_fill_gradientn(limits = c(-5, 5), colors = c("blue", "white", "red"), na.value = "grey") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  + 
  ggtitle("Bulk RNAseq")

# boxplot for selected proteins at baseline------------------------------------------
# iMED cohort 2015
protein <- "CD83"

plotDat_boxplot <- metadat_iMED_2015 %>% 
  full_join(inputDat_iMED_2015 %>% t() %>% 
              as.data.frame %>% select(protein) %>% rownames_to_column("SampleName")) %>%
  mutate(H1N1_abFC = ifelse(H1N1_reclassify == "LH" |H1N1_reclassify == "HH", "R", "NR")) %>%
  mutate(H1N1_abBaseline = ifelse(H1N1_T1 > 40, "high", "low")) %>%
  mutate(H1N1_abBaseline = factor(H1N1_abBaseline, levels = c("low", "high")))


plotDat_boxplot %>% count(H1N1_reclassify)
plotDat_boxplot %>% count(H3N2_reclassify)
plotDat_boxplot %>% count(B_reclassify)

# compare re-classfiy groups
compare_responder <- list( c("NR", "Other"), c("Other", "TR"), c("NR", "TR") )
compare_reClass <- list( c("LL", "LH"), c("LL", "HL"), c("LL", "HH"))
compare_abFC <- list(c("R", "NR"))
compare_abT1 <- list(c("low", "high"))

# based on reclassification 
cowplot::plot_grid(
  plotDat_boxplot %>%
    ggboxplot(x = "H1N1_reclassify", y = protein,
              paletter = "jco", add = "jitter") + 
    stat_compare_means(comparisons = compare_reClass, method = "t.test"),
  plotDat_boxplot %>%
    ggboxplot(x = "H3N2_reclassify", y = protein,
              paletter = "jco", add = "jitter") + 
    stat_compare_means(comparisons = compare_reClass, method = "t.test"),
  plotDat_boxplot %>%
    ggboxplot(x = "B_reclassify", y = protein,
              paletter = "jco", add = "jitter") + 
    stat_compare_means(comparisons = compare_reClass, method = "t.test"),
  nrow = 1
)

# based on abFC: NR vs R
ggboxplot(plotDat_boxplot, x = "H1N1_abFC", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_abFC, method = "t.test")

# based on abT1: low Ab vs high Ab
ggboxplot(plotDat_boxplot, x = "H1N1_abBaseline", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_abT1, method = "t.test")

# based on responders group: NR, other, TR
ggboxplot(plotDat_boxplot, x = "responder", y = protein,
          paletter = "jco", add = "jitter") + facet_wrap(~season, nrow = 1) + 
  stat_compare_means(comparisons = compare_responder, method = "t.test")

# boxplot for selected proteins  overtime ------------------------------------------

iMED_transcrip_overtime <- transcriptomeList[[2]] %>%
  filter(SampleTime %in% c("T1", "T4"))

iMED_transcripDat <- transcriptomeList[[1]] %>% as.data.frame %>% 
  select(iMED_transcrip_overtime$SampleName)

metadat_iMED_2015_overtime <- metadata_healthy %>% select(-time) %>%
  filter(season == "2015", cohort == "iMED") %>%
  right_join(iMED_transcrip_overtime, by = c("probandID" = "patientID")) %>% 
  rename("time" = "SampleTime")


## plot protein expression overtime -----------------
protein <- "CD83"
plotDat_boxplot_overtime <- metadat_iMED_2015_overtime %>% 
  full_join(iMED_transcripDat %>% t() %>% 
              as.data.frame %>% select(protein) %>% rownames_to_column("SampleName")) 
# plot fo H1N1
stat.test <- plotDat_boxplot_overtime %>% 
  group_by(season, H1N1_reclassify) %>%  
  t_test(CD83 ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(H1N1_reclassify == "LL", 0.9, 
                       ifelse(H1N1_reclassify == "LH", 1.9,
                              ifelse(H1N1_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(H1N1_reclassify == "LL", 1.1, 
                       ifelse(H1N1_reclassify == "LH", 2.1,
                              ifelse(H1N1_reclassify == "HL", 3.1, 4.1))))

bxp_H1N1 <- ggboxplot(plotDat_boxplot_overtime, x = "H1N1_reclassify", y = protein, color = "time",
                 paletter = "jco", add = "jitter") + stat_pvalue_manual(stat.test, label = "p.signif")

# plot fo H3N2
stat.test <- plotDat_boxplot_overtime %>% 
  group_by(season, H3N2_reclassify) %>%  
  t_test(CD83 ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(H3N2_reclassify == "LL", 0.9, 
                       ifelse(H3N2_reclassify == "LH", 1.9,
                              ifelse(H3N2_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(H3N2_reclassify == "LL", 1.1, 
                       ifelse(H3N2_reclassify == "LH", 2.1,
                              ifelse(H3N2_reclassify == "HL", 3.1, 4.1))))

bxp_H3N2 <- ggboxplot(plotDat_boxplot_overtime, x = "H3N2_reclassify", y = protein, color = "time",
                 paletter = "jco", add = "jitter") + stat_pvalue_manual(stat.test, label = "p.signif")

# plot fo B
stat.test <- plotDat_boxplot_overtime %>% 
  group_by(season, B_reclassify) %>%  
  t_test(CD83 ~ time) %>% 
  add_xy_position() %>% add_significance() %>%
  mutate(xmin = ifelse(B_reclassify == "LL", 0.9, 
                       ifelse(B_reclassify == "LH", 1.9,
                              ifelse(B_reclassify == "HL", 2.9, 3.9))),
         xmax = ifelse(B_reclassify == "LL", 1.1, 
                       ifelse(B_reclassify == "LH", 2.1,
                              ifelse(B_reclassify == "HL", 3.1, 4.1))))

bxp_B <- ggboxplot(plotDat_boxplot_overtime, x = "B_reclassify", y = protein, color = "time",
                 paletter = "jco", add = "jitter") + stat_pvalue_manual(stat.test, label = "p.signif")

# all plot together
cowplot::plot_grid(bxp_H1N1, bxp_H3N2, bxp_B, nrow = 1)

# cor(CD83, other proteins) at transcriptome level -------------------------
bulkRNAseq_baseline <- iMED_transcripDat %>% t() %>% as.data.frame

CD83_relatedPro <- c("CD1A", "CD1B", "CD1C", "CD1D", "CD1E", "CD40", "CD80", "CD86", "TNF", "CCR7") # based on stringDB

outcome <- matrix(NA, nrow = length(CD83_relatedPro), ncol = 2) %>% as.data.frame
names(outcome) <- c("p.value", "cor.estimate")
rownames(outcome) <- CD83_relatedPro

for (i in 1:length(CD83_relatedPro)) {
  outcome$p.value[i] <- cor.test(bulkRNAseq_baseline$CD83, bulkRNAseq_baseline[[CD83_relatedPro[i]]])$p.value
  outcome$cor.estimate[i] <- cor.test(bulkRNAseq_baseline$CD83, bulkRNAseq_baseline[[CD83_relatedPro[i]]])$estimate
}

outcome %>% filter(p.value < 0.05)

